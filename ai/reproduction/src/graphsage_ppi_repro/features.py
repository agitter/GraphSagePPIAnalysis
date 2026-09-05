"""Reconstruct the 50 GraphSAGE PPI feature columns from MSigDB GMT files."""

from __future__ import annotations

import csv
import hashlib
import os
import struct
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .provenance import sha256_file, write_json_atomic, write_tsv_atomic
from .topology import read_node_mapping

FEATURE_SPEC_FIELDS = (
    "column_index_0based",
    "collection",
    "source_file",
    "source_row_1based",
    "canonical_name",
    "historical_alias",
    "membership_sha256",
    "source_unique_entrez_count",
    "expected_positive_rows",
    "all_zero_after_projection",
    "evidence_status",
    "evidence_reference",
)


class FeatureError(RuntimeError):
    """Raised when MSigDB selection or projection violates the specification."""


@dataclass(frozen=True)
class GeneSet:
    collection: str
    source_file: str
    source_row_1based: int
    name: str
    description: str
    members: frozenset[int]
    raw_member_fields: int

    @property
    def membership_sha256(self) -> str:
        # A newline-terminated sorted representation is stable across Python
        # versions and independent of hash-table iteration.
        payload = "".join(f"{member}\n" for member in sorted(self.members)).encode("ascii")
        return hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class ExpectedFeature:
    column_index_0based: int
    collection: str
    source_file: str
    source_row_1based: int
    canonical_name: str
    historical_alias: str
    membership_sha256: str
    source_unique_entrez_count: int
    expected_positive_rows: int
    all_zero_after_projection: bool
    evidence_status: str
    evidence_reference: str


def parse_gmt(path: Path, *, collection: str) -> list[GeneSet]:
    """Parse an Entrez-ID GMT file without changing deposited row order."""

    records: list[GeneSet] = []
    with path.open(encoding="utf-8") as handle:
        for source_row_1based, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise FeatureError(f"Malformed GMT row at {path}:{source_row_1based}")
            members: set[int] = set()
            for field in fields[2:]:
                try:
                    members.add(int(field))
                except ValueError as exc:
                    raise FeatureError(
                        f"Non-integer Entrez field {field!r} at {path}:{source_row_1based}"
                    ) from exc
            records.append(
                GeneSet(
                    collection=collection,
                    source_file=path.name,
                    source_row_1based=source_row_1based,
                    name=fields[0],
                    description=fields[1],
                    members=frozenset(members),
                    raw_member_fields=len(fields) - 2,
                )
            )
    if not records:
        raise FeatureError(f"No GMT records found in {path}")
    return records


def select_gene_sets(
    collections: Sequence[tuple[str, Path]],
    *,
    minimum_source_members: int,
    maximum_columns: int,
) -> tuple[list[GeneSet], dict[str, int]]:
    """Apply the recovered C1-then-C3 source-order selection rule."""

    selected: list[GeneSet] = []
    qualifying_counts: dict[str, int] = {}
    for collection, path in collections:
        rows = parse_gmt(path, collection=collection)
        qualifying = [row for row in rows if len(row.members) >= minimum_source_members]
        qualifying_counts[collection] = len(qualifying)
        for row in qualifying:
            if len(selected) == maximum_columns:
                break
            selected.append(row)
        if len(selected) == maximum_columns:
            break
    if len(selected) != maximum_columns:
        raise FeatureError(
            f"Selection produced {len(selected)} columns; expected {maximum_columns}"
        )
    return selected, qualifying_counts


def read_expected_features(path: Path) -> list[ExpectedFeature]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FEATURE_SPEC_FIELDS:
            raise FeatureError(
                f"Unexpected columns in {path}: {reader.fieldnames!r}; "
                f"expected {list(FEATURE_SPEC_FIELDS)!r}"
            )
        rows = [
            ExpectedFeature(
                column_index_0based=int(row["column_index_0based"]),
                collection=row["collection"],
                source_file=row["source_file"],
                source_row_1based=int(row["source_row_1based"]),
                canonical_name=row["canonical_name"],
                historical_alias=row["historical_alias"],
                membership_sha256=row["membership_sha256"],
                source_unique_entrez_count=int(row["source_unique_entrez_count"]),
                expected_positive_rows=int(row["expected_positive_rows"]),
                all_zero_after_projection=row["all_zero_after_projection"] == "1",
                evidence_status=row["evidence_status"],
                evidence_reference=row["evidence_reference"],
            )
            for row in reader
        ]
    if [row.column_index_0based for row in rows] != list(range(len(rows))):
        raise FeatureError(
            "feature_columns.tsv indices are not consecutive and file-ordered"
        )
    return rows


def _npy_header(array: np.ndarray) -> bytes:
    """Build the NumPy 1.0 header style used by the 2017 deposited file.

    Modern NumPy aligns headers to 64 bytes, while the deposited GraphSAGE file
    used the older 16-byte alignment.  The data array is identical either way;
    writing the historical header also permits byte-level validation of the
    individual ``ppi-feats.npy`` file.
    """

    if array.dtype.hasobject:
        raise FeatureError("Object arrays are not supported")
    dtype = array.dtype
    if dtype.byteorder == "=":
        dtype = dtype.newbyteorder("<" if np.little_endian else ">")
    descr = np.lib.format.dtype_to_descr(dtype)
    shape = tuple(int(value) for value in array.shape)
    header_text = (
        "{'descr': "
        + repr(descr)
        + ", 'fortran_order': "
        + ("True" if np.isfortran(array) and not array.flags.c_contiguous else "False")
        + ", 'shape': "
        + repr(shape)
        + ", }"
    )
    prefix_size = len(np.lib.format.magic(1, 0)) + 2
    padding = (-(prefix_size + len(header_text.encode("latin1")) + 1)) % 16
    header = (header_text + " " * padding + "\n").encode("latin1")
    if len(header) > 65_535:
        raise FeatureError("NumPy 1.0 header exceeds uint16 length")
    return np.lib.format.magic(1, 0) + struct.pack("<H", len(header)) + header


def write_npy_v1_16byte_aligned(path: Path, array: np.ndarray) -> None:
    """Atomically write a non-object NumPy array with the historical header."""

    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(file_descriptor)
    temporary = Path(temporary_name)
    try:
        contiguous = np.ascontiguousarray(array)
        with temporary.open("wb") as handle:
            handle.write(_npy_header(contiguous))
            handle.write(contiguous.tobytes(order="C"))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _validate_selected_against_spec(
    selected: Sequence[GeneSet],
    expected: Sequence[ExpectedFeature],
    positive_rows: Sequence[int],
) -> None:
    if len(selected) != len(expected):
        raise FeatureError(
            f"Derived {len(selected)} features but feature_columns.tsv has {len(expected)}"
        )
    for index, (observed, wanted) in enumerate(zip(selected, expected, strict=True)):
        comparisons = {
            "column index": index == wanted.column_index_0based,
            "collection": observed.collection == wanted.collection,
            "source file": observed.source_file == wanted.source_file,
            "source row": observed.source_row_1based == wanted.source_row_1based,
            "canonical name": observed.name == wanted.canonical_name,
            "membership SHA-256": observed.membership_sha256 == wanted.membership_sha256,
            "source member count": (
                len(observed.members) == wanted.source_unique_entrez_count
            ),
            "projected positive rows": (
                positive_rows[index] == wanted.expected_positive_rows
            ),
            "all-zero status": (
                (positive_rows[index] == 0) == wanted.all_zero_after_projection
            ),
        }
        failed = [name for name, passed in comparisons.items() if not passed]
        if failed:
            raise FeatureError(
                f"Derived feature column {index} disagrees with the accepted "
                "specification: "
                + ", ".join(failed)
            )


def reconstruct_features(
    *,
    mapping_path: Path,
    c1_path: Path,
    c3_path: Path,
    feature_spec_path: Path,
    matrix_output: Path,
    selected_output: Path,
    summary_output: Path,
    minimum_source_members: int = 200,
    maximum_columns: int = 50,
) -> dict[str, object]:
    """Derive, project, and verify the complete 56,944 x 50 feature matrix."""

    mapping = read_node_mapping(mapping_path)
    node_ids = [int(row["graphsage_node_id"]) for row in mapping]
    if node_ids != list(range(len(mapping))):
        raise FeatureError("Node mapping is not in consecutive GraphSAGE row order")
    genes = np.fromiter(
        (int(row["entrez_gene_id"]) for row in mapping), dtype=np.int64, count=len(mapping)
    )
    selected, qualifying_counts = select_gene_sets(
        (("C1", c1_path), ("C3", c3_path)),
        minimum_source_members=minimum_source_members,
        maximum_columns=maximum_columns,
    )

    matrix = np.zeros((len(genes), len(selected)), dtype=np.float64)
    for column, gene_set in enumerate(selected):
        matrix[:, column] = np.fromiter(
            (int(gene) in gene_set.members for gene in genes),
            dtype=np.float64,
            count=len(genes),
        )
    positive_rows = [int(matrix[:, column].sum()) for column in range(matrix.shape[1])]
    expected = read_expected_features(feature_spec_path)
    _validate_selected_against_spec(selected, expected, positive_rows)

    selected_rows = []
    for column, (gene_set, expected_row) in enumerate(zip(selected, expected, strict=True)):
        selected_rows.append(
            {
                "column_index_0based": column,
                "collection": gene_set.collection,
                "source_file": gene_set.source_file,
                "source_row_1based": gene_set.source_row_1based,
                "canonical_name": gene_set.name,
                "historical_alias": expected_row.historical_alias,
                "description": gene_set.description,
                "source_unique_entrez_count": len(gene_set.members),
                "raw_member_fields": gene_set.raw_member_fields,
                "membership_sha256": gene_set.membership_sha256,
                "projected_positive_rows": positive_rows[column],
                "all_zero_after_projection": int(positive_rows[column] == 0),
            }
        )
    write_tsv_atomic(selected_output, tuple(selected_rows[0]), selected_rows)
    write_npy_v1_16byte_aligned(matrix_output, matrix)

    float64_data_hash = hashlib.sha256(matrix.tobytes(order="C")).hexdigest()
    uint8_data_hash = hashlib.sha256(matrix.astype(np.uint8).tobytes(order="C")).hexdigest()
    collection_counts = {
        collection: sum(row.collection == collection for row in selected)
        for collection in ("C1", "C3")
    }
    summary: dict[str, object] = {
        "schema_version": 1,
        "algorithm": {
            "collection_order": ["C1", "C3"],
            "minimum_source_members": minimum_source_members,
            "threshold_operator": ">=",
            "preserve_gmt_source_order": True,
            "global_column_cap": maximum_columns,
            "project_to_graphsage_rows_after_selection": True,
        },
        "inputs": {
            "node_mapping": str(mapping_path.resolve()),
            "node_mapping_sha256": sha256_file(mapping_path),
            "c1_gmt": str(c1_path.resolve()),
            "c1_gmt_sha256": sha256_file(c1_path),
            "c3_gmt": str(c3_path.resolve()),
            "c3_gmt_sha256": sha256_file(c3_path),
            "feature_columns_spec": str(feature_spec_path.resolve()),
            "feature_columns_spec_sha256": sha256_file(feature_spec_path),
        },
        "counts": {
            "rows": int(matrix.shape[0]),
            "columns": int(matrix.shape[1]),
            "cells": int(matrix.size),
            "selected_by_collection": collection_counts,
            "qualifying_by_collection_before_global_cap": qualifying_counts,
            "all_zero_columns_0based": [
                index for index, count in enumerate(positive_rows) if count == 0
            ],
        },
        "hashes": {
            "float64_c_order_data_sha256": float64_data_hash,
            "uint8_c_order_data_sha256": uint8_data_hash,
            "npy_file_sha256": sha256_file(matrix_output),
        },
        "outputs": {
            "matrix": str(matrix_output.resolve()),
            "selected_features": str(selected_output.resolve()),
            "selected_features_sha256": sha256_file(selected_output),
        },
    }
    write_json_atomic(summary_output, summary)
    return summary

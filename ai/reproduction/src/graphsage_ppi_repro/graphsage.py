"""Assemble the deterministic files consumed by supervised GraphSAGE.

The scientific reconstruction stages deliberately emit transparent intermediate
files: an explicit node-to-gene table, an explicit edge table, and NumPy feature
and label matrices.  This module converts those intermediates into the four
files expected by the original GraphSAGE loader.

Three files can be reproduced byte-for-byte because their serialization is
fully determined by recovered behavior: ``ppi-id_map.json``,
``ppi-class_map.json``, and ``ppi-feats.npy``.  JSON graph link order is not a
scientific property, so ``ppi-G.json`` uses a documented canonical ordering of
undirected edges and is validated structurally against the release.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import TextIO

import numpy as np

from .legacy_order import ordered_string_keys
from .provenance import sha256_file, write_json_atomic, write_text_atomic
from .topology import read_node_mapping

GRAPH_FILENAME = "ppi-G.json"
ID_MAP_FILENAME = "ppi-id_map.json"
CLASS_MAP_FILENAME = "ppi-class_map.json"
FEATURE_FILENAME = "ppi-feats.npy"
CHECKSUM_FILENAME = "SHA256SUMS"
VALID_SPLITS = frozenset({"train", "validation", "test"})


class GraphSAGEError(RuntimeError):
    """Raised when reconstructed intermediates cannot form a valid dataset."""


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _temporary_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    return Path(name)


def _replace_after_write(temporary: Path, destination: Path) -> None:
    try:
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _read_canonical_edges(path: Path, row_count: int) -> list[tuple[int, int]]:
    """Read, validate, canonicalize, and sort unique undirected edges."""

    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"source_node_id", "target_node_id"}
        if not required.issubset(reader.fieldnames or ()):
            raise GraphSAGEError(
                f"Edge table {path} lacks required columns {sorted(required)}"
            )

        edges: list[tuple[int, int]] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                source = int(row["source_node_id"])
                target = int(row["target_node_id"])
            except (KeyError, TypeError, ValueError) as exc:
                raise GraphSAGEError(
                    f"Invalid edge identifiers at {path}:{line_number}"
                ) from exc
            if not 0 <= source < row_count or not 0 <= target < row_count:
                raise GraphSAGEError(
                    f"Out-of-range edge ({source}, {target}) at {path}:{line_number}"
                )
            edges.append((source, target) if source <= target else (target, source))

    unique_edges = set(edges)
    if len(unique_edges) != len(edges):
        raise GraphSAGEError(
            f"Edge table {path} contains duplicate undirected edge records"
        )
    if not edges:
        raise GraphSAGEError(f"Edge table {path} is empty")
    return sorted(unique_edges)


def _validate_mapping(
    rows: Sequence[dict[str, str]],
) -> tuple[list[str], list[int]]:
    """Validate row alignment and return split and graph assignments."""

    if not rows:
        raise GraphSAGEError("Node mapping is empty")
    splits: list[str] = []
    graph_indices: list[int] = []
    split_by_graph: dict[int, str] = {}
    previous_graph_index = 0

    for expected_id, row in enumerate(rows):
        try:
            node_id = int(row["graphsage_node_id"])
            feature_row = int(row["feature_label_row_index"])
            graph_index = int(row["graph_index_1based"])
            split = row["split"]
        except (KeyError, TypeError, ValueError) as exc:
            raise GraphSAGEError(
                f"Invalid node-mapping row at zero-based position {expected_id}"
            ) from exc
        if node_id != expected_id or feature_row != expected_id:
            raise GraphSAGEError(
                "Node mapping must be consecutive and aligned with feature/label "
                f"rows; position {expected_id} contains node={node_id}, "
                f"row={feature_row}"
            )
        if graph_index <= 0 or graph_index < previous_graph_index:
            raise GraphSAGEError(
                "Graph indices must be positive and grouped in ascending order; "
                f"node {node_id} has graph index {graph_index}"
            )
        if split not in VALID_SPLITS:
            raise GraphSAGEError(f"Unknown split {split!r} for GraphSAGE node {node_id}")
        prior_split = split_by_graph.setdefault(graph_index, split)
        if prior_split != split:
            raise GraphSAGEError(
                f"Graph {graph_index} contains both {prior_split!r} and {split!r} rows"
            )

        previous_graph_index = graph_index
        splits.append(split)
        graph_indices.append(graph_index)

    observed_graphs = sorted(split_by_graph)
    expected_graphs = list(range(1, len(observed_graphs) + 1))
    if observed_graphs != expected_graphs:
        raise GraphSAGEError(
            f"Graph indices must be consecutive and one-based; observed {observed_graphs}"
        )
    return splits, graph_indices


def _write_graph_json(
    path: Path,
    splits: Sequence[str],
    edges: Iterable[tuple[int, int]],
) -> None:
    """Write deterministic NetworkX node-link JSON without a trailing newline."""

    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write('{"directed": false, "graph": {}, "nodes": [')
            for node_id, split in enumerate(splits):
                if node_id:
                    handle.write(", ")
                node = {
                    "test": split == "test",
                    "id": node_id,
                    "val": split == "validation",
                }
                json.dump(node, handle)
            handle.write('], "links": [')
            for edge_index, (source, target) in enumerate(edges):
                if edge_index:
                    handle.write(", ")
                json.dump({"source": source, "target": target}, handle)
            handle.write('], "multigraph": false}')
        _replace_after_write(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_identity_id_map(path: Path, row_count: int) -> None:
    """Write the released identity-map representation without a trailing newline."""

    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("{")
            for node_id in range(row_count):
                if node_id:
                    handle.write(", ")
                json.dump(str(node_id), handle)
                handle.write(f": {node_id}")
            handle.write("}")
        _replace_after_write(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_legacy_class_map(
    path: Path,
    labels: np.ndarray,
    *,
    word_size_bits: int,
) -> list[str]:
    """Write labels in recovered CPython-2 string-dictionary key order."""

    row_count = int(labels.shape[0])
    ordered_keys = ordered_string_keys(
        (str(node_id) for node_id in range(row_count)),
        word_size_bits=word_size_bits,
    )
    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("{")
            for output_index, key in enumerate(ordered_keys):
                if output_index:
                    handle.write(", ")
                json.dump(key, handle)
                handle.write(": ")
                json.dump(labels[int(key)].tolist(), handle)
            handle.write("}")
        _replace_after_write(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return ordered_keys


def _copy_atomic(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise GraphSAGEError(f"Required reconstructed file is missing: {source}")
    temporary = _temporary_path(destination)
    try:
        shutil.copyfile(source, temporary)
        _replace_after_write(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _hash_lines(lines: Iterator[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("ascii"))
    return digest.hexdigest()


def _write_checksums(path: Path, files: Sequence[Path]) -> None:
    """Write paths relative to the checksum file so ``sha256sum -c`` works."""

    lines = [
        f"{sha256_file(item)}  {item.relative_to(path.parent).as_posix()}" for item in files
    ]
    write_text_atomic(path, "\n".join(lines) + "\n")


def assemble_graphsage(
    *,
    mapping_path: Path,
    edge_path: Path,
    feature_matrix_path: Path,
    label_matrix_path: Path,
    output_directory: Path,
    summary_output: Path,
    word_size_bits: int = 64,
) -> dict[str, object]:
    """Assemble the four deterministic files used by supervised GraphSAGE."""

    mapping = read_node_mapping(mapping_path)
    splits, graph_indices = _validate_mapping(mapping)
    row_count = len(mapping)
    edges = _read_canonical_edges(edge_path, row_count)
    cross_graph_edges = [
        (source, target)
        for source, target in edges
        if graph_indices[source] != graph_indices[target]
    ]
    if cross_graph_edges:
        raise GraphSAGEError(
            "Edge table contains cross-graph records; first observed "
            f"{cross_graph_edges[0]}"
        )

    features = np.load(feature_matrix_path, allow_pickle=False, mmap_mode="r")
    labels = np.load(label_matrix_path, allow_pickle=False, mmap_mode="r")
    if features.ndim != 2 or features.shape[0] != row_count:
        raise GraphSAGEError(
            f"Feature shape {features.shape} does not align with {row_count} rows"
        )
    if labels.ndim != 2 or labels.shape[0] != row_count:
        raise GraphSAGEError(
            f"Label shape {labels.shape} does not align with {row_count} rows"
        )
    if features.dtype != np.dtype(np.float64):
        raise GraphSAGEError(
            f"Expected float64 GraphSAGE features, observed {features.dtype}"
        )
    if labels.dtype != np.dtype(np.uint8):
        raise GraphSAGEError(f"Expected uint8 labels, observed {labels.dtype}")
    if not np.all((features == 0.0) | (features == 1.0)):
        raise GraphSAGEError("GraphSAGE feature matrix is not binary")
    if not np.all((labels == 0) | (labels == 1)):
        raise GraphSAGEError("GraphSAGE label matrix is not binary")

    output_directory.mkdir(parents=True, exist_ok=True)
    graph_path = output_directory / GRAPH_FILENAME
    id_map_path = output_directory / ID_MAP_FILENAME
    class_map_path = output_directory / CLASS_MAP_FILENAME
    feature_path = output_directory / FEATURE_FILENAME
    checksum_path = output_directory.parent / CHECKSUM_FILENAME

    _write_graph_json(graph_path, splits, edges)
    _write_identity_id_map(id_map_path, row_count)
    class_map_keys = _write_legacy_class_map(
        class_map_path,
        labels,
        word_size_bits=word_size_bits,
    )
    _copy_atomic(feature_matrix_path, feature_path)
    generated_files = [graph_path, id_map_path, class_map_path, feature_path]
    _write_checksums(checksum_path, generated_files)

    split_counts = {split: splits.count(split) for split in sorted(set(splits))}
    output_hashes = {item.name: sha256_file(item) for item in generated_files}
    summary: dict[str, object] = {
        "schema_version": 1,
        "scope": "deterministic files consumed by supervised GraphSAGE",
        "algorithm": {
            "graph_format": "NetworkX node-link JSON",
            "graph_node_order": "ascending GraphSAGE node ID",
            "graph_link_order": "lexicographically sorted canonical undirected pairs",
            "historical_graph_link_order_reproduced": False,
            "id_map_order": "ascending numeric node ID",
            "class_map_order": (
                f"{word_size_bits}-bit unrandomized CPython 2.7 string-dict order"
            ),
            "json_whitespace": "Python json default separators",
            "json_trailing_newline": False,
            "feature_handling": "byte-preserving copy of reconstructed ppi-feats.npy",
            "unsupervised_walks_included": False,
        },
        "inputs": {
            "node_mapping": str(mapping_path.resolve()),
            "node_mapping_sha256": sha256_file(mapping_path),
            "edge_table": str(edge_path.resolve()),
            "edge_table_sha256": sha256_file(edge_path),
            "feature_matrix": str(feature_matrix_path.resolve()),
            "feature_matrix_sha256": sha256_file(feature_matrix_path),
            "label_matrix": str(label_matrix_path.resolve()),
            "label_matrix_sha256": sha256_file(label_matrix_path),
        },
        "counts": {
            "graphs": len(set(graph_indices)),
            "rows": row_count,
            "edge_records": len(edges),
            "feature_columns": int(features.shape[1]),
            "label_columns": int(labels.shape[1]),
            "split_rows": split_counts,
        },
        "content_hashes": {
            "canonical_edges_sha256": _hash_lines(
                f"{source}\t{target}\n" for source, target in edges
            ),
            "split_flags_sha256": _hash_lines(
                f"{node_id}\t{split}\n" for node_id, split in enumerate(splits)
            ),
            "class_map_key_order_sha256": _hash_lines(f"{key}\n" for key in class_map_keys),
        },
        "outputs": {
            "directory": str(output_directory.resolve()),
            "files": {
                item.name: {
                    "path": str(item.resolve()),
                    "size_bytes": item.stat().st_size,
                    "sha256": output_hashes[item.name],
                }
                for item in generated_files
            },
            "checksums": str(checksum_path.resolve()),
            "checksums_sha256": sha256_file(checksum_path),
        },
    }
    write_json_atomic(summary_output, summary)
    return summary

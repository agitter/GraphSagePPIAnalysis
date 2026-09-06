"""Tests for source-order MSigDB selection and feature projection."""

from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path

import numpy as np

from graphsage_ppi_repro.features import (
    parse_gmt,
    reconstruct_features,
    select_gene_sets,
    write_npy_v1_16byte_aligned,
)


def _membership_hash(members: set[int]) -> str:
    payload = "".join(f"{member}\n" for member in sorted(members)).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def _write_mapping(path: Path) -> None:
    fields = [
        "graphsage_node_id",
        "feature_label_row_index",
        "graph_index_1based",
        "tissue",
        "split",
        "local_node_index_0based",
        "entrez_gene_id",
        "python2_dict_table_slot",
        "python2_dict_table_size",
    ]
    genes = [1, 2, 3, 1]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for node_id, gene in enumerate(genes):
            writer.writerow(
                {
                    "graphsage_node_id": node_id,
                    "feature_label_row_index": node_id,
                    "graph_index_1based": 1,
                    "tissue": "example",
                    "split": "train",
                    "local_node_index_0based": node_id,
                    "entrez_gene_id": gene,
                    "python2_dict_table_slot": node_id,
                    "python2_dict_table_size": 8,
                }
            )


def _write_feature_spec(
    path: Path, rows: list[tuple[str, int, str, set[int], int]]
) -> None:
    fields = [
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
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, (
            collection,
            source_row,
            name,
            members,
            positive_rows,
        ) in enumerate(rows):
            writer.writerow(
                {
                    "column_index_0based": index,
                    "collection": collection,
                    "source_file": f"{collection.lower()}.gmt",
                    "source_row_1based": source_row,
                    "canonical_name": name,
                    "historical_alias": "",
                    "membership_sha256": _membership_hash(members),
                    "source_unique_entrez_count": len(members),
                    "expected_positive_rows": positive_rows,
                    "all_zero_after_projection": int(positive_rows == 0),
                    "evidence_status": "synthetic",
                    "evidence_reference": "tests/test_features.py",
                }
            )


def test_parse_gmt_deduplicates_members_without_reordering_rows(tmp_path: Path) -> None:
    path = tmp_path / "c1.gmt"
    path.write_text("first\tdesc\t3\t1\t1\nsecond\tdesc\t9\t8\n", encoding="utf-8")
    rows = parse_gmt(path, collection="C1")
    assert [row.name for row in rows] == ["first", "second"]
    assert rows[0].members == frozenset({1, 3})
    assert rows[0].raw_member_fields == 3


def test_source_order_threshold_and_global_cap(tmp_path: Path) -> None:
    c1 = tmp_path / "c1.gmt"
    c3 = tmp_path / "c3.gmt"
    c1.write_text(
        "small\tdesc\t1\nfirst\tdesc\t1\t2\nsecond\tdesc\t2\t3\n",
        encoding="utf-8",
    )
    c3.write_text(
        "third\tdesc\t1\t4\nfourth\tdesc\t2\t4\n",
        encoding="utf-8",
    )
    selected, qualifying = select_gene_sets(
        (("C1", c1), ("C3", c3)),
        minimum_source_members=2,
        maximum_columns=3,
    )
    assert [row.name for row in selected] == ["first", "second", "third"]
    assert qualifying == {"C1": 2, "C3": 2}


def test_full_feature_reconstruction_on_tiny_inputs(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.tsv.gz"
    _write_mapping(mapping)
    c1 = tmp_path / "c1.gmt"
    c3 = tmp_path / "c3.gmt"
    c1.write_text(
        "c1_first\tdesc\t1\t2\nc1_second\tdesc\t2\t3\n",
        encoding="utf-8",
    )
    c3.write_text("c3_first\tdesc\t1\t4\n", encoding="utf-8")
    feature_spec = tmp_path / "feature_columns.tsv"
    _write_feature_spec(
        feature_spec,
        [
            ("C1", 1, "c1_first", {1, 2}, 3),
            ("C1", 2, "c1_second", {2, 3}, 2),
            ("C3", 1, "c3_first", {1, 4}, 2),
        ],
    )

    matrix_path = tmp_path / "features.npy"
    selected_path = tmp_path / "selected.tsv"
    summary_path = tmp_path / "summary.json"
    reconstruct_features(
        mapping_path=mapping,
        c1_path=c1,
        c3_path=c3,
        feature_spec_path=feature_spec,
        matrix_output=matrix_path,
        selected_output=selected_path,
        summary_output=summary_path,
        minimum_source_members=2,
        maximum_columns=3,
    )
    observed = np.load(matrix_path, allow_pickle=False)
    expected = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 0], [1, 0, 1]], dtype=np.float64)
    assert np.array_equal(observed, expected)


def test_historical_npy_writer_is_loadable_and_deterministic(tmp_path: Path) -> None:
    array = np.arange(12, dtype=np.float64).reshape(3, 4)
    first = tmp_path / "first.npy"
    second = tmp_path / "second.npy"
    write_npy_v1_16byte_aligned(first, array)
    write_npy_v1_16byte_aligned(second, array)
    assert first.read_bytes() == second.read_bytes()
    assert len(first.read_bytes()[:80]) == 80
    assert np.array_equal(np.load(first, allow_pickle=False), array)

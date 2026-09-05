"""Tests for source-line parsing and deterministic graph-row reconstruction."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import tarfile
from pathlib import Path

from graphsage_ppi_repro.legacy_order import ordered_string_keys
from graphsage_ppi_repro.topology import reconstruct_topology


def _write_network_archive(path: Path, member: str, payload: bytes) -> None:
    with tarfile.open(path, "w:gz") as archive:
        info = tarfile.TarInfo(member)
        info.size = len(payload)
        info.mtime = 0
        archive.addfile(info, io.BytesIO(payload))


def _read_tsv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_reconstruct_small_network(tmp_path: Path) -> None:
    member = "bio-tissue-networks/example.edgelist"
    payload = b"10 20\n30 10\n20 40\n40 40\n"
    archive = tmp_path / "networks.tar.gz"
    _write_network_archive(archive, member, payload)

    graph_spec = tmp_path / "selected_graphs.tsv"
    graph_spec.write_text(
        "\t".join(
            [
                "graph_index_1based",
                "tissue",
                "split",
                "ohmnet_member",
                "ohmnet_member_sha256",
                "node_count",
                "edge_count",
                "evidence_status",
                "evidence_reference",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "1",
                "example",
                "train",
                member,
                hashlib.sha256(payload).hexdigest(),
                "4",
                "4",
                "synthetic",
                "tests/test_topology.py",
            ]
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )

    mapping = tmp_path / "node_to_entrez.tsv.gz"
    edges = tmp_path / "edges.tsv.gz"
    summary = tmp_path / "summary.json"
    result = reconstruct_topology(
        archive_path=archive,
        graph_spec_path=graph_spec,
        mapping_output=mapping,
        edge_output=edges,
        summary_output=summary,
    )

    expected_genes = ordered_string_keys(["10", "20", "30", "10", "20", "40", "40", "40"])
    mapping_rows = _read_tsv_gz(mapping)
    assert [row["entrez_gene_id"] for row in mapping_rows] == expected_genes
    assert [int(row["graphsage_node_id"]) for row in mapping_rows] == list(range(4))

    node_by_gene = {
        row["entrez_gene_id"]: int(row["graphsage_node_id"]) for row in mapping_rows
    }
    edge_rows = _read_tsv_gz(edges)
    assert [
        (int(row["source_node_id"]), int(row["target_node_id"])) for row in edge_rows
    ] == [
        (node_by_gene["10"], node_by_gene["20"]),
        (node_by_gene["30"], node_by_gene["10"]),
        (node_by_gene["20"], node_by_gene["40"]),
        (node_by_gene["40"], node_by_gene["40"]),
    ]
    assert result["counts"] == {
        "graphs": 1,
        "rows": 4,
        "distinct_entrez_gene_ids": 4,
        "edge_records": 4,
        "split_graph_counts": {"train": 1},
        "split_row_counts": {"train": 4},
    }


def test_reconstruction_outputs_are_byte_deterministic(tmp_path: Path) -> None:
    member = "bio-tissue-networks/example.edgelist"
    payload = b"1 2\n2 3\n"
    archive = tmp_path / "networks.tar.gz"
    _write_network_archive(archive, member, payload)
    graph_spec = tmp_path / "selected_graphs.tsv"
    graph_spec.write_text(
        "graph_index_1based\ttissue\tsplit\tohmnet_member\tohmnet_member_sha256\t"
        "node_count\tedge_count\tevidence_status\tevidence_reference\n"
        f"1\texample\ttest\t{member}\t{hashlib.sha256(payload).hexdigest()}\t"
        "3\t2\tsynthetic\ttest\n",
        encoding="utf-8",
        newline="\n",
    )

    output_hashes: list[tuple[str, str]] = []
    for run in (1, 2):
        mapping = tmp_path / f"mapping-{run}.tsv.gz"
        edges = tmp_path / f"edges-{run}.tsv.gz"
        reconstruct_topology(
            archive_path=archive,
            graph_spec_path=graph_spec,
            mapping_output=mapping,
            edge_output=edges,
            summary_output=tmp_path / f"summary-{run}.json",
        )
        output_hashes.append(
            (
                hashlib.sha256(mapping.read_bytes()).hexdigest(),
                hashlib.sha256(edges.read_bytes()).hexdigest(),
            )
        )
    assert output_hashes[0] == output_hashes[1]

"""Tests for deterministic GraphSAGE file assembly and validation."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from graphsage_ppi_repro.graphsage import GraphSAGEError, assemble_graphsage
from graphsage_ppi_repro.validate import ValidationError, validate_graphsage_dataset


def _write_mapping(path: Path) -> None:
    fields = (
        "graphsage_node_id",
        "feature_label_row_index",
        "graph_index_1based",
        "tissue",
        "split",
        "local_node_index_0based",
        "entrez_gene_id",
    )
    rows = (
        (0, 0, 1, "training", "train", 0, 10),
        (1, 1, 1, "training", "train", 1, 20),
        (2, 2, 2, "validation", "validation", 0, 30),
        (3, 3, 3, "testing", "test", 0, 40),
    )
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


def _write_edges(path: Path, edges: list[tuple[int, int]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("source_node_id", "target_node_id"))
        writer.writerows(edges)


def _write_matrices(tmp_path: Path) -> tuple[Path, Path]:
    features = tmp_path / "features.npy"
    labels = tmp_path / "labels.npy"
    np.save(
        features,
        np.asarray(
            [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]],
            dtype=np.float64,
        ),
        allow_pickle=False,
    )
    np.save(
        labels,
        np.asarray([[1, 0], [0, 1], [1, 1], [0, 0]], dtype=np.uint8),
        allow_pickle=False,
    )
    return features, labels


def _assemble(tmp_path: Path, name: str = "output") -> tuple[Path, dict[str, object]]:
    mapping = tmp_path / "mapping.tsv.gz"
    edges = tmp_path / "edges.tsv.gz"
    if not mapping.exists():
        _write_mapping(mapping)
    if not edges.exists():
        _write_edges(edges, [(1, 0), (3, 3), (2, 2), (1, 1)])
    features, labels = _write_matrices(tmp_path)
    output_root = tmp_path / name
    summary = assemble_graphsage(
        mapping_path=mapping,
        edge_path=edges,
        feature_matrix_path=features,
        label_matrix_path=labels,
        output_directory=output_root / "ppi",
        summary_output=tmp_path / f"{name}-summary.json",
    )
    return output_root, summary


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_assemble_graphsage_writes_deterministic_supervised_files(
    tmp_path: Path,
) -> None:
    first, first_summary = _assemble(tmp_path, "first")
    second, second_summary = _assemble(tmp_path, "second")

    expected_files = (
        "ppi-G.json",
        "ppi-id_map.json",
        "ppi-class_map.json",
        "ppi-feats.npy",
    )
    for filename in expected_files:
        assert (first / "ppi" / filename).read_bytes() == (
            second / "ppi" / filename
        ).read_bytes()

    graph = json.loads((first / "ppi" / "ppi-G.json").read_text())
    assert graph == {
        "directed": False,
        "graph": {},
        "nodes": [
            {"test": False, "id": 0, "val": False},
            {"test": False, "id": 1, "val": False},
            {"test": False, "id": 2, "val": True},
            {"test": True, "id": 3, "val": False},
        ],
        "links": [
            {"source": 0, "target": 1},
            {"source": 1, "target": 1},
            {"source": 2, "target": 2},
            {"source": 3, "target": 3},
        ],
        "multigraph": False,
    }
    assert json.loads((first / "ppi" / "ppi-id_map.json").read_text()) == {
        str(index): index for index in range(4)
    }
    assert json.loads((first / "ppi" / "ppi-class_map.json").read_text()) == {
        "0": [1, 0],
        "1": [0, 1],
        "2": [1, 1],
        "3": [0, 0],
    }
    assert (first / "ppi" / "ppi-feats.npy").read_bytes() == (
        tmp_path / "features.npy"
    ).read_bytes()
    assert not (first / "ppi" / "ppi-walks.txt").exists()

    checksum_lines = (first / "SHA256SUMS").read_text().splitlines()
    assert checksum_lines == [
        f"{_sha256(first / 'ppi' / filename)}  ppi/{filename}"
        for filename in expected_files
    ]
    assert first_summary["counts"] == second_summary["counts"]
    assert first_summary["counts"] == {
        "graphs": 3,
        "rows": 4,
        "edge_records": 4,
        "feature_columns": 2,
        "label_columns": 2,
        "split_rows": {"test": 1, "train": 2, "validation": 1},
    }


def test_full_validation_accepts_structural_graph_equivalence(tmp_path: Path) -> None:
    output_root, _summary = _assemble(tmp_path)
    graph_path = output_root / "ppi" / "ppi-G.json"
    id_map_path = output_root / "ppi" / "ppi-id_map.json"
    class_map_path = output_root / "ppi" / "ppi-class_map.json"
    feature_path = output_root / "ppi" / "ppi-feats.npy"

    graph = json.loads(graph_path.read_text())
    graph["links"] = list(reversed(graph["links"]))
    reference = tmp_path / "reference.zip"
    with zipfile.ZipFile(reference, "w") as archive:
        archive.writestr("ppi/ppi-G.json", json.dumps(graph))
        archive.writestr("ppi/ppi-id_map.json", id_map_path.read_bytes())
        archive.writestr("ppi/ppi-class_map.json", class_map_path.read_bytes())
        archive.writestr("ppi/ppi-feats.npy", feature_path.read_bytes())

    result = validate_graphsage_dataset(
        reference_archive=reference,
        graph_path=graph_path,
        id_map_path=id_map_path,
        class_map_path=class_map_path,
        feature_path=feature_path,
        json_output=tmp_path / "validation.json",
        markdown_output=tmp_path / "validation.md",
    )
    assert result["all_checks_pass"]
    assert not result["byte_equality"]["ppi-G.json"]
    assert result["checks"]["graph_edge_multiset_exact"]
    assert all(
        result["byte_equality"][filename]
        for filename in ("ppi-id_map.json", "ppi-class_map.json", "ppi-feats.npy")
    )


def test_validation_rejects_changed_class_map(tmp_path: Path) -> None:
    output_root, _summary = _assemble(tmp_path)
    graph_path = output_root / "ppi" / "ppi-G.json"
    id_map_path = output_root / "ppi" / "ppi-id_map.json"
    class_map_path = output_root / "ppi" / "ppi-class_map.json"
    feature_path = output_root / "ppi" / "ppi-feats.npy"

    changed = json.loads(class_map_path.read_text())
    changed["0"][0] = 0
    reference = tmp_path / "reference.zip"
    with zipfile.ZipFile(reference, "w") as archive:
        archive.writestr("ppi/ppi-G.json", graph_path.read_bytes())
        archive.writestr("ppi/ppi-id_map.json", id_map_path.read_bytes())
        archive.writestr("ppi/ppi-class_map.json", json.dumps(changed))
        archive.writestr("ppi/ppi-feats.npy", feature_path.read_bytes())

    with pytest.raises(ValidationError, match="class_map"):
        validate_graphsage_dataset(
            reference_archive=reference,
            graph_path=graph_path,
            id_map_path=id_map_path,
            class_map_path=class_map_path,
            feature_path=feature_path,
            json_output=tmp_path / "validation.json",
            markdown_output=tmp_path / "validation.md",
        )


def test_assembly_rejects_duplicate_or_cross_graph_edges(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.tsv.gz"
    _write_mapping(mapping)
    features, labels = _write_matrices(tmp_path)

    duplicate_edges = tmp_path / "duplicate-edges.tsv.gz"
    _write_edges(duplicate_edges, [(0, 1), (1, 0)])
    with pytest.raises(GraphSAGEError, match="duplicate"):
        assemble_graphsage(
            mapping_path=mapping,
            edge_path=duplicate_edges,
            feature_matrix_path=features,
            label_matrix_path=labels,
            output_directory=tmp_path / "duplicate-output" / "ppi",
            summary_output=tmp_path / "duplicate-summary.json",
        )

    cross_graph_edges = tmp_path / "cross-graph-edges.tsv.gz"
    _write_edges(cross_graph_edges, [(0, 2)])
    with pytest.raises(GraphSAGEError, match="cross-graph"):
        assemble_graphsage(
            mapping_path=mapping,
            edge_path=cross_graph_edges,
            feature_matrix_path=features,
            label_matrix_path=labels,
            output_directory=tmp_path / "cross-output" / "ppi",
            summary_output=tmp_path / "cross-summary.json",
        )

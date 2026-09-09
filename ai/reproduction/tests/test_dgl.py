from __future__ import annotations

import csv
import json
import zipfile
from pathlib import Path

import numpy as np
import pytest

from graphsage_ppi_repro.dgl import assemble_dgl
from graphsage_ppi_repro.validate import ValidationError, validate_dgl_dataset


def _write_graph_spec(path: Path) -> None:
    fields = (
        "graph_index_1based",
        "tissue",
        "split",
        "ohmnet_member",
        "ohmnet_member_sha256",
        "node_count",
        "edge_count",
        "evidence_status",
        "evidence_reference",
    )
    rows = [
        (1, "train_a", "train", 4, 3),
        (2, "train_b", "train", 5, 3),
        (3, "valid_a", "validation", 2, 1),
        (4, "test_a", "test", 2, 1),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for graph_index, tissue, split, nodes, edges in rows:
            writer.writerow(
                {
                    "graph_index_1based": graph_index,
                    "tissue": tissue,
                    "split": split,
                    "ohmnet_member": f"networks/{tissue}.edgelist",
                    "ohmnet_member_sha256": "0" * 64,
                    "node_count": nodes,
                    "edge_count": edges,
                    "evidence_status": "synthetic",
                    "evidence_reference": "test",
                }
            )


def _write_graphsage_fixture(root: Path) -> dict[str, Path]:
    root.mkdir(parents=True)
    splits = ["train"] * 9 + ["validation"] * 2 + ["test"] * 2
    nodes = [
        {
            "test": split == "test",
            "id": index,
            "val": split == "validation",
        }
        for index, split in enumerate(splits)
    ]
    edges = [(0, 1), (1, 1), (2, 2), (4, 5), (5, 6), (7, 8), (9, 10), (11, 12)]
    graph = {
        "directed": False,
        "graph": {},
        "nodes": nodes,
        "links": [{"source": source, "target": target} for source, target in edges],
        "multigraph": False,
    }
    graph_path = root / "ppi-G.json"
    graph_path.write_text(json.dumps(graph), encoding="utf-8")

    id_map_path = root / "ppi-id_map.json"
    id_map_path.write_text(
        json.dumps({str(index): index for index in range(len(nodes))}),
        encoding="utf-8",
    )
    labels = np.asarray(
        [[index % 2, (index // 2) % 2, (index // 3) % 2] for index in range(len(nodes))],
        dtype=np.int64,
    )
    class_map_path = root / "ppi-class_map.json"
    class_map_path.write_text(
        json.dumps({str(index): labels[index].tolist() for index in range(len(nodes))}),
        encoding="utf-8",
    )
    features = np.asarray(
        [[index % 2, 0.0, (index // 2) % 2] for index in range(len(nodes))],
        dtype=np.float64,
    )
    feature_path = root / "ppi-feats.npy"
    np.save(feature_path, features, allow_pickle=False)
    return {
        "graph": graph_path,
        "id_map": id_map_path,
        "class_map": class_map_path,
        "features": feature_path,
        "labels": labels,
        "raw_features": features,
    }


def _build_synthetic_dgl(tmp_path: Path) -> tuple[Path, dict[str, object], dict[str, Path]]:
    inputs = _write_graphsage_fixture(tmp_path / "graphsage")
    graph_spec = tmp_path / "selected_graphs.tsv"
    _write_graph_spec(graph_spec)
    output = tmp_path / "dgl" / "ppi"
    summary = assemble_dgl(
        graph_path=inputs["graph"],
        id_map_path=inputs["id_map"],
        class_map_path=inputs["class_map"],
        feature_path=inputs["features"],
        graph_spec_path=graph_spec,
        output_directory=output,
        summary_output=tmp_path / "build" / "dgl-summary.json",
    )
    return output, summary, inputs


def _archive_directory(directory: Path, destination: Path) -> None:
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(directory.iterdir()):
            archive.write(path, path.name)


def test_dgl_assembly_reproduces_component_grouping_and_scaling(tmp_path: Path) -> None:
    output, summary, inputs = _build_synthetic_dgl(tmp_path)

    train_graph_ids = np.load(output / "train_graph_id.npy", allow_pickle=False)
    assert train_graph_ids.tolist() == [1, 1, 1, 1, 1, 1, 2, 2, 2]
    expected_order = np.asarray([0, 1, 2, 3, 7, 8, 4, 5, 6])
    train_labels = np.load(output / "train_labels.npy", allow_pickle=False)
    assert train_labels.dtype == np.dtype(np.int64)
    assert np.array_equal(train_labels, inputs["labels"][expected_order])

    raw_train = inputs["raw_features"][:9]
    mean = raw_train.mean(axis=0)
    scale = np.sqrt(raw_train.var(axis=0))
    scale[scale == 0.0] = 1.0
    expected_features = (raw_train - mean) / scale
    train_features = np.load(output / "train_feats.npy", allow_pickle=False)
    assert train_features.dtype == np.dtype(np.float64)
    assert np.array_equal(train_features, expected_features[expected_order])
    assert summary["scaler"]["zero_variance_columns_0based"] == [1]

    assignments = summary["component_assignments"]
    second = assignments[1]
    assert second["largest_component_nodes"] == 3
    assert second["non_largest_component_nodes"] == 2
    assert second["largest_component_graph_id"] == 2
    assert second["non_largest_component_graph_id"] == 1


def test_dgl_graphs_have_bidirectional_edges_and_one_loop_per_node(tmp_path: Path) -> None:
    output, _, _ = _build_synthetic_dgl(tmp_path)
    graph = json.loads((output / "train_graph.json").read_text(encoding="utf-8"))
    edges = [(link["source"], link["target"]) for link in graph["links"]]
    edge_set = set(edges)

    assert graph["directed"] is True
    assert graph["multigraph"] is False
    assert edges == sorted(edges)
    assert len(edges) == len(edge_set)
    assert sum(source == target for source, target in edges) == 9
    assert all((node, node) in edge_set for node in range(9))
    assert (0, 1) in edge_set and (1, 0) in edge_set
    assert (6, 7) in edge_set and (7, 6) in edge_set


def test_dgl_discrete_arrays_and_checksums_are_deterministic(tmp_path: Path) -> None:
    first, _, _ = _build_synthetic_dgl(tmp_path / "first")
    second, _, _ = _build_synthetic_dgl(tmp_path / "second")
    for split in ("train", "valid", "test"):
        for suffix in ("graph_id.npy", "labels.npy", "graph.json", "feats.npy"):
            assert (first / f"{split}_{suffix}").read_bytes() == (
                second / f"{split}_{suffix}"
            ).read_bytes()
    assert (first.parent / "SHA256SUMS").read_bytes() == (
        second.parent / "SHA256SUMS"
    ).read_bytes()


def test_dgl_validation_uses_fixed_numeric_tolerance(tmp_path: Path) -> None:
    output, _, _ = _build_synthetic_dgl(tmp_path)
    reference = tmp_path / "reference.zip"
    _archive_directory(output, reference)

    train_features_path = output / "train_feats.npy"
    train_features = np.load(train_features_path, allow_pickle=False)
    train_features[0, 0] += 5.0e-13
    np.save(train_features_path, train_features, allow_pickle=False)

    result = validate_dgl_dataset(
        reference_archive=reference,
        reconstructed_directory=output,
        json_output=tmp_path / "validation.json",
        markdown_output=tmp_path / "validation.md",
        feature_atol=1.0e-12,
        feature_rtol=0.0,
    )
    assert result["all_checks_pass"]
    assert not result["splits"]["train"]["byte_equality"]["feats"]

    with pytest.raises(ValidationError, match="feature_values_within_tolerance"):
        validate_dgl_dataset(
            reference_archive=reference,
            reconstructed_directory=output,
            json_output=tmp_path / "strict-validation.json",
            markdown_output=tmp_path / "strict-validation.md",
            feature_atol=1.0e-14,
            feature_rtol=0.0,
        )


def test_dgl_reconstruction_does_not_open_reference_archive() -> None:
    reproduction_root = Path(__file__).resolve().parents[1]
    module = (reproduction_root / "src/graphsage_ppi_repro/dgl.py").read_text(
        encoding="utf-8"
    )
    snakefile = (reproduction_root / "Snakefile").read_text(encoding="utf-8")

    assert "dgl_ppi.zip" not in module
    assemble_block = snakefile.split("rule assemble_dgl:", 1)[1].split("\n\nrule ", 1)[0]
    assert "DGL_REFERENCE" not in assemble_block
    verify_block = snakefile.split("rule verify_dgl_reference:", 1)[1].split(
        "\n\nrule ", 1
    )[0]
    assert "archive=str(DGL_REFERENCE)" in verify_block
    assert "output:" in verify_block
    output_block = verify_block.split("output:", 1)[1].split("params:", 1)[0]
    assert "archive=str(DGL_REFERENCE)" not in output_block

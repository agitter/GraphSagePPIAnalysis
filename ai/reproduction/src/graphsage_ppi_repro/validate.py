"""Independent comparison of reconstructed artifacts with released targets.

No reconstruction module imports this file.  Keeping all target access here is
an auditable non-circularity boundary: ``pixi run reproduce`` can complete when
reference archives are absent, while ``pixi run validate`` requires them.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO

import numpy as np

from .provenance import sha256_file, write_json_atomic, write_text_atomic
from .topology import read_node_mapping


class ValidationError(RuntimeError):
    """Raised when reconstructed data disagree with the released reference."""


def _canonical_edge(source: int, target: int) -> tuple[int, int]:
    return (source, target) if source <= target else (target, source)


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _read_edges(path: Path) -> list[tuple[int, int]]:
    with _open_text(path) as handle:
        return [
            (int(row["source_node_id"]), int(row["target_node_id"]))
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def _member_sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _graphwise_edge_multisets(
    edges: list[tuple[int, int]], node_to_graph: list[int]
) -> tuple[dict[int, Counter[tuple[int, int]]], int]:
    by_graph: dict[int, Counter[tuple[int, int]]] = defaultdict(Counter)
    cross_graph = 0
    for source, target in edges:
        source_graph = node_to_graph[source]
        target_graph = node_to_graph[target]
        if source_graph != target_graph:
            cross_graph += 1
            continue
        by_graph[source_graph][_canonical_edge(source, target)] += 1
    return dict(by_graph), cross_graph


def _write_markdown(path: Path, result: dict[str, object]) -> None:
    checks = result["checks"]
    if not isinstance(checks, dict):
        raise TypeError("Validation result checks must be a dictionary")

    lines = [
        "# GraphSAGE PPI reconstruction validation - topology and features",
        "",
        f"Overall result: **{'PASS' if result['all_checks_pass'] else 'FAIL'}**",
        "",
        "This implementation milestone validates the complete biological row order,",
        "all 818,716 logical undirected edge records, the deposited split flags, and",
        "all 2,847,200 feature cells. GO labels, final GraphSAGE serialization, and",
        "the DGL transformation remain later implementation stages.",
        "",
        "| Check | Result |",
        "|---|---:|",
    ]
    for name, passed in checks.items():
        lines.append(f"| `{name}` | {'PASS' if passed else 'FAIL'} |")
    lines.extend(
        [
            "",
            "## Comparison levels",
            "",
            "- Row order: consecutive target IDs and exact split flags.",
            "- Topology: exact graph-wise undirected edge multisets.",
            "- Features: exact dtype, shape, values, and deposited NPY bytes.",
            "- Reference access: confined to `graphsage_ppi_repro.validate`.",
            "",
        ]
    )
    write_text_atomic(path, "\n".join(lines))


def validate_topology_and_features(
    *,
    graphsage_reference_zip: Path,
    mapping_path: Path,
    edge_path: Path,
    feature_matrix_path: Path,
    output_json: Path,
    output_markdown: Path,
) -> dict[str, object]:
    """Compare milestone-one reconstruction results to the released archive."""

    mapping = read_node_mapping(mapping_path)
    reconstructed_edges = _read_edges(edge_path)
    reconstructed_features = np.load(feature_matrix_path, allow_pickle=False)

    with zipfile.ZipFile(graphsage_reference_zip) as archive:
        graph_bytes = archive.read("ppi/ppi-G.json")
        id_map_bytes = archive.read("ppi/ppi-id_map.json")
        feature_bytes = archive.read("ppi/ppi-feats.npy")
    graph = json.loads(graph_bytes)
    id_map = json.loads(id_map_bytes)
    target_features = np.load(io.BytesIO(feature_bytes), allow_pickle=False)
    target_edges = [(int(edge["source"]), int(edge["target"])) for edge in graph["links"]]

    expected_node_ids = list(range(len(mapping)))
    mapping_node_ids = [int(row["graphsage_node_id"]) for row in mapping]
    target_nodes = sorted(graph["nodes"], key=lambda node: int(node["id"]))
    target_node_ids = [int(node["id"]) for node in target_nodes]
    mapping_splits = [row["split"] for row in mapping]
    target_splits = [
        (
            "test"
            if bool(node.get("test"))
            else "validation"
            if bool(node.get("val"))
            else "train"
        )
        for node in target_nodes
    ]
    node_to_graph = [int(row["graph_index_1based"]) for row in mapping]

    reconstructed_by_graph, reconstructed_cross_graph = _graphwise_edge_multisets(
        reconstructed_edges, node_to_graph
    )
    target_by_graph, target_cross_graph = _graphwise_edge_multisets(
        target_edges, node_to_graph
    )

    checks = {
        "mapping_rows_consecutive": mapping_node_ids == expected_node_ids,
        "reference_node_ids_consecutive": target_node_ids == expected_node_ids,
        "reference_id_map_is_identity": len(id_map) == len(mapping)
        and all(int(id_map[str(index)]) == index for index in expected_node_ids),
        "split_flags_exact": mapping_splits == target_splits,
        "reconstructed_has_no_cross_graph_edges": reconstructed_cross_graph == 0,
        "reference_has_no_cross_graph_edges": target_cross_graph == 0,
        "edge_record_count_exact": len(reconstructed_edges) == len(target_edges),
        "graphwise_undirected_edge_multisets_exact": (
            reconstructed_by_graph == target_by_graph
        ),
        "feature_shape_exact": reconstructed_features.shape == target_features.shape,
        "feature_dtype_exact": reconstructed_features.dtype == target_features.dtype,
        "feature_array_exact": np.array_equal(reconstructed_features, target_features),
        "feature_npy_bytes_exact": feature_matrix_path.read_bytes() == feature_bytes,
    }

    mismatching_graphs = sorted(
        graph_index
        for graph_index in set(reconstructed_by_graph) | set(target_by_graph)
        if reconstructed_by_graph.get(graph_index, Counter())
        != target_by_graph.get(graph_index, Counter())
    )
    result: dict[str, object] = {
        "schema_version": 1,
        "scope": "complete topology, row order, split flags, and 50 input features",
        "reference": {
            "archive": str(graphsage_reference_zip.resolve()),
            "archive_sha256": sha256_file(graphsage_reference_zip),
            "ppi_G_json_sha256": _member_sha256(graph_bytes),
            "ppi_id_map_json_sha256": _member_sha256(id_map_bytes),
            "ppi_feats_npy_sha256": _member_sha256(feature_bytes),
        },
        "reconstructed": {
            "node_mapping": str(mapping_path.resolve()),
            "node_mapping_sha256": sha256_file(mapping_path),
            "edge_table": str(edge_path.resolve()),
            "edge_table_sha256": sha256_file(edge_path),
            "feature_matrix": str(feature_matrix_path.resolve()),
            "feature_matrix_sha256": sha256_file(feature_matrix_path),
        },
        "counts": {
            "rows": len(mapping),
            "graphs": len(set(node_to_graph)),
            "edge_records": len(reconstructed_edges),
            "feature_shape": [int(value) for value in reconstructed_features.shape],
            "feature_cells": int(reconstructed_features.size),
        },
        "diagnostics": {
            "mismatching_graph_indices_1based": mismatching_graphs,
            "reconstructed_cross_graph_edges": reconstructed_cross_graph,
            "reference_cross_graph_edges": target_cross_graph,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    write_json_atomic(output_json, result)
    _write_markdown(output_markdown, result)
    if not result["all_checks_pass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise ValidationError("Validation failed: " + ", ".join(failed))
    return result

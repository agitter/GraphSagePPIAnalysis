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


def _write_label_markdown(path: Path, result: dict[str, object]) -> None:
    checks = result["checks"]
    if not isinstance(checks, dict):
        raise TypeError("Label-validation checks must be a dictionary")
    lines = [
        "# GraphSAGE PPI reconstruction validation - GO labels",
        "",
        f"Overall result: **{'PASS' if result['all_checks_pass'] else 'FAIL'}**",
        "",
        "The reconstructed 56,944 x 121 label matrix and numeric-node-ordered",
        "class-map semantics were compared with the released `ppi-class_map.json`.",
        "The reference archive is not read during reconstruction.",
        "",
        "| Check | Result |",
        "|---|---:|",
    ]
    for name, passed in checks.items():
        lines.append(f"| `{name}` | {'PASS' if passed else 'FAIL'} |")
    lines.extend(
        [
            "",
            f"- Cells compared: {result['cells_compared']:,}",
            f"- Differing cells: {result['differing_cells']}",
            f"- Differing rows: {result['differing_rows']}",
            "",
        ]
    )
    write_text_atomic(path, "\n".join(lines))


def validate_labels_against_graphsage(
    *,
    reference_archive: Path,
    reconstructed_matrix: Path,
    reconstructed_class_map: Path,
    json_output: Path,
    markdown_output: Path,
    expected_rows: int = 56_944,
    expected_columns: int = 121,
) -> dict[str, object]:
    """Compare reconstructed labels with the released GraphSAGE class map."""

    reconstructed = np.load(reconstructed_matrix, allow_pickle=False)
    class_map = json.loads(reconstructed_class_map.read_text(encoding="utf-8"))
    if not isinstance(class_map, dict):
        raise ValidationError("Reconstructed class map must be a JSON object")

    with zipfile.ZipFile(reference_archive) as archive:
        member = next(
            (name for name in archive.namelist() if name.endswith("ppi-class_map.json")),
            None,
        )
        if member is None:
            raise ValidationError(
                f"No ppi-class_map.json member found in {reference_archive}"
            )
        target_bytes = archive.read(member)
    target = json.loads(target_bytes)
    if not isinstance(target, dict):
        raise ValidationError("Released class map must be a JSON object")

    expected_keys = [str(index) for index in range(len(target))]
    target_keys_exact = set(target) == set(expected_keys)
    reconstructed_keys_exact = list(class_map) == expected_keys
    if target_keys_exact:
        target_matrix = np.asarray([target[key] for key in expected_keys], dtype=np.uint8)
    else:
        target_matrix = np.empty((0, 0), dtype=np.uint8)
    if reconstructed_keys_exact:
        class_map_matrix = np.asarray(
            [class_map[key] for key in expected_keys], dtype=np.uint8
        )
    else:
        class_map_matrix = np.empty((0, 0), dtype=np.uint8)

    shape_exact = reconstructed.shape == target_matrix.shape
    class_map_shape_exact = class_map_matrix.shape == target_matrix.shape
    if shape_exact:
        difference = reconstructed != target_matrix
        differing_cells = int(np.count_nonzero(difference))
        differing_rows = int(np.count_nonzero(np.any(difference, axis=1)))
        mismatch_positions = np.argwhere(difference)
        first_mismatch = (
            [int(value) for value in mismatch_positions[0]]
            if mismatch_positions.size
            else None
        )
    else:
        differing_cells = None
        differing_rows = None
        first_mismatch = None

    checks = {
        "reference_keys_complete": target_keys_exact,
        "reconstructed_keys_numeric_order": reconstructed_keys_exact,
        "matrix_dtype_uint8": reconstructed.dtype == np.dtype(np.uint8),
        "matrix_shape_exact": shape_exact,
        "matrix_values_exact": shape_exact and differing_cells == 0,
        "class_map_shape_exact": class_map_shape_exact,
        "class_map_values_exact": class_map_shape_exact
        and np.array_equal(class_map_matrix, target_matrix),
        "label_column_count_exact": (
            shape_exact and reconstructed.shape[1] == expected_columns
        ),
        "label_row_count_exact": shape_exact and reconstructed.shape[0] == expected_rows,
    }
    result: dict[str, object] = {
        "schema_version": 1,
        "scope": (
            f"complete {expected_rows:,}-row by {expected_columns}-column GO label matrix"
        ),
        "reference": {
            "archive": str(reference_archive.resolve()),
            "archive_sha256": sha256_file(reference_archive),
            "member": member,
            "member_sha256": _member_sha256(target_bytes),
        },
        "reconstructed": {
            "matrix": str(reconstructed_matrix.resolve()),
            "matrix_sha256": sha256_file(reconstructed_matrix),
            "class_map": str(reconstructed_class_map.resolve()),
            "class_map_sha256": sha256_file(reconstructed_class_map),
        },
        "reconstructed_shape": [int(value) for value in reconstructed.shape],
        "reference_shape": [int(value) for value in target_matrix.shape],
        "reconstructed_dtype": str(reconstructed.dtype),
        "positive_cells_reconstructed": int(reconstructed.sum()),
        "positive_cells_reference": int(target_matrix.sum()),
        "cells_compared": int(target_matrix.size) if shape_exact else 0,
        "differing_cells": differing_cells,
        "differing_rows": differing_rows,
        "first_mismatch_row_column": first_mismatch,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    write_json_atomic(json_output, result)
    _write_label_markdown(markdown_output, result)
    if not result["all_checks_pass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise ValidationError("Label validation failed: " + ", ".join(failed))
    return result


def _edge_multiset(graph: dict[str, object]) -> Counter[tuple[int, int]]:
    links = graph.get("links")
    if not isinstance(links, list):
        raise ValidationError("Graph node-link JSON must contain a links list")
    edges: Counter[tuple[int, int]] = Counter()
    for index, link in enumerate(links):
        if not isinstance(link, dict):
            raise ValidationError(f"Graph link {index} is not a JSON object")
        try:
            source = int(link["source"])
            target = int(link["target"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError(f"Invalid graph link at index {index}") from exc
        edges[_canonical_edge(source, target)] += 1
    return edges


def _ordered_graph_edges(graph: dict[str, object]) -> list[tuple[int, int]]:
    links = graph.get("links")
    if not isinstance(links, list):
        raise ValidationError("Graph node-link JSON must contain a links list")
    return [
        (int(link["source"]), int(link["target"]))
        for link in links
        if isinstance(link, dict)
    ]


def _numeric_class_matrix(
    value: object,
    *,
    description: str,
) -> tuple[list[str], np.ndarray]:
    if not isinstance(value, dict):
        raise ValidationError(f"{description} must be a JSON object")
    expected_keys = [str(index) for index in range(len(value))]
    if set(value) != set(expected_keys):
        raise ValidationError(
            f"{description} does not contain consecutive numeric string keys"
        )
    try:
        matrix = np.asarray([value[key] for key in expected_keys])
    except (TypeError, ValueError) as exc:
        raise ValidationError(
            f"{description} does not contain a rectangular matrix"
        ) from exc
    if matrix.ndim != 2:
        raise ValidationError(f"{description} must contain one label vector per node")
    if not np.issubdtype(matrix.dtype, np.integer):
        raise ValidationError(f"{description} labels must be integers")
    if not np.all((matrix == 0) | (matrix == 1)):
        raise ValidationError(f"{description} labels must be binary")
    return expected_keys, matrix.astype(np.uint8, copy=False)


def _write_graphsage_markdown(path: Path, result: dict[str, object]) -> None:
    checks = result["checks"]
    byte_equality = result["byte_equality"]
    if not isinstance(checks, dict) or not isinstance(byte_equality, dict):
        raise TypeError("GraphSAGE validation checks must be dictionaries")

    lines = [
        "# Complete supervised GraphSAGE PPI artifact validation",
        "",
        f"Overall result: **{'PASS' if result['all_checks_pass'] else 'FAIL'}**",
        "",
        "The four reconstructed files consumed by supervised GraphSAGE were",
        "compared with the corresponding members of the released archive.",
        "`ppi-walks.txt` is optional unsupervised input and remains out of scope.",
        "",
        "| Required comparison | Result |",
        "|---|---:|",
    ]
    for name, passed in checks.items():
        lines.append(f"| `{name}` | {'PASS' if passed else 'FAIL'} |")
    lines.extend(
        [
            "",
            "## Byte-level comparisons",
            "",
            "| File | Byte-identical to release |",
            "|---|---:|",
        ]
    )
    for name, passed in byte_equality.items():
        lines.append(f"| `{name}` | {'YES' if passed else 'NO'} |")
    lines.extend(
        [
            "",
            "The ID map, class map, and NumPy feature file are required to be",
            "byte-identical. The graph JSON is required to be structurally exact;",
            "the reconstruction deliberately sorts undirected links into a simple",
            "canonical order instead of reproducing NetworkX/Python-2 dictionary",
            "iteration order. This preserves the logical graph but can change the",
            "exact order supplied to stochastic neighbor-sampling code.",
            "",
            f"- Nodes compared: {result['counts']['nodes']:,}",
            f"- Undirected edge records compared: {result['counts']['edges']:,}",
            f"- Feature cells compared: {result['counts']['feature_cells']:,}",
            f"- Label cells compared: {result['counts']['label_cells']:,}",
            "",
        ]
    )
    write_text_atomic(path, "\n".join(lines))


def validate_graphsage_dataset(
    *,
    reference_archive: Path,
    graph_path: Path,
    id_map_path: Path,
    class_map_path: Path,
    feature_path: Path,
    json_output: Path,
    markdown_output: Path,
) -> dict[str, object]:
    """Validate the complete deterministic supervised GraphSAGE data contract."""

    reconstructed_bytes = {
        "ppi-G.json": graph_path.read_bytes(),
        "ppi-id_map.json": id_map_path.read_bytes(),
        "ppi-class_map.json": class_map_path.read_bytes(),
        "ppi-feats.npy": feature_path.read_bytes(),
    }
    try:
        with zipfile.ZipFile(reference_archive) as archive:
            reference_bytes = {
                name: archive.read(f"ppi/{name}") for name in reconstructed_bytes
            }
    except (KeyError, zipfile.BadZipFile) as exc:
        raise ValidationError(
            f"Invalid or incomplete GraphSAGE reference archive: {reference_archive}"
        ) from exc

    try:
        reconstructed_graph = json.loads(reconstructed_bytes["ppi-G.json"])
        reference_graph = json.loads(reference_bytes["ppi-G.json"])
        reconstructed_id_map = json.loads(reconstructed_bytes["ppi-id_map.json"])
        reference_id_map = json.loads(reference_bytes["ppi-id_map.json"])
        reconstructed_class_map = json.loads(reconstructed_bytes["ppi-class_map.json"])
        reference_class_map = json.loads(reference_bytes["ppi-class_map.json"])
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValidationError(f"Invalid reconstructed or reference JSON: {exc}") from exc

    if not isinstance(reconstructed_graph, dict) or not isinstance(reference_graph, dict):
        raise ValidationError("Graph files must contain JSON objects")
    if not isinstance(reconstructed_id_map, dict) or not isinstance(reference_id_map, dict):
        raise ValidationError("ID-map files must contain JSON objects")

    reconstructed_keys, reconstructed_labels = _numeric_class_matrix(
        reconstructed_class_map,
        description="Reconstructed class map",
    )
    reference_keys, reference_labels = _numeric_class_matrix(
        reference_class_map,
        description="Reference class map",
    )
    reconstructed_features = np.load(
        io.BytesIO(reconstructed_bytes["ppi-feats.npy"]), allow_pickle=False
    )
    reference_features = np.load(
        io.BytesIO(reference_bytes["ppi-feats.npy"]), allow_pickle=False
    )

    reconstructed_nodes = reconstructed_graph.get("nodes")
    reference_nodes = reference_graph.get("nodes")
    if not isinstance(reconstructed_nodes, list) or not isinstance(reference_nodes, list):
        raise ValidationError("Graph files must contain node lists")

    expected_node_ids = list(range(len(reconstructed_nodes)))
    reconstructed_node_ids = [
        int(node.get("id", -1)) if isinstance(node, dict) else -1
        for node in reconstructed_nodes
    ]
    reconstructed_edges = _edge_multiset(reconstructed_graph)
    reference_edges = _edge_multiset(reference_graph)
    ordered_reconstructed_edges = _ordered_graph_edges(reconstructed_graph)
    expected_ids = [str(index) for index in range(len(reference_nodes))]
    id_map_semantics_exact = (
        set(reconstructed_id_map) == set(expected_ids)
        and all(reconstructed_id_map[key] == int(key) for key in expected_ids)
        and reconstructed_id_map == reference_id_map
    )
    cross_file_rows = {
        len(reconstructed_nodes),
        len(reconstructed_id_map),
        int(reconstructed_features.shape[0]),
        int(reconstructed_labels.shape[0]),
    }
    graph_metadata_fields = ("directed", "graph", "multigraph")
    byte_equality = {
        name: reconstructed_bytes[name] == reference_bytes[name]
        for name in reconstructed_bytes
    }
    checks = {
        "graph_schema_exact": set(reconstructed_graph)
        == {"directed", "graph", "nodes", "links", "multigraph"},
        "graph_node_ids_consecutive": reconstructed_node_ids == expected_node_ids,
        "graph_metadata_exact": all(
            reconstructed_graph.get(field) == reference_graph.get(field)
            for field in graph_metadata_fields
        ),
        "graph_nodes_exact": reconstructed_nodes == reference_nodes,
        "graph_edge_multiset_exact": reconstructed_edges == reference_edges,
        "graph_edges_unique": all(count == 1 for count in reconstructed_edges.values()),
        "graph_edges_canonical_and_sorted": ordered_reconstructed_edges
        == sorted(reconstructed_edges),
        "id_map_semantics_exact": id_map_semantics_exact,
        "id_map_bytes_exact": byte_equality["ppi-id_map.json"],
        "class_map_key_set_exact": reconstructed_keys == reference_keys,
        "class_map_values_exact": np.array_equal(reconstructed_labels, reference_labels),
        "class_map_bytes_exact": byte_equality["ppi-class_map.json"],
        "feature_shape_exact": reconstructed_features.shape == reference_features.shape,
        "feature_dtype_exact": reconstructed_features.dtype == reference_features.dtype,
        "feature_values_exact": np.array_equal(reconstructed_features, reference_features),
        "feature_npy_bytes_exact": byte_equality["ppi-feats.npy"],
        "cross_file_row_counts_agree": len(cross_file_rows) == 1,
        "walk_file_absent": not (graph_path.parent / "ppi-walks.txt").exists(),
    }

    result: dict[str, object] = {
        "schema_version": 1,
        "scope": "complete deterministic files consumed by supervised GraphSAGE",
        "reference": {
            "archive": str(reference_archive.resolve()),
            "archive_sha256": sha256_file(reference_archive),
            "members": {
                name: {
                    "size_bytes": len(data),
                    "sha256": _member_sha256(data),
                }
                for name, data in reference_bytes.items()
            },
        },
        "reconstructed": {
            name: {
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for name, path in {
                "ppi-G.json": graph_path,
                "ppi-id_map.json": id_map_path,
                "ppi-class_map.json": class_map_path,
                "ppi-feats.npy": feature_path,
            }.items()
        },
        "counts": {
            "nodes": len(reconstructed_nodes),
            "edges": sum(reconstructed_edges.values()),
            "feature_cells": int(reconstructed_features.size),
            "label_cells": int(reconstructed_labels.size),
        },
        "comparison_policy": {
            "ppi-G.json": "structural graph equality; canonical link order",
            "ppi-id_map.json": "byte equality",
            "ppi-class_map.json": "byte equality",
            "ppi-feats.npy": "array and byte equality",
        },
        "byte_equality": byte_equality,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    write_json_atomic(json_output, result)
    _write_graphsage_markdown(markdown_output, result)
    if not result["all_checks_pass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise ValidationError("Complete GraphSAGE validation failed: " + ", ".join(failed))
    return result

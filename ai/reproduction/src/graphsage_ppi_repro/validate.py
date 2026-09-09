"""Independent comparison of reconstructed artifacts with released targets.

No reconstruction module imports this file. Keeping all target access here is
an auditable non-circularity boundary: ``pixi run reproduce`` can complete when
reference archives are absent, while ``pixi run validate`` requires them.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np

from .provenance import sha256_file, write_json_atomic, write_text_atomic


class ValidationError(RuntimeError):
    """Raised when reconstructed data disagree with the released reference."""


def _canonical_edge(source: int, target: int) -> tuple[int, int]:
    return (source, target) if source <= target else (target, source)


def _member_sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


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


def _inspect_directed_dgl_graph(
    data: bytes,
    *,
    description: str,
) -> dict[str, object]:
    """Parse one DGL node-link graph into compact validation structures."""

    try:
        graph = json.loads(data)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValidationError(f"Invalid {description} JSON: {exc}") from exc
    if not isinstance(graph, dict):
        raise ValidationError(f"{description} must be a JSON object")
    nodes = graph.get("nodes")
    links = graph.get("links")
    if not isinstance(nodes, list) or not isinstance(links, list):
        raise ValidationError(f"{description} must contain node and link lists")

    node_ids: list[int] = []
    for index, node in enumerate(nodes):
        if not isinstance(node, dict):
            raise ValidationError(f"{description} node {index} is not an object")
        try:
            node_ids.append(int(node["id"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError(f"{description} node {index} has an invalid ID") from exc

    edge_set: set[tuple[int, int]] = set()
    edges_unique = True
    edges_canonical_and_sorted = True
    previous_edge: tuple[int, int] | None = None
    self_loop_count = 0
    for index, link in enumerate(links):
        if not isinstance(link, dict):
            raise ValidationError(f"{description} link {index} is not an object")
        try:
            source = int(link["source"])
            target = int(link["target"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError(
                f"{description} link {index} has invalid endpoints"
            ) from exc
        if not 0 <= source < len(nodes) or not 0 <= target < len(nodes):
            raise ValidationError(
                f"{description} link {index} has an out-of-range endpoint"
            )
        edge = (source, target)
        if edge in edge_set:
            edges_unique = False
        edge_set.add(edge)
        if previous_edge is not None and edge < previous_edge:
            edges_canonical_and_sorted = False
        previous_edge = edge
        if source == target:
            self_loop_count += 1

    return {
        "schema": set(graph),
        "directed": graph.get("directed"),
        "multigraph": graph.get("multigraph"),
        "graph_metadata": graph.get("graph"),
        "node_ids": node_ids,
        "node_count": len(nodes),
        "edge_set": edge_set,
        "edges_unique": edges_unique,
        "edges_canonical_and_sorted": edges_canonical_and_sorted,
        "self_loop_count": self_loop_count,
    }


def _load_npy_bytes(data: bytes, description: str) -> np.ndarray:
    try:
        return np.load(io.BytesIO(data), allow_pickle=False)
    except (OSError, ValueError) as exc:
        raise ValidationError(f"Could not read {description}: {exc}") from exc


def _write_dgl_markdown(path: Path, result: dict[str, object]) -> None:
    splits = result["splits"]
    if not isinstance(splits, dict):
        raise TypeError("DGL validation splits must be a dictionary")
    lines = [
        "# DGL PPI derivative validation",
        "",
        f"Overall result: **{'PASS' if result['all_checks_pass'] else 'FAIL'}**",
        "",
        "The reconstructed DGL-format JSON and NumPy files were compared with",
        "the independently acquired DGL PPI archive. Graph IDs and labels require",
        "exact equality. Features require float64 equality within the predeclared",
        "absolute tolerance. Directed graph links require exact set equality; their",
        "serialization order is canonicalized by the reproduction.",
        "",
        "| Split | Rows | Graphs | Directed edges | Max feature difference | Result |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for split_name in OUTPUT_SPLITS_FOR_VALIDATION:
        split = splits[split_name]
        if not isinstance(split, dict):
            raise TypeError(f"DGL split {split_name} must be a dictionary")
        difference = split["feature_max_abs_difference"]
        difference_text = "n/a" if difference is None else f"{difference:.3g}"
        lines.append(
            f"| `{split_name}` | {split['counts']['rows']:,} | "
            f"{split['counts']['graphs']:,} | {split['counts']['directed_edges']:,} | "
            f"{difference_text} | "
            f"{'PASS' if split['all_checks_pass'] else 'FAIL'} |"
        )

    lines.extend(["", "## Required comparisons", ""])
    for split_name in OUTPUT_SPLITS_FOR_VALIDATION:
        split = splits[split_name]
        lines.extend(
            [
                f"### {split_name}",
                "",
                "| Check | Result |",
                "|---|---:|",
            ]
        )
        for name, passed in split["checks"].items():
            lines.append(f"| `{name}` | {'PASS' if passed else 'FAIL'} |")
        lines.extend(
            [
                "",
                f"- Feature tolerance: atol={result['feature_tolerance']['atol']}, "
                f"rtol={result['feature_tolerance']['rtol']}",
                f"- Feature arrays byte-identical: "
                f"{'YES' if split['byte_equality']['feats'] else 'NO'}",
                f"- Graph JSON byte-identical: "
                f"{'YES' if split['byte_equality']['graph'] else 'NO'}",
                "",
            ]
        )
    write_text_atomic(path, "\n".join(lines))


OUTPUT_SPLITS_FOR_VALIDATION = ("train", "valid", "test")


def validate_dgl_dataset(
    *,
    reference_archive: Path,
    reconstructed_directory: Path,
    json_output: Path,
    markdown_output: Path,
    feature_atol: float = 1.0e-12,
    feature_rtol: float = 0.0,
) -> dict[str, object]:
    """Validate the reconstructed DGL interchange files at meaningful levels."""

    if feature_atol < 0.0 or feature_rtol < 0.0:
        raise ValidationError("DGL feature tolerances must be non-negative")
    try:
        archive = zipfile.ZipFile(reference_archive)
    except zipfile.BadZipFile as exc:
        raise ValidationError(
            f"Invalid DGL reference archive: {reference_archive}"
        ) from exc

    split_results: dict[str, object] = {}
    reconstructed_file_metadata: dict[str, object] = {}
    reference_member_metadata: dict[str, object] = {}
    try:
        available = set(archive.namelist())
        for split_name in OUTPUT_SPLITS_FOR_VALIDATION:
            filenames = {
                "graph": f"{split_name}_graph.json",
                "feats": f"{split_name}_feats.npy",
                "labels": f"{split_name}_labels.npy",
                "graph_id": f"{split_name}_graph_id.npy",
            }
            missing_reference = sorted(set(filenames.values()) - available)
            if missing_reference:
                raise ValidationError(
                    "DGL reference archive is missing: " + ", ".join(missing_reference)
                )

            reconstructed_bytes: dict[str, bytes] = {}
            reference_bytes: dict[str, bytes] = {}
            for kind, filename in filenames.items():
                path = reconstructed_directory / filename
                if not path.is_file():
                    raise ValidationError(f"Missing reconstructed DGL file: {path}")
                reconstructed_bytes[kind] = path.read_bytes()
                reference_bytes[kind] = archive.read(filename)
                reconstructed_file_metadata[filename] = {
                    "path": str(path.resolve()),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                reference_member_metadata[filename] = {
                    "size_bytes": len(reference_bytes[kind]),
                    "sha256": _member_sha256(reference_bytes[kind]),
                }

            reference_graph = _inspect_directed_dgl_graph(
                reference_bytes["graph"],
                description=f"reference {split_name} graph",
            )
            reconstructed_graph = _inspect_directed_dgl_graph(
                reconstructed_bytes["graph"],
                description=f"reconstructed {split_name} graph",
            )
            reference_features = _load_npy_bytes(
                reference_bytes["feats"], f"reference {split_name} features"
            )
            reconstructed_features = _load_npy_bytes(
                reconstructed_bytes["feats"],
                f"reconstructed {split_name} features",
            )
            reference_labels = _load_npy_bytes(
                reference_bytes["labels"], f"reference {split_name} labels"
            )
            reconstructed_labels = _load_npy_bytes(
                reconstructed_bytes["labels"],
                f"reconstructed {split_name} labels",
            )
            reference_graph_ids = _load_npy_bytes(
                reference_bytes["graph_id"],
                f"reference {split_name} graph IDs",
            )
            reconstructed_graph_ids = _load_npy_bytes(
                reconstructed_bytes["graph_id"],
                f"reconstructed {split_name} graph IDs",
            )

            feature_shape_exact = reconstructed_features.shape == reference_features.shape
            if feature_shape_exact and reconstructed_features.size:
                feature_max_difference = float(
                    np.max(np.abs(reconstructed_features - reference_features))
                )
                features_within_tolerance = bool(
                    np.allclose(
                        reconstructed_features,
                        reference_features,
                        atol=feature_atol,
                        rtol=feature_rtol,
                    )
                )
            else:
                feature_max_difference = None
                features_within_tolerance = False

            row_count = int(reconstructed_graph["node_count"])
            expected_node_ids = list(range(row_count))
            reconstructed_edge_set = reconstructed_graph["edge_set"]
            reference_edge_set = reference_graph["edge_set"]
            cross_graph_edges = 0
            if reconstructed_graph_ids.shape == (row_count,):
                for source, target in reconstructed_edge_set:
                    if reconstructed_graph_ids[source] != reconstructed_graph_ids[target]:
                        cross_graph_edges += 1

            checks = {
                "graph_schema_exact": reconstructed_graph["schema"]
                == {"directed", "multigraph", "graph", "nodes", "links"},
                "graph_is_directed_non_multigraph": (
                    reconstructed_graph["directed"] is True
                    and reconstructed_graph["multigraph"] is False
                ),
                "graph_metadata_empty": reconstructed_graph["graph_metadata"] == {},
                "graph_node_ids_consecutive": (
                    reconstructed_graph["node_ids"] == expected_node_ids
                ),
                "graph_nodes_match_reference": (
                    reconstructed_graph["node_ids"] == reference_graph["node_ids"]
                ),
                "graph_edge_set_exact": reconstructed_edge_set == reference_edge_set,
                "graph_edges_unique": bool(reconstructed_graph["edges_unique"]),
                "graph_edges_canonical_and_sorted": bool(
                    reconstructed_graph["edges_canonical_and_sorted"]
                ),
                "one_self_loop_per_node": (
                    reconstructed_graph["self_loop_count"] == row_count
                    and all(
                        (node, node) in reconstructed_edge_set for node in range(row_count)
                    )
                ),
                "no_cross_graph_id_edges": cross_graph_edges == 0,
                "graph_id_shape_exact": (
                    reconstructed_graph_ids.shape
                    == reference_graph_ids.shape
                    == (row_count,)
                ),
                "graph_id_dtype_exact": (
                    reconstructed_graph_ids.dtype
                    == reference_graph_ids.dtype
                    == np.dtype(np.int64)
                ),
                "graph_id_values_exact": np.array_equal(
                    reconstructed_graph_ids, reference_graph_ids
                ),
                "graph_id_npy_bytes_exact": (
                    reconstructed_bytes["graph_id"] == reference_bytes["graph_id"]
                ),
                "label_shape_exact": (
                    reconstructed_labels.shape == reference_labels.shape
                    and reconstructed_labels.shape[0] == row_count
                ),
                "label_dtype_exact": (
                    reconstructed_labels.dtype
                    == reference_labels.dtype
                    == np.dtype(np.int64)
                ),
                "label_values_exact": np.array_equal(
                    reconstructed_labels, reference_labels
                ),
                "label_npy_bytes_exact": (
                    reconstructed_bytes["labels"] == reference_bytes["labels"]
                ),
                "feature_shape_exact": feature_shape_exact
                and reconstructed_features.shape[0] == row_count,
                "feature_dtype_exact": (
                    reconstructed_features.dtype
                    == reference_features.dtype
                    == np.dtype(np.float64)
                ),
                "feature_values_within_tolerance": features_within_tolerance,
                "cross_file_row_counts_agree": (
                    reconstructed_features.shape[0]
                    == reconstructed_labels.shape[0]
                    == reconstructed_graph_ids.shape[0]
                    == row_count
                ),
            }
            values = np.unique(reconstructed_graph_ids)
            split_results[split_name] = {
                "counts": {
                    "rows": row_count,
                    "graphs": int(len(values)),
                    "directed_edges": len(reconstructed_edge_set),
                    "self_loops": int(reconstructed_graph["self_loop_count"]),
                    "cross_graph_id_edges": cross_graph_edges,
                },
                "graph_ids": [int(value) for value in values],
                "feature_max_abs_difference": feature_max_difference,
                "byte_equality": {
                    kind: reconstructed_bytes[kind] == reference_bytes[kind]
                    for kind in filenames
                },
                "checks": checks,
                "all_checks_pass": all(checks.values()),
            }
    finally:
        archive.close()

    result: dict[str, object] = {
        "schema_version": 1,
        "scope": "DGL PPI derivative of the reconstructed GraphSAGE data",
        "reference": {
            "archive": str(reference_archive.resolve()),
            "archive_sha256": sha256_file(reference_archive),
            "members": reference_member_metadata,
        },
        "reconstructed": {
            "directory": str(reconstructed_directory.resolve()),
            "files": reconstructed_file_metadata,
        },
        "comparison_policy": {
            "graph_json": "directed edge-set and node-level structural equality",
            "graph_id": "array and byte equality",
            "labels": "array and byte equality",
            "features": "float64 numeric equality within fixed tolerance",
            "archive_container": "not compared",
        },
        "feature_tolerance": {"atol": feature_atol, "rtol": feature_rtol},
        "splits": split_results,
        "all_checks_pass": all(
            bool(split["all_checks_pass"]) for split in split_results.values()
        ),
    }
    write_json_atomic(json_output, result)
    _write_dgl_markdown(markdown_output, result)
    if not result["all_checks_pass"]:
        failed = [
            f"{split_name}:{name}"
            for split_name, split in split_results.items()
            for name, passed in split["checks"].items()
            if not passed
        ]
        raise ValidationError("DGL validation failed: " + ", ".join(failed))
    return result

"""Derive the DGL PPI files from reconstructed GraphSAGE artifacts.

The released DGL dataset does more than split the 56,944 GraphSAGE rows into
three files.  It gives each tissue's largest connected component its own graph
ID, assigns every smaller component in a split to that split's first graph ID,
and then groups rows by graph ID.  Features are standardized from statistics
fit on the original GraphSAGE training rows.  Finally, undirected links become
directed arcs and every node receives exactly one self-loop.

This module implements that deterministic transformation without importing
DGL.  The released interchange format is ordinary node-link JSON plus NumPy
arrays, so avoiding the framework dependency keeps the reproduction small and
makes each operation directly inspectable.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .provenance import sha256_file, write_json_atomic, write_text_atomic
from .topology import GraphSpec, read_graph_specs

SPEC_TO_OUTPUT_SPLIT = {
    "train": "train",
    "validation": "valid",
    "test": "test",
}
CHECKSUM_FILENAME = "SHA256SUMS"


class DGLReconstructionError(RuntimeError):
    """Raised when GraphSAGE artifacts cannot form the expected DGL data."""


@dataclass(frozen=True)
class SplitPlan:
    """One contiguous GraphSAGE split and its DGL output naming."""

    specification_name: str
    output_name: str
    graph_indices: tuple[int, ...]
    start_row: int
    end_row: int

    @property
    def first_graph_id(self) -> int:
        return self.graph_indices[0]


@dataclass(frozen=True)
class GraphSAGEInputs:
    """Validated GraphSAGE arrays and graph structure used by the transform."""

    features: np.ndarray
    labels: np.ndarray
    nodes: tuple[dict[str, object], ...]
    edges: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class ComponentAssignment:
    """How one original tissue block contributes to DGL graph IDs."""

    graph_index: int
    tissue: str
    split: str
    full_nodes: int
    component_count: int
    largest_component_nodes: int
    non_largest_component_nodes: int
    largest_component_graph_id: int
    non_largest_component_graph_id: int


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


def _write_npy_atomic(path: Path, array: np.ndarray) -> None:
    temporary = _temporary_path(path)
    try:
        with temporary.open("wb") as handle:
            np.save(handle, np.ascontiguousarray(array), allow_pickle=False)
        _replace_after_write(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _hash_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("ascii"))
    return digest.hexdigest()


def _array_data_sha256(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    return hashlib.sha256(contiguous.tobytes(order="C")).hexdigest()


def _read_json_object(path: Path, description: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DGLReconstructionError(
            f"Could not read {description} from {path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise DGLReconstructionError(f"{description} must be a JSON object: {path}")
    return value


def _read_identity_map(path: Path, row_count: int) -> None:
    value = _read_json_object(path, "GraphSAGE ID map")
    expected_keys = {str(index) for index in range(row_count)}
    if set(value) != expected_keys:
        raise DGLReconstructionError(
            "GraphSAGE ID map does not contain exactly the consecutive row keys"
        )
    for index in range(row_count):
        if value[str(index)] != index:
            raise DGLReconstructionError(
                f"GraphSAGE ID map is not the identity at row {index}"
            )


def _read_class_map(path: Path, row_count: int) -> np.ndarray:
    value = _read_json_object(path, "GraphSAGE class map")
    expected_keys = {str(index) for index in range(row_count)}
    if set(value) != expected_keys:
        raise DGLReconstructionError(
            "GraphSAGE class map does not contain exactly the consecutive row keys"
        )
    try:
        labels = np.asarray(
            [value[str(index)] for index in range(row_count)], dtype=np.int64
        )
    except (TypeError, ValueError) as exc:
        raise DGLReconstructionError(
            "GraphSAGE class map does not contain a rectangular integer matrix"
        ) from exc
    if labels.ndim != 2 or not np.all((labels == 0) | (labels == 1)):
        raise DGLReconstructionError("GraphSAGE class-map labels must be a binary matrix")
    return labels


def _read_graphsage_inputs(
    *,
    graph_path: Path,
    id_map_path: Path,
    class_map_path: Path,
    feature_path: Path,
) -> GraphSAGEInputs:
    graph = _read_json_object(graph_path, "GraphSAGE graph")
    if graph.get("directed") is not False or graph.get("multigraph") is not False:
        raise DGLReconstructionError(
            "GraphSAGE input must be an undirected, non-multigraph node-link graph"
        )
    if graph.get("graph") != {}:
        raise DGLReconstructionError("GraphSAGE graph metadata must be empty")

    raw_nodes = graph.get("nodes")
    raw_links = graph.get("links")
    if not isinstance(raw_nodes, list) or not isinstance(raw_links, list):
        raise DGLReconstructionError("GraphSAGE graph must contain node and link lists")
    row_count = len(raw_nodes)
    nodes: list[dict[str, object]] = []
    for expected_id, raw_node in enumerate(raw_nodes):
        if not isinstance(raw_node, dict) or raw_node.get("id") != expected_id:
            raise DGLReconstructionError(
                f"GraphSAGE nodes must have consecutive IDs; failed at {expected_id}"
            )
        nodes.append(raw_node)

    edges: list[tuple[int, int]] = []
    seen_edges: set[tuple[int, int]] = set()
    for index, raw_link in enumerate(raw_links):
        if not isinstance(raw_link, dict):
            raise DGLReconstructionError(f"GraphSAGE link {index} is not an object")
        try:
            source = int(raw_link["source"])
            target = int(raw_link["target"])
        except (KeyError, TypeError, ValueError) as exc:
            raise DGLReconstructionError(
                f"GraphSAGE link {index} has invalid endpoints"
            ) from exc
        if not 0 <= source < row_count or not 0 <= target < row_count:
            raise DGLReconstructionError(
                f"GraphSAGE link {index} has an out-of-range endpoint"
            )
        edge = (source, target) if source <= target else (target, source)
        if edge in seen_edges:
            raise DGLReconstructionError(
                f"GraphSAGE graph contains duplicate undirected edge {edge}"
            )
        seen_edges.add(edge)
        edges.append(edge)

    _read_identity_map(id_map_path, row_count)
    labels = _read_class_map(class_map_path, row_count)
    features = np.load(feature_path, allow_pickle=False)
    if features.ndim != 2 or features.shape[0] != row_count:
        raise DGLReconstructionError(
            f"GraphSAGE feature shape {features.shape} does not align with {row_count} rows"
        )
    if features.dtype != np.dtype(np.float64):
        raise DGLReconstructionError(
            f"Expected float64 GraphSAGE features, observed {features.dtype}"
        )
    if not np.all((features == 0.0) | (features == 1.0)):
        raise DGLReconstructionError("GraphSAGE features must be binary")
    if labels.shape[0] != row_count:
        raise DGLReconstructionError("GraphSAGE feature and label row counts disagree")

    return GraphSAGEInputs(
        features=features,
        labels=labels,
        nodes=tuple(nodes),
        edges=tuple(edges),
    )


def _block_bounds(graph_specs: Sequence[GraphSpec]) -> np.ndarray:
    sizes = np.asarray([spec.node_count for spec in graph_specs], dtype=np.int64)
    return np.concatenate((np.asarray([0], dtype=np.int64), np.cumsum(sizes)))


def _split_plans(graph_specs: Sequence[GraphSpec], bounds: np.ndarray) -> list[SplitPlan]:
    plans: list[SplitPlan] = []
    for specification_name in ("train", "validation", "test"):
        positions = [
            index
            for index, spec in enumerate(graph_specs)
            if spec.split == specification_name
        ]
        if not positions:
            raise DGLReconstructionError(
                f"selected_graphs.tsv has no {specification_name!r} graph"
            )
        expected_positions = list(range(positions[0], positions[-1] + 1))
        if positions != expected_positions:
            raise DGLReconstructionError(
                f"{specification_name!r} graph blocks are not contiguous"
            )
        graph_indices = tuple(graph_specs[index].graph_index_1based for index in positions)
        plans.append(
            SplitPlan(
                specification_name=specification_name,
                output_name=SPEC_TO_OUTPUT_SPLIT[specification_name],
                graph_indices=graph_indices,
                start_row=int(bounds[positions[0]]),
                end_row=int(bounds[positions[-1] + 1]),
            )
        )
    flattened = [index for plan in plans for index in plan.graph_indices]
    expected = [spec.graph_index_1based for spec in graph_specs]
    if flattened != expected:
        raise DGLReconstructionError(
            "selected_graphs.tsv must order training, validation, then test blocks"
        )
    return plans


def _validate_block_metadata(
    inputs: GraphSAGEInputs,
    graph_specs: Sequence[GraphSpec],
    bounds: np.ndarray,
) -> None:
    if int(bounds[-1]) != len(inputs.nodes):
        raise DGLReconstructionError(
            f"Selected graph sizes sum to {int(bounds[-1])}, but GraphSAGE has "
            f"{len(inputs.nodes)} nodes"
        )
    block_for_row = np.empty(len(inputs.nodes), dtype=np.int64)
    for block, spec in enumerate(graph_specs):
        start, end = int(bounds[block]), int(bounds[block + 1])
        block_for_row[start:end] = block
        expected_val = spec.split == "validation"
        expected_test = spec.split == "test"
        for row in range(start, end):
            node = inputs.nodes[row]
            if node.get("val") is not expected_val or node.get("test") is not expected_test:
                raise DGLReconstructionError(
                    f"GraphSAGE split flags disagree with selected_graphs.tsv at row {row}"
                )
    for source, target in inputs.edges:
        if block_for_row[source] != block_for_row[target]:
            raise DGLReconstructionError(
                f"GraphSAGE edge ({source}, {target}) crosses tissue blocks"
            )


def _adjacency(row_count: int, edges: Sequence[tuple[int, int]]) -> list[list[int]]:
    adjacency: list[list[int]] = [[] for _ in range(row_count)]
    for source, target in edges:
        if source == target:
            continue
        adjacency[source].append(target)
        adjacency[target].append(source)
    for neighbors in adjacency:
        neighbors.sort()
    return adjacency


def _connected_components(
    start: int,
    end: int,
    adjacency: Sequence[Sequence[int]],
) -> list[list[int]]:
    seen = np.zeros(end - start, dtype=np.bool_)
    components: list[list[int]] = []
    for root in range(start, end):
        offset = root - start
        if seen[offset]:
            continue
        seen[offset] = True
        stack = [root]
        component: list[int] = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in adjacency[node]:
                neighbor_offset = neighbor - start
                if neighbor_offset < 0 or neighbor_offset >= end - start:
                    raise DGLReconstructionError(
                        f"Adjacency from row {node} leaves its tissue block"
                    )
                if not seen[neighbor_offset]:
                    seen[neighbor_offset] = True
                    stack.append(neighbor)
        components.append(component)
    return components


def _assign_graph_ids(
    graph_specs: Sequence[GraphSpec],
    plans: Sequence[SplitPlan],
    bounds: np.ndarray,
    adjacency: Sequence[Sequence[int]],
) -> tuple[np.ndarray, list[ComponentAssignment]]:
    graph_ids = np.empty(int(bounds[-1]), dtype=np.int64)
    assignments: list[ComponentAssignment] = []
    plan_by_split = {plan.specification_name: plan for plan in plans}

    for plan in plans:
        graph_ids[plan.start_row : plan.end_row] = plan.first_graph_id

    for block, spec in enumerate(graph_specs):
        start, end = int(bounds[block]), int(bounds[block + 1])
        components = _connected_components(start, end, adjacency)
        if not components:
            raise DGLReconstructionError(f"No component found for graph {spec.tissue}")
        largest_index = max(
            range(len(components)), key=lambda index: len(components[index])
        )
        largest = components[largest_index]
        graph_ids[np.asarray(largest, dtype=np.int64)] = spec.graph_index_1based
        first_graph_id = plan_by_split[spec.split].first_graph_id
        assignments.append(
            ComponentAssignment(
                graph_index=spec.graph_index_1based,
                tissue=spec.tissue,
                split=SPEC_TO_OUTPUT_SPLIT[spec.split],
                full_nodes=end - start,
                component_count=len(components),
                largest_component_nodes=len(largest),
                non_largest_component_nodes=(end - start) - len(largest),
                largest_component_graph_id=spec.graph_index_1based,
                non_largest_component_graph_id=first_graph_id,
            )
        )
    return graph_ids, assignments


def _row_order(graph_ids: np.ndarray, plan: SplitPlan) -> np.ndarray:
    pieces = [np.flatnonzero(graph_ids == graph_id) for graph_id in plan.graph_indices]
    order = np.concatenate(pieces).astype(np.int64, copy=False)
    expected = np.arange(plan.start_row, plan.end_row, dtype=np.int64)
    if not np.array_equal(np.sort(order), expected):
        raise DGLReconstructionError(
            f"DGL row order for {plan.output_name} is not a permutation of its split"
        )
    return order


def _standardize_features(
    features: np.ndarray,
    train_start: int,
    train_end: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    training = np.asarray(features[train_start:train_end], dtype=np.float64)
    mean = training.mean(axis=0, dtype=np.float64)
    variance = training.var(axis=0, dtype=np.float64)
    scale = np.sqrt(variance)
    scale[scale == 0.0] = 1.0
    standardized = (np.asarray(features, dtype=np.float64) - mean) / scale
    if not np.all(np.isfinite(standardized)):
        raise DGLReconstructionError("Feature standardization produced non-finite values")
    return standardized, mean, variance, scale


def _directed_edges_for_split(
    edges: Sequence[tuple[int, int]],
    order: np.ndarray,
    row_count: int,
) -> tuple[list[tuple[int, int]], int, int]:
    local_index = np.full(row_count, -1, dtype=np.int64)
    local_index[order] = np.arange(len(order), dtype=np.int64)
    directed: list[tuple[int, int]] = []
    has_loop = np.zeros(len(order), dtype=np.bool_)
    undirected_count = 0
    original_loop_count = 0

    for source, target in edges:
        local_source = int(local_index[source])
        local_target = int(local_index[target])
        if local_source < 0 and local_target < 0:
            continue
        if local_source < 0 or local_target < 0:
            raise DGLReconstructionError(
                f"GraphSAGE edge ({source}, {target}) crosses DGL split boundaries"
            )
        undirected_count += 1
        if local_source == local_target:
            original_loop_count += 1
            has_loop[local_source] = True
            directed.append((local_source, local_source))
        else:
            directed.append((local_source, local_target))
            directed.append((local_target, local_source))

    directed.extend((index, index) for index in np.flatnonzero(~has_loop).tolist())
    directed.sort()
    if any(left == right for left, right in zip(directed, directed[1:], strict=False)):
        raise DGLReconstructionError("DGL directed edge construction produced duplicates")
    return directed, undirected_count, original_loop_count


def _write_graph_json(
    path: Path,
    node_count: int,
    edges: Iterable[tuple[int, int]],
) -> None:
    temporary = _temporary_path(path)
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write('{"directed": true, "multigraph": false, "graph": {}, "nodes": [')
            for node_id in range(node_count):
                if node_id:
                    handle.write(", ")
                json.dump({"id": node_id}, handle)
            handle.write('], "links": [')
            for edge_index, (source, target) in enumerate(edges):
                if edge_index:
                    handle.write(", ")
                json.dump({"source": source, "target": target}, handle)
            handle.write("]}")
        _replace_after_write(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_checksums(path: Path, files: Sequence[Path]) -> None:
    lines = [
        f"{sha256_file(item)}  {item.relative_to(path.parent).as_posix()}" for item in files
    ]
    write_text_atomic(path, "\n".join(lines) + "\n")


def _assignment_dict(assignment: ComponentAssignment) -> dict[str, object]:
    return {
        "graph_index": assignment.graph_index,
        "tissue": assignment.tissue,
        "split": assignment.split,
        "full_nodes": assignment.full_nodes,
        "component_count": assignment.component_count,
        "largest_component_nodes": assignment.largest_component_nodes,
        "non_largest_component_nodes": assignment.non_largest_component_nodes,
        "largest_component_graph_id": assignment.largest_component_graph_id,
        "non_largest_component_graph_id": assignment.non_largest_component_graph_id,
    }


def assemble_dgl(
    *,
    graph_path: Path,
    id_map_path: Path,
    class_map_path: Path,
    feature_path: Path,
    graph_spec_path: Path,
    output_directory: Path,
    summary_output: Path,
) -> dict[str, object]:
    """Build the 12 released-format DGL files from reconstructed GraphSAGE data."""

    inputs = _read_graphsage_inputs(
        graph_path=graph_path,
        id_map_path=id_map_path,
        class_map_path=class_map_path,
        feature_path=feature_path,
    )
    graph_specs = read_graph_specs(graph_spec_path)
    bounds = _block_bounds(graph_specs)
    plans = _split_plans(graph_specs, bounds)
    _validate_block_metadata(inputs, graph_specs, bounds)
    adjacency = _adjacency(len(inputs.nodes), inputs.edges)
    graph_ids, assignments = _assign_graph_ids(graph_specs, plans, bounds, adjacency)

    training_plan = next(plan for plan in plans if plan.output_name == "train")
    standardized, mean, variance, scale = _standardize_features(
        inputs.features,
        training_plan.start_row,
        training_plan.end_row,
    )

    output_directory.mkdir(parents=True, exist_ok=True)
    generated_files: list[Path] = []
    split_summaries: dict[str, object] = {}
    for plan in plans:
        order = _row_order(graph_ids, plan)
        directed_edges, original_edges, original_loops = _directed_edges_for_split(
            inputs.edges,
            order,
            len(inputs.nodes),
        )
        graph_output = output_directory / f"{plan.output_name}_graph.json"
        feature_output = output_directory / f"{plan.output_name}_feats.npy"
        label_output = output_directory / f"{plan.output_name}_labels.npy"
        graph_id_output = output_directory / f"{plan.output_name}_graph_id.npy"

        split_features = np.ascontiguousarray(standardized[order], dtype=np.float64)
        split_labels = np.ascontiguousarray(inputs.labels[order], dtype=np.int64)
        split_graph_ids = np.ascontiguousarray(graph_ids[order], dtype=np.int64)
        _write_graph_json(graph_output, len(order), directed_edges)
        _write_npy_atomic(feature_output, split_features)
        _write_npy_atomic(label_output, split_labels)
        _write_npy_atomic(graph_id_output, split_graph_ids)
        generated_files.extend(
            (graph_output, feature_output, label_output, graph_id_output)
        )

        values, counts = np.unique(split_graph_ids, return_counts=True)
        split_summaries[plan.output_name] = {
            "graph_ids": [int(value) for value in values],
            "graph_id_counts": {
                str(int(value)): int(count)
                for value, count in zip(values, counts, strict=True)
            },
            "rows": int(len(order)),
            "feature_shape": [int(value) for value in split_features.shape],
            "feature_dtype": str(split_features.dtype),
            "label_shape": [int(value) for value in split_labels.shape],
            "label_dtype": str(split_labels.dtype),
            "graph_id_dtype": str(split_graph_ids.dtype),
            "directed_edges": len(directed_edges),
            "self_loops": int(len(order)),
            "original_undirected_edges_including_loops": original_edges,
            "original_self_loops": original_loops,
            "row_order_sha256": _array_data_sha256(order),
            "graph_id_data_sha256": _array_data_sha256(split_graph_ids),
            "label_data_sha256": _array_data_sha256(split_labels),
            "feature_data_sha256": _array_data_sha256(split_features),
            "directed_edges_sha256": _hash_lines(
                f"{source}\t{target}\n" for source, target in directed_edges
            ),
        }

    checksum_path = output_directory.parent / CHECKSUM_FILENAME
    _write_checksums(checksum_path, generated_files)
    output_files = {
        path.name: {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in generated_files
    }
    summary: dict[str, object] = {
        "schema_version": 1,
        "scope": "DGL PPI derivative built only from reconstructed GraphSAGE artifacts",
        "algorithm": {
            "component_rule": (
                "assign each tissue largest connected component to its graph index; "
                "assign all other components to the first graph ID of that split"
            ),
            "row_order": "ascending graph ID, then stable original GraphSAGE row order",
            "feature_scaling": (
                "StandardScaler-compatible population mean and variance fit on all "
                "original GraphSAGE training rows in float64"
            ),
            "edge_conversion": (
                "expand every non-loop undirected edge in both directions and retain "
                "exactly one self-loop per node"
            ),
            "graph_link_order": "lexicographically sorted directed pairs",
            "reference_target_used": False,
        },
        "inputs": {
            "graph": str(graph_path.resolve()),
            "graph_sha256": sha256_file(graph_path),
            "id_map": str(id_map_path.resolve()),
            "id_map_sha256": sha256_file(id_map_path),
            "class_map": str(class_map_path.resolve()),
            "class_map_sha256": sha256_file(class_map_path),
            "features": str(feature_path.resolve()),
            "features_sha256": sha256_file(feature_path),
            "graph_specification": str(graph_spec_path.resolve()),
            "graph_specification_sha256": sha256_file(graph_spec_path),
        },
        "scaler": {
            "fit_split": "train",
            "fit_rows": training_plan.end_row - training_plan.start_row,
            "computation_dtype": "float64",
            "population_variance_ddof": 0,
            "zero_variance_columns_0based": [
                int(index) for index in np.flatnonzero(variance == 0.0)
            ],
            "mean": mean.tolist(),
            "variance": variance.tolist(),
            "scale": scale.tolist(),
        },
        "component_assignments": [
            _assignment_dict(assignment) for assignment in assignments
        ],
        "splits": split_summaries,
        "outputs": {
            "directory": str(output_directory.resolve()),
            "files": output_files,
            "checksums": str(checksum_path.resolve()),
            "checksums_sha256": sha256_file(checksum_path),
        },
    }
    write_json_atomic(summary_output, summary)
    return summary

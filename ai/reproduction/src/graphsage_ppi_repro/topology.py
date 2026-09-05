"""Reconstruct the 24 GraphSAGE PPI graph blocks and biological row identities."""

from __future__ import annotations

import csv
import gzip
import hashlib
import tarfile
from collections import Counter
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from .legacy_order import Python2InsertionDict, python2_string_hash
from .provenance import sha256_file, write_json_atomic, write_tsv_atomic

GRAPH_SPEC_FIELDS = (
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
VALID_SPLITS = frozenset({"train", "validation", "test"})


class TopologyError(RuntimeError):
    """Raised when source topology disagrees with the frozen specification."""


@dataclass(frozen=True)
class GraphSpec:
    """One accepted GraphSAGE graph block and its OhmNet source member."""

    graph_index_1based: int
    tissue: str
    split: str
    ohmnet_member: str
    ohmnet_member_sha256: str
    node_count: int
    edge_count: int
    evidence_status: str
    evidence_reference: str


@dataclass(frozen=True)
class ParsedNetwork:
    """One selected OhmNet layer after legacy node-order reconstruction."""

    spec: GraphSpec
    raw_edges: tuple[tuple[str, str], ...]
    ordered_nodes: tuple[tuple[int, str], ...]
    dict_table_size: int
    member_sha256: str


def read_graph_specs(path: Path) -> list[GraphSpec]:
    """Read the explicit data-level graph identity, order, and split table."""

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != GRAPH_SPEC_FIELDS:
            raise TopologyError(
                f"Unexpected columns in {path}: {reader.fieldnames!r}; "
                f"expected {list(GRAPH_SPEC_FIELDS)!r}"
            )

        graphs: list[GraphSpec] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                graph = GraphSpec(
                    graph_index_1based=int(row["graph_index_1based"]),
                    tissue=row["tissue"].strip(),
                    split=row["split"].strip(),
                    ohmnet_member=row["ohmnet_member"].strip(),
                    ohmnet_member_sha256=row["ohmnet_member_sha256"].strip().lower(),
                    node_count=int(row["node_count"]),
                    edge_count=int(row["edge_count"]),
                    evidence_status=row["evidence_status"].strip(),
                    evidence_reference=row["evidence_reference"].strip(),
                )
            except (KeyError, ValueError) as exc:
                raise TopologyError(
                    f"Invalid graph row {line_number} in {path}: {exc}"
                ) from exc

            location = f"{path}:{line_number}"
            if graph.split not in VALID_SPLITS:
                raise TopologyError(f"Invalid split {graph.split!r} at {location}")
            if not graph.tissue or not graph.ohmnet_member:
                raise TopologyError(f"Blank tissue or archive member at {location}")
            if graph.node_count < 1 or graph.edge_count < 1:
                raise TopologyError(f"Non-positive graph dimensions at {location}")
            if len(graph.ohmnet_member_sha256) != 64 or any(
                character not in "0123456789abcdef"
                for character in graph.ohmnet_member_sha256
            ):
                raise TopologyError(f"Invalid member SHA-256 at {location}")
            graphs.append(graph)

    expected_indices = list(range(1, len(graphs) + 1))
    observed_indices = [graph.graph_index_1based for graph in graphs]
    if observed_indices != expected_indices:
        raise TopologyError(
            "Graph indices must be consecutive and file-ordered; "
            f"observed {observed_indices}"
        )
    if len({graph.tissue for graph in graphs}) != len(graphs):
        raise TopologyError("selected_graphs.tsv contains duplicate tissue names")
    if len({graph.ohmnet_member for graph in graphs}) != len(graphs):
        raise TopologyError("selected_graphs.tsv contains duplicate archive members")
    return graphs


def _parse_edgelist(data: bytes, *, member_name: str) -> tuple[tuple[str, str], ...]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TopologyError(f"OhmNet member is not ASCII: {member_name}") from exc

    edges: list[tuple[str, str]] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 2:
            raise TopologyError(f"Malformed edge at {member_name}:{line_number}")
        source, target = fields[0], fields[1]
        try:
            int(source)
            int(target)
        except ValueError as exc:
            raise TopologyError(
                f"Non-integer Entrez identifier at {member_name}:{line_number}"
            ) from exc
        edges.append((source, target))
    if not edges:
        raise TopologyError(f"No edges found in {member_name}")
    return tuple(edges)


def _canonical_entrez_edge(source: str, target: str) -> tuple[int, int]:
    left, right = int(source), int(target)
    return (left, right) if left <= right else (right, left)


def parse_selected_networks(
    archive_path: Path,
    graph_specs: Sequence[GraphSpec],
    *,
    word_size_bits: int = 64,
) -> list[ParsedNetwork]:
    """Read selected OhmNet members and reconstruct their legacy node order."""

    parsed: list[ParsedNetwork] = []
    with tarfile.open(archive_path, mode="r:*") as archive:
        available = {member.name for member in archive.getmembers() if member.isfile()}
        for spec in graph_specs:
            if spec.ohmnet_member not in available:
                raise TopologyError(
                    f"Archive {archive_path} lacks selected member {spec.ohmnet_member}"
                )
            extracted = archive.extractfile(spec.ohmnet_member)
            if extracted is None:
                raise TopologyError(
                    f"Could not read {spec.ohmnet_member} from {archive_path}"
                )
            data = extracted.read()
            member_hash = hashlib.sha256(data).hexdigest()
            if member_hash != spec.ohmnet_member_sha256:
                raise TopologyError(
                    f"Decompressed SHA-256 mismatch for {spec.ohmnet_member}: "
                    f"expected {spec.ohmnet_member_sha256}, observed {member_hash}"
                )

            edges = _parse_edgelist(data, member_name=spec.ohmnet_member)
            canonical_edges = {
                _canonical_entrez_edge(source, target) for source, target in edges
            }
            if len(canonical_edges) != len(edges):
                raise TopologyError(
                    f"{spec.ohmnet_member} contains duplicate undirected edge records"
                )
            if len(edges) != spec.edge_count:
                raise TopologyError(
                    f"Edge count mismatch for {spec.tissue}: expected {spec.edge_count}, "
                    f"observed {len(edges)}"
                )

            node_table: Python2InsertionDict[str] = Python2InsertionDict(
                python2_string_hash, word_size_bits=word_size_bits
            )
            for source, target in edges:
                node_table.insert(source)
                node_table.insert(target)
            ordered_nodes = tuple(node_table.slots_and_keys())
            if len(ordered_nodes) != spec.node_count:
                raise TopologyError(
                    f"Node count mismatch for {spec.tissue}: expected {spec.node_count}, "
                    f"observed {len(ordered_nodes)}"
                )

            parsed.append(
                ParsedNetwork(
                    spec=spec,
                    raw_edges=edges,
                    ordered_nodes=ordered_nodes,
                    dict_table_size=node_table.table_size,
                    member_sha256=member_hash,
                )
            )
    return parsed


def _open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def read_node_mapping(path: Path) -> list[dict[str, str]]:
    """Read a reconstructed node mapping in plain or gzip-compressed TSV form."""

    with _open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _hash_lines(lines: Iterator[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("ascii"))
    return digest.hexdigest()


def reconstruct_topology(
    *,
    archive_path: Path,
    graph_spec_path: Path,
    mapping_output: Path,
    edge_output: Path,
    summary_output: Path,
    word_size_bits: int = 64,
) -> dict[str, object]:
    """Reconstruct all rows and logical undirected edges from OhmNet."""

    specs = read_graph_specs(graph_spec_path)
    networks = parse_selected_networks(
        archive_path, specs, word_size_bits=word_size_bits
    )

    mapping_rows: list[dict[str, object]] = []
    graph_rows: list[dict[str, object]] = []
    local_maps: list[dict[str, int]] = []
    next_global_id = 0
    all_gene_ids: set[int] = set()

    for network in networks:
        start = next_global_id
        local_by_gene: dict[str, int] = {}
        for local_index, (slot, gene_text) in enumerate(network.ordered_nodes):
            global_id = start + local_index
            local_by_gene[gene_text] = global_id
            gene_id = int(gene_text)
            all_gene_ids.add(gene_id)
            mapping_rows.append(
                {
                    "graphsage_node_id": global_id,
                    "feature_label_row_index": global_id,
                    "graph_index_1based": network.spec.graph_index_1based,
                    "tissue": network.spec.tissue,
                    "split": network.spec.split,
                    "local_node_index_0based": local_index,
                    "entrez_gene_id": gene_id,
                    "python2_dict_table_slot": slot,
                    "python2_dict_table_size": network.dict_table_size,
                }
            )
        next_global_id += len(network.ordered_nodes)
        local_maps.append(local_by_gene)
        graph_rows.append(
            {
                "graph_index_1based": network.spec.graph_index_1based,
                "tissue": network.spec.tissue,
                "split": network.spec.split,
                "row_start_inclusive": start,
                "row_end_exclusive": next_global_id,
                "node_count": len(network.ordered_nodes),
                "edge_count": len(network.raw_edges),
                "python2_dict_table_size": network.dict_table_size,
                "ohmnet_member": network.spec.ohmnet_member,
                "ohmnet_member_sha256": network.member_sha256,
            }
        )

    mapping_fields = (
        "graphsage_node_id",
        "feature_label_row_index",
        "graph_index_1based",
        "tissue",
        "split",
        "local_node_index_0based",
        "entrez_gene_id",
        "python2_dict_table_slot",
        "python2_dict_table_size",
    )
    write_tsv_atomic(mapping_output, mapping_fields, mapping_rows)

    edge_fields = (
        "graph_index_1based",
        "tissue",
        "source_line_1based",
        "source_node_id",
        "target_node_id",
        "source_entrez_gene_id",
        "target_entrez_gene_id",
    )
    edge_sequence_digest = hashlib.sha256()

    def edge_rows() -> Iterator[dict[str, object]]:
        for network, local_by_gene in zip(networks, local_maps, strict=True):
            for source_line_1based, (source_gene, target_gene) in enumerate(
                network.raw_edges, start=1
            ):
                source_node_id = local_by_gene[source_gene]
                target_node_id = local_by_gene[target_gene]
                edge_sequence_digest.update(
                    f"{source_node_id}\t{target_node_id}\n".encode("ascii")
                )
                yield {
                    "graph_index_1based": network.spec.graph_index_1based,
                    "tissue": network.spec.tissue,
                    "source_line_1based": source_line_1based,
                    "source_node_id": source_node_id,
                    "target_node_id": target_node_id,
                    "source_entrez_gene_id": int(source_gene),
                    "target_entrez_gene_id": int(target_gene),
                }

    write_tsv_atomic(edge_output, edge_fields, edge_rows())

    split_graph_counts = Counter(network.spec.split for network in networks)
    split_row_counts = Counter(str(row["split"]) for row in mapping_rows)
    row_gene_hash = _hash_lines(
        iter(f"{row['entrez_gene_id']}\n" for row in mapping_rows)
    )
    node_gene_pair_hash = _hash_lines(
        iter(
            f"{row['graphsage_node_id']}\t{row['entrez_gene_id']}\n"
            for row in mapping_rows
        )
    )
    edge_count = sum(len(network.raw_edges) for network in networks)

    summary: dict[str, object] = {
        "schema_version": 1,
        "algorithm": {
            "read_edgelists_in_source_line_order": True,
            "entrez_identifiers_kept_as_strings_during_insertion": True,
            "dictionary_model": (
                f"{word_size_bits}-bit unrandomized CPython 2.7 insertion-only dict"
            ),
            "global_rows_are_concatenated_in_selected_graph_order": True,
        },
        "inputs": {
            "ohmnet_archive": str(archive_path.resolve()),
            "ohmnet_archive_sha256": sha256_file(archive_path),
            "selected_graphs": str(graph_spec_path.resolve()),
            "selected_graphs_sha256": sha256_file(graph_spec_path),
        },
        "counts": {
            "graphs": len(networks),
            "rows": len(mapping_rows),
            "distinct_entrez_gene_ids": len(all_gene_ids),
            "edge_records": edge_count,
            "split_graph_counts": dict(sorted(split_graph_counts.items())),
            "split_row_counts": dict(sorted(split_row_counts.items())),
        },
        "content_hashes": {
            "entrez_row_sequence_sha256": row_gene_hash,
            "node_id_entrez_pairs_sha256": node_gene_pair_hash,
            "source_order_node_edge_pairs_sha256": edge_sequence_digest.hexdigest(),
        },
        "graphs": graph_rows,
        "outputs": {
            "node_mapping": str(mapping_output.resolve()),
            "node_mapping_sha256": sha256_file(mapping_output),
            "edges": str(edge_output.resolve()),
            "edges_sha256": sha256_file(edge_output),
        },
    }
    write_json_atomic(summary_output, summary)
    return summary

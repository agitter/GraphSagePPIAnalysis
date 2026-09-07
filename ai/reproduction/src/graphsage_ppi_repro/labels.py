"""Reconstruct the 121 GraphSAGE PPI GO-label columns from dated sources.

The implementation follows one global policy.  It does not inspect the released
GraphSAGE class map while reconstructing labels.  The deposited target is read
only by :mod:`graphsage_ppi_repro.validate` after reconstruction is complete.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

import numpy as np

from .identifiers import IdentifierError, build_identifier_mapping
from .provenance import (
    sha256_file,
    write_json_atomic,
    write_tsv_atomic,
)
from .topology import read_node_mapping

LABEL_COLUMN_FIELDS = (
    "column_index",
    "go_id",
    "namespace",
    "duplicate_vector_group",
    "identity_status",
    "evidence_reference",
)
VALID_NAMESPACES = frozenset(
    {"biological_process", "cellular_component", "molecular_function"}
)
ASPECT_RELATIONS = {
    "P": "involved_in",
    "C": "part_of",
    "F": "enables",
}
GAF_COLUMN_COUNT = 17


class LabelError(RuntimeError):
    """Raised when GO sources or reconstructed labels violate the specification."""


@dataclass(frozen=True)
class LabelColumn:
    """One released GraphSAGE label column."""

    column_index: int
    go_id: str
    namespace: str
    duplicate_vector_group: str
    identity_status: str
    evidence_reference: str


@dataclass(frozen=True)
class OntologyTerm:
    """The OBO fields needed for label propagation and reporting."""

    go_id: str
    name: str
    namespace: str
    alt_ids: tuple[str, ...]
    is_a_parents: tuple[str, ...]
    is_obsolete: bool


@dataclass
class Ontology:
    """A minimal GO ontology with alternate-ID and ``is_a`` traversal."""

    terms: dict[str, OntologyTerm]
    alt_to_primary: dict[str, str]
    _ancestor_cache: dict[str, frozenset[str]] = field(default_factory=dict)

    def canonical(self, go_id: str) -> str:
        """Return the primary GO ID when *go_id* is an alternate ID."""

        return self.alt_to_primary.get(go_id, go_id)

    def ancestors(self, go_id: str) -> frozenset[str]:
        """Return the canonical term and every transitive ``is_a`` ancestor."""

        canonical = self.canonical(go_id)
        return self._ancestors(canonical, visiting=set())

    def _ancestors(self, go_id: str, *, visiting: set[str]) -> frozenset[str]:
        cached = self._ancestor_cache.get(go_id)
        if cached is not None:
            return cached
        if go_id in visiting:
            raise LabelError(f"Cycle detected in GO is_a ancestry at {go_id}")
        visiting.add(go_id)
        found = {go_id}
        term = self.terms.get(go_id)
        if term is not None:
            for parent in term.is_a_parents:
                found.update(self._ancestors(self.canonical(parent), visiting=visiting))
        visiting.remove(go_id)
        result = frozenset(found)
        self._ancestor_cache[go_id] = result
        return result


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_label_columns(path: Path) -> list[LabelColumn]:
    """Read and validate the evidence-backed released column specification."""

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != LABEL_COLUMN_FIELDS:
            raise LabelError(
                f"Unexpected columns in {path}: {reader.fieldnames!r}; "
                f"expected {list(LABEL_COLUMN_FIELDS)!r}"
            )
        columns: list[LabelColumn] = []
        for line_number, row in enumerate(reader, start=2):
            try:
                column = LabelColumn(
                    column_index=int(row["column_index"]),
                    go_id=row["go_id"].strip(),
                    namespace=row["namespace"].strip(),
                    duplicate_vector_group=row["duplicate_vector_group"].strip(),
                    identity_status=row["identity_status"].strip(),
                    evidence_reference=row["evidence_reference"].strip(),
                )
            except (KeyError, ValueError) as exc:
                raise LabelError(
                    f"Invalid label column at {path}:{line_number}: {exc}"
                ) from exc
            location = f"{path}:{line_number}"
            if not column.go_id.startswith("GO:"):
                raise LabelError(f"Invalid GO ID at {location}: {column.go_id!r}")
            if column.namespace not in VALID_NAMESPACES:
                raise LabelError(f"Invalid namespace {column.namespace!r} at {location}")
            if not column.identity_status:
                raise LabelError(f"Missing identity status at {location}")
            columns.append(column)

    expected_indices = list(range(len(columns)))
    observed_indices = [column.column_index for column in columns]
    if observed_indices != expected_indices:
        raise LabelError(
            "Label column indices must be consecutive and file-ordered; "
            f"observed {observed_indices}"
        )
    go_ids = [column.go_id for column in columns]
    if len(set(go_ids)) != len(go_ids):
        raise LabelError(f"Duplicate primary GO IDs in {path}")
    return columns


def _finalize_obo_term(
    values: Mapping[str, object],
    *,
    terms: dict[str, OntologyTerm],
    alt_to_primary: dict[str, str],
    path: Path,
) -> None:
    if not values or "id" not in values:
        return
    go_id = str(values["id"])
    alt_ids = tuple(str(value) for value in values.get("alt_id", ()))
    parents = tuple(str(value).split(maxsplit=1)[0] for value in values.get("is_a", ()))
    term = OntologyTerm(
        go_id=go_id,
        name=str(values.get("name", "")),
        namespace=str(values.get("namespace", "")),
        alt_ids=alt_ids,
        is_a_parents=parents,
        is_obsolete=str(values.get("is_obsolete", "false")) == "true",
    )
    if go_id in terms:
        raise LabelError(f"Duplicate ontology term {go_id} in {path}")
    terms[go_id] = term
    for alt_id in alt_ids:
        previous = alt_to_primary.get(alt_id)
        if previous is not None and previous != go_id:
            raise LabelError(
                f"Alternate GO ID {alt_id} maps to both {previous} and {go_id}"
            )
        alt_to_primary[alt_id] = go_id


def read_ontology(path: Path) -> Ontology:
    """Parse GO terms, alternate IDs, and direct ``is_a`` parent edges."""

    terms: dict[str, OntologyTerm] = {}
    alt_to_primary: dict[str, str] = {}
    stanza: dict[str, object] | None = None

    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line == "[Term]":
                if stanza is not None:
                    _finalize_obo_term(
                        stanza,
                        terms=terms,
                        alt_to_primary=alt_to_primary,
                        path=path,
                    )
                stanza = {}
                continue
            if line.startswith("["):
                if stanza is not None:
                    _finalize_obo_term(
                        stanza,
                        terms=terms,
                        alt_to_primary=alt_to_primary,
                        path=path,
                    )
                stanza = None
                continue
            if stanza is None:
                continue
            if not line:
                _finalize_obo_term(
                    stanza,
                    terms=terms,
                    alt_to_primary=alt_to_primary,
                    path=path,
                )
                stanza = None
                continue
            if ": " not in line:
                continue
            key, value = line.split(": ", 1)
            if key in {"alt_id", "is_a"}:
                sequence = stanza.setdefault(key, [])
                if not isinstance(sequence, list):
                    raise LabelError(f"Malformed repeated OBO field {key} in {path}")
                sequence.append(value)
            else:
                stanza[key] = value

    if stanza is not None:
        _finalize_obo_term(
            stanza,
            terms=terms,
            alt_to_primary=alt_to_primary,
            path=path,
        )
    if not terms:
        raise LabelError(f"No [Term] stanzas found in {path}")
    return Ontology(terms=terms, alt_to_primary=alt_to_primary)


def _normalized_relation(qualifier: str, aspect: str) -> tuple[str | None, bool]:
    values = [value for value in qualifier.split("|") if value]
    if "NOT" in values:
        return None, True
    relations = [value for value in values if value != "NOT"]
    if relations:
        return "|".join(relations), False
    try:
        return ASPECT_RELATIONS[aspect], False
    except KeyError as exc:
        raise LabelError(f"Unknown GAF aspect {aspect!r}") from exc


def build_term_gene_memberships(
    *,
    gaf_path: Path,
    accession_to_genes: Mapping[str, frozenset[int]],
    ontology: Ontology,
    evidence_codes: frozenset[str],
    allowed_relations: frozenset[str],
) -> tuple[dict[str, set[int]], dict[str, int]]:
    """Apply the accepted GAF filters and build propagated term memberships."""

    term_to_genes: dict[str, set[int]] = defaultdict(set)
    counters: Counter[str] = Counter()
    with _open_text(gaf_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("!") or not raw_line.strip():
                continue
            counters["rows_total"] += 1
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != GAF_COLUMN_COUNT:
                raise LabelError(
                    f"Expected {GAF_COLUMN_COUNT} GAF fields at "
                    f"{gaf_path}:{line_number}; observed {len(fields)}"
                )
            (
                database,
                accession,
                _symbol,
                qualifier,
                go_id,
                _reference,
                evidence,
                _with_from,
                aspect,
                _name,
                _synonyms,
                _object_type,
                _taxon,
                _date,
                _assigned_by,
                _annotation_extension,
                _gene_product_form_id,
            ) = fields
            if database != "UniProtKB":
                counters["skipped_database"] += 1
                continue
            relation, negated = _normalized_relation(qualifier, aspect)
            if negated:
                counters["skipped_not"] += 1
                continue
            if evidence not in evidence_codes:
                counters["skipped_evidence"] += 1
                continue
            if relation not in allowed_relations:
                counters["skipped_relation"] += 1
                continue
            genes = accession_to_genes.get(accession, frozenset())
            if not genes:
                counters["skipped_unmapped_accession"] += 1
                continue
            canonical = ontology.canonical(go_id)
            if canonical not in ontology.terms:
                counters["accepted_unknown_ontology_term"] += 1
            for propagated_term in ontology.ancestors(canonical):
                term_to_genes[propagated_term].update(genes)
            counters["accepted_rows"] += 1
    return dict(term_to_genes), dict(sorted(counters.items()))


def _rank_terms(
    term_to_genes: Mapping[str, set[int]],
) -> list[tuple[int, str]]:
    return sorted(
        ((len(genes), go_id) for go_id, genes in term_to_genes.items()),
        key=lambda item: (-item[0], item[1]),
    )


def _write_npy_atomic(path: Path, matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("wb") as handle:
            np.save(handle, matrix, allow_pickle=False)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_class_map_atomic(path: Path, matrix: np.ndarray) -> None:
    """Write a numeric-node-ordered JSON class map without a large Python dict."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("{")
            for row_index, row in enumerate(matrix):
                if row_index:
                    handle.write(",")
                json.dump(str(row_index), handle)
                handle.write(":")
                json.dump(row.tolist(), handle, separators=(",", ":"))
            handle.write("}\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _sha256_lines(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _validate_columns_against_ontology(
    columns: Sequence[LabelColumn], ontology: Ontology
) -> None:
    for column in columns:
        canonical = ontology.canonical(column.go_id)
        if canonical != column.go_id:
            raise LabelError(
                f"Label column {column.column_index} uses alternate GO ID "
                f"{column.go_id}; canonical ID is {canonical}"
            )
        term = ontology.terms.get(column.go_id)
        if term is None:
            raise LabelError(f"Label column GO ID absent from ontology: {column.go_id}")
        if term.namespace != column.namespace:
            raise LabelError(
                f"Namespace mismatch for {column.go_id}: specification has "
                f"{column.namespace}, ontology has {term.namespace}"
            )


def _validate_duplicate_groups(
    columns: Sequence[LabelColumn], matrix: np.ndarray
) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = defaultdict(list)
    for column in columns:
        if column.duplicate_vector_group:
            groups[column.duplicate_vector_group].append(column.column_index)
    for group, indices in groups.items():
        expected_go_ids = set(group.split("|"))
        observed_go_ids = {columns[index].go_id for index in indices}
        if observed_go_ids != expected_go_ids:
            raise LabelError(
                f"Duplicate-vector group {group!r} contains columns {observed_go_ids}"
            )
        first = matrix[:, indices[0]]
        if any(not np.array_equal(first, matrix[:, index]) for index in indices[1:]):
            raise LabelError(f"Expected duplicate membership vectors differ: {group}")
    return dict(sorted(groups.items()))


def reconstruct_labels(
    *,
    mapping_path: Path,
    gaf_path: Path,
    gpi_path: Path,
    gp2protein_path: Path,
    ontology_path: Path,
    label_columns_path: Path,
    identifier_decisions_path: Path,
    matrix_output: Path,
    class_map_output: Path,
    selected_labels_output: Path,
    summary_output: Path,
    evidence_codes: Iterable[str],
    allowed_relations: Iterable[str],
    selected_term_count: int = 121,
) -> dict[str, object]:
    """Reconstruct all row labels and write deterministic intermediate artifacts."""

    if selected_term_count < 1:
        raise LabelError("selected_term_count must be positive")
    evidence = frozenset(evidence_codes)
    relations = frozenset(allowed_relations)
    if not evidence or not relations:
        raise LabelError("Evidence-code and relation sets must be nonempty")

    mapping_rows = read_node_mapping(mapping_path)
    node_ids = [int(row["graphsage_node_id"]) for row in mapping_rows]
    expected_node_ids = list(range(len(mapping_rows)))
    if node_ids != expected_node_ids:
        raise LabelError("Node mapping must be ordered by consecutive GraphSAGE node ID")
    gene_ids = [int(row["entrez_gene_id"]) for row in mapping_rows]
    graph_genes = set(gene_ids)

    columns = read_label_columns(label_columns_path)
    if len(columns) != selected_term_count:
        raise LabelError(
            f"Expected {selected_term_count} label columns; observed {len(columns)}"
        )
    ontology = read_ontology(ontology_path)
    _validate_columns_against_ontology(columns, ontology)

    try:
        _products, accession_to_genes, decisions, identifier_statistics = (
            build_identifier_mapping(
                gpi_path=gpi_path,
                gp2protein_path=gp2protein_path,
                decisions_path=identifier_decisions_path,
            )
        )
    except IdentifierError as exc:
        raise LabelError(str(exc)) from exc

    term_to_genes, gaf_statistics = build_term_gene_memberships(
        gaf_path=gaf_path,
        accession_to_genes=accession_to_genes,
        ontology=ontology,
        evidence_codes=evidence,
        allowed_relations=relations,
    )
    ranked_terms = _rank_terms(term_to_genes)
    if len(ranked_terms) < selected_term_count:
        raise LabelError(
            f"Only {len(ranked_terms)} GO terms had mapped annotations; "
            f"cannot select {selected_term_count}"
        )
    derived_terms = [go_id for _count, go_id in ranked_terms[:selected_term_count]]
    specified_terms = [column.go_id for column in columns]
    if set(derived_terms) != set(specified_terms):
        missing = sorted(set(specified_terms) - set(derived_terms))
        unexpected = sorted(set(derived_terms) - set(specified_terms))
        raise LabelError(
            "Top-term selection disagrees with label_columns.tsv; "
            f"missing={missing}, unexpected={unexpected}"
        )
    rank_by_term = {go_id: rank for rank, (_count, go_id) in enumerate(ranked_terms, 1)}

    matrix = np.zeros((len(gene_ids), len(columns)), dtype=np.uint8)
    for column in columns:
        member_genes = term_to_genes.get(column.go_id, set())
        matrix[:, column.column_index] = np.fromiter(
            (gene_id in member_genes for gene_id in gene_ids),
            dtype=np.uint8,
            count=len(gene_ids),
        )

    duplicate_groups = _validate_duplicate_groups(columns, matrix)
    repeated_gene_vectors: dict[int, bytes] = {}
    repeated_gene_conflicts: list[int] = []
    for row_index, gene_id in enumerate(gene_ids):
        vector = matrix[row_index].tobytes()
        previous = repeated_gene_vectors.setdefault(gene_id, vector)
        if previous != vector:
            repeated_gene_conflicts.append(gene_id)
    if repeated_gene_conflicts:
        raise LabelError(
            "Repeated GeneIDs received inconsistent label vectors: "
            f"{sorted(set(repeated_gene_conflicts))[:20]}"
        )

    _write_npy_atomic(matrix_output, matrix)
    _write_class_map_atomic(class_map_output, matrix)

    selected_rows: list[dict[str, object]] = []
    namespace_counts: Counter[str] = Counter()
    for column in columns:
        term = ontology.terms[column.go_id]
        namespace_counts[column.namespace] += 1
        selected_rows.append(
            {
                "column_index": column.column_index,
                "prevalence_rank": rank_by_term[column.go_id],
                "go_id": column.go_id,
                "go_name": term.name,
                "namespace": column.namespace,
                "historical_gene_count": len(term_to_genes.get(column.go_id, set())),
                "graph_gene_count": len(
                    term_to_genes.get(column.go_id, set()) & graph_genes
                ),
                "positive_rows": int(matrix[:, column.column_index].sum()),
                "duplicate_vector_group": column.duplicate_vector_group,
                "identity_status": column.identity_status,
                "evidence_reference": column.evidence_reference,
            }
        )
    write_tsv_atomic(
        selected_labels_output,
        (
            "column_index",
            "prevalence_rank",
            "go_id",
            "go_name",
            "namespace",
            "historical_gene_count",
            "graph_gene_count",
            "positive_rows",
            "duplicate_vector_group",
            "identity_status",
            "evidence_reference",
        ),
        selected_rows,
    )

    mapped_graph_genes = graph_genes & set().union(*accession_to_genes.values())
    data_hash = hashlib.sha256(matrix.tobytes(order="C")).hexdigest()
    distinct_vectors = len(
        {matrix[:, column_index].tobytes() for column_index in range(matrix.shape[1])}
    )
    boundary_start = max(0, selected_term_count - 2)
    boundary_end = min(len(ranked_terms), selected_term_count + 3)
    prevalence_boundary = [
        {
            "rank": rank,
            "go_id": go_id,
            "historical_gene_count": count,
        }
        for rank, (count, go_id) in enumerate(
            ranked_terms[boundary_start:boundary_end], start=boundary_start + 1
        )
    ]

    summary: dict[str, object] = {
        "schema_version": 1,
        "algorithm": {
            "evidence_codes": sorted(evidence),
            "exclude_not": True,
            "ordinary_relations": sorted(relations),
            "canonicalize_alternate_go_ids": True,
            "ontology_propagation": ["is_a"],
            "term_selection": f"top_{selected_term_count}_by_distinct_geneid_prevalence",
            "released_column_order": str(label_columns_path.resolve()),
        },
        "inputs": {
            "node_mapping": str(mapping_path.resolve()),
            "node_mapping_sha256": sha256_file(mapping_path),
            "gaf": str(gaf_path.resolve()),
            "gaf_sha256": sha256_file(gaf_path),
            "gpi": str(gpi_path.resolve()),
            "gpi_sha256": sha256_file(gpi_path),
            "gp2protein": str(gp2protein_path.resolve()),
            "gp2protein_sha256": sha256_file(gp2protein_path),
            "ontology": str(ontology_path.resolve()),
            "ontology_sha256": sha256_file(ontology_path),
            "label_columns": str(label_columns_path.resolve()),
            "label_columns_sha256": sha256_file(label_columns_path),
            "identifier_decisions": str(identifier_decisions_path.resolve()),
            "identifier_decisions_sha256": sha256_file(identifier_decisions_path),
        },
        "counts": {
            "rows": int(matrix.shape[0]),
            "columns": int(matrix.shape[1]),
            "cells": int(matrix.size),
            "positive_cells": int(matrix.sum()),
            "distinct_graph_gene_ids": len(graph_genes),
            "mapped_graph_gene_ids": len(mapped_graph_genes),
            "unmapped_graph_gene_ids": sorted(graph_genes - mapped_graph_genes),
            "selected_by_namespace": dict(sorted(namespace_counts.items())),
            "distinct_column_membership_vectors": distinct_vectors,
            "duplicate_vector_groups": duplicate_groups,
            "repeated_gene_vector_conflicts": 0,
            "ontology_terms": len(ontology.terms),
            "ontology_alternate_ids": len(ontology.alt_to_primary),
        },
        "identifier_mapping": {
            **identifier_statistics,
            "decision_ids": [decision.decision_id for decision in decisions],
        },
        "gaf_processing": gaf_statistics,
        "term_selection": {
            "derived_set_matches_column_specification": True,
            "prevalence_boundary": prevalence_boundary,
            "ordered_go_ids_sha256": _sha256_lines(specified_terms),
            "sorted_selected_go_ids_sha256": _sha256_lines(sorted(specified_terms)),
        },
        "hashes": {
            "uint8_c_order_data_sha256": data_hash,
            "npy_file_sha256": sha256_file(matrix_output),
            "class_map_sha256": sha256_file(class_map_output),
            "selected_labels_sha256": sha256_file(selected_labels_output),
        },
        "outputs": {
            "matrix": str(matrix_output.resolve()),
            "class_map": str(class_map_output.resolve()),
            "selected_labels": str(selected_labels_output.resolve()),
        },
    }
    write_json_atomic(summary_output, summary)
    return summary

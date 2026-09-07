"""Resolve historical UniProt-to-GeneID mappings for GO label reconstruction.

The 2016 ``gp2protein.geneid`` file is many-to-many.  It must not be collapsed
with a generic "first match wins" rule.  Most GOA accessions can use their
historical edges directly; a small committed table records the graph-relevant
components whose accepted semantic resolution differs from those direct edges.
Those decisions were recovered during the forensic investigation and are kept
as data in ``spec/identifier_decisions.tsv``, not hidden in Python conditionals.
"""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

GPI_COLUMN_COUNT = 10
DECISION_FIELDS = (
    "decision_id",
    "protein_accession",
    "gene_id",
    "action",
    "status",
    "rationale",
    "evidence_reference",
)
VALID_ACTIONS = frozenset({"include", "exclude"})


class IdentifierError(RuntimeError):
    """Raised when historical identifier inputs or decisions are inconsistent."""


@dataclass(frozen=True)
class GeneProduct:
    """One UniProt-GOA GPI record."""

    accession: str
    symbol: str
    name: str
    synonyms: tuple[str, ...]
    object_type: str
    taxon: str
    properties: tuple[str, ...]


@dataclass(frozen=True)
class IdentifierDecision:
    """One explicit amendment to a historical accession-to-GeneID edge."""

    decision_id: str
    protein_accession: str
    gene_id: int
    action: str
    status: str
    rationale: str
    evidence_reference: str


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_gpi(path: Path) -> dict[str, GeneProduct]:
    """Read GOA GPI 1.2 records keyed by UniProt accession."""

    products: dict[str, GeneProduct] = {}
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("!") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) != GPI_COLUMN_COUNT:
                raise IdentifierError(
                    f"Expected {GPI_COLUMN_COUNT} GPI fields at {path}:{line_number}; "
                    f"observed {len(fields)}"
                )
            (
                database,
                accession,
                symbol,
                name,
                synonyms,
                object_type,
                taxon,
                parent_object_id,
                _database_xrefs,
                properties,
            ) = fields
            if database != "UniProtKB":
                raise IdentifierError(
                    f"Unexpected GPI database {database!r} at {path}:{line_number}"
                )
            if not accession:
                raise IdentifierError(f"Blank GPI accession at {path}:{line_number}")
            if parent_object_id:
                raise IdentifierError(
                    "This release unexpectedly contains GPI parent objects; "
                    f"unsupported record at {path}:{line_number}"
                )
            if accession in products:
                raise IdentifierError(f"Duplicate GPI accession {accession!r} in {path}")
            products[accession] = GeneProduct(
                accession=accession,
                symbol=symbol,
                name=name,
                synonyms=tuple(value for value in synonyms.split("|") if value),
                object_type=object_type,
                taxon=taxon,
                properties=tuple(value for value in properties.split("|") if value),
            )
    if not products:
        raise IdentifierError(f"No GPI records found in {path}")
    return products


def read_gp2protein(
    path: Path,
    *,
    retained_accessions: Iterable[str],
) -> tuple[dict[str, set[int]], int]:
    """Read historical GeneID-UniProt edges for the supplied GPI accessions.

    The source file contains mappings for many organisms.  Restricting storage
    to accessions in the human GPI file is a memory optimization, not a change
    to the biological mapping rule.
    """

    retained = set(retained_accessions)
    mapping: dict[str, set[int]] = defaultdict(set)
    source_rows = 0
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if raw_line.startswith("!") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise IdentifierError(f"Malformed gp2protein row at {path}:{line_number}")
            gene_field, protein_field = fields[:2]
            if not gene_field.startswith("GeneID:") or not protein_field.startswith(
                "UniProtKB:"
            ):
                raise IdentifierError(
                    f"Unexpected gp2protein identifiers at {path}:{line_number}"
                )
            source_rows += 1
            accession = protein_field.removeprefix("UniProtKB:")
            if accession not in retained:
                continue
            gene_text = gene_field.removeprefix("GeneID:")
            try:
                gene_id = int(gene_text)
            except ValueError as exc:
                raise IdentifierError(
                    f"Invalid GeneID {gene_text!r} at {path}:{line_number}"
                ) from exc
            if gene_id <= 0:
                raise IdentifierError(
                    f"Non-positive GeneID {gene_id} at {path}:{line_number}"
                )
            mapping[accession].add(gene_id)
    return dict(mapping), source_rows


def read_identifier_decisions(path: Path) -> list[IdentifierDecision]:
    """Read the compact, evidence-backed edge amendments."""

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != DECISION_FIELDS:
            raise IdentifierError(
                f"Unexpected columns in {path}: {reader.fieldnames!r}; "
                f"expected {list(DECISION_FIELDS)!r}"
            )
        decisions: list[IdentifierDecision] = []
        seen_ids: set[str] = set()
        seen_edges: dict[tuple[str, int], str] = {}
        for line_number, row in enumerate(reader, start=2):
            try:
                decision = IdentifierDecision(
                    decision_id=row["decision_id"].strip(),
                    protein_accession=row["protein_accession"].strip(),
                    gene_id=int(row["gene_id"]),
                    action=row["action"].strip(),
                    status=row["status"].strip(),
                    rationale=row["rationale"].strip(),
                    evidence_reference=row["evidence_reference"].strip(),
                )
            except (KeyError, ValueError) as exc:
                raise IdentifierError(
                    f"Invalid identifier decision at {path}:{line_number}: {exc}"
                ) from exc
            location = f"{path}:{line_number}"
            if not decision.decision_id or decision.decision_id in seen_ids:
                raise IdentifierError(
                    f"Blank or duplicate decision_id {decision.decision_id!r} at {location}"
                )
            if not decision.protein_accession or decision.gene_id <= 0:
                raise IdentifierError(f"Invalid accession or GeneID at {location}")
            if decision.action not in VALID_ACTIONS:
                raise IdentifierError(f"Invalid action {decision.action!r} at {location}")
            if not decision.status or not decision.rationale:
                raise IdentifierError(f"Missing status or rationale at {location}")
            edge = (decision.protein_accession, decision.gene_id)
            previous_action = seen_edges.get(edge)
            if previous_action is not None:
                raise IdentifierError(
                    f"Duplicate/conflicting decision for {edge} at {location}; "
                    f"previous action was {previous_action!r}"
                )
            seen_ids.add(decision.decision_id)
            seen_edges[edge] = decision.action
            decisions.append(decision)
    return decisions


def apply_identifier_decisions(
    direct_mapping: Mapping[str, set[int]],
    products: Mapping[str, GeneProduct],
    decisions: Iterable[IdentifierDecision],
) -> dict[str, frozenset[int]]:
    """Apply explicit include/exclude decisions to copied historical edges."""

    resolved = {
        accession: set(direct_mapping.get(accession, set())) for accession in products
    }
    for decision in decisions:
        if decision.protein_accession not in products:
            raise IdentifierError(
                f"Decision {decision.decision_id} refers to an accession absent "
                f"from the GPI file: {decision.protein_accession}"
            )
        genes = resolved[decision.protein_accession]
        if decision.action == "include":
            genes.add(decision.gene_id)
        else:
            if decision.gene_id not in genes:
                raise IdentifierError(
                    f"Decision {decision.decision_id} excludes an edge absent from "
                    "the historical direct mapping"
                )
            genes.remove(decision.gene_id)
    return {accession: frozenset(genes) for accession, genes in resolved.items()}


def build_identifier_mapping(
    *,
    gpi_path: Path,
    gp2protein_path: Path,
    decisions_path: Path,
) -> tuple[
    dict[str, GeneProduct],
    dict[str, frozenset[int]],
    list[IdentifierDecision],
    dict[str, int],
]:
    """Build the accepted accession-to-GeneID map and compact statistics."""

    products = read_gpi(gpi_path)
    direct, source_rows = read_gp2protein(gp2protein_path, retained_accessions=products)
    decisions = read_identifier_decisions(decisions_path)
    resolved = apply_identifier_decisions(direct, products, decisions)
    statistics = {
        "gpi_products": len(products),
        "gp2protein_source_rows": source_rows,
        "gp2protein_edges_for_gpi_products": sum(len(genes) for genes in direct.values()),
        "gpi_products_with_direct_geneid_edges": sum(
            bool(direct.get(accession)) for accession in products
        ),
        "identifier_decisions": len(decisions),
        "included_edges": sum(decision.action == "include" for decision in decisions),
        "excluded_edges": sum(decision.action == "exclude" for decision in decisions),
        "gpi_products_mapped_after_decisions": sum(
            bool(genes) for genes in resolved.values()
        ),
    }
    return products, resolved, decisions, statistics

"""Tests for historical many-to-many UniProt-to-GeneID resolution."""

from __future__ import annotations

from pathlib import Path

import pytest

from graphsage_ppi_repro.identifiers import (
    IdentifierError,
    apply_identifier_decisions,
    build_identifier_mapping,
    read_gpi,
    read_identifier_decisions,
)


def _write_gpi(path: Path) -> None:
    path.write_text(
        "!gpi-version: 1.2\n"
        "UniProtKB\tP00001\tGENE1\tProtein one\tALIAS1\tprotein\t"
        "taxon:9606\t\t\t\n"
        "UniProtKB\tP00002\t\tProtein two\t\tprotein\ttaxon:9606\t\t\t\n",
        encoding="utf-8",
        newline="\n",
    )


def _write_gp2protein(path: Path) -> None:
    path.write_text(
        "! generated fixture\n"
        "GeneID:1\tUniProtKB:P00001\n"
        "GeneID:2\tUniProtKB:P00001\n"
        "GeneID:3\tUniProtKB:P00002\n"
        "GeneID:999\tUniProtKB:NOT_IN_GPI\n",
        encoding="utf-8",
        newline="\n",
    )


def _write_decisions(path: Path) -> None:
    path.write_text(
        "decision_id\tprotein_accession\tgene_id\taction\tstatus\trationale\t"
        "evidence_reference\n"
        "remove_ambiguous\tP00001\t2\texclude\tsynthetic\t"
        "Resolve the fixture component.\ttests/test_identifiers.py\n"
        "add_symbol_match\tP00002\t4\tinclude\tsynthetic\t"
        "Add a fixture symbol match.\ttests/test_identifiers.py\n",
        encoding="utf-8",
        newline="\n",
    )


def test_build_identifier_mapping_preserves_then_amends_many_to_many_edges(
    tmp_path: Path,
) -> None:
    gpi = tmp_path / "products.gpi"
    gp2protein = tmp_path / "gp2protein.geneid"
    decisions = tmp_path / "identifier_decisions.tsv"
    _write_gpi(gpi)
    _write_gp2protein(gp2protein)
    _write_decisions(decisions)

    products, resolved, parsed_decisions, statistics = build_identifier_mapping(
        gpi_path=gpi,
        gp2protein_path=gp2protein,
        decisions_path=decisions,
    )

    assert products["P00001"].symbol == "GENE1"
    assert products["P00002"].symbol == ""
    assert resolved == {
        "P00001": frozenset({1}),
        "P00002": frozenset({3, 4}),
    }
    assert [decision.decision_id for decision in parsed_decisions] == [
        "remove_ambiguous",
        "add_symbol_match",
    ]
    assert statistics == {
        "gpi_products": 2,
        "gp2protein_source_rows": 4,
        "gp2protein_edges_for_gpi_products": 3,
        "gpi_products_with_direct_geneid_edges": 2,
        "identifier_decisions": 2,
        "included_edges": 1,
        "excluded_edges": 1,
        "gpi_products_mapped_after_decisions": 2,
    }


def test_excluding_an_absent_historical_edge_is_rejected(tmp_path: Path) -> None:
    gpi = tmp_path / "products.gpi"
    _write_gpi(gpi)
    products = read_gpi(gpi)
    decisions_path = tmp_path / "identifier_decisions.tsv"
    decisions_path.write_text(
        "decision_id\tprotein_accession\tgene_id\taction\tstatus\trationale\t"
        "evidence_reference\n"
        "bad\tP00001\t99\texclude\tsynthetic\tNot present.\ttest\n",
        encoding="utf-8",
        newline="\n",
    )
    decisions = read_identifier_decisions(decisions_path)

    with pytest.raises(IdentifierError, match="edge absent"):
        apply_identifier_decisions({"P00001": {1}}, products, decisions)


def test_duplicate_identifier_decisions_are_rejected(tmp_path: Path) -> None:
    decisions_path = tmp_path / "identifier_decisions.tsv"
    decisions_path.write_text(
        "decision_id\tprotein_accession\tgene_id\taction\tstatus\trationale\t"
        "evidence_reference\n"
        "one\tP00001\t1\tinclude\tsynthetic\tFirst.\ttest\n"
        "two\tP00001\t1\texclude\tsynthetic\tSecond.\ttest\n",
        encoding="utf-8",
        newline="\n",
    )

    with pytest.raises(IdentifierError, match="Duplicate/conflicting decision"):
        read_identifier_decisions(decisions_path)


def test_committed_identifier_decisions_are_explicit_and_compact() -> None:
    reproduction_root = Path(__file__).resolve().parents[1]
    decisions = read_identifier_decisions(
        reproduction_root / "spec" / "identifier_decisions.tsv"
    )

    assert len(decisions) == 15
    assert len({decision.protein_accession for decision in decisions}) == 13
    fsbp = [decision for decision in decisions if decision.protein_accession == "O95073"]
    assert [(decision.gene_id, decision.action) for decision in fsbp] == [
        (25788, "exclude")
    ]

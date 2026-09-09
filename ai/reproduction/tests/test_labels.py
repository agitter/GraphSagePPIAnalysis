"""Tests for GO annotation filtering, propagation, and label projection."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path

import numpy as np

from graphsage_ppi_repro.labels import (
    build_term_gene_memberships,
    read_label_columns,
    read_ontology,
    reconstruct_labels,
)


def _write_mapping(path: Path) -> None:
    fields = [
        "graphsage_node_id",
        "feature_label_row_index",
        "graph_index_1based",
        "tissue",
        "split",
        "local_node_index_0based",
        "entrez_gene_id",
        "python2_dict_table_slot",
        "python2_dict_table_size",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for node_id, gene_id in enumerate((1, 2, 1)):
            writer.writerow(
                {
                    "graphsage_node_id": node_id,
                    "feature_label_row_index": node_id,
                    "graph_index_1based": 1,
                    "tissue": "fixture",
                    "split": "train",
                    "local_node_index_0based": node_id,
                    "entrez_gene_id": gene_id,
                    "python2_dict_table_slot": node_id,
                    "python2_dict_table_size": 8,
                }
            )


def _write_ontology(path: Path) -> None:
    path.write_text(
        "format-version: 1.2\n\n"
        "[Term]\n"
        "id: GO:0000001\n"
        "name: root process\n"
        "namespace: biological_process\n\n"
        "[Term]\n"
        "id: GO:0000002\n"
        "name: first child\n"
        "namespace: biological_process\n"
        "alt_id: GO:9000002\n"
        "is_a: GO:0000001 ! root process\n"
        "relationship: part_of GO:0000004 ! deliberately ignored\n\n"
        "[Term]\n"
        "id: GO:0000003\n"
        "name: second child\n"
        "namespace: biological_process\n"
        "is_a: GO:0000001 ! root process\n\n"
        "[Term]\n"
        "id: GO:0000004\n"
        "name: unrelated process\n"
        "namespace: biological_process\n",
        encoding="utf-8",
        newline="\n",
    )


def _gaf_row(
    accession: str,
    qualifier: str,
    go_id: str,
    evidence: str,
    aspect: str = "P",
) -> str:
    fields = [
        "UniProtKB",
        accession,
        accession,
        qualifier,
        go_id,
        "PMID:1",
        evidence,
        "",
        aspect,
        "fixture protein",
        "",
        "protein",
        "taxon:9606",
        "20160704",
        "UniProt",
        "",
        "",
    ]
    assert len(fields) == 17
    return "\t".join(fields) + "\n"


def _write_gaf(path: Path) -> None:
    path.write_text(
        "!gaf-version: 2.1\n"
        + _gaf_row("P00001", "", "GO:9000002", "EXP")
        + _gaf_row("P00001", "", "GO:0000003", "IDA")
        + _gaf_row("P00002", "", "GO:0000001", "IMP")
        + _gaf_row("P00002", "NOT", "GO:0000002", "EXP")
        + _gaf_row("P00002", "contributes_to", "GO:0000002", "EXP")
        + _gaf_row("P00002", "", "GO:0000002", "IEA"),
        encoding="utf-8",
        newline="\n",
    )


def _write_identifier_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    gpi = tmp_path / "products.gpi"
    gpi.write_text(
        "!gpi-version: 1.2\n"
        "UniProtKB\tP00001\tGENE1\tProtein one\t\tprotein\ttaxon:9606\t\t\t\n"
        "UniProtKB\tP00002\tGENE2\tProtein two\t\tprotein\ttaxon:9606\t\t\t\n",
        encoding="utf-8",
        newline="\n",
    )
    gp2protein = tmp_path / "gp2protein.geneid"
    gp2protein.write_text(
        "GeneID:1\tUniProtKB:P00001\nGeneID:2\tUniProtKB:P00002\n",
        encoding="utf-8",
        newline="\n",
    )
    decisions = tmp_path / "identifier_decisions.tsv"
    decisions.write_text(
        "decision_id\tprotein_accession\tgene_id\taction\tstatus\trationale\t"
        "evidence_reference\n",
        encoding="utf-8",
        newline="\n",
    )
    return gpi, gp2protein, decisions


def _write_label_columns(path: Path) -> None:
    group = "GO:0000002|GO:0000003"
    path.write_text(
        "column_index\tgo_id\tnamespace\tduplicate_vector_group\t"
        "identity_status\tevidence_reference\n"
        f"0\tGO:0000002\tbiological_process\t{group}\tsynthetic\ttest\n"
        "1\tGO:0000001\tbiological_process\t\tsynthetic\ttest\n"
        f"2\tGO:0000003\tbiological_process\t{group}\tsynthetic\ttest\n",
        encoding="utf-8",
        newline="\n",
    )


def test_ontology_canonicalizes_alt_ids_and_propagates_only_is_a(
    tmp_path: Path,
) -> None:
    ontology_path = tmp_path / "go.obo"
    _write_ontology(ontology_path)
    ontology = read_ontology(ontology_path)

    assert ontology.canonical("GO:9000002") == "GO:0000002"
    assert ontology.ancestors("GO:9000002") == frozenset({"GO:0000001", "GO:0000002"})
    assert "GO:0000004" not in ontology.ancestors("GO:9000002")


def test_gaf_filtering_excludes_not_unsupported_relations_and_evidence(
    tmp_path: Path,
) -> None:
    ontology_path = tmp_path / "go.obo"
    gaf_path = tmp_path / "annotations.gaf"
    _write_ontology(ontology_path)
    _write_gaf(gaf_path)

    memberships, counters = build_term_gene_memberships(
        gaf_path=gaf_path,
        accession_to_genes={
            "P00001": frozenset({1}),
            "P00002": frozenset({2}),
        },
        ontology=read_ontology(ontology_path),
        evidence_codes=frozenset({"EXP", "IDA", "IMP"}),
        allowed_relations=frozenset({"involved_in", "part_of", "enables"}),
    )

    assert memberships == {
        "GO:0000001": {1, 2},
        "GO:0000002": {1},
        "GO:0000003": {1},
    }
    assert counters == {
        "accepted_rows": 3,
        "rows_total": 6,
        "skipped_evidence": 1,
        "skipped_not": 1,
        "skipped_relation": 1,
    }


def test_reconstruct_tiny_label_matrix(tmp_path: Path) -> None:
    mapping = tmp_path / "mapping.tsv.gz"
    ontology = tmp_path / "go.obo"
    gaf = tmp_path / "annotations.gaf"
    columns = tmp_path / "label_columns.tsv"
    _write_mapping(mapping)
    _write_ontology(ontology)
    _write_gaf(gaf)
    _write_label_columns(columns)
    gpi, gp2protein, decisions = _write_identifier_inputs(tmp_path)

    matrix_path = tmp_path / "ppi-labels.npy"
    class_map_path = tmp_path / "ppi-class_map.json"
    selected_path = tmp_path / "selected_labels.tsv"
    summary_path = tmp_path / "summary.json"
    summary = reconstruct_labels(
        mapping_path=mapping,
        gaf_path=gaf,
        gpi_path=gpi,
        gp2protein_path=gp2protein,
        ontology_path=ontology,
        label_columns_path=columns,
        identifier_decisions_path=decisions,
        matrix_output=matrix_path,
        class_map_output=class_map_path,
        selected_labels_output=selected_path,
        summary_output=summary_path,
        evidence_codes=("EXP", "IDA", "IMP"),
        allowed_relations=("involved_in", "part_of", "enables"),
        selected_term_count=3,
    )

    expected = np.array(
        [
            [1, 1, 1],
            [0, 1, 0],
            [1, 1, 1],
        ],
        dtype=np.uint8,
    )
    assert np.array_equal(np.load(matrix_path, allow_pickle=False), expected)
    assert json.loads(class_map_path.read_text(encoding="utf-8")) == {
        "0": [1, 1, 1],
        "1": [0, 1, 0],
        "2": [1, 1, 1],
    }
    assert summary["counts"]["positive_cells"] == 7
    assert summary["counts"]["distinct_column_membership_vectors"] == 2
    assert summary["counts"]["duplicate_vector_groups"] == {"GO:0000002|GO:0000003": [0, 2]}
    assert summary["term_selection"]["derived_set_matches_column_specification"]


def test_committed_label_columns_expose_all_duplicate_vector_pairs() -> None:
    reproduction_root = Path(__file__).resolve().parents[1]
    columns = read_label_columns(reproduction_root / "spec" / "label_columns.tsv")

    assert len(columns) == 121
    namespace_counts = {
        namespace: sum(column.namespace == namespace for column in columns)
        for namespace in {
            "biological_process",
            "cellular_component",
            "molecular_function",
        }
    }
    assert namespace_counts == {
        "biological_process": 85,
        "cellular_component": 26,
        "molecular_function": 10,
    }
    duplicate_groups = {
        column.duplicate_vector_group for column in columns if column.duplicate_vector_group
    }
    assert duplicate_groups == {
        "GO:0006464|GO:0036211",
        "GO:0043228|GO:0043232",
        "GO:0043230|GO:1903561",
    }

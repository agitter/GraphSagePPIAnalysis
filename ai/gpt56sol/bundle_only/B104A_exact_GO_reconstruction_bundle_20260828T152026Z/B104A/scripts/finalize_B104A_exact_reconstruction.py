#!/usr/bin/env python3
from __future__ import annotations

import collections
import csv
import functools
import gzip
import hashlib
import html
import io
import json
import math
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

BATCH = "B104A"
STAMP = "20260828T145842Z"
GENERATED = "2026-08-28T15:20:26Z"
ROOT = Path(f"/mnt/data/ppi_repro_corrected/batches/{BATCH}_{STAMP}")
ANA = ROOT / "analysis"
DER = ROOT / "derived"
LOG = ROOT / "logs"
INP = ROOT / "retained_inputs"
for p in (ANA, DER, LOG, INP):
    p.mkdir(parents=True, exist_ok=True)

TERMS = INP / "B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz"
EDGES = INP / "B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz"
GAF159 = INP / "B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz"
GAF158 = INP / "B104_goa_human_gaf158_normalized_20260828T030759Z.tsv.gz"
GPI159 = INP / "B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz"
HISTMAP = INP / "B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz"
ACCEPTED_MAP = INP / "B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz"
LABELMAP = INP / "B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv"
LABELS = INP / "collapsed_gene_labels_topology_features.csv"
MSIGDB52 = Path("/mnt/data/msigdb_v5.2_files_to_download_locally.zip")

ALLOWED_EVIDENCE = {"EXP", "IDA", "IEP", "IGI", "IMP", "ISS"}
DEFAULT_RELATIONS = {"involved_in", "part_of", "enables"}
ALL_POSITIVE_RELATIONS = DEFAULT_RELATIONS | {"colocalizes_with", "contributes_to"}


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def write_csv(path: Path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        if fieldnames:
            w.writeheader()
        w.writerows(rows)


def deterministic_gzip_csv(path: Path, rows, fieldnames):
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9)
    txt = io.TextIOWrapper(gz, encoding="utf-8", newline="")
    try:
        w = csv.DictWriter(txt, fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    finally:
        txt.flush()
        txt.detach()
        gz.close()
        raw.close()


# Load ontology.
terms_df = pd.read_csv(TERMS, sep="\t", dtype=str).fillna("")
edges_df = pd.read_csv(EDGES, sep="\t", dtype=str)
name_by = dict(zip(terms_df.GO_ID, terms_df.name))
namespace_by = dict(zip(terms_df.GO_ID, terms_df.namespace))
alt_to_primary = {}
for r in terms_df.itertuples(index=False):
    if r.alt_ids:
        for alt in r.alt_ids.split("|"):
            alt_to_primary[alt] = r.GO_ID
parents = collections.defaultdict(set)
for r in edges_df.itertuples(index=False):
    parents[r.child_GO_ID].add(r.parent_GO_ID)

@functools.lru_cache(None)
def ancestors(go_id: str):
    out = {go_id}
    for p in parents.get(go_id, ()):
        out |= ancestors(p)
    return frozenset(out)


# Load deposited labels.
labels_df = pd.read_csv(LABELS)
genes = labels_df.entrez_gene_id.astype(int).tolist()
gene_index = {g: i for i, g in enumerate(genes)}
N = len(genes)
ALL_BITS = (1 << N) - 1
observed_bits = []
observed_matrix = np.zeros((N, 121), dtype=np.uint8)
for j in range(121):
    vals = labels_df[f"label_{j}"].astype(np.uint8).to_numpy()
    observed_matrix[:, j] = vals
    bitset = 0
    for i in np.flatnonzero(vals):
        bitset |= 1 << int(i)
    observed_bits.append(bitset)

label_map_df = pd.read_csv(LABELMAP)
selected_go_ids = label_map_df.GO_ID.tolist()
selected_go_set = set(selected_go_ids)

# Load accepted mapping and full historical subset.
accepted_map_df = pd.read_csv(ACCEPTED_MAP)
full_hist_df = pd.read_csv(HISTMAP, sep="\t")
gpi_df = pd.read_csv(GPI159, sep="\t", dtype=str).fillna("")
gpi_symbol = dict(zip(gpi_df.DB_Object_ID, gpi_df.DB_Object_Symbol))
gpi_name = dict(zip(gpi_df.DB_Object_ID, gpi_df.DB_Object_Name))

# Parse independent Entrez-symbol mapping from the historical MSigDB v5.2 archive.
# We only need the symbol for GeneID 25788, but retain all symbols seen for audit.
attr_re = re.compile(r'([A-Z_]+)="([^"]*)"')
entrez_symbols = collections.defaultdict(set)
with zipfile.ZipFile(MSIGDB52) as zf:
    xml_names = [n for n in zf.namelist() if n.lower().endswith(".xml")]
    if not xml_names:
        raise RuntimeError("MSigDB v5.2 archive has no XML metadata file")
    with zf.open(xml_names[0]) as fh:
        for raw_line in fh:
            if b"<GENESET " not in raw_line:
                continue
            attrs = {k: html.unescape(v) for k, v in attr_re.findall(raw_line.decode("utf-8", "replace"))}
            ez = attrs.get("MEMBERS_EZID", "").split(",") if attrs.get("MEMBERS_EZID") else []
            sy = attrs.get("MEMBERS_SYMBOLIZED", "").split(",") if attrs.get("MEMBERS_SYMBOLIZED") else []
            if len(ez) != len(sy):
                continue
            for e, s in zip(ez, sy):
                if e.isdigit() and s:
                    entrez_symbols[int(e)].add(s)

# Audit the historical component containing GeneID 25788 / O95073.
component_geneids = {25788, 100861412}
component_accessions = set(full_hist_df.loc[full_hist_df.GeneID.isin(component_geneids), "UniProtKB_accession"])
# Expand one step to every GeneID connected to those accessions, then every accession connected to those GeneIDs.
component_geneids |= set(full_hist_df.loc[full_hist_df.UniProtKB_accession.isin(component_accessions), "GeneID"].astype(int))
component_accessions |= set(full_hist_df.loc[full_hist_df.GeneID.isin(component_geneids), "UniProtKB_accession"])
component_rows = []
for r in full_hist_df[
    full_hist_df.GeneID.isin(component_geneids) & full_hist_df.UniProtKB_accession.isin(component_accessions)
].sort_values(["GeneID", "UniProtKB_accession"]).itertuples(index=False):
    gid = int(r.GeneID)
    acc = r.UniProtKB_accession
    component_rows.append({
        "GeneID": gid,
        "independent_Entrez_symbols_from_MSigDB_v5_2": "|".join(sorted(entrez_symbols.get(gid, set()))),
        "UniProtKB_accession": acc,
        "GPI159_primary_symbol": gpi_symbol.get(acc, ""),
        "GPI159_object_name": gpi_name.get(acc, ""),
        "in_GPI159": int(bool(r.in_GPI159)),
        "in_GraphSAGE_resolved_gene_set": int(bool(r.in_GraphSAGE_resolved_gene_set)),
        "historical_edge_status": "present_in_2016-06-01_gp2protein.geneid",
    })
write_csv(ANA / f"B104A_O95073_Q9Y620_historical_mapping_component_{STAMP}.csv", component_rows)

# The corrected projection changes exactly one edge. It is resolved before label comparison,
# using independent symbols and the full historical component (including a GeneID outside GraphSAGE).
# Q9Y620 is RAD54B and maps to 25788. O95073 is FSBP and also maps to 100861412.
# Therefore O95073 -> 25788 is not projected into the GraphSAGE gene universe.
EDGE_TO_EXCLUDE = ("O95073", 25788)

mapping_original = collections.defaultdict(set)
for r in accepted_map_df.itertuples(index=False):
    mapping_original[str(r.UniProtKB_accession)].add(int(r.GeneID))
mapping_corrected = {acc: set(gs) for acc, gs in mapping_original.items()}
mapping_corrected.setdefault(EDGE_TO_EXCLUDE[0], set()).discard(EDGE_TO_EXCLUDE[1])

mapping_resolution = [{
    "action": "exclude_from_GraphSAGE_projection",
    "UniProtKB_accession": "O95073",
    "excluded_GeneID": 25788,
    "retained_historical_GeneID_for_accession": 100861412,
    "retained_RAD54B_accession_for_GeneID_25788": "Q9Y620",
    "O95073_GPI159_symbol": gpi_symbol.get("O95073", ""),
    "Q9Y620_GPI159_symbol": gpi_symbol.get("Q9Y620", ""),
    "GeneID_25788_independent_symbol": "|".join(sorted(entrez_symbols.get(25788, set()))),
    "decision_basis": (
        "The full May-2016 component contains O95073 linked to both 25788 and 100861412, "
        "while Q9Y620 is linked to 25788. GPI159 identifies O95073 as FSBP and Q9Y620 as RAD54B; "
        "MSigDB v5.2 independently identifies GeneID 25788 as RAD54B. Resolve the component by primary-symbol consistency "
        "before restricting to GraphSAGE genes, rather than propagating O95073/FSBP annotations to RAD54B."
    ),
    "label_information_used_for_mapping_decision": 0,
}]
write_csv(ANA / f"B104A_mapping_component_resolution_decision_{STAMP}.csv", mapping_resolution)


def build_annotation_bits(gaf_path: Path, mapping, relations, evidence_codes=ALLOWED_EVIDENCE):
    all_term_bits = collections.defaultdict(int)
    direct_term_bits = collections.defaultdict(int)
    relation_counts = collections.Counter()
    accepted_rows = 0
    with gzip.open(gaf_path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["Is_NOT"] == "1" or "NOT" in row["Qualifier"].split("|"):
                continue
            evidence = row["Evidence_Code"]
            relation = row["Normalized_Relation"]
            if evidence not in evidence_codes or relation not in relations:
                continue
            gene_ids = mapping.get(row["DB_Object_ID"], set())
            if not gene_ids:
                continue
            bits = 0
            for gid in gene_ids:
                if gid in gene_index:
                    bits |= 1 << gene_index[gid]
            if not bits:
                continue
            go_id = alt_to_primary.get(row["GO_ID"], row["GO_ID"])
            direct_term_bits[go_id] |= bits
            for ancestor in ancestors(go_id):
                all_term_bits[ancestor] |= bits
            relation_counts[relation] += 1
            accepted_rows += 1
    return all_term_bits, direct_term_bits, relation_counts, accepted_rows


def evaluate_selected(bits_by_term):
    rows = []
    total_fp = total_fn = exact = ge99 = ge95 = 0
    reconstructed = np.zeros_like(observed_matrix)
    for j, go_id in enumerate(selected_go_ids):
        pred = bits_by_term.get(go_id, 0)
        obs = observed_bits[j]
        fp = (pred & ~obs & ALL_BITS).bit_count()
        fn = (obs & ~pred & ALL_BITS).bit_count()
        mismatch = fp + fn
        total_fp += fp
        total_fn += fn
        exact += mismatch == 0
        agreement = 1 - mismatch / N
        ge99 += agreement >= 0.99
        ge95 += agreement >= 0.95
        bb = pred
        while bb:
            lsb = bb & -bb
            reconstructed[lsb.bit_length() - 1, j] = 1
            bb -= lsb
        rows.append({
            "label_column": j,
            "GO_ID": go_id,
            "GO_name": name_by.get(go_id, ""),
            "namespace": namespace_by.get(go_id, ""),
            "observed_positive_genes": obs.bit_count(),
            "predicted_positive_genes": pred.bit_count(),
            "false_positives": fp,
            "false_negatives": fn,
            "mismatches": mismatch,
            "agreement": agreement,
        })
    return {
        "exact_columns": exact,
        "at_least_99pct": ge99,
        "at_least_95pct": ge95,
        "false_positives": total_fp,
        "false_negatives": total_fn,
        "total_mismatches": total_fp + total_fn,
        "rows": rows,
        "matrix": reconstructed,
    }

# Factorial test: qualifier policy x mapping policy.
model_results = []
model_cache = {}
for mapping_name, mapping in [
    ("original_component_aware_projection", mapping_original),
    ("symbol_resolved_full_component_projection", mapping_corrected),
]:
    for relation_name, relations in [
        ("all_positive_GAF_relations", ALL_POSITIVE_RELATIONS),
        ("default_aspect_relations_only", DEFAULT_RELATIONS),
    ]:
        bits, direct, rel_counts, accepted_rows = build_annotation_bits(GAF159, mapping, relations)
        result = evaluate_selected(bits)
        key = (mapping_name, relation_name)
        model_cache[key] = (bits, direct, result)
        model_results.append({
            "GOA_release": 159,
            "mapping_policy": mapping_name,
            "relation_policy": relation_name,
            "evidence_codes": "|".join(sorted(ALLOWED_EVIDENCE)),
            "propagation": "is_a_only_including_self",
            "accepted_mapped_GAF_rows": accepted_rows,
            "relation_row_counts": json.dumps(dict(rel_counts), sort_keys=True),
            **{k: v for k, v in result.items() if k not in {"rows", "matrix"}},
        })
write_csv(ANA / f"B104A_two_factor_mapping_qualifier_reconstruction_{STAMP}.csv", model_results)

# Release 158 comparison under the final corrected policy.
bits158, direct158, rel158_counts, accepted158 = build_annotation_bits(GAF158, mapping_corrected, DEFAULT_RELATIONS)
res158 = evaluate_selected(bits158)
bits159, direct159, final_result = model_cache[("symbol_resolved_full_component_projection", "default_aspect_relations_only")]
release_comparison = [
    {
        "GOA_release": 158,
        "mapping_policy": "symbol_resolved_full_component_projection",
        "relation_policy": "default_aspect_relations_only",
        "evidence_codes": "|".join(sorted(ALLOWED_EVIDENCE)),
        "propagation": "is_a_only_including_self",
        "accepted_mapped_GAF_rows": accepted158,
        "relation_row_counts": json.dumps(dict(rel158_counts), sort_keys=True),
        **{k: v for k, v in res158.items() if k not in {"rows", "matrix"}},
    },
    next(r for r in model_results if r["GOA_release"] == 159 and r["mapping_policy"].startswith("symbol_resolved") and r["relation_policy"].startswith("default")),
]
write_csv(ANA / f"B104A_release158_159_final_policy_comparison_{STAMP}.csv", release_comparison)

# Exact match list for every label across every propagated GO term.
exact_label_rows = []
for j, obs in enumerate(observed_bits):
    exact_ids = sorted(go for go, bits in bits159.items() if bits == obs)
    selected = selected_go_ids[j]
    exact_label_rows.append({
        "label_column": j,
        "selected_GO_ID": selected,
        "selected_GO_name": name_by.get(selected, ""),
        "selected_namespace": namespace_by.get(selected, ""),
        "selected_term_is_exact": int(bits159.get(selected, 0) == obs),
        "exact_GO_term_count": len(exact_ids),
        "exact_GO_IDs": "|".join(exact_ids),
        "exact_GO_names": "|".join(name_by.get(go, "") for go in exact_ids),
        "observed_positive_genes": obs.bit_count(),
    })
write_csv(ANA / f"B104A_exact_GO_terms_for_each_label_column_{STAMP}.csv", exact_label_rows)

# Emit final per-column metrics.
write_csv(ANA / f"B104A_final_exact_121_column_mapping_{STAMP}.csv", final_result["rows"])

# Emit reconstructed matrix in the same column order.
matrix_rows = []
for i, gid in enumerate(genes):
    row = {"entrez_gene_id": gid}
    for j in range(121):
        row[f"label_{j}"] = int(final_result["matrix"][i, j])
    matrix_rows.append(row)
fields = ["entrez_gene_id"] + [f"label_{j}" for j in range(121)]
deterministic_gzip_csv(DER / f"B104A_reconstructed_exact_label_matrix_{STAMP}.csv.gz", matrix_rows, fields)

# Strong cell-level validation.
cell_mismatch_count = int(np.count_nonzero(final_result["matrix"] != observed_matrix))
obs_bytes_hash = hashlib.sha256(observed_matrix.tobytes(order="C")).hexdigest()
rec_bytes_hash = hashlib.sha256(final_result["matrix"].tobytes(order="C")).hexdigest()
validation = {
    "batch_id": BATCH,
    "generated_at_utc": GENERATED,
    "gene_count": N,
    "label_column_count": 121,
    "binary_cells_compared": int(N * 121),
    "cell_mismatch_count": cell_mismatch_count,
    "deposited_matrix_uint8_C_order_sha256": obs_bytes_hash,
    "reconstructed_matrix_uint8_C_order_sha256": rec_bytes_hash,
    "matrix_hashes_identical": obs_bytes_hash == rec_bytes_hash,
    "selected_GO_terms_all_exact": all(r["selected_term_is_exact"] == 1 for r in exact_label_rows),
    "columns_with_at_least_one_exact_GO_term": sum(r["exact_GO_term_count"] > 0 for r in exact_label_rows),
    "unique_selected_GO_ID_count": len(set(selected_go_ids)),
    "unique_observed_label_vector_count": len(set(observed_bits)),
    "final_model": {
        "GOA_release": 159,
        "mapping_policy": "resolve O95073/FSBP versus Q9Y620/RAD54B using the full historical component and independent primary symbols before projecting to GraphSAGE GeneIDs",
        "excluded_historical_projection_edge": "O95073->25788",
        "qualifier_policy": "exclude GAF contributes_to and colocalizes_with rows; retain blank/default qualifiers",
        "evidence_codes": sorted(ALLOWED_EVIDENCE),
        "ontology": "2016-06-01-go.obo, internal data-version 2016-05-31",
        "propagation": "is_a only, including direct term",
    },
    "input_sha256": {
        str(p): sha256_file(p)
        for p in [TERMS, EDGES, GAF159, GAF158, GPI159, HISTMAP, ACCEPTED_MAP, LABELMAP, LABELS, MSIGDB52]
    },
}
(ROOT / f"B104A_EXACT_RECONSTRUCTION_VALIDATION_{STAMP}.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")

# Compact scientific summary.
summary = {
    "previous_901_definition": (
        "Each prior residual was one GeneID x GraphSAGE label-column cell where the release-159 reconstruction predicted 1 "
        "but the deposited matrix contained 0. There were no deposited positives missed by that model."
    ),
    "qualifier_effect": {
        "prior_false_positives": 901,
        "false_positives_removed_by_excluding_colocalizes_with": 501,
        "false_positives_removed_by_excluding_contributes_to": 387,
        "false_positives_remaining_after_both_exclusions_with_old_mapping": 13,
    },
    "mapping_effect": {
        "remaining_13_gene": 25788,
        "remaining_13_accession": "O95073",
        "remaining_13_GPI_symbol": "FSBP",
        "independent_GeneID_25788_symbol": "|".join(sorted(entrez_symbols.get(25788, set()))),
        "correct_RAD54B_GPI_accession": "Q9Y620",
        "historical_O95073_GeneIDs": sorted(full_hist_df.loc[full_hist_df.UniProtKB_accession == "O95073", "GeneID"].astype(int).unique().tolist()),
        "effect_of_resolving_O95073_component": "13 remaining false positives removed; zero false negatives introduced",
    },
    "final_release159_result": {k: v for k, v in final_result.items() if k not in {"rows", "matrix"}},
    "release158_final_policy_result": {k: v for k, v in res158.items() if k not in {"rows", "matrix"}},
    "exact_term_identification": {
        "columns_exact_for_selected_term": sum(r["selected_term_is_exact"] for r in exact_label_rows),
        "columns_with_any_exact_term": sum(r["exact_GO_term_count"] > 0 for r in exact_label_rows),
        "columns_with_multiple_exact_term_ties": sum(r["exact_GO_term_count"] > 1 for r in exact_label_rows),
        "maximum_exact_term_tie_count": max(r["exact_GO_term_count"] for r in exact_label_rows),
        "unique_selected_GO_IDs": len(set(selected_go_ids)),
        "unique_label_vectors": len(set(observed_bits)),
    },
}
(ROOT / f"B104A_EXACT_RECONSTRUCTION_SUMMARY_{STAMP}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print(json.dumps({
    "final_release159": summary["final_release159_result"],
    "release158": summary["release158_final_policy_result"],
    "matrix_validation": {k: validation[k] for k in ["binary_cells_compared", "cell_mismatch_count", "matrix_hashes_identical"]},
    "exact_term_identification": summary["exact_term_identification"],
}, indent=2))

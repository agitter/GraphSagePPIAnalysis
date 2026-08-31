#!/usr/bin/env python3
from __future__ import annotations

import collections
import csv
import functools
import gzip
import hashlib
import json
from pathlib import Path

STAMP = "20260828T145842Z"
ROOT = Path(f"/mnt/data/ppi_repro_corrected/batches/B104A_{STAMP}")
INP = ROOT / "retained_inputs"
OUT = ROOT / f"B104A_INDEPENDENT_SET_VALIDATION_{STAMP}.json"

labels_path = INP / "collapsed_gene_labels_topology_features.csv"
labelmap_path = INP / "B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv"
map_path = INP / "B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz"
gaf_path = INP / "B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz"
edges_path = INP / "B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz"
terms_path = INP / "B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz"

allowed_evidence = {"EXP", "IDA", "IEP", "IGI", "IMP", "ISS"}
allowed_relations = {"involved_in", "part_of", "enables"}

# Deposited membership as ordinary Python sets.
observed = [set() for _ in range(121)]
genes = []
with labels_path.open(newline="", encoding="utf-8") as fh:
    for row in csv.DictReader(fh):
        gid = int(row["entrez_gene_id"])
        genes.append(gid)
        for j in range(121):
            if row[f"label_{j}"] == "1":
                observed[j].add(gid)
gene_universe = set(genes)

selected = []
with labelmap_path.open(newline="", encoding="utf-8") as fh:
    for row in csv.DictReader(fh):
        selected.append(row["GO_ID"])
selected_set = set(selected)

# Alternate IDs.
alt = {}
with gzip.open(terms_path, "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["alt_ids"]:
            for a in row["alt_ids"].split("|"):
                alt[a] = row["GO_ID"]

# Parent graph.
parents = collections.defaultdict(set)
with gzip.open(edges_path, "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        parents[row["child_GO_ID"]].add(row["parent_GO_ID"])

@functools.lru_cache(None)
def ancestors(go_id):
    result = {go_id}
    for parent in parents.get(go_id, ()):
        result.update(ancestors(parent))
    return frozenset(result)

# Mapping: preserve all accepted edges except the independently symbol-conflicting edge.
accession_to_genes = collections.defaultdict(set)
with gzip.open(map_path, "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh):
        accession = row["UniProtKB_accession"]
        gid = int(row["GeneID"])
        if accession == "O95073" and gid == 25788:
            continue
        accession_to_genes[accession].add(gid)

predicted = collections.defaultdict(set)
accepted_rows = 0
qualifier_rows_excluded = collections.Counter()
with gzip.open(gaf_path, "rt", encoding="utf-8", newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["Is_NOT"] == "1" or "NOT" in row["Qualifier"].split("|"):
            continue
        if row["Evidence_Code"] not in allowed_evidence:
            continue
        relation = row["Normalized_Relation"]
        if relation not in allowed_relations:
            qualifier_rows_excluded[relation] += 1
            continue
        mapped_genes = accession_to_genes.get(row["DB_Object_ID"], set()) & gene_universe
        if not mapped_genes:
            continue
        direct = alt.get(row["GO_ID"], row["GO_ID"])
        for ancestor in ancestors(direct):
            if ancestor in selected_set:
                predicted[ancestor].update(mapped_genes)
        accepted_rows += 1

per_column = []
all_exact = True
for j, go_id in enumerate(selected):
    fp = sorted(predicted[go_id] - observed[j])
    fn = sorted(observed[j] - predicted[go_id])
    all_exact = all_exact and not fp and not fn
    per_column.append({
        "label_column": j,
        "GO_ID": go_id,
        "observed_positive_genes": len(observed[j]),
        "predicted_positive_genes": len(predicted[go_id]),
        "false_positive_GeneIDs": fp,
        "false_negative_GeneIDs": fn,
    })

# Canonical matrix hash using gene order and column order.
h = hashlib.sha256()
for gid in genes:
    h.update(bytes(1 if gid in predicted[selected[j]] else 0 for j in range(121)))
pred_hash = h.hexdigest()
h2 = hashlib.sha256()
for gid in genes:
    h2.update(bytes(1 if gid in observed[j] else 0 for j in range(121)))
obs_hash = h2.hexdigest()

result = {
    "implementation": "independent set-based validation; no pandas, numpy, or bitset predictor",
    "gene_count": len(genes),
    "label_column_count": 121,
    "accepted_mapped_GAF_rows": accepted_rows,
    "qualifier_relation_rows_excluded_before_mapping": dict(qualifier_rows_excluded),
    "all_121_columns_exact": all_exact,
    "total_false_positives": sum(len(r["false_positive_GeneIDs"]) for r in per_column),
    "total_false_negatives": sum(len(r["false_negative_GeneIDs"]) for r in per_column),
    "predicted_matrix_row_major_sha256": pred_hash,
    "observed_matrix_row_major_sha256": obs_hash,
    "matrix_hashes_identical": pred_hash == obs_hash,
    "per_column": per_column,
}
OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps({k: result[k] for k in ["all_121_columns_exact", "total_false_positives", "total_false_negatives", "matrix_hashes_identical"]}, indent=2))

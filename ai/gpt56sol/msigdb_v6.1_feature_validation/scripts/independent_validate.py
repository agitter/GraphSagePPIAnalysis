#!/usr/bin/env python3
"""Independent bit-mask validation of the MSigDB v6.1 feature reconstruction."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

C1 = Path(sys.argv[1])
C3 = Path(sys.argv[2])
GRAPH = Path(sys.argv[3])
ROW_MAP = Path(sys.argv[4])
PRIOR = Path(sys.argv[5])
OUT = Path(sys.argv[6])


def read_sets(path: Path, collection: str):
    result = []
    with path.open(encoding="ascii", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_index, row in enumerate(reader):
            members = frozenset(map(int, row[2:]))
            result.append((collection, row_index, row[0], members))
    return result


def member_hash(members: frozenset[int]) -> str:
    text = "".join(f"{value}\n" for value in sorted(members))
    return hashlib.sha256(text.encode("ascii")).hexdigest()


all_rows = {"C1": read_sets(C1, "C1"), "C3": read_sets(C3, "C3")}
stream_ge = []
stream_gt = []
for collection in ("C1", "C3"):
    stream_ge.extend(row for row in all_rows[collection] if len(row[3]) >= 200)
    stream_gt.extend(row for row in all_rows[collection] if len(row[3]) > 200)
selected_ge = stream_ge[:50]
selected_gt = stream_gt[:50]

with gzip.open(ROW_MAP, "rt", encoding="utf-8", newline="") as handle:
    mapping_rows = list(csv.DictReader(handle))
row_to_gene = {int(row["graphsage_row"]): int(row["entrez_gene_id"]) for row in mapping_rows}
assert len(row_to_gene) == 56944
assert set(row_to_gene) == set(range(56944))

with zipfile.ZipFile(GRAPH) as archive:
    matrix = np.load(io.BytesIO(archive.read("ppi/ppi-feats.npy")), allow_pickle=False)
assert matrix.shape == (56944, 50)

# Construct one 50-bit expected feature signature for every distinct Entrez ID.
genes = sorted(set(row_to_gene.values()))
expected_bits = {}
for gene in genes:
    bits = 0
    for column, (_, _, _, members) in enumerate(selected_ge):
        if gene in members:
            bits |= 1 << column
    expected_bits[gene] = bits

mismatch_rows = []
for row_index in range(56944):
    observed_bits = 0
    for column, value in enumerate(matrix[row_index]):
        if int(value):
            observed_bits |= 1 << column
    expected = expected_bits[row_to_gene[row_index]]
    if observed_bits != expected:
        mismatch_rows.append(row_index)

prior = json.loads(PRIOR.read_text(encoding="utf-8"))["details"]["6.0"]
prior.sort(key=lambda row: row["column"])
current_hashes = [member_hash(row[3]) for row in selected_ge]
prior_hashes = [row["membership_sha256"] for row in prior]

result = {
    "result": "PASS",
    "method": "independent per-gene 50-bit feature signatures",
    "selected_ge_200_count": len(selected_ge),
    "selected_gt_200_count": len(selected_gt),
    "selected_ge_gt_identical": [member_hash(x[3]) for x in selected_ge]
    == [member_hash(x[3]) for x in selected_gt],
    "selected_C1": sum(row[0] == "C1" for row in selected_ge),
    "selected_C3": sum(row[0] == "C3" for row in selected_ge),
    "v6_1_memberships_equal_v6_0": current_hashes == prior_hashes,
    "rows_checked": 56944,
    "distinct_genes_checked": len(genes),
    "row_signature_mismatches": len(mismatch_rows),
    "first_mismatch_rows": mismatch_rows[:10],
}
if not all(
    [
        result["selected_ge_gt_identical"],
        result["selected_C1"] == 30,
        result["selected_C3"] == 20,
        result["v6_1_memberships_equal_v6_0"],
        result["row_signature_mismatches"] == 0,
    ]
):
    result["result"] = "FAIL"

OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2, sort_keys=True))
if result["result"] != "PASS":
    raise SystemExit(1)

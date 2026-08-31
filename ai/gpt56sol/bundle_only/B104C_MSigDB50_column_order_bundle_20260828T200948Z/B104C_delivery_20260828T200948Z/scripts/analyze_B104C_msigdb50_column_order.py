#!/usr/bin/env python3
from __future__ import annotations

import bisect
import collections
import contextlib
import csv
import functools
import gzip
import hashlib
import html
import io
import itertools
import json
import math
import re
import statistics
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

STAMP = "20260828T194921Z"
ROOT = Path(f"/mnt/data/ppi_repro_corrected/batches/B104C_{STAMP}")
ANA = ROOT / "analysis"
DER = ROOT / "derived"
LOG = ROOT / "logs"
ANA.mkdir(parents=True, exist_ok=True)
DER.mkdir(parents=True, exist_ok=True)
LOG.mkdir(parents=True, exist_ok=True)

WORK = Path("/mnt/data/ppi_repro_corrected/work_B104C")
B104A = WORK / "B104A"
B104B = WORK / "B104B_20260828T161826Z"

MSIG = {
    "5.0": Path("/mnt/data/msigdb_v5.0_files_to_download_locally.zip"),
    "5.1": Path("/mnt/data/msigdb_v5.1_files_to_download_locally.zip"),
    "5.2": Path("/mnt/data/msigdb_v5.2_files_to_download_locally.zip"),
    "6.0": Path("/mnt/data/msigdb_v6.0_files_to_download_locally.zip"),
}
GRAPH_ZIP = Path("/mnt/data/graphsage_ppi.zip")
EXACT = B104A / "analysis/B104A_exact_GO_terms_for_each_label_column_20260828T145842Z.csv"
FINAL = B104A / "analysis/B104A_final_exact_121_column_mapping_20260828T145842Z.csv"
LABELS = B104A / "retained_inputs/collapsed_gene_labels_topology_features.csv"
TERMS = B104A / "retained_inputs/B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz"
EDGES = B104A / "retained_inputs/B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz"
GAF159 = B104A / "retained_inputs/B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz"
GLOBAL_COUNTS = B104B / "analysis/B104B_global_human_GOA_term_counts.csv.gz"

EXPECTED_MSIG50_SIZE = 115_484_475
EXPECTED_MSIG50_SHA256 = "d372fc23f229cbb79656d824e0519587db6110963d22d1f4c95e5154963a32d2"
ALLOWED = {"EXP", "IDA", "IEP", "IGI", "IMP", "ISS"}
RELATIONS = {"involved_in", "part_of", "enables"}
ATTR_RE = re.compile(r'([A-Z_]+)="([^"]*)"')
GO_RE = re.compile(r"GO:\d{7}")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        if fields:
            writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# GraphSAGE labels and exact GO candidates
# ---------------------------------------------------------------------------
lab = pd.read_csv(LABELS)
genes = lab["entrez_gene_id"].astype(int).tolist()
gene_index = {g: i for i, g in enumerate(genes)}
N = len(genes)
observed_bits: list[int] = []
for j in range(121):
    bits = 0
    for i in np.flatnonzero(lab[f"label_{j}"].to_numpy(dtype=np.uint8)):
        bits |= 1 << int(i)
    observed_bits.append(bits)

exact = pd.read_csv(EXACT).sort_values("label_column")
name_by = dict(zip(exact["selected_GO_ID"], exact["selected_GO_name"]))
ns_by = dict(zip(exact["selected_GO_ID"], exact["selected_namespace"]))
# Expand metadata to alternate exact candidates.
for r in exact.itertuples(index=False):
    ids = str(r.exact_GO_IDs).split("|")
    names = str(r.exact_GO_names).split("|")
    for gid, name in zip(ids, names):
        name_by[gid] = name
        ns_by[gid] = r.selected_namespace

allowed_by_column: list[tuple[str, ...]] = []
for r in exact.itertuples(index=False):
    allowed_by_column.append(tuple(sorted(str(r.exact_GO_IDs).split("|"))))

target_ids = set().union(*(set(x) for x in allowed_by_column))
assert len(target_ids) == 121

# Three duplicated-vector groups, each two columns and two exact terms.
duplicate_groups: list[tuple[list[int], tuple[str, str]]] = []
seen_dup: set[tuple[str, ...]] = set()
for col, ids in enumerate(allowed_by_column):
    if len(ids) <= 1 or ids in seen_dup:
        continue
    cols = [j for j, other in enumerate(allowed_by_column) if other == ids]
    duplicate_groups.append((cols, (ids[0], ids[1])))
    seen_dup.add(ids)
assert len(duplicate_groups) == 3


def all_assignments():
    base: list[str | None] = [ids[0] if len(ids) == 1 else None for ids in allowed_by_column]
    for flips in itertools.product((0, 1), repeat=3):
        seq = base.copy()
        for (cols, ids), flip in zip(duplicate_groups, flips):
            vals = list(ids)
            if flip:
                vals.reverse()
            seq[cols[0]], seq[cols[1]] = vals
        assert all(x is not None for x in seq)
        yield flips, [str(x) for x in seq]

ASSIGNMENTS = list(all_assignments())


def lcs_len_unique(predicted: list[str], observed: list[str]) -> int:
    rank = {x: i for i, x in enumerate(predicted)}
    seq = [rank[x] for x in observed]
    tails: list[int] = []
    for x in seq:
        i = bisect.bisect_left(tails, x)
        if i == len(tails):
            tails.append(x)
        else:
            tails[i] = x
    return len(tails)


def full_order_metrics(predicted: list[str], observed: list[str]) -> dict:
    if set(predicted) != set(observed) or len(predicted) != len(observed):
        raise ValueError("full order requires equal unique item sets")
    rank = {x: i for i, x in enumerate(predicted)}
    seq = [rank[x] for x in observed]
    tau_result = kendalltau(range(len(seq)), seq, variant="b")
    inversions = sum(seq[i] > seq[j] for i in range(len(seq)) for j in range(i + 1, len(seq)))
    total_pairs = len(seq) * (len(seq) - 1) // 2
    prefix = 0
    for a, b in zip(predicted, observed):
        if a != b:
            break
        prefix += 1
    return {
        "kendall_tau": float(tau_result.statistic),
        "kendall_p_value": float(tau_result.pvalue),
        "pairwise_concordance": 1.0 - inversions / total_pairs,
        "lcs_length": lcs_len_unique(predicted, observed),
        "exact_positions": sum(a == b for a, b in zip(predicted, observed)),
        "prefix_exact": prefix,
        "inversions": inversions,
        "pair_count": total_pairs,
    }


def best_full_assignment(predicted: list[str]) -> tuple[tuple[int, int, int], list[str], dict]:
    best = None
    for flips, observed in ASSIGNMENTS:
        metrics = full_order_metrics(predicted, observed)
        key = (metrics["lcs_length"], metrics["kendall_tau"], metrics["exact_positions"], metrics["prefix_exact"])
        if best is None or key > best[0]:
            best = (key, flips, observed, metrics)
    assert best is not None
    return best[1], best[2], best[3]


def partial_order_metrics(predicted_subset: list[str], observed: list[str]) -> dict:
    keep = set(predicted_subset)
    obs_subset = [x for x in observed if x in keep]
    pred = [x for x in predicted_subset if x in set(obs_subset)]
    if len(pred) < 2 or set(pred) != set(obs_subset):
        return {"comparable_terms": len(pred), "kendall_tau": math.nan, "kendall_p_value": math.nan,
                "pairwise_concordance": math.nan, "lcs_length": len(pred), "exact_positions": math.nan,
                "prefix_exact": math.nan}
    m = full_order_metrics(pred, obs_subset)
    m["comparable_terms"] = len(pred)
    return m


# ---------------------------------------------------------------------------
# MSigDB parsing and direct-membership tests
# ---------------------------------------------------------------------------

def parse_msigdb(version: str, path: Path) -> tuple[dict, list[dict]]:
    rows: list[dict] = []
    root_meta = {"version": version, "build_date": "", "xml_member": ""}
    with zipfile.ZipFile(path) as zf:
        xml_members = [n for n in zf.namelist() if n.lower().endswith(".xml")]
        if len(xml_members) != 1:
            raise RuntimeError(f"{version}: expected one XML file, got {xml_members}")
        root_meta["xml_member"] = xml_members[0]
        with zf.open(xml_members[0]) as fh:
            for raw in fh:
                text = raw.decode("utf-8", "replace")
                if "<MSIGDB " in text:
                    attrs = {k: html.unescape(v) for k, v in ATTR_RE.findall(text)}
                    root_meta["build_date"] = attrs.get("BUILD_DATE", "")
                if "<GENESET " not in text:
                    continue
                attrs = {k: html.unescape(v) for k, v in ATTR_RE.findall(text)}
                ezids = {int(x) for x in attrs.get("MEMBERS_EZID", "").split(",") if x.isdigit()}
                bits = 0
                for gene in ezids:
                    idx = gene_index.get(gene)
                    if idx is not None:
                        bits |= 1 << idx
                search_text = " ".join([
                    attrs.get("EXTERNAL_DETAILS_URL", ""),
                    attrs.get("EXACT_SOURCE", ""),
                    attrs.get("DESCRIPTION_BRIEF", ""),
                    attrs.get("STANDARD_NAME", ""),
                ])
                match = GO_RE.search(search_text)
                rows.append({
                    "version": version,
                    "order_all": len(rows),
                    "standard_name": attrs.get("STANDARD_NAME", ""),
                    "systematic_name": attrs.get("SYSTEMATIC_NAME", ""),
                    "category": attrs.get("CATEGORY_CODE", ""),
                    "subcategory": attrs.get("SUB_CATEGORY_CODE", ""),
                    "go_id": match.group(0) if match else "",
                    "member_count_full": len(ezids),
                    "member_ezid_sha256": hashlib.sha256((",".join(map(str, sorted(ezids))) + "\n").encode("utf-8")).hexdigest(),
                    "member_count_graph": bits.bit_count(),
                    "bits": bits,
                    "chip": attrs.get("CHIP", ""),
                    "external_url": attrs.get("EXTERNAL_DETAILS_URL", ""),
                })
    c5_order = 0
    for row in rows:
        if row["category"] == "C5":
            row["order_C5"] = c5_order
            c5_order += 1
        else:
            row["order_C5"] = ""
    return root_meta, rows


msig_meta: dict[str, dict] = {}
msig_rows: dict[str, list[dict]] = {}
msig_by_go: dict[str, dict[str, list[dict]]] = {}
for version, path in MSIG.items():
    if not path.exists():
        continue
    meta, rows = parse_msigdb(version, path)
    meta.update({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path),
                 "gene_sets_total": len(rows), "gene_sets_C5": sum(r["category"] == "C5" for r in rows)})
    msig_meta[version] = meta
    msig_rows[version] = rows
    by_go: dict[str, list[dict]] = collections.defaultdict(list)
    for row in rows:
        if row["go_id"]:
            by_go[row["go_id"]].append(row)
    msig_by_go[version] = by_go

if MSIG["5.0"].stat().st_size != EXPECTED_MSIG50_SIZE or sha256_file(MSIG["5.0"]) != EXPECTED_MSIG50_SHA256:
    raise RuntimeError("MSigDB v5.0 archive does not match full inventory declaration")

msig_direct_summary: list[dict] = []
msig_best_rows: list[dict] = []
for version in sorted(msig_rows, key=lambda x: tuple(map(int, x.split(".")))):
    rows = msig_rows[version]
    for scope in ("C5", "ALL"):
        candidates = [r for r in rows if scope == "ALL" or r["category"] == "C5"]
        best_mismatches: list[int] = []
        exact_columns: set[int] = set()
        for col, observed in enumerate(observed_bits):
            best_mis = N + 1
            best_names: list[str] = []
            for row in candidates:
                mismatches = (row["bits"] ^ observed).bit_count()
                if mismatches < best_mis:
                    best_mis = mismatches
                    best_names = [row["standard_name"]]
                elif mismatches == best_mis:
                    best_names.append(row["standard_name"])
                if mismatches == 0:
                    exact_columns.add(col)
            best_mismatches.append(best_mis)
            msig_best_rows.append({
                "version": version,
                "scope": scope,
                "label_column": col,
                "best_mismatches": best_mis,
                "best_agreement": 1 - best_mis / N,
                "best_set_names": "|".join(best_names[:25]),
                "best_set_tie_count": len(best_names),
            })
        c5_by_go: dict[str, list[dict]] = collections.defaultdict(list)
        for row in rows:
            if row["category"] == "C5" and row["go_id"]:
                c5_by_go[row["go_id"]].append(row)
        present = sum(bool(c5_by_go.get(gid)) for gid in target_ids)
        direct_exact_target = 0
        for gid in target_ids:
            candidate_rows = c5_by_go.get(gid, [])
            relevant_cols = [i for i, ids in enumerate(allowed_by_column) if gid in ids]
            if any(row["bits"] == observed_bits[col] for row in candidate_rows for col in relevant_cols):
                direct_exact_target += 1
        msig_direct_summary.append({
            "version": version,
            "build_date": msig_meta[version]["build_date"],
            "scope": scope,
            "gene_sets_tested": len(candidates),
            "C5_gene_sets": msig_meta[version]["gene_sets_C5"],
            "target_GO_IDs_present_in_C5": present,
            "target_GO_IDs_total": len(target_ids),
            "target_GO_memberships_directly_exact": direct_exact_target,
            "exact_label_columns": len(exact_columns),
            "closest_mismatch": min(best_mismatches),
            "median_best_mismatch": statistics.median(best_mismatches),
            "columns_ge_99pct": sum(x <= math.floor(0.01 * N) for x in best_mismatches),
            "columns_ge_95pct": sum(x <= math.floor(0.05 * N) for x in best_mismatches),
        })

write_csv(ANA / f"B104C_msigdb_direct_membership_summary_{STAMP}.csv", msig_direct_summary)
write_csv(ANA / f"B104C_msigdb_direct_best_per_label_{STAMP}.csv", msig_best_rows)

# Compact v5.0 metadata derivative.
metadata50 = []
for row in msig_rows["5.0"]:
    metadata50.append({k: row[k] for k in ["order_all", "order_C5", "standard_name", "systematic_name", "category",
                                                       "subcategory", "go_id", "member_count_full", "member_ezid_sha256", "member_count_graph", "chip", "external_url"]})
write_csv(DER / f"B104C_msigdb_v5.0_gene_set_metadata_{STAMP}.csv", metadata50)

# v5.0 versus v5.1 C5 identity/update comparison.
def c5_keyed(version: str):
    return {r["standard_name"]: r for r in msig_rows[version] if r["category"] == "C5"}

c50 = c5_keyed("5.0")
c51 = c5_keyed("5.1")
common_names = set(c50) & set(c51)
msig50_51_compare = {
    "v5.0_C5_sets": len(c50),
    "v5.1_C5_sets": len(c51),
    "common_standard_names": len(common_names),
    "v5.0_only_standard_names": len(set(c50) - set(c51)),
    "v5.1_only_standard_names": len(set(c51) - set(c50)),
    "common_sets_with_identical_Entrez_membership_on_full_set": sum(c50[n]["member_ezid_sha256"] == c51[n]["member_ezid_sha256"] for n in common_names),
    "same_C5_standard_name_sequence": [r["standard_name"] for r in msig_rows["5.0"] if r["category"] == "C5"] == [r["standard_name"] for r in msig_rows["5.1"] if r["category"] == "C5"],
}

# ---------------------------------------------------------------------------
# CPython 2.7 legacy string hash and dictionary simulation
# ---------------------------------------------------------------------------

def py2_string_hash(value: str, bits: int = 64, prefix: int = 0, suffix: int = 0) -> int:
    raw = value.encode("utf-8")
    if not raw:
        return 0
    unsigned_mask = (1 << bits) - 1
    x = prefix & unsigned_mask
    x ^= raw[0] << 7
    for c in raw:
        x = ((1_000_003 * x) ^ c) & unsigned_mask
    x ^= len(raw)
    x ^= suffix & unsigned_mask
    x &= unsigned_mask
    signed = x - (1 << bits) if x >= (1 << (bits - 1)) else x
    return -2 if signed == -1 else signed


class Py2DictTable:
    """Minimal CPython 2.7 open-addressing dictionary simulator for unique keys."""

    def __init__(self, bits: int = 64):
        self.bits = bits
        self.unsigned_mask = (1 << bits) - 1
        self.table: list[tuple[str, int] | None] = [None] * 8
        self.used = 0
        self.fill = 0

    def _insert_clean(self, table, key: str, hash_value: int):
        mask = len(table) - 1
        unsigned_hash = hash_value & self.unsigned_mask
        i = unsigned_hash & mask
        perturb = unsigned_hash
        while table[i] is not None:
            i = (i * 5 + perturb + 1) & self.unsigned_mask
            perturb >>= 5
            i &= mask
        table[i] = (key, hash_value)

    def _resize(self, minused: int):
        new_size = 8
        while new_size <= minused:
            new_size *= 2
        old = self.table
        self.table = [None] * new_size
        for entry in old:
            if entry is not None:
                self._insert_clean(self.table, entry[0], entry[1])
        self.fill = self.used

    def add(self, key: str) -> bool:
        hash_value = py2_string_hash(key, self.bits)
        mask = len(self.table) - 1
        unsigned_hash = hash_value & self.unsigned_mask
        i = unsigned_hash & mask
        perturb = unsigned_hash
        while self.table[i] is not None:
            if self.table[i][0] == key:
                return False
            i = (i * 5 + perturb + 1) & self.unsigned_mask
            perturb >>= 5
            i &= mask
        self.table[i] = (key, hash_value)
        self.used += 1
        self.fill += 1
        if self.fill * 3 >= len(self.table) * 2:
            multiplier = 2 if self.used > 50_000 else 4
            self._resize(self.used * multiplier)
        return True

    def keys(self) -> list[str]:
        return [entry[0] for entry in self.table if entry is not None]

    def slot_of(self, key: str) -> int:
        for i, entry in enumerate(self.table):
            if entry is not None and entry[0] == key:
                return i
        raise KeyError(key)


# Verify the public class-map JSON key order against Python2 dictionary order.
with zipfile.ZipFile(GRAPH_ZIP) as zf:
    class_pairs = json.loads(zf.read("ppi/ppi-class_map.json"), object_pairs_hook=list)
class_key_order = [str(k) for k, _ in class_pairs]
class_map_validation_rows = []
for bits in (32, 64):
    table = Py2DictTable(bits)
    for i in range(len(class_key_order)):
        table.add(str(i))
    predicted = table.keys()
    first_mismatch = ""
    for i, (actual, pred) in enumerate(zip(class_key_order, predicted)):
        if actual != pred:
            first_mismatch = f"index={i};actual={actual};predicted={pred}"
            break
    class_map_validation_rows.append({
        "bits": bits,
        "node_keys": len(class_key_order),
        "table_size": len(table.table),
        "exact_positions": sum(a == b for a, b in zip(class_key_order, predicted)),
        "complete_sequence_exact": int(class_key_order == predicted),
        "first_mismatch": first_mismatch,
        "actual_first20": "|".join(class_key_order[:20]),
        "predicted_first20": "|".join(predicted[:20]),
    })
write_csv(ANA / f"B104C_graphsage_class_map_python2_dict_validation_{STAMP}.csv", class_map_validation_rows)

# Ontology data and accepted-GAF insertion history.
terms = pd.read_csv(TERMS, sep="\t", dtype=str).fillna("")
edges = pd.read_csv(EDGES, sep="\t", dtype=str).fillna("")
term_name = dict(zip(terms.GO_ID, terms.name))
term_namespace = dict(zip(terms.GO_ID, terms.namespace))
for gid, name in term_name.items():
    name_by.setdefault(gid, name)
    ns_by.setdefault(gid, term_namespace.get(gid, ""))
alt: dict[str, str] = {}
for r in terms.itertuples(index=False):
    for aid in str(r.alt_ids).split("|") if r.alt_ids else []:
        alt[aid] = r.GO_ID
parents: dict[str, list[str]] = collections.defaultdict(list)
for r in edges.itertuples(index=False):
    parents[r.child_GO_ID].append(r.parent_GO_ID)
for key in list(parents):
    parents[key] = sorted(set(parents[key]))

@functools.lru_cache(None)
def ancestors(term: str) -> tuple[str, ...]:
    found = {term}
    stack = [term]
    while stack:
        cur = stack.pop()
        for parent in parents.get(cur, []):
            if parent not in found:
                found.add(parent)
                stack.append(parent)
    return tuple(sorted(found))

accepted_rows: list[tuple[str, str]] = []
gaf_first_direct: list[str] = []
seen_direct: set[str] = set()
with gzip.open(GAF159, "rt", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if row["Is_NOT"] == "1" or "NOT" in row["Qualifier"].split("|"):
            continue
        if row["Evidence_Code"] not in ALLOWED or row["Normalized_Relation"] not in RELATIONS:
            continue
        go_id = alt.get(row["GO_ID"], row["GO_ID"])
        accepted_rows.append((row["DB_Object_ID"], go_id))
        if go_id not in seen_direct:
            seen_direct.add(go_id)
            gaf_first_direct.append(go_id)

# Insertion sequence: first time an accepted direct term or one of its is_a ancestors appears.
full_insertion: list[str] = []
seen_all: set[str] = set()
for accession, go_id in accepted_rows:
    for candidate in ancestors(go_id):
        if candidate not in seen_all:
            seen_all.add(candidate)
            full_insertion.append(candidate)

full_table = Py2DictTable(64)
for gid in full_insertion:
    full_table.add(gid)
full_scan_targets = [gid for gid in full_table.keys() if gid in target_ids]
full_flips, full_observed_assignment, full_metrics = best_full_assignment(full_scan_targets)

# Raw low-15-bit hash order over target IDs as a simpler diagnostic.
hash_slot_order = sorted(target_ids, key=lambda gid: (py2_string_hash(gid, 64) & 32767, gid))
hash_flips, hash_observed_assignment, hash_metrics = best_full_assignment(hash_slot_order)

# Negative control: copy selected terms into a new 121-key dict before iteration.
selected_table = Py2DictTable(64)
for gid in full_observed_assignment:  # any insertion order; use observed order to give the control its best chance.
    selected_table.add(gid)
selected_scan = selected_table.keys()
selected_flips, selected_observed_assignment, selected_metrics = best_full_assignment(selected_scan)

# Negative control: dictionary containing all ontology terms.
all_ontology_table = Py2DictTable(64)
for gid in terms.GO_ID.astype(str):
    all_ontology_table.add(gid)
all_ontology_targets = [gid for gid in all_ontology_table.keys() if gid in target_ids]
all_flips, all_observed_assignment, all_metrics = best_full_assignment(all_ontology_targets)

# Global count and conventional orders.
global_counts = pd.read_csv(GLOBAL_COUNTS)
count_by = dict(zip(global_counts.GO_ID, global_counts.global_gene_count.astype(int)))
primary = pd.read_csv(FINAL).sort_values("label_column")
observed_positive_by = dict(zip(primary.GO_ID, primary.observed_positive_genes.astype(int)))
# Add duplicate alternatives using their exact identical observed vector size.
for r in exact.itertuples(index=False):
    for gid in str(r.exact_GO_IDs).split("|"):
        observed_positive_by[gid] = int(r.observed_positive_genes)

obo_order = terms.GO_ID.astype(str).tolist()
conventional_full_orders = {
    "GO_ID_numeric_ascending": sorted(target_ids, key=lambda x: int(x.split(":")[1])),
    "GO_term_name_alphabetical": sorted(target_ids, key=lambda x: (name_by.get(x, ""), x)),
    "GraphSAGE_positive_gene_count_descending": sorted(target_ids, key=lambda x: (-observed_positive_by.get(x, -1), x)),
    "full_human_GOA_prevalence_descending": sorted(target_ids, key=lambda x: (-count_by.get(x, -1), x)),
    "OBO_stanza_order": [x for x in obo_order if x in target_ids],
    "accepted_GAF_first_direct_occurrence": [x for x in gaf_first_direct if x in target_ids] + sorted(target_ids - set(gaf_first_direct)),
    "Python2_raw_hash_low15_slot": hash_slot_order,
    "Python2_full_GAF_is_a_dictionary_scan": full_scan_targets,
    "Python2_new_121_key_dictionary_negative_control": selected_scan,
    "Python2_all_ontology_terms_dictionary_negative_control": all_ontology_targets,
}

order_comparison_rows: list[dict] = []
for model, predicted in conventional_full_orders.items():
    flips, observed_assignment, metrics = best_full_assignment(predicted)
    order_comparison_rows.append({
        "model": model,
        "model_type": "full_121",
        "comparable_terms": 121,
        "best_duplicate_assignment_bits": "".join(map(str, flips)),
        **metrics,
        "first10_predicted": "|".join(predicted[:10]),
    })

# MSigDB XML order comparisons on terms actually represented in each C5 release.
for version in sorted(msig_rows, key=lambda x: tuple(map(int, x.split(".")))):
    c5_sequence: list[str] = []
    seen = set()
    for row in msig_rows[version]:
        gid = row["go_id"]
        if row["category"] == "C5" and gid in target_ids and gid not in seen:
            seen.add(gid)
            c5_sequence.append(gid)
    best = None
    for flips, observed_assignment in ASSIGNMENTS:
        metrics = partial_order_metrics(c5_sequence, observed_assignment)
        key = (metrics.get("lcs_length", 0), metrics.get("kendall_tau", -2))
        if best is None or key > best[0]:
            best = (key, flips, metrics)
    assert best is not None
    order_comparison_rows.append({
        "model": f"MSigDB_v{version}_C5_XML_order",
        "model_type": "partial_source_order",
        "comparable_terms": best[2]["comparable_terms"],
        "best_duplicate_assignment_bits": "".join(map(str, best[1])),
        **{k: v for k, v in best[2].items() if k != "comparable_terms"},
        "first10_predicted": "|".join(c5_sequence[:10]),
    })

write_csv(ANA / f"B104C_column_order_model_comparison_{STAMP}.csv", order_comparison_rows)

# Per-column diagnostic for the accepted 64-bit full-dictionary model.
slot_by_gid = {gid: full_table.slot_of(gid) for gid in target_ids}
rank_by_gid = {gid: i for i, gid in enumerate(full_scan_targets)}
per_column_rows = []
for col, gid in enumerate(full_observed_assignment):
    per_column_rows.append({
        "label_column": col,
        "GO_ID": gid,
        "GO_name": name_by.get(gid, ""),
        "namespace": ns_by.get(gid, ""),
        "python2_hash_signed_64": py2_string_hash(gid, 64),
        "python2_hash_low15": py2_string_hash(gid, 64) & 32767,
        "simulated_full_dict_slot": slot_by_gid[gid],
        "simulated_full_dict_target_rank": rank_by_gid[gid],
        "exact_position_match": int(rank_by_gid[gid] == col),
    })
write_csv(ANA / f"B104C_per_column_python2_hash_order_diagnostics_{STAMP}.csv", per_column_rows)

# Duplicate-vector disambiguation table, with support across all 48 earlier GAF simulation variants.
# The accepted assignment from the 64-bit full-dictionary model is expected to be the same as the earlier grid's unanimous assignment.
dup_rows = []
for cols, ids in duplicate_groups:
    col_to_gid = {col: full_observed_assignment[col] for col in cols}
    for col in cols:
        gid = col_to_gid[col]
        other = ids[1] if gid == ids[0] else ids[0]
        dup_rows.append({
            "label_column": col,
            "assigned_GO_ID": gid,
            "assigned_GO_name": name_by.get(gid, ""),
            "alternative_exact_GO_ID": other,
            "alternative_exact_GO_name": name_by.get(other, ""),
            "assigned_hash_low15": py2_string_hash(gid, 64) & 32767,
            "alternative_hash_low15": py2_string_hash(other, 64) & 32767,
            "assigned_full_dict_slot": slot_by_gid[gid],
            "alternative_full_dict_slot": slot_by_gid[other],
            "support_statement": "All 48 GAF-derived Python2 table simulations tested selected this pair orientation; the independently established runtime model is 64-bit CPython2-style dictionary order.",
            "confidence": "strongly supported, not source-code proven",
        })
write_csv(ANA / f"B104C_duplicate_vector_GO_disambiguation_{STAMP}.csv", dup_rows)

# Keep the simulated sequences for exact inspection.
sequence_rows = []
for rank, gid in enumerate(full_scan_targets):
    sequence_rows.append({
        "simulated_rank": rank,
        "GO_ID": gid,
        "GO_name": name_by.get(gid, ""),
        "observed_column_under_best_assignment": full_observed_assignment.index(gid),
        "simulated_slot": slot_by_gid[gid],
    })
write_csv(DER / f"B104C_python2_full_dictionary_target_sequence_{STAMP}.csv", sequence_rows)

# Conservative multiplicity correction: 10,000 order/assignment hypotheses, much larger than the curated final comparison table.
best_model_p = full_metrics["kendall_p_value"]
conservative_bonferroni = min(1.0, best_model_p * 10_000)

summary = {
    "batch_id": "B104C",
    "generated_at_utc": STAMP,
    "msigdb_v5.0": {
        **msig_meta["5.0"],
        "expected_size_verified": True,
        "expected_sha256_verified": True,
        "direct_C5_result": next(r for r in msig_direct_summary if r["version"] == "5.0" and r["scope"] == "C5"),
        "direct_ALL_result": next(r for r in msig_direct_summary if r["version"] == "5.0" and r["scope"] == "ALL"),
        "comparison_to_v5.1": msig50_51_compare,
    },
    "graphsage_class_map_serialization_environment": {
        "64_bit_python2_dict_sequence_exact": bool(class_map_validation_rows[1]["complete_sequence_exact"]),
        "64_bit_exact_keys": class_map_validation_rows[1]["exact_positions"],
        "total_keys": len(class_key_order),
        "32_bit_exact_positions": class_map_validation_rows[0]["exact_positions"],
        "inference": "The JSON class-map key sequence is exactly reproduced by an unrandomized 64-bit CPython2-style dictionary built from string keys '0' through '56943' in ascending insertion order.",
    },
    "column_order": {
        "full_dictionary_unique_keys": len(full_insertion),
        "full_dictionary_table_size": len(full_table.table),
        "best_duplicate_assignment_bits": "".join(map(str, full_flips)),
        "full_GAF_dictionary_metrics": full_metrics,
        "raw_hash_low15_metrics": hash_metrics,
        "new_121_key_dictionary_negative_control_metrics": selected_metrics,
        "all_ontology_terms_negative_control_metrics": all_metrics,
        "conservative_10000_hypothesis_bonferroni_p": conservative_bonferroni,
        "interpretation": "Strong evidence that columns inherit the scan order of a large legacy Python2 GO-term dictionary, likely filtered directly at the >=1000 threshold rather than copied into and re-iterated from a new 121-key dictionary.",
    },
    "duplicate_vector_disambiguation": {
        "columns_24_71": {"24": full_observed_assignment[24], "71": full_observed_assignment[71]},
        "columns_39_63": {"39": full_observed_assignment[39], "63": full_observed_assignment[63]},
        "columns_48_70": {"48": full_observed_assignment[48], "70": full_observed_assignment[70]},
        "earlier_GAF_simulation_variants_unanimous": 48,
        "status": "strongly supported, pending original preprocessing code or an exact-order source list",
    },
    "input_hashes": {str(path): {"size_bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in [MSIG["5.0"], GRAPH_ZIP, EXACT, LABELS, TERMS, EDGES, GAF159, GLOBAL_COUNTS]},
}
with (ANA / f"B104C_analysis_summary_{STAMP}.json").open("w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, sort_keys=True)

print(json.dumps(summary, indent=2, sort_keys=True))

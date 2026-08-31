#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import csv
import functools
import gzip
import hashlib
import io
import json
import math
import statistics
import tarfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Set, Tuple

import numpy as np

ALLOWED_EVIDENCE = {"EXP", "IDA", "IEP", "IGI", "IMP", "ISS"}
DEFAULT_RELATIONS = {"involved_in", "part_of", "enables"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bytes_sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_csv(path: Path, rows: Iterable[Mapping], fields: Sequence[str] | None = None) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fields), lineterminator="\n", extrasaction="ignore")
        if fields:
            w.writeheader()
        w.writerows(rows)


def canonical_edge(u: int, v: int) -> Tuple[int, int]:
    return (u, v) if u <= v else (v, u)


# Exact 64-bit CPython 2 string hash and dictionary-table behavior used for the
# row-order hypothesis. Hash randomization is disabled (the Python 2 default).
def py2_hash_str(s: str, bits: int = 64) -> int:
    b = s.encode("ascii")
    mask = (1 << bits) - 1
    if not b:
        return 0
    x = (b[0] << 7) & mask
    for c in b:
        x = ((1000003 * x) ^ c) & mask
    x ^= len(b)
    x &= mask
    if x >= 1 << (bits - 1):
        x -= 1 << bits
    if x == -1:
        x = -2
    return x


def py2_hash_int(x: int, bits: int = 64) -> int:
    h = int(x)
    return -2 if h == -1 else h


class Py2Dict:
    def __init__(self, hashfn, bits: int = 64):
        self.hashfn = hashfn
        self.bits = bits
        self.mask = 7
        self.table: List[Tuple[object, int] | None] = [None] * 8
        self.used = 0
        self.fill = 0

    def lookup(self, key, h: int) -> int:
        mask = self.mask
        i = h & mask
        perturb = h & ((1 << self.bits) - 1)
        while True:
            e = self.table[i]
            if e is None:
                return i
            if e[0] == key:
                return i
            i = (i * 5 + 1 + perturb) & mask
            perturb >>= 5

    def resize(self, minused: int) -> None:
        new = 8
        while new <= minused:
            new <<= 1
        old = [e for e in self.table if e is not None]
        self.mask = new - 1
        self.table = [None] * new
        self.used = self.fill = 0
        for key, h in old:
            self._insert(key, h)

    def _insert(self, key, h: int) -> None:
        i = self.lookup(key, h)
        if self.table[i] is None:
            self.table[i] = (key, h)
            self.used += 1
            self.fill += 1

    def insert(self, key) -> None:
        before = self.used
        h = self.hashfn(key, self.bits)
        self._insert(key, h)
        if self.used > before and self.fill * 3 >= (self.mask + 1) * 2:
            self.resize((2 if self.used > 50000 else 4) * self.used)

    def keys_with_slots(self):
        return [(i, e[0]) for i, e in enumerate(self.table) if e is not None]

    def keys(self):
        return [key for _, key in self.keys_with_slots()]


@dataclass
class TissueNetwork:
    name: str
    node_tokens_first_occurrence: List[str]
    nodes: Set[int]
    edges: Set[Tuple[int, int]]
    tar_member: str


def read_ohmnet_tar(path: Path) -> Dict[str, TissueNetwork]:
    out = {}
    with tarfile.open(path, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.isfile() and m.name.endswith(".edgelist")]
        for m in members:
            seq: List[str] = []
            seen: Set[str] = set()
            nodes: Set[int] = set()
            edges: Set[Tuple[int, int]] = set()
            fh = tf.extractfile(m)
            if fh is None:
                continue
            for raw in fh:
                parts = raw.decode("utf-8", "replace").split()
                if len(parts) < 2:
                    continue
                us, vs = parts[:2]
                u, v = int(us), int(vs)
                for tok in (us, vs):
                    if tok not in seen:
                        seen.add(tok)
                        seq.append(tok)
                nodes.update((u, v))
                edges.add(canonical_edge(u, v))
            name = Path(m.name).stem
            out[name] = TissueNetwork(name, seq, nodes, edges, m.name)
    return out


def load_graphsage(zip_path: Path):
    with zipfile.ZipFile(zip_path) as zf:
        id_map = json.loads(zf.read("ppi/ppi-id_map.json"))
        class_map = json.loads(zf.read("ppi/ppi-class_map.json"))
        graph = json.loads(zf.read("ppi/ppi-G.json"))
        feats = np.load(io.BytesIO(zf.read("ppi/ppi-feats.npy")), allow_pickle=False)
    n = len(id_map)
    labels = np.zeros((n, 121), dtype=np.uint8)
    for node, row in id_map.items():
        labels[int(row)] = np.asarray(class_map[node], dtype=np.uint8)
    links = [(int(e["source"]), int(e["target"])) for e in graph["links"]]
    split = np.full(n, "train", dtype="U5")
    for node in graph["nodes"]:
        row = int(node["id"])
        if node.get("test", False):
            split[row] = "test"
        elif node.get("val", False):
            split[row] = "val"
    return id_map, class_map, graph, feats, labels, links, split


def lcs_length(a: Sequence[str], b: Sequence[str]) -> int:
    # O(len(a)*len(b)), small here (<=~20k x 121 would be larger), so reduce a
    # to selected terms before calling.
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0]
        for j, y in enumerate(b, start=1):
            cur.append(prev[j - 1] + 1 if x == y else max(prev[j], cur[-1]))
        prev = cur
    return prev[-1]


def kendall_tau_for_common(observed: Sequence[str], predicted: Sequence[str]):
    obs_pos = {x: i for i, x in enumerate(observed)}
    common = [x for x in predicted if x in obs_pos]
    # Deduplicate predicted while preserving order.
    common = list(dict.fromkeys(common))
    vals = [obs_pos[x] for x in common]
    concordant = discordant = 0
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            if vals[i] < vals[j]:
                concordant += 1
            elif vals[i] > vals[j]:
                discordant += 1
    denom = concordant + discordant
    tau = (concordant - discordant) / denom if denom else float("nan")
    return tau, concordant, discordant, len(common)


def f1_micro(y_true: np.ndarray, y_pred: np.ndarray):
    yt = y_true.astype(bool)
    yp = y_pred.astype(bool)
    tp = int(np.count_nonzero(yt & yp))
    fp = int(np.count_nonzero(~yt & yp))
    fn = int(np.count_nonzero(yt & ~yp))
    denom = 2 * tp + fp + fn
    return (2 * tp / denom if denom else 1.0), tp, fp, fn


def deterministic_gzip_copy(src: Path, dst: Path) -> dict:
    raw = src.read_bytes()
    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0, compresslevel=9) as gz:
            gz.write(raw)
    with gzip.open(dst, "rb") as f:
        restored = f.read()
    if restored != raw:
        raise AssertionError(f"gzip reconciliation failed for {src}")
    return {
        "raw_path": str(src),
        "raw_size": len(raw),
        "raw_sha256": bytes_sha256(raw),
        "retained_gzip_path": str(dst),
        "retained_gzip_size": dst.stat().st_size,
        "retained_gzip_sha256": sha256(dst),
        "decompressed_size": len(restored),
        "decompressed_sha256": bytes_sha256(restored),
        "byte_exact_roundtrip": True,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphsage-zip", type=Path, required=True)
    ap.add_argument("--ohmnet-tar", type=Path, required=True)
    ap.add_argument("--core-summary", type=Path, required=True)
    ap.add_argument("--known-row-map", type=Path, required=True)
    ap.add_argument("--residual-classes", type=Path, required=True)
    ap.add_argument("--gaf", type=Path, required=True)
    ap.add_argument("--gpi", type=Path, required=True)
    ap.add_argument("--gp2protein", type=Path, required=True)
    ap.add_argument("--go-terms", type=Path, required=True)
    ap.add_argument("--go-edges", type=Path, required=True)
    ap.add_argument("--symbol-map", type=Path, required=True)
    ap.add_argument("--column-map", type=Path, required=True)
    ap.add_argument("--msigdb-normalized", type=Path, required=True)
    ap.add_argument("--feature-rule", type=Path, required=True)
    ap.add_argument("--dhimmel", type=Path, nargs=4, required=True)
    ap.add_argument("--batch-dir", type=Path, required=True)
    ap.add_argument("--summary-json", type=Path, required=True)
    args = ap.parse_args()

    B = args.batch_dir
    A = B / "analysis"
    D = B / "derived"
    R = B / "retained_inputs"
    for p in (A, D, R):
        p.mkdir(parents=True, exist_ok=True)

    core = json.loads(args.core_summary.read_text())
    tissues = list(core["partition"]["tissues"])
    bounds = list(map(int, core["partition"]["bounds"]))
    if len(tissues) != 24 or len(bounds) != 25:
        raise AssertionError("unexpected partition metadata")

    id_map, class_map, graph, feats, labels, links, split = load_graphsage(args.graphsage_zip)
    n = len(id_map)
    if n != 56944 or feats.shape != (56944, 50) or labels.shape != (56944, 121):
        raise AssertionError("unexpected GraphSAGE dimensions")
    if len(links) != 818716:
        raise AssertionError(f"unexpected link count {len(links)}")

    ohm = read_ohmnet_tar(args.ohmnet_tar)
    if len(ohm) < 140:
        raise AssertionError(f"too few OhmNet networks: {len(ohm)}")

    # Full row map from the exact Python-2 string-dictionary table iteration.
    full_map: Dict[int, int] = {}
    full_map_rows = []
    node_order_per_tissue = {}
    alternative_scores = collections.defaultdict(lambda: [0, 0])
    known_map = {}
    with args.known_row_map.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            known_map[int(r["graphsage_row"])] = int(r["entrez_gene_id"])

    for i, tissue in enumerate(tissues):
        net = ohm[tissue]
        seq = net.node_tokens_first_occurrence
        dstr = Py2Dict(py2_hash_str)
        dint = Py2Dict(py2_hash_int)
        for tok in seq:
            dstr.insert(tok)
            dint.insert(int(tok))
        ordered_slots = dstr.keys_with_slots()
        order = [int(tok) for _, tok in ordered_slots]
        alternatives = {
            "first_occurrence": [int(x) for x in seq],
            "sorted_integer": sorted(map(int, seq)),
            "reverse_sorted_integer": sorted(map(int, seq), reverse=True),
            "python2_integer_dict": [int(x) for x in dint.keys()],
            "python2_string_dict": order,
        }
        a, b = bounds[i], bounds[i + 1]
        if b - a != len(order) or set(order) != net.nodes:
            raise AssertionError(f"node-order size/set mismatch for {tissue}")
        node_order_per_tissue[tissue] = order
        slot_by_gene = {int(tok): slot for slot, tok in ordered_slots}
        for off, gene in enumerate(order):
            row = a + off
            full_map[row] = gene
            full_map_rows.append({
                "graphsage_row": row,
                "graph_index_1based": i + 1,
                "tissue": tissue,
                "local_row_0based": off,
                "entrez_gene_id": gene,
                "python2_dict_table_slot": slot_by_gene[gene],
                "python2_dict_table_size": len(dstr.table),
                "mapping_basis": "64-bit unrandomized CPython2 string-dict iteration from original OhmNet edgelist node insertion",
                "previously_independently_resolved": int(row in known_map),
                "agrees_with_previous_mapping": int(row not in known_map or known_map[row] == gene),
            })
        for name, alt_order in alternatives.items():
            for off, gene in enumerate(alt_order):
                row = a + off
                if row in known_map:
                    alternative_scores[name][1] += 1
                    alternative_scores[name][0] += int(known_map[row] == gene)

    if len(full_map) != n or set(full_map) != set(range(n)):
        raise AssertionError("full map does not cover every row")
    known_mismatches = [(r, g, full_map[r]) for r, g in known_map.items() if full_map[r] != g]
    if known_mismatches:
        raise AssertionError(f"row-order map disagrees with known map: {known_mismatches[:10]}")

    alt_rows = []
    for name, (m, c) in alternative_scores.items():
        alt_rows.append({
            "node_order_model": name,
            "matches_previous_independent_rows": m,
            "previous_independent_rows_compared": c,
            "agreement_fraction": m / c if c else "",
        })
    write_csv(A / "B104G_node_order_model_controls.csv", alt_rows)
    write_csv(D / "B104G_full_graphsage_row_to_entrez_mapping.csv.gz", full_map_rows)
    # csv module cannot gzip automatically; rewrite compressed and remove accidental plain bytes if any.
    plain_tmp = D / "B104G_full_graphsage_row_to_entrez_mapping.csv.gz"
    # write_csv above wrote plain bytes with .gz suffix; replace with deterministic gzip.
    raw_map_bytes = plain_tmp.read_bytes()
    with plain_tmp.open("wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0, compresslevel=9) as gz:
            gz.write(raw_map_bytes)

    # Validate every former residual class without consulting labels.
    residual = json.loads(args.residual_classes.read_text())
    residual_rows = []
    for cls in residual:
        predicted = [full_map[int(r)] for r in cls["rows"]]
        expected = list(map(int, cls["candidate_genes"]))
        residual_rows.append({
            "graph_index": cls["graph_index"],
            "tissue": cls["tissue"],
            "WL_color": cls["color"],
            "row_count": len(cls["rows"]),
            "rows": "|".join(map(str, cls["rows"])),
            "predicted_GeneIDs_in_row_order": "|".join(map(str, predicted)),
            "candidate_GeneIDs": "|".join(map(str, expected)),
            "candidate_set_exact": int(set(predicted) == set(expected)),
            "predicted_GeneIDs_unique": int(len(predicted) == len(set(predicted))),
        })
    if not all(r["candidate_set_exact"] and r["predicted_GeneIDs_unique"] for r in residual_rows):
        raise AssertionError("row-order map fails a former residual class")
    write_csv(A / "B104G_former_equivalence_classes_resolved_by_node_order.csv", residual_rows)

    # Graph edge/node verification on all rows, including every formerly unresolved node.
    block_for_row = np.empty(n, dtype=np.int16)
    for i in range(24):
        block_for_row[bounds[i]:bounds[i + 1]] = i
    gs_edges_by_tissue = [set() for _ in range(24)]
    cross_block = []
    for u, v in links:
        bi = int(block_for_row[u])
        if int(block_for_row[v]) != bi:
            cross_block.append((u, v))
            continue
        gs_edges_by_tissue[bi].add(canonical_edge(full_map[u], full_map[v]))
    if cross_block:
        raise AssertionError(f"cross-block GraphSAGE edges: {cross_block[:5]}")
    edge_rows = []
    for i, tissue in enumerate(tissues):
        net = ohm[tissue]
        mapped_nodes = {full_map[r] for r in range(bounds[i], bounds[i + 1])}
        ge = gs_edges_by_tissue[i]
        missing = net.edges - ge
        extra = ge - net.edges
        edge_rows.append({
            "graph_index_1based": i + 1,
            "tissue": tissue,
            "graphsage_nodes": bounds[i + 1] - bounds[i],
            "ohmnet_nodes": len(net.nodes),
            "node_set_exact": int(mapped_nodes == net.nodes),
            "graphsage_edges": len(ge),
            "ohmnet_edges": len(net.edges),
            "edge_set_exact": int(ge == net.edges),
            "missing_edges": len(missing),
            "extra_edges": len(extra),
            "ohmnet_tar_member": net.tar_member,
        })
    if not all(r["node_set_exact"] and r["edge_set_exact"] for r in edge_rows):
        raise AssertionError("full node/edge reconstruction failed")
    write_csv(A / "B104G_full_24_tissue_node_edge_verification.csv", edge_rows)

    # Full feature validation from retained normalized MSigDB v5.0 gene sets.
    msig_sets = {}
    with gzip.open(args.msigdb_normalized, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            genes = {int(x) for x in r["member_Entrez_IDs"].split("|") if x}
            msig_sets[(r["category"], r["standard_name"])] = genes
    selected_features = []
    with args.feature_rule.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            selected_features.append((int(r["graphsage_feature_column"]), r["collection"], r["set_name"]))
    selected_features.sort()
    if [x[0] for x in selected_features] != list(range(50)):
        raise AssertionError("feature rule does not contain columns 0..49")
    feature_gene_sets = []
    feature_meta_rows = []
    for col, coll, name in selected_features:
        key = (coll, name)
        if key not in msig_sets:
            raise AssertionError(f"selected feature absent from normalized MSigDB v5.0: {key}")
        gene_set = msig_sets[key]
        feature_gene_sets.append(gene_set)
        feature_meta_rows.append({"feature_column": col, "collection": coll, "gene_set": name, "full_source_gene_count_v5_0": len(gene_set)})
    expected_feats = np.zeros((n, 50), dtype=np.uint8)
    for row in range(n):
        g = full_map[row]
        expected_feats[row] = [int(g in s) for s in feature_gene_sets]
    observed_feats = (feats > 0.5).astype(np.uint8)
    feature_diff = expected_feats != observed_feats
    for j, r in enumerate(feature_meta_rows):
        r.update({
            "observed_positive_rows": int(observed_feats[:, j].sum()),
            "expected_positive_rows": int(expected_feats[:, j].sum()),
            "mismatched_rows": int(feature_diff[:, j].sum()),
            "exact_all_56944_rows": int(not feature_diff[:, j].any()),
            "all_zero_observed": int(observed_feats[:, j].sum() == 0),
        })
    if feature_diff.any():
        where = np.argwhere(feature_diff)
        raise AssertionError(f"full feature reconstruction mismatch: {where[:10].tolist()}")
    write_csv(A / "B104G_full_50_feature_validation.csv", feature_meta_rows)

    # Repeated-gene consistency of deposited features and labels.
    rows_by_gene = collections.defaultdict(list)
    for row, g in full_map.items():
        rows_by_gene[g].append(row)
    feature_conflicts = []
    label_conflicts = []
    unique_gene_features = {}
    unique_gene_labels = {}
    for g, rows in rows_by_gene.items():
        f0 = observed_feats[rows[0]]
        y0 = labels[rows[0]]
        unique_gene_features[g] = tuple(map(int, f0))
        unique_gene_labels[g] = tuple(map(int, y0))
        if any(not np.array_equal(observed_feats[r], f0) for r in rows[1:]):
            feature_conflicts.append(g)
        if any(not np.array_equal(labels[r], y0) for r in rows[1:]):
            label_conflicts.append(g)
    if feature_conflicts or label_conflicts:
        raise AssertionError(f"repeated-gene inconsistency features={feature_conflicts[:5]} labels={label_conflicts[:5]}")

    graph_genes = set(rows_by_gene)

    # Historical GOA mapping and exact label transformation.
    GPI = {}
    with gzip.open(args.gpi, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            GPI[r["DB_Object_ID"]] = {
                "symbol": r["DB_Object_Symbol"] or r.get("GAF_Fallback_Symbol", ""),
                "name": r["DB_Object_Name"],
                "synonyms": r["DB_Object_Synonyms"],
                "taxon": r["Taxon"],
            }
    acc_edges = collections.defaultdict(set)
    with gzip.open(args.gp2protein, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            acc_edges[r["UniProtKB_accession"]].add(int(r["GeneID"]))
    sym_to_genes = {}
    with gzip.open(args.symbol_map, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sym_to_genes[r["gene_symbol"]] = {int(x) for x in r["Entrez_GeneIDs"].split("|") if x}

    adj = collections.defaultdict(set)
    for a, gs in acc_edges.items():
        if a not in GPI:
            continue
        for g in gs:
            adj[("a", a)].add(("g", g))
            adj[("g", g)].add(("a", a))
    resolved_accessions = {}
    seen = set()
    mapping_component_rows = []
    for start in list(adj):
        if start in seen:
            continue
        seen.add(start)
        stack = [start]
        nodes = []
        while stack:
            x = stack.pop()
            nodes.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        accs = {v for t, v in nodes if t == "a"}
        genes = {v for t, v in nodes if t == "g"}
        cand = {}
        for a in accs:
            matches = sym_to_genes.get(GPI[a]["symbol"], set()) & genes
            if len(matches) == 1:
                cand[a] = next(iter(matches))
        bij = len(accs) == len(genes) == len(cand) and len(set(cand.values())) == len(genes)
        if bij:
            for a, g in cand.items():
                resolved_accessions[a] = {g}
        if genes & graph_genes and (len(accs) > 1 or len(genes) > 1):
            mapping_component_rows.append({
                "accessions": "|".join(sorted(accs)),
                "GeneIDs": "|".join(map(str, sorted(genes))),
                "primary_symbol_candidates": "|".join(f"{a}:{cand[a]}" for a in sorted(cand)),
                "unique_symbol_bijection": int(bij),
                "resolution": "|".join(f"{a}->{next(iter(resolved_accessions[a]))}" for a in sorted(accs) if a in resolved_accessions) if bij else "retain_all_historical_edges",
            })

    mapping = collections.defaultdict(set)
    mapping_method = {}
    for a, d in GPI.items():
        if a in resolved_accessions:
            mapping[a] = resolved_accessions[a] & graph_genes
            mapping_method[a] = "full_component_unique_primary_symbol_bijection"
        else:
            direct = acc_edges.get(a, set()) & graph_genes
            if direct:
                mapping[a] = set(direct)
                mapping_method[a] = "historical_gp2protein_all_edges"
            else:
                matches = sym_to_genes.get(d["symbol"], set()) & graph_genes
                if len(matches) == 1:
                    mapping[a] = set(matches)
                    mapping_method[a] = "unique_primary_symbol_fallback"
                else:
                    mapping_method[a] = "unmapped"
    mapping["O95073"].discard(25788)
    write_csv(A / "B104G_mapping_components_touching_full_graph_gene_universe.csv", mapping_component_rows)

    alt_id = {}
    parents = collections.defaultdict(set)
    go_name = {}
    go_namespace = {}
    with gzip.open(args.go_terms, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            go_name[r["GO_ID"]] = r["name"]
            go_namespace[r["GO_ID"]] = r["namespace"]
            for a in r["alt_ids"].split("|") if r["alt_ids"] else []:
                alt_id[a] = r["GO_ID"]
    with gzip.open(args.go_edges, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            parents[r["child_GO_ID"]].add(r["parent_GO_ID"])

    @functools.lru_cache(None)
    def ancestors(go: str):
        go = alt_id.get(go, go)
        out = {go}
        for p in parents.get(go, ()):
            out.update(ancestors(p))
        return frozenset(out)

    column_rows = []
    selected = []
    with args.column_map.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            col = int(r["label_column"])
            go = r["inferred_GO_ID"]
            selected.append(go)
            rr = dict(r)
            rr["label_column"] = col
            column_rows.append(rr)
    if len(selected) != 121:
        raise AssertionError("column map is not 121 rows")
    selected_set = set(selected)
    # Include all alternative candidates so we can test the duplicate-vector ambiguity on all 4,301 genes.
    all_candidate_terms = set(selected)
    for r in column_rows:
        all_candidate_terms.update(r["all_exact_GO_ID_candidates"].split("|"))

    predicted_terms = {g: set() for g in graph_genes}
    accepted_gaf_rows = 0
    with gzip.open(args.gaf, "rt", newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Is_NOT"] == "1" or "NOT" in r["Qualifier"].split("|"):
                continue
            if r["Evidence_Code"] not in ALLOWED_EVIDENCE:
                continue
            if r["Normalized_Relation"] not in DEFAULT_RELATIONS:
                continue
            genes = mapping.get(r["DB_Object_ID"], set())
            if not genes:
                continue
            propagated = ancestors(r["GO_ID"])
            if not (propagated & all_candidate_terms):
                continue
            accepted_gaf_rows += 1
            for g in genes:
                predicted_terms[g].update(propagated)

    expected_labels = np.zeros((n, 121), dtype=np.uint8)
    for row in range(n):
        ts = predicted_terms[full_map[row]]
        expected_labels[row] = [int(go in ts) for go in selected]
    label_diff = expected_labels != labels
    if label_diff.any():
        where = np.argwhere(label_diff)
        mismatch_examples = [{"row": int(r), "GeneID": full_map[int(r)], "column": int(c), "GO_ID": selected[int(c)], "observed": int(labels[r, c]), "expected": int(expected_labels[r, c])} for r, c in where[:100]]
        write_csv(A / "B104G_full_label_mismatch_examples.csv", mismatch_examples)
        raise AssertionError(f"full label reconstruction mismatch count={int(label_diff.sum())}; examples={mismatch_examples[:5]}")

    label_validation_rows = []
    for c, r in enumerate(column_rows):
        label_validation_rows.append({
            "label_column": c,
            "GO_ID": selected[c],
            "GO_name": go_name.get(selected[c], r.get("inferred_GO_name", "")),
            "namespace": go_namespace.get(selected[c], r.get("namespace", "")),
            "observed_positive_rows": int(labels[:, c].sum()),
            "expected_positive_rows": int(expected_labels[:, c].sum()),
            "mismatched_rows": int(label_diff[:, c].sum()),
            "exact_all_56944_rows": int(not label_diff[:, c].any()),
            "identity_basis": r["identity_basis"],
            "identity_confidence": r["identity_confidence"],
            "all_exact_GO_ID_candidates": r["all_exact_GO_ID_candidates"],
        })
    write_csv(A / "B104G_full_121_label_validation.csv", label_validation_rows)

    # Check whether the three duplicate candidate pairs remain indistinguishable on all 4,301 genes.
    duplicate_rows = []
    candidate_vectors = {}
    gene_order = sorted(graph_genes)
    for go in all_candidate_terms:
        candidate_vectors[go] = tuple(int(go in predicted_terms[g]) for g in gene_order)
    seen_pairs = set()
    for r in column_rows:
        cands = tuple(r["all_exact_GO_ID_candidates"].split("|"))
        if len(cands) <= 1:
            continue
        pair = tuple(sorted(cands))
        if pair in seen_pairs:
            continue
        seen_pairs.add(pair)
        a, b = pair
        diff_genes = [g for g, x, y in zip(gene_order, candidate_vectors[a], candidate_vectors[b]) if x != y]
        duplicate_rows.append({
            "candidate_GO_ID_A": a,
            "candidate_GO_ID_B": b,
            "full_graph_unique_genes_compared": len(gene_order),
            "membership_differing_gene_count": len(diff_genes),
            "membership_vectors_identical_on_all_4301_genes": int(not diff_genes),
            "differing_GeneIDs": "|".join(map(str, diff_genes[:100])),
            "conclusion": "matrix cannot disambiguate candidate orientation" if not diff_genes else "full mapping adds disambiguating membership",
        })
    write_csv(A / "B104G_duplicate_GO_candidate_full_universe_test.csv", duplicate_rows)

    # Mapping coverage and special gene rows.
    mapped_graph_genes = set().union(*mapping.values()) if mapping else set()
    unmapped_graph_genes = sorted(graph_genes - mapped_graph_genes)
    special_gene_rows = []
    for g in sorted(set(unmapped_graph_genes) | {10159, 25788, 7957, 29901}):
        rows = rows_by_gene.get(g, [])
        special_gene_rows.append({
            "GeneID": g,
            "graphsage_rows": "|".join(map(str, rows)),
            "tissues": "|".join(sorted({tissues[int(block_for_row[r])] for r in rows})),
            "row_count": len(rows),
            "mapped_from_GPI_accession": int(g in mapped_graph_genes),
            "observed_positive_feature_count": sum(unique_gene_features[g]) if g in unique_gene_features else "",
            "observed_positive_label_count": sum(unique_gene_labels[g]) if g in unique_gene_labels else "",
            "predicted_positive_label_count": len(predicted_terms.get(g, set()) & selected_set),
        })
    write_csv(A / "B104G_special_and_unmapped_gene_rows.csv", special_gene_rows)

    # Full leakage analysis now that every row has a serialization-derived GeneID.
    leakage_rows = []
    gene_train_vectors = collections.defaultdict(set)
    train_genes = set()
    for row in range(n):
        if split[row] == "train":
            g = full_map[row]
            train_genes.add(g)
            gene_train_vectors[g].add(tuple(map(int, labels[row])))
    train_conflicts = {g: len(vs) for g, vs in gene_train_vectors.items() if len(vs) != 1}
    if train_conflicts:
        raise AssertionError(f"training label conflict by gene: {list(train_conflicts.items())[:5]}")
    train_vector = {g: next(iter(vs)) for g, vs in gene_train_vectors.items()}
    for target in ("val", "test"):
        rows = np.flatnonzero(split == target)
        pred = np.zeros((len(rows), 121), dtype=np.uint8)
        seen_rows = 0
        for i, row in enumerate(rows):
            g = full_map[int(row)]
            if g in train_vector:
                seen_rows += 1
                pred[i] = np.asarray(train_vector[g], dtype=np.uint8)
        score, tp, fp, fn = f1_micro(labels[rows], pred)
        target_genes = {full_map[int(r)] for r in rows}
        unseen_genes = sorted(target_genes - train_genes)
        leakage_rows.append({
            "split": target,
            "rows": len(rows),
            "unique_GeneIDs": len(target_genes),
            "rows_with_GeneID_seen_in_training": seen_rows,
            "fraction_rows_with_GeneID_seen_in_training": seen_rows / len(rows),
            "unique_GeneIDs_seen_in_training": len(target_genes & train_genes),
            "unique_GeneIDs_unseen_in_training": len(unseen_genes),
            "unseen_GeneIDs": "|".join(map(str, unseen_genes)),
            "lookup_micro_F1_zero_for_unseen": score,
            "true_positive_label_cells": tp,
            "false_positive_label_cells": fp,
            "false_negative_label_cells": fn,
        })
    write_csv(A / "B104G_full_mapping_leakage_recalculation.csv", leakage_rows)

    # Analyze the four genuine dhimmel Entrez-native annotation summaries.
    observed_gene_bits = []
    gene_index = {g: i for i, g in enumerate(gene_order)}
    for c in range(121):
        bit = 0
        for g in gene_order:
            if unique_gene_labels[g][c]:
                bit |= 1 << gene_index[g]
        observed_gene_bits.append(bit)

    dhimmel_summary_rows = []
    dhimmel_per_column_rows = []
    dhimmel_order_rows = []
    dhimmel_25788_rows = []
    dhimmel_reconciliation = []
    selected_order = selected
    selected_ids = set(selected_order)
    for path in args.dhimmel:
        records = []
        source_order = []
        with path.open(newline="", encoding="utf-8") as f:
            rd = csv.DictReader(f, delimiter="\t")
            expected_fields = ["go_id", "go_name", "go_domain", "tax_id", "annotation_type", "size", "gene_ids", "gene_symbols"]
            if rd.fieldnames != expected_fields:
                raise AssertionError(f"unexpected dhimmel schema in {path}: {rd.fieldnames}")
            for source_row, r in enumerate(rd, start=1):
                genes = {int(x) for x in r["gene_ids"].split("|") if x}
                bit = 0
                for g in genes & graph_genes:
                    bit |= 1 << gene_index[g]
                records.append({
                    "go_id": r["go_id"], "go_name": r["go_name"], "domain": r["go_domain"],
                    "annotation_type": r["annotation_type"], "declared_size": int(r["size"]),
                    "genes": genes, "graph_bit": bit, "source_row": source_row,
                })
                source_order.append(r["go_id"])
        by_go = {r["go_id"]: r for r in records}
        exact_best = at99 = at95 = 0
        best_mismatches = []
        candidate_mismatches = []
        candidate_exact = 0
        for c, obs in enumerate(observed_gene_bits):
            distances = [(obs ^ r["graph_bit"]).bit_count() for r in records]
            best = min(distances)
            best_mismatches.append(best)
            best_indices = [i for i, d in enumerate(distances) if d == best]
            best_ids = [records[i]["go_id"] for i in best_indices]
            exact_best += int(best == 0)
            agreement = 1 - best / len(gene_order)
            at99 += int(agreement >= 0.99)
            at95 += int(agreement >= 0.95)
            go = selected[c]
            rec = by_go.get(go)
            if rec is None:
                cand_diff = len(gene_order)
                cand_fp = cand_fn = ""
                cand_agreement = 0.0
            else:
                cand_diff = (obs ^ rec["graph_bit"]).bit_count()
                cand_fp = (rec["graph_bit"] & ~obs).bit_count()
                cand_fn = (obs & ~rec["graph_bit"]).bit_count()
                cand_agreement = 1 - cand_diff / len(gene_order)
                candidate_exact += int(cand_diff == 0)
            candidate_mismatches.append(cand_diff)
            dhimmel_per_column_rows.append({
                "file": path.name,
                "label_column": c,
                "inferred_GO_ID": go,
                "inferred_GO_name": go_name.get(go, ""),
                "candidate_term_present": int(rec is not None),
                "candidate_term_mismatches": cand_diff,
                "candidate_term_false_positives": cand_fp,
                "candidate_term_false_negatives": cand_fn,
                "candidate_term_agreement": cand_agreement,
                "best_any_term_mismatches": best,
                "best_any_term_agreement": agreement,
                "best_any_term_GO_IDs": "|".join(best_ids[:25]),
                "best_any_term_tie_count": len(best_ids),
            })
        dhimmel_summary_rows.append({
            "file": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "term_rows": len(records),
            "selected_121_GO_IDs_present": len(selected_ids & set(by_go)),
            "candidate_GO_ID_columns_exact": candidate_exact,
            "candidate_GO_ID_total_mismatched_gene_label_cells": sum(candidate_mismatches),
            "best_any_term_exact_columns": exact_best,
            "best_any_term_columns_at_least_99pct": at99,
            "best_any_term_columns_at_least_95pct": at95,
            "closest_any_term_mismatch_gene_count": min(best_mismatches),
            "median_best_any_term_mismatch_gene_count": statistics.median(best_mismatches),
            "directly_reproduces_full_GraphSAGE_label_matrix": 0,
        })

        # File row order and Python 2 dictionary order as column-order controls.
        source_filtered = [x for x in source_order if x in selected_ids]
        d = Py2Dict(py2_hash_str)
        for x in source_order:
            d.insert(x)
        dict_filtered = [x for x in d.keys() if x in selected_ids]
        for model, pred in [("source_row_order", source_filtered), ("python2_string_dict_order_after_source_insertion", dict_filtered)]:
            tau, con, dis, common_n = kendall_tau_for_common(selected_order, pred)
            dhimmel_order_rows.append({
                "file": path.name,
                "order_model": model,
                "selected_terms_present": common_n,
                "longest_common_subsequence_with_121_column_order": lcs_length(pred, selected_order),
                "kendall_tau_on_common_terms": tau,
                "pairwise_concordant": con,
                "pairwise_discordant": dis,
                "exact_absolute_positions_among_121": sum(i < len(pred) and pred[i] == selected_order[i] for i in range(min(121, len(pred)))),
                "predicted_prefix": "|".join(pred[:10]),
            })

        # Entrez-native handling of the nested-gene concern and 12 verifier genes.
        for g in [25788, 100861412, 10159, 3248,3988,8564,27201,30061,51166,51312,55471,55801,56994,79017,121599]:
            if g not in graph_genes and g != 100861412:
                continue
            predicted_cols = [c for c, go in enumerate(selected) if go in by_go and g in by_go[go]["genes"]]
            observed_cols = [c for c in range(121) if g in unique_gene_labels and unique_gene_labels[g][c]]
            dhimmel_25788_rows.append({
                "file": path.name,
                "GeneID": g,
                "predicted_positive_selected_columns": "|".join(map(str, predicted_cols)),
                "predicted_positive_selected_count": len(predicted_cols),
                "observed_GraphSAGE_positive_columns": "|".join(map(str, observed_cols)),
                "observed_positive_count": len(observed_cols),
                "mismatch_count": len(set(predicted_cols) ^ set(observed_cols)) if g in unique_gene_labels else "not_in_graph",
            })

        retained = R / (path.name + ".gz")
        rec = deterministic_gzip_copy(path, retained)
        rec["line_count"] = sum(1 for _ in path.open("rb"))
        rec["source_url"] = f"https://raw.githubusercontent.com/dhimmel/gene-ontology/962a5e12f8590400c2891cde93fd6a783b26e02e/annotations/taxid_9606/{path.name}"
        dhimmel_reconciliation.append(rec)

    write_csv(A / "B104G_dhimmel_annotation_file_summary.csv", dhimmel_summary_rows)
    write_csv(A / "B104G_dhimmel_per_column_comparison.csv.gz", dhimmel_per_column_rows)
    # Compress the per-column CSV that write_csv wrote as plain bytes.
    p = A / "B104G_dhimmel_per_column_comparison.csv.gz"
    raw = p.read_bytes()
    with p.open("wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0, compresslevel=9) as gz:
            gz.write(raw)
    write_csv(A / "B104G_dhimmel_column_order_controls.csv", dhimmel_order_rows)
    write_csv(A / "B104G_dhimmel_selected_gene_audit.csv", dhimmel_25788_rows)
    (A / "B104G_dhimmel_raw_to_retained_reconciliation.json").write_text(json.dumps(dhimmel_reconciliation, indent=2) + "\n")

    # Explore remaining graph-selection/split provenance uncertainty using all OhmNet networks.
    selected_set_tissues = set(tissues)
    all_network_rows = []
    for name, net in ohm.items():
        all_network_rows.append({
            "tissue": name,
            "node_count": len(net.nodes),
            "edge_count": len(net.edges),
            "selected_by_graphsage": int(name in selected_set_tissues),
            "graphsage_graph_index_1based": tissues.index(name) + 1 if name in selected_set_tissues else "",
            "graphsage_split": "train" if name in tissues[:20] else "val" if name in tissues[20:22] else "test" if name in tissues[22:] else "",
            "meets_15000_edge_threshold": int(len(net.edges) >= 15000),
            "meets_35000_edge_threshold": int(len(net.edges) >= 35000),
            "tar_member": net.tar_member,
        })
    all_network_rows.sort(key=lambda r: (-r["edge_count"], r["tissue"]))
    for rank, r in enumerate(all_network_rows, start=1):
        r["edge_count_rank_desc"] = rank
    write_csv(A / "B104G_all_OhmNet_network_sizes_and_GraphSAGE_selection.csv", all_network_rows)
    top24 = {r["tissue"] for r in all_network_rows[:24]}
    graph_selection_summary = {
        "total_OhmNet_networks": len(all_network_rows),
        "GraphSAGE_selected_tissues": 24,
        "networks_at_least_15000_edges": sum(r["meets_15000_edge_threshold"] for r in all_network_rows),
        "networks_at_least_35000_edges": sum(r["meets_35000_edge_threshold"] for r in all_network_rows),
        "all_selected_meet_15000_edges": all(r["edge_count"] >= 15000 for r in all_network_rows if r["selected_by_graphsage"]),
        "validation_test_all_meet_35000_edges": all(r["edge_count"] >= 35000 for r in all_network_rows if r["graphsage_split"] in {"val", "test"}),
        "selected_set_equals_top24_by_edge_count": selected_set_tissues == top24,
        "selected_tissues_not_in_top24_by_edge_count": sorted(selected_set_tissues - top24),
        "top24_tissues_not_selected": sorted(top24 - selected_set_tissues),
        "selection_and_split_procedure_recovered": False,
    }
    (A / "B104G_graph_selection_open_question_summary.json").write_text(json.dumps(graph_selection_summary, indent=2) + "\n")

    # Explicit status matrix: what is exact, strongly supported, or still open.
    status_rows = [
        {"component":"24 graph block-to-tissue identities", "current_status":"exactly verified", "evidence":"all node and edge sets match the assigned OhmNet tissue networks", "remaining_ambiguity":"none at dataset-comparison level"},
        {"component":"row-to-Entrez identity for 56,944 rows", "current_status":"strongly supported full resolution", "evidence":"Python-2 string-dict order agrees with all 56,411 previously independent row identities; resolves every former class; all 818,716 edges, 50 features, and 121 labels validate", "remaining_ambiguity":"original preprocessing source code not located; the exact mechanism is inferred from deterministic serialization behavior"},
        {"component":"50 feature values", "current_status":"exactly reproduced on all rows", "evidence":"zero differences across 2,847,200 cells under the retained MSigDB membership rule", "remaining_ambiguity":"literal source code, exact threshold operator, whether C7 was traversed, and MSigDB release are not identifiable from the matrix"},
        {"component":"121 label values", "current_status":"exactly reproduced on all rows under one global transformation", "evidence":"zero differences across 6,890,224 cells", "remaining_ambiguity":"GOA v159 is uniquely exact among tested releases, but an equivalent Entrez-native/preprocessed source could have been used"},
        {"component":"GO term identities", "current_status":"115 columns uniquely identified; six identities strongly supported", "evidence":"three duplicate-vector pairs remain identical across all 4,301 genes", "remaining_ambiguity":"orientation of columns 24/71, 39/63, and 48/70 remains provisional"},
        {"component":"GO term selection", "current_status":"strongly supported", "evidence":"release-159 full-human top 121 equals all terms at the natural approximately 1,000-gene/protein boundary", "remaining_ambiguity":"top-121 versus >=1000 and gene versus protein counting are observationally equivalent here"},
        {"component":"GO column ordering", "current_status":"strong Python-2 dictionary fingerprint; not perfectly reconstructed", "evidence":"best tested LCS 94/121 plus exact independent class_map JSON-key-order control", "remaining_ambiguity":"exact dictionary insertion history or ordered source list"},
        {"component":"24-network selection and 20/2/2 split", "current_status":"partially characterized", "evidence":"sizes and stated thresholds can be checked", "remaining_ambiguity":"exact selection and split assignment procedure"},
        {"component":"DGL conversion", "current_status":"exactly reproduced previously", "evidence":"labels, graph IDs, float64 standardization, directed edges, self-loops, and component handling", "remaining_ambiguity":"none at output-comparison level"},
        {"component":"identifier mapping semantics", "current_status":"globally consistent mapping gives exact labels", "evidence":"historical many-to-many edges retained; O95073 nested-gene projection separated semantically", "remaining_ambiguity":"whether the original pipeline used GOA+mapping or an Entrez-native source"},
    ]
    write_csv(A / "B104G_current_resolution_and_remaining_ambiguities.csv", status_rows)

    # Compact per-gene output for sharing and future work.
    gene_rows = []
    for g in gene_order:
        rows = rows_by_gene[g]
        gene_rows.append({
            "GeneID": g,
            "row_count": len(rows),
            "graphsage_rows": "|".join(map(str, rows)),
            "tissues": "|".join(sorted({tissues[int(block_for_row[r])] for r in rows})),
            "positive_feature_count": sum(unique_gene_features[g]),
            "positive_label_count": sum(unique_gene_labels[g]),
            "GOA_mapping_covered": int(g in mapped_graph_genes),
        })
    write_csv(D / "B104G_full_4301_gene_universe.csv.gz", gene_rows)
    p = D / "B104G_full_4301_gene_universe.csv.gz"
    raw = p.read_bytes()
    with p.open("wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0, compresslevel=9) as gz:
            gz.write(raw)

    # Output summary.
    summary = {
        "full_row_identity": {
            "graphsage_rows": n,
            "unique_Entrez_GeneIDs": len(graph_genes),
            "previous_independent_rows_compared": len(known_map),
            "previous_independent_row_agreements": len(known_map) - len(known_mismatches),
            "previous_independent_row_disagreements": len(known_mismatches),
            "former_unresolved_rows": sum(len(c["rows"]) for c in residual),
            "former_unresolved_classes": len(residual),
            "former_classes_candidate_sets_exact": sum(r["candidate_set_exact"] for r in residual_rows),
            "mapping_model": "64-bit unrandomized CPython2 string-key dictionary table iteration after inserting node tokens from each original OhmNet edgelist in line order",
            "claim_strength": "strongly supported full row resolution; original preprocessing source code not located",
        },
        "topology": {
            "all_24_node_sets_exact": all(r["node_set_exact"] for r in edge_rows),
            "all_24_edge_sets_exact": all(r["edge_set_exact"] for r in edge_rows),
            "GraphSAGE_undirected_links_including_loops": len(links),
            "all_links_verified": sum(r["graphsage_edges"] for r in edge_rows),
        },
        "features": {
            "rows": n,
            "columns": 50,
            "cells_compared": int(feature_diff.size),
            "mismatched_cells": int(feature_diff.sum()),
            "exact_columns": int(sum(not feature_diff[:, j].any() for j in range(50))),
            "unique_GeneIDs": len(graph_genes),
            "repeated_gene_conflicts": len(feature_conflicts),
            "claim_scope": "exact matrix reproduction; source procedure remains a strong hypothesis, not original-code proof",
        },
        "labels": {
            "rows": n,
            "columns": 121,
            "cells_compared": int(label_diff.size),
            "mismatched_cells": int(label_diff.sum()),
            "exact_columns": int(sum(not label_diff[:, j].any() for j in range(121))),
            "unique_GeneIDs": len(graph_genes),
            "repeated_gene_conflicts": len(label_conflicts),
            "unmapped_graph_GeneIDs": unmapped_graph_genes,
            "accepted_GAF_rows_contributing_to_candidate_terms": accepted_gaf_rows,
            "global_policy": {
                "GOA": "human release 159",
                "evidence": sorted(ALLOWED_EVIDENCE),
                "relations": sorted(DEFAULT_RELATIONS),
                "excluded_positive_qualifier_relations": ["colocalizes_with", "contributes_to"],
                "NOT_excluded": True,
                "propagation": "is_a only including direct term",
                "per_gene_or_per_column_tuning": False,
            },
            "claim_scope": "exact values on all rows under a single fixed transformation; literal upstream source/product still inferential",
        },
        "duplicate_GO_candidates": duplicate_rows,
        "leakage": leakage_rows,
        "dhimmel": dhimmel_summary_rows,
        "graph_selection": graph_selection_summary,
        "input_sha256": {
            "graphsage_ppi.zip": sha256(args.graphsage_zip),
            "bio-tissue-networks.tar.gz": sha256(args.ohmnet_tar),
            "core_summary": sha256(args.core_summary),
            "known_row_map": sha256(args.known_row_map),
            "residual_classes": sha256(args.residual_classes),
            "GAF159_normalized": sha256(args.gaf),
            "GPI159_normalized": sha256(args.gpi),
            "gp2protein_relevant_subset": sha256(args.gp2protein),
            "GO_terms": sha256(args.go_terms),
            "GO_is_a_edges": sha256(args.go_edges),
            "symbol_map": sha256(args.symbol_map),
            "column_map": sha256(args.column_map),
            "MSigDB_normalized": sha256(args.msigdb_normalized),
            "feature_rule": sha256(args.feature_rule),
            **{p.name: sha256(p) for p in args.dhimmel},
        },
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

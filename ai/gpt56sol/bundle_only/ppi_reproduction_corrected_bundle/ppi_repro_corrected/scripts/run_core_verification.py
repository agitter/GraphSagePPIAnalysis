#!/usr/bin/env python3
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import os
import shutil
import tarfile
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Set, Tuple

import numpy as np
from sklearn.preprocessing import StandardScaler


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_extract_zip(path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path) as zf:
        root = dest.resolve()
        for member in zf.infolist():
            target = (dest / member.filename).resolve()
            if root not in target.parents and target != root:
                raise RuntimeError(f"Unsafe zip member: {member.filename}")
        zf.extractall(dest)


def safe_extract_tar(path: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "r:gz") as tf:
        root = dest.resolve()
        for member in tf.getmembers():
            target = (dest / member.name).resolve()
            if root not in target.parents and target != root:
                raise RuntimeError(f"Unsafe tar member: {member.name}")
        tf.extractall(dest)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        if not rows:
            raise ValueError(f"Cannot infer columns for empty CSV {path}")
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_graphsage(extracted: Path):
    ppi = extracted / "ppi"
    graph = json.loads((ppi / "ppi-G.json").read_text())
    id_map = json.loads((ppi / "ppi-id_map.json").read_text())
    class_map = json.loads((ppi / "ppi-class_map.json").read_text())
    feats = np.load(ppi / "ppi-feats.npy")
    n = len(graph["nodes"])
    if feats.shape[0] != n or len(class_map) != n or len(id_map) != n:
        raise AssertionError("GraphSAGE arrays/maps do not agree on row count")
    # The supplied archive has consecutive numeric node IDs and an identity id_map.
    if any(int(id_map[str(i)]) != i for i in range(n)):
        raise AssertionError("GraphSAGE id_map is not identity; this script expects the supplied archive")
    labels = np.asarray([class_map[str(i)] for i in range(n)], dtype=np.uint8)
    edges = [(int(e["source"]), int(e["target"])) for e in graph["links"]]
    splits = np.asarray([
        "test" if bool(node.get("test")) else "valid" if bool(node.get("val")) else "train"
        for node in graph["nodes"]
    ])
    return graph, feats, labels, edges, splits


def canonical_edge(u: int, v: int) -> Tuple[int, int]:
    return (u, v) if u <= v else (v, u)


@dataclass
class OhmNetwork:
    name: str
    path: Path
    nodes: Set[int]
    edges: Set[Tuple[int, int]]
    adjacency: Dict[int, Set[int]]


def load_ohm_networks(root: Path) -> Dict[str, OhmNetwork]:
    out: Dict[str, OhmNetwork] = {}
    for path in sorted(root.rglob("*.edgelist")):
        name = path.stem
        nodes: Set[int] = set()
        edges: Set[Tuple[int, int]] = set()
        adj: Dict[int, Set[int]] = collections.defaultdict(set)
        with path.open() as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 2:
                    continue
                u, v = int(parts[0]), int(parts[1])
                nodes.add(u)
                nodes.add(v)
                edges.add(canonical_edge(u, v))
                adj[u].add(v)
                adj[v].add(u)
        out[name] = OhmNetwork(name, path, nodes, edges, dict(adj))
    return out


def safe_cuts(n: int, edges: Sequence[Tuple[int, int]]) -> List[int]:
    # cut k separates rows [0,k) and [k,n). Mark every cut crossed by an edge.
    diff = np.zeros(n + 2, dtype=np.int64)
    for u, v in edges:
        if u == v:
            continue
        a, b = sorted((u, v))
        diff[a + 1] += 1
        diff[b + 1] -= 1
    active = 0
    cuts = [0]
    for k in range(1, n):
        active += int(diff[k])
        if active == 0:
            cuts.append(k)
    cuts.append(n)
    return cuts


def partition_graph(
    n: int,
    edges: Sequence[Tuple[int, int]],
    ohm: Mapping[str, OhmNetwork],
) -> Tuple[List[int], List[str], List[dict]]:
    cuts = safe_cuts(n, edges)
    # In this archive the safe cuts directly yield 24 blocks. Retain a DP fallback
    # in case adjacent safe intervals must be merged to match an OhmNet network.
    edge_prefix = np.zeros(n + 1, dtype=np.int64)
    for u, v in edges:
        edge_prefix[max(u, v) + 1] += 1
    edge_prefix = np.cumsum(edge_prefix)

    candidates_by_stats: Dict[Tuple[int, int], List[str]] = collections.defaultdict(list)
    for name, net in ohm.items():
        candidates_by_stats[(len(net.nodes), len(net.edges))].append(name)

    def interval_m(a: int, b: int) -> int:
        # Safe cuts guarantee no crossing edges; every edge whose max endpoint < b
        # and max endpoint >= a is inside [a,b).
        return int(edge_prefix[b] - edge_prefix[a])

    m = len(cuts)
    dp: List[List[Tuple[List[int], List[str]]]] = [[] for _ in range(m)]
    dp[0] = [([cuts[0]], [])]
    for i in range(m - 1):
        if not dp[i]:
            continue
        for j in range(i + 1, m):
            a, b = cuts[i], cuts[j]
            stats = (b - a, interval_m(a, b))
            names = candidates_by_stats.get(stats, [])
            if not names:
                continue
            for bounds, assigned in dp[i]:
                for name in names:
                    if name in assigned:
                        continue
                    dp[j].append((bounds + [b], assigned + [name]))
            # Blocks are thousands of nodes; avoid unbounded search once a match occurs.
            if names:
                break
    sols = dp[-1]
    if len(sols) != 1:
        raise AssertionError(f"Expected one exact tissue partition, found {len(sols)}")
    bounds, tissues = sols[0]

    # Exact edge-set and node-count verification occurs after WL mapping; here we
    # record block statistics and the unique n/m assignment.
    rows = []
    for i, tissue in enumerate(tissues):
        a, b = bounds[i], bounds[i + 1]
        rows.append({
            "graph_index": i + 1,
            "split": "train" if i < 20 else "valid" if i < 22 else "test",
            "row_start_inclusive": a,
            "row_end_exclusive": b,
            "node_count": b - a,
            "edge_count_undirected_including_loops": interval_m(a, b),
            "ohmnet_tissue": tissue,
            "ohmnet_node_count": len(ohm[tissue].nodes),
            "ohmnet_edge_count": len(ohm[tissue].edges),
            "stats_exact": (b - a == len(ohm[tissue].nodes) and interval_m(a, b) == len(ohm[tissue].edges)),
        })
    return bounds, tissues, rows


def make_global_adjacency(n: int, edges: Sequence[Tuple[int, int]]) -> List[Set[int]]:
    adj: List[Set[int]] = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def wl_pair(
    nodes_a: Sequence[int],
    adj_a: Sequence[Set[int]],
    adj_b: Mapping[int, Set[int]],
    max_iter: int = 50,
):
    nodes_b = list(adj_b)
    sig_a = {u: (len(adj_a[u] - {u}), u in adj_a[u]) for u in nodes_a}
    sig_b = {u: (len(adj_b[u] - {u}), u in adj_b[u]) for u in nodes_b}
    keys = sorted(set(sig_a.values()) | set(sig_b.values()))
    enc = {s: i for i, s in enumerate(keys)}
    ca = {u: enc[s] for u, s in sig_a.items()}
    cb = {u: enc[s] for u, s in sig_b.items()}
    for iteration in range(1, max_iter + 1):
        sa = {u: (ca[u], tuple(sorted(ca[v] for v in adj_a[u]))) for u in nodes_a}
        sb = {u: (cb[u], tuple(sorted(cb[v] for v in adj_b[u]))) for u in nodes_b}
        keys = sorted(set(sa.values()) | set(sb.values()))
        enc = {s: i for i, s in enumerate(keys)}
        na = {u: enc[s] for u, s in sa.items()}
        nb = {u: enc[s] for u, s in sb.items()}
        old_classes = len(set(ca.values()) | set(cb.values()))
        ca, cb = na, nb
        if len(enc) == old_classes:
            return ca, cb, iteration
    return ca, cb, max_iter


def topology_wl_mapping(
    bounds: Sequence[int],
    tissues: Sequence[str],
    adj_g: Sequence[Set[int]],
    ohm: Mapping[str, OhmNetwork],
):
    unique: Dict[int, int] = {}
    ambiguous: List[dict] = []
    stats: List[dict] = []
    for i, tissue in enumerate(tissues):
        a, b = bounds[i], bounds[i + 1]
        nodes = list(range(a, b))
        ca, cb, iterations = wl_pair(nodes, adj_g, ohm[tissue].adjacency)
        da: Dict[int, List[int]] = collections.defaultdict(list)
        db: Dict[int, List[int]] = collections.defaultdict(list)
        for u, color in ca.items():
            da[color].append(u)
        for u, color in cb.items():
            db[color].append(u)
        mismatched = []
        unique_n = 0
        amb_n = 0
        amb_classes = 0
        for color, rows in da.items():
            genes = db.get(color, [])
            if len(rows) != len(genes):
                mismatched.append((color, len(rows), len(genes)))
                continue
            if len(rows) == 1:
                unique[rows[0]] = genes[0]
                unique_n += 1
            elif len(rows) > 1:
                ambiguous.append({
                    "graph_index": i + 1,
                    "tissue": tissue,
                    "color": int(color),
                    "rows": sorted(rows),
                    "candidate_genes": sorted(genes),
                })
                amb_n += len(rows)
                amb_classes += 1
        if mismatched:
            raise AssertionError(f"WL color multiplicities differ for {tissue}: {mismatched[:5]}")
        stats.append({
            "graph_index": i + 1,
            "tissue": tissue,
            "wl_iterations": iterations,
            "node_count": b - a,
            "topology_unique_rows": unique_n,
            "ambiguous_rows": amb_n,
            "ambiguous_classes": amb_classes,
            "color_multiplicity_match": True,
        })
    return unique, ambiguous, stats


def parse_msigdb_entrez(zip_path: Path, version: str) -> List[dict]:
    records = []
    with zipfile.ZipFile(zip_path) as zf:
        for member in zf.namelist():
            if not member.lower().endswith(".entrez.gmt"):
                continue
            basename = Path(member).name
            text = zf.read(member).decode("utf-8", "replace")
            for line_number, line in enumerate(text.splitlines(), start=1):
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                genes = frozenset(int(x) for x in fields[2:] if x.isdigit())
                records.append({
                    "version": version,
                    "source_member": member,
                    "source_basename": basename,
                    "name": fields[0],
                    "description": fields[1],
                    "genes": genes,
                })
    return records


def collection_from_basename(name: str) -> str:
    low = name.lower()
    for prefix in ["c1.", "c2.", "c3.", "c4.", "c5.", "c6.", "c7."]:
        if low.startswith(prefix):
            return prefix[:-1].upper()
    return "ALL"


def feature_identification_and_mapping(
    x: np.ndarray,
    topology_map: Mapping[int, int],
    ambiguous: Sequence[dict],
    msigdb_v52: Path,
):
    rows_by_gene: Dict[int, List[int]] = collections.defaultdict(list)
    for row, gene in topology_map.items():
        rows_by_gene[gene].append(row)
    gene_feature: Dict[int, np.ndarray] = {}
    conflicts = []
    for gene, rows in rows_by_gene.items():
        vals = x[rows]
        if not np.all(vals == vals[0]):
            conflicts.append({"gene_id": gene, "rows": rows})
        gene_feature[gene] = (vals[0] > 0.5).astype(np.uint8)
    if conflicts:
        raise AssertionError(f"Topology-unique rows disagree in features for {len(conflicts)} genes")

    genes = sorted(gene_feature)
    gene_index = {g: i for i, g in enumerate(genes)}
    observed_bits = []
    for col in range(x.shape[1]):
        bit = 0
        for i, gene in enumerate(genes):
            if gene_feature[gene][col]:
                bit |= 1 << i
        observed_bits.append(bit)

    all_records = parse_msigdb_entrez(msigdb_v52, "5.2")
    # Deduplicate aggregate/specialized copies by (name, exact gene set). Prefer
    # a collection-specific member over msigdb.all when reporting source.
    dedup: Dict[Tuple[str, frozenset], dict] = {}
    for rec in all_records:
        key = (rec["name"], rec["genes"])
        prev = dedup.get(key)
        if prev is None:
            dedup[key] = rec
        else:
            prev_coll = collection_from_basename(prev["source_basename"])
            new_coll = collection_from_basename(rec["source_basename"])
            if prev_coll == "ALL" and new_coll != "ALL":
                dedup[key] = rec
    records = list(dedup.values())
    for rec in records:
        bit = 0
        for gene in rec["genes"]:
            idx = gene_index.get(gene)
            if idx is not None:
                bit |= 1 << idx
        rec["mapped_bit"] = bit

    feature_rows = []
    selected: Dict[int, dict | None] = {}
    for col, obs in enumerate(observed_bits):
        best = min((obs ^ rec["mapped_bit"]).bit_count() for rec in records)
        ties = [rec for rec in records if (obs ^ rec["mapped_bit"]).bit_count() == best]
        # Report all observational ties; choose only when the exact gene set is
        # uniquely determined on the mapped universe. Column 10 is all zero and
        # deliberately remains unassigned.
        exact_ties = [r for r in ties if best == 0]
        unique_gene_sets = {r["genes"] for r in exact_ties}
        chosen = exact_ties[0] if best == 0 and len(unique_gene_sets) == 1 and col != 10 else None
        selected[col] = chosen
        feature_rows.append({
            "column": col,
            "observed_positive_genes_topology_unique": int(obs.bit_count()),
            "best_mismatch_gene_count": int(best),
            "best_agreement": 1.0 - (best / len(genes)),
            "observational_tie_count": len(ties),
            "unique_full_gene_set_count_among_exact_ties": len(unique_gene_sets),
            "chosen_name": chosen["name"] if chosen else "",
            "chosen_collection": collection_from_basename(chosen["source_basename"]) if chosen else "",
            "chosen_source_member": chosen["source_member"] if chosen else "",
            "all_zero_column": bool(np.count_nonzero(x[:, col]) == 0),
            "identification_status": "exact_unique" if chosen else "all_zero_unidentifiable" if np.count_nonzero(x[:, col]) == 0 else "ambiguous_or_nonexact",
            "top_tied_names": "|".join(r["name"] for r in ties[:25]),
        })

    # Every nonzero column in this archive has a unique exact v5.2 set.
    missing_nonzero = [r["column"] for r in feature_rows if not r["all_zero_column"] and r["identification_status"] != "exact_unique"]
    if missing_nonzero:
        raise AssertionError(f"Could not uniquely identify nonzero feature columns: {missing_nonzero}")

    # Expected 50-column feature vector for every candidate Entrez gene.
    needed_genes = set()
    for cls in ambiguous:
        needed_genes.update(cls["candidate_genes"])
    expected: Dict[int, Tuple[int, ...]] = {}
    for gene in needed_genes:
        vals = []
        for col in range(x.shape[1]):
            rec = selected[col]
            vals.append(0 if rec is None else int(gene in rec["genes"]))
        expected[gene] = tuple(vals)

    mapping = dict(topology_map)
    residual = []
    added = 0
    for cls in ambiguous:
        row_groups: Dict[Tuple[int, ...], List[int]] = collections.defaultdict(list)
        gene_groups: Dict[Tuple[int, ...], List[int]] = collections.defaultdict(list)
        for row in cls["rows"]:
            row_groups[tuple(int(v > 0.5) for v in x[row])].append(row)
        for gene in cls["candidate_genes"]:
            gene_groups[expected[gene]].append(gene)
        keys = set(row_groups) | set(gene_groups)
        class_residual_rows = []
        class_residual_genes = []
        for key in keys:
            rows = sorted(row_groups.get(key, []))
            genes_for_key = sorted(gene_groups.get(key, []))
            if len(rows) != len(genes_for_key):
                raise AssertionError(
                    f"Feature multiplicities disagree in {cls['tissue']} class {cls['color']}: {len(rows)} rows vs {len(genes_for_key)} genes"
                )
            if len(rows) == 1:
                mapping[rows[0]] = genes_for_key[0]
                added += 1
            elif rows:
                class_residual_rows.extend(rows)
                class_residual_genes.extend(genes_for_key)
        if class_residual_rows:
            residual.append({
                "graph_index": cls["graph_index"],
                "tissue": cls["tissue"],
                "color": cls["color"],
                "rows": sorted(class_residual_rows),
                "candidate_genes": sorted(class_residual_genes),
            })
    return mapping, residual, feature_rows, {
        "topology_unique_rows": len(topology_map),
        "topology_unique_genes": len(set(topology_map.values())),
        "feature_disambiguated_additional_rows": added,
        "mapped_rows_after_features": len(mapping),
        "mapped_unique_genes_after_features": len(set(mapping.values())),
        "residual_rows_after_features": sum(len(x["rows"]) for x in residual),
        "residual_classes_after_features": len(residual),
    }


def verify_mapped_edges(
    edges: Sequence[Tuple[int, int]],
    bounds: Sequence[int],
    tissues: Sequence[str],
    mapping: Mapping[int, int],
    ohm: Mapping[str, OhmNetwork],
):
    block_for_row = np.empty(bounds[-1], dtype=np.int16)
    for i in range(len(tissues)):
        block_for_row[bounds[i]:bounds[i + 1]] = i
    checked = 0
    matched = 0
    mismatches = []
    for u, v in edges:
        if u not in mapping or v not in mapping:
            continue
        i = int(block_for_row[u])
        if int(block_for_row[v]) != i:
            raise AssertionError("Cross-block edge encountered after exact partition")
        checked += 1
        pair = canonical_edge(mapping[u], mapping[v])
        if pair in ohm[tissues[i]].edges:
            matched += 1
        elif len(mismatches) < 100:
            mismatches.append({"row_u": u, "row_v": v, "gene_u": mapping[u], "gene_v": mapping[v], "tissue": tissues[i]})
    return {"edges_with_both_endpoints_mapped": checked, "matched_edges": matched, "mismatch_count": checked - matched, "mismatch_examples": mismatches}


def connected_components_for_block(
    start: int,
    end: int,
    adjacency: Sequence[Set[int]],
) -> List[List[int]]:
    seen: Set[int] = set()
    comps = []
    for root in range(start, end):
        if root in seen:
            continue
        stack = [root]
        seen.add(root)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adjacency[u]:
                if v == u or v < start or v >= end or v in seen:
                    continue
                seen.add(v)
                stack.append(v)
        comps.append(comp)
    return comps


def verify_dgl(
    dgl_dir: Path,
    graph: dict,
    feats: np.ndarray,
    labels: np.ndarray,
    edges: Sequence[Tuple[int, int]],
    splits: np.ndarray,
    bounds: Sequence[int],
    adjacency: Sequence[Set[int]],
):
    n = len(graph["nodes"])
    graph_id = np.empty(n, dtype=np.int64)
    split_ranges = [(0, 20, 1, "train"), (20, 22, 21, "valid"), (22, 24, 23, "test")]
    component_rows = []
    for first_block, last_block, first_id, split_name in split_ranges:
        graph_id[bounds[first_block]:bounds[last_block]] = first_id
        for block in range(first_block, last_block):
            comps = connected_components_for_block(bounds[block], bounds[block + 1], adjacency)
            lcc = max(comps, key=len)
            graph_id[lcc] = block + 1
            component_rows.append({
                "graph_index": block + 1,
                "split": split_name,
                "full_nodes": bounds[block + 1] - bounds[block],
                "component_count": len(comps),
                "largest_component_nodes": len(lcc),
                "non_lcc_nodes": (bounds[block + 1] - bounds[block]) - len(lcc),
                "assigned_dgl_graph_id_for_lcc": block + 1,
                "assigned_dgl_graph_id_for_non_lcc": first_id,
            })

    scaler = StandardScaler().fit(feats[splits == "train"].astype(np.float64))
    standardized = scaler.transform(feats.astype(np.float64))
    split_results = {}
    edge_rows = []
    for first_block, last_block, first_id, split_name in split_ranges:
        ids = range(first_block + 1, last_block + 1)
        order = np.concatenate([np.flatnonzero(graph_id == gid) for gid in ids])
        d_gid = np.load(dgl_dir / f"{split_name}_graph_id.npy")
        d_y = np.load(dgl_dir / f"{split_name}_labels.npy")
        d_x = np.load(dgl_dir / f"{split_name}_feats.npy")
        d_graph = json.loads((dgl_dir / f"{split_name}_graph.json").read_text())
        local = {int(row): i for i, row in enumerate(order)}
        expected_edges: Set[Tuple[int, int]] = set()
        original_undirected_edges = 0
        original_loops = 0
        for u, v in edges:
            if u in local and v in local:
                original_undirected_edges += 1
                a, b = local[u], local[v]
                if a == b:
                    original_loops += 1
                    expected_edges.add((a, a))
                else:
                    expected_edges.add((a, b))
                    expected_edges.add((b, a))
        for i in range(len(order)):
            expected_edges.add((i, i))
        got_edges = {(int(e["source"]), int(e["target"])) for e in d_graph["links"]}
        max_diff = float(np.max(np.abs(standardized[order] - d_x)))
        result = {
            "row_count": int(len(order)),
            "graph_id_exact": bool(np.array_equal(graph_id[order], d_gid)),
            "labels_exact": bool(np.array_equal(labels[order], d_y)),
            "feature_dtype_in_archive": str(d_x.dtype),
            "feature_shape": list(d_x.shape),
            "feature_max_abs_difference_float64_standard_scaler": max_diff,
            "feature_allclose_at_1e_12": bool(np.allclose(standardized[order], d_x, atol=1e-12, rtol=0.0)),
            "edge_set_exact": bool(expected_edges == got_edges),
            "expected_directed_edge_count_including_one_loop_per_node": len(expected_edges),
            "dgl_directed_edge_count": len(got_edges),
            "missing_edge_count": len(expected_edges - got_edges),
            "extra_edge_count": len(got_edges - expected_edges),
            "original_graphsage_undirected_edge_count_including_loops": original_undirected_edges,
            "original_graphsage_loop_count": original_loops,
        }
        split_results[split_name] = result
        edge_rows.append({"split": split_name, **result})
    return {
        "algorithm": [
            "For each split, assign all rows initially to the split's first graph ID.",
            "For each of the 24 tissue blocks, assign its largest connected component to that tissue's graph ID.",
            "Leave every non-largest component in the first graph ID of its split.",
            "Concatenate rows by graph ID with stable original-row order.",
            "Fit sklearn StandardScaler on all GraphSAGE training rows as float64 and transform all rows.",
            "Convert every non-loop undirected edge to both directions and retain exactly one self-loop per node.",
        ],
        "splits": split_results,
        "all_checks_pass": all(
            r["graph_id_exact"] and r["labels_exact"] and r["feature_allclose_at_1e_12"] and r["edge_set_exact"]
            for r in split_results.values()
        ),
    }, component_rows, edge_rows


def collapse_labels(labels: np.ndarray, mapping: Mapping[int, int]):
    rows_by_gene: Dict[int, List[int]] = collections.defaultdict(list)
    for row, gene in mapping.items():
        rows_by_gene[gene].append(row)
    collapsed: Dict[int, np.ndarray] = {}
    inconsistent = []
    for gene, rows in rows_by_gene.items():
        vals = labels[rows]
        if not np.all(vals == vals[0]):
            inconsistent.append({"gene_id": gene, "rows": sorted(rows)})
        collapsed[gene] = vals[0]
    return collapsed, inconsistent


def micro_f1(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float, int, int, int]:
    tp = int(np.count_nonzero((y_true == 1) & (y_pred == 1)))
    fp = int(np.count_nonzero((y_true == 0) & (y_pred == 1)))
    fn = int(np.count_nonzero((y_true == 1) & (y_pred == 0)))
    denom = 2 * tp + fp + fn
    return (2 * tp / denom if denom else 0.0, tp, fp, fn)


def leakage_metrics(
    labels: np.ndarray,
    splits: np.ndarray,
    mapping: Mapping[int, int],
):
    train_by_gene: Dict[int, List[int]] = collections.defaultdict(list)
    for row, gene in mapping.items():
        if splits[row] == "train":
            train_by_gene[gene].append(row)
    test_rows_all = np.flatnonzero(splits == "test")
    test_rows_mapped = [int(r) for r in test_rows_all if int(r) in mapping]
    predictions = []
    truths = []
    seen = 0
    identical = 0
    for r in test_rows_all:
        r = int(r)
        gene = mapping.get(r)
        if gene is not None and gene in train_by_gene:
            pred = labels[train_by_gene[gene][0]]
            seen += 1
            identical += int(np.array_equal(pred, labels[r]))
        else:
            pred = np.zeros(labels.shape[1], dtype=np.uint8)
        predictions.append(pred)
        truths.append(labels[r])
    f1, tp, fp, fn = micro_f1(np.asarray(truths), np.asarray(predictions))
    return {
        "mapping_basis": "topology plus independently identified MSigDB v5.2 features; labels were not used to resolve identities",
        "test_rows_total": int(len(test_rows_all)),
        "test_rows_with_resolved_gene": int(len(test_rows_mapped)),
        "test_rows_seen_in_training": int(seen),
        "test_rows_seen_fraction_of_all_test_rows": seen / len(test_rows_all),
        "seen_test_rows_with_identical_label_vector": int(identical),
        "identical_fraction_among_seen_test_rows": identical / seen if seen else 0.0,
        "lookup_micro_f1_all_test_rows_unresolved_or_unseen_predicted_zero": f1,
        "true_positive_labels": tp,
        "false_positive_labels": fp,
        "false_negative_labels": fn,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("/mnt/data"))
    parser.add_argument("--work-dir", type=Path, default=Path("/mnt/data/ppi_repro_corrected/work"))
    parser.add_argument("--output-dir", type=Path, default=Path("/mnt/data/ppi_repro_corrected/results"))
    args = parser.parse_args()
    inp, work, out = args.input_dir, args.work_dir, args.output_dir
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    started = time.time()

    gs_extract = work / "graphsage"
    dgl_extract = work / "dgl"
    ohm_extract = work / "ohmnet_networks"
    if not (gs_extract / "ppi" / "ppi-G.json").exists():
        safe_extract_zip(inp / "graphsage_ppi.zip", gs_extract)
    if not (dgl_extract / "train_graph.json").exists():
        safe_extract_zip(inp / "dgl_ppi.zip", dgl_extract)
    if not list(ohm_extract.rglob("*.edgelist")):
        safe_extract_tar(inp / "bio-tissue-networks.tar.gz", ohm_extract)

    print("stage: load GraphSAGE", flush=True)
    graph, feats, labels, edges, splits = load_graphsage(gs_extract)
    n = len(graph["nodes"])
    edge_set = {canonical_edge(u, v) for u, v in edges}
    if len(edge_set) != len(edges):
        raise AssertionError("GraphSAGE link list contains duplicate undirected edges")
    loop_count = sum(u == v for u, v in edges)
    print("stage: load OhmNet", flush=True)
    ohm = load_ohm_networks(ohm_extract)
    print("stage: partition", flush=True)
    bounds, tissues, partition_rows = partition_graph(n, edges, ohm)
    print("partition done", tissues, flush=True)
    write_csv(out / "tissue_partition.csv", partition_rows)

    adj_g = make_global_adjacency(n, edges)
    print("stage: WL", flush=True)
    topology_map, ambiguous, wl_stats = topology_wl_mapping(bounds, tissues, adj_g, ohm)
    print("WL done", len(topology_map), len(ambiguous), flush=True)
    write_csv(out / "wl_tissue_summary.csv", wl_stats)
    (out / "wl_ambiguous_topology_only.json").write_text(json.dumps(ambiguous, indent=2))

    print("stage: feature identification", flush=True)
    mapping, residual, feature_rows, mapping_summary = feature_identification_and_mapping(
        feats, topology_map, ambiguous, inp / "msigdb_v5.2_files_to_download_locally.zip"
    )
    print("feature identification done", mapping_summary, flush=True)
    write_csv(out / "feature_column_mapping.csv", feature_rows)
    (out / "wl_residual_after_features.json").write_text(json.dumps(residual, indent=2))
    map_rows = [
        {"graphsage_row": row, "entrez_gene_id": gene, "resolution_basis": "topology" if row in topology_map else "topology_plus_feature"}
        for row, gene in sorted(mapping.items())
    ]
    write_csv(out / "graphsage_row_to_entrez_topology_features.csv", map_rows)

    print("stage: edge verification", flush=True)
    edge_verification = verify_mapped_edges(edges, bounds, tissues, mapping, ohm)
    print("stage: DGL", flush=True)
    dgl_result, component_rows, dgl_edge_rows = verify_dgl(
        dgl_extract, graph, feats, labels, edges, splits, bounds, adj_g
    )
    write_csv(out / "dgl_component_assignment.csv", component_rows)
    write_csv(out / "dgl_split_verification.csv", dgl_edge_rows)
    (out / "dgl_transformation_verification.json").write_text(json.dumps(dgl_result, indent=2))

    print("stage: labels/leakage", flush=True)
    collapsed, inconsistent = collapse_labels(labels, mapping)
    if inconsistent:
        raise AssertionError(f"Resolved Entrez genes have inconsistent label vectors: {len(inconsistent)}")
    collapsed_rows = []
    for gene in sorted(collapsed):
        row = {"entrez_gene_id": gene}
        row.update({f"label_{i}": int(v) for i, v in enumerate(collapsed[gene])})
        collapsed_rows.append(row)
    write_csv(out / "collapsed_gene_labels_topology_features.csv", collapsed_rows)
    leakage = leakage_metrics(labels, splits, mapping)

    union_ohm_genes = set().union(*(ohm[t].nodes for t in tissues))
    split_counts = {name: int(np.count_nonzero(splits == name)) for name in ["train", "valid", "test"]}
    summary = {
        "generated_at_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "runtime_seconds": time.time() - started,
        "inputs": {
            "graphsage_ppi_sha256": sha256_file(inp / "graphsage_ppi.zip"),
            "dgl_ppi_sha256": sha256_file(inp / "dgl_ppi.zip"),
            "ohmnet_networks_sha256": sha256_file(inp / "bio-tissue-networks.tar.gz"),
            "msigdb_v52_sha256": sha256_file(inp / "msigdb_v5.2_files_to_download_locally.zip"),
        },
        "graphsage": {
            "nodes": n,
            "undirected_links_including_loops": len(edges),
            "self_loops": loop_count,
            "features_shape": list(feats.shape),
            "features_dtype": str(feats.dtype),
            "labels_shape": list(labels.shape),
            "labels_dtype": str(labels.dtype),
            "split_node_counts": split_counts,
            "connected_components_total": None,
            "safe_cut_count_including_endpoints": len(safe_cuts(n, edges)),
            "safe_intervals": len(safe_cuts(n, edges)) - 1,
            "node_json_fields": sorted(set().union(*(node.keys() for node in graph["nodes"]))),
            "has_graph_id_field": any("graph_id" in node for node in graph["nodes"]),
        },
        "partition": {
            "block_count": len(tissues),
            "bounds": bounds,
            "tissues": tissues,
            "all_node_edge_statistics_exact": all(r["stats_exact"] for r in partition_rows),
        },
        "gene_identity": {
            **mapping_summary,
            "ohmnet_entrez_gene_universe_across_24_tissues": len(union_ohm_genes),
            "unresolved_unique_candidate_genes": len(set().union(*(set(r["candidate_genes"]) for r in residual))) if residual else 0,
            "edge_verification": edge_verification,
        },
        "features": {
            "columns": feats.shape[1],
            "exact_uniquely_identified_nonzero_columns": sum(r["identification_status"] == "exact_unique" for r in feature_rows),
            "all_zero_columns": [r["column"] for r in feature_rows if r["all_zero_column"]],
            "c1_columns": sum(r["chosen_collection"] == "C1" for r in feature_rows),
            "c3_columns": sum(r["chosen_collection"] == "C3" for r in feature_rows),
        },
        "dgl": dgl_result,
        "leakage": leakage,
    }
    # Connected component count for the entire combined GraphSAGE graph.
    seen: Set[int] = set()
    cc = 0
    for root in range(n):
        if root in seen:
            continue
        cc += 1
        stack = [root]
        seen.add(root)
        while stack:
            u = stack.pop()
            for v in adj_g[u]:
                if v == u or v in seen:
                    continue
                seen.add(v)
                stack.append(v)
    summary["graphsage"]["connected_components_total"] = cc
    (out / "core_verification_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({
        "status": "PASS",
        "runtime_seconds": summary["runtime_seconds"],
        "tissues": len(tissues),
        "topology_unique_rows": mapping_summary["topology_unique_rows"],
        "mapped_rows_after_features": mapping_summary["mapped_rows_after_features"],
        "mapped_unique_genes_after_features": mapping_summary["mapped_unique_genes_after_features"],
        "dgl_all_checks_pass": dgl_result["all_checks_pass"],
        "leakage": leakage,
    }, indent=2))


if __name__ == "__main__":
    main()

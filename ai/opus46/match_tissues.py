"""Attempt to recover gene identity by matching GraphSAGE's 24 graphs
to OhmNet tissue edgelists. Tests cheap hypotheses first."""
import json, glob, os, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
OH = "/home/claude/ohm/bio-tissue-networks/bio-tissue-networks/"

gid = np.load(D + "train_graph_id.npy")
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")
N = G.number_of_nodes()

# per-graph subgraphs of the GraphSAGE release (self-loops KEPT, to match OhmNet)
sage = {}
for g in np.unique(gid):
    nodes = np.where(gid == g)[0]
    sub = G.subgraph(nodes.tolist())
    sage[int(g)] = (sub.number_of_nodes(), sub.number_of_edges(),
                    nx.number_of_selfloops(sub), sub)

print("GraphSAGE per-graph (nodes, edges, selfloops):")
for g in sorted(sage):
    n, e, s, _ = sage[g]
    print(f"  g{g:<3d} n={n:5d} e={e:7d} selfloops={s:5d}")

# OhmNet
ohm = {}
for f in sorted(glob.glob(OH + "*.edgelist")):
    name = os.path.basename(f)[:-9]
    H = nx.read_edgelist(f, nodetype=int)
    ohm[name] = (H.number_of_nodes(), H.number_of_edges(),
                 nx.number_of_selfloops(H), H)

print(f"\nOhmNet: {len(ohm)} tissue networks")

print("\n--- exact (nodes, edges) signature match ---")
by_sig = collections.defaultdict(list)
for name, (n, e, s, _) in ohm.items():
    by_sig[(n, e)].append(name)

matches = {}
for g in sorted(sage):
    n, e, s, _ = sage[g]
    cand = by_sig.get((n, e), [])
    print(f"  g{g:<3d} ({n},{e}) -> {cand if cand else 'NO EXACT MATCH'}")
    if len(cand) == 1:
        matches[g] = cand[0]

print(f"\nunique exact matches: {len(matches)}/24")

if len(matches) < 24:
    print("\n--- fallback: nearest OhmNet net by (nodes, edges) ---")
    for g in sorted(sage):
        if g in matches:
            continue
        n, e, s, _ = sage[g]
        best = sorted(ohm.items(),
                      key=lambda kv: abs(kv[1][0] - n) + abs(kv[1][1] - e) / 100)[:3]
        print(f"  g{g:<3d} sage=({n},{e}) closest: " +
              ", ".join(f"{k}=({v[0]},{v[1]})" for k, v in best))

json.dump(matches, open("/home/claude/matches.json", "w"))

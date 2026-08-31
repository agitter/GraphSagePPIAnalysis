"""Given matched (sage graph, ohmnet tissue) pairs, try to recover the actual
node-index -> Entrez gene ID mapping.

Cheap hypothesis first: node indices within a graph block follow sorted Entrez ID
order. Verify by comparing per-node degree sequences IN ORDER (not sorted)."""
import json, glob, os, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
OH = "/home/claude/ohm/bio-tissue-networks/bio-tissue-networks/"
matches = json.load(open("/home/claude/matches.json"))

gid = np.load(D + "train_graph_id.npy")
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")


def pruned(name):
    H = nx.read_edgelist(OH + name + ".edgelist", nodetype=int)
    drop = [n for n in H if H.degree(n) - 2 * H.has_edge(n, n) == 0]
    H.remove_nodes_from(drop)
    return H


print("Testing: does sage node order == sorted(Entrez ID) order?\n")
recovered = {}
for g, name in sorted(matches.items(), key=lambda kv: int(kv[0])):
    g = int(g)
    idx = np.sort(np.where(gid == g)[0])
    sub = G.subgraph(idx.tolist())
    H = pruned(name)
    genes = sorted(H.nodes())
    dsage = np.array([sub.degree(i) for i in idx])
    dohm = np.array([H.degree(x) for x in genes])
    exact = np.array_equal(dsage, dohm)
    print(f"  g{g:<3d} {name:24s} in-order degree match: "
          f"{'YES' if exact else 'no'}  ({100*(dsage==dohm).mean():5.1f}% positionwise)")
    if exact:
        recovered[g] = dict(zip(idx.tolist(), genes))

print(f"\nrecovered by sorted-order hypothesis: {len(recovered)}/{len(matches)}")

# Stronger check on a recovered graph: does the full edge set match?
if recovered:
    g = sorted(recovered)[0]
    name = matches[str(g)]
    m = recovered[g]
    H = pruned(name)
    sub = G.subgraph(list(m.keys()))
    mapped = {frozenset((m[u], m[v])) for u, v in sub.edges()}
    truth = {frozenset((u, v)) for u, v in H.edges()}
    print(f"\nEDGE-SET verification on g{g} ({name}):")
    print(f"  sage edges mapped : {len(mapped)}")
    print(f"  ohmnet edges      : {len(truth)}")
    print(f"  identical         : {mapped == truth}")
    if mapped != truth:
        print(f"  jaccard           : {len(mapped&truth)/len(mapped|truth):.4f}")

json.dump({str(k): {str(a): b for a, b in v.items()}
           for k, v in recovered.items()},
          open("/home/claude/node2gene.json", "w"))

"""Hypothesis: sage_graph == ohmnet_graph minus nodes whose ONLY incident edge
is a self-loop (that removal deletes 1 node and 1 edge each -> dn == de).
Test it, then match exactly, then verify by degree-sequence isomorphism invariant."""
import json, glob, os, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
OH = "/home/claude/ohm/bio-tissue-networks/bio-tissue-networks/"

gid = np.load(D + "train_graph_id.npy")
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")

sage = {}
for g in np.unique(gid):
    sub = G.subgraph(np.where(gid == g)[0].tolist()).copy()
    sage[int(g)] = sub

ohm = {}
for f in sorted(glob.glob(OH + "*.edgelist")):
    name = os.path.basename(f)[:-9]
    H = nx.read_edgelist(f, nodetype=int)
    H2 = H.copy()
    # drop nodes with no non-self-loop edge
    drop = [n for n in H2 if H2.degree(n) - 2 * (H2.has_edge(n, n)) == 0]
    H2.remove_nodes_from(drop)
    ohm[name] = (H, H2, len(drop))

print("transform check (a few tissues): raw -> pruned")
for k in list(ohm)[:3]:
    H, H2, d = ohm[k]
    print(f"  {k:24s} ({H.number_of_nodes()},{H.number_of_edges()}) -> "
          f"({H2.number_of_nodes()},{H2.number_of_edges()}) dropped={d}")

sig = collections.defaultdict(list)
for k, (H, H2, d) in ohm.items():
    sig[(H2.number_of_nodes(), H2.number_of_edges())].append(k)

print("\n--- exact match after transform ---")
matches, ambiguous = {}, {}
for g in sorted(sage):
    s = sage[g]
    key = (s.number_of_nodes(), s.number_of_edges())
    cand = sig.get(key, [])
    if len(cand) == 1:
        matches[g] = cand[0]
    elif len(cand) > 1:
        ambiguous[g] = cand
    print(f"  g{g:<3d} {key} -> {cand if cand else 'NONE'}")

print(f"\nunique: {len(matches)}  ambiguous: {len(ambiguous)}  unmatched: "
      f"{24-len(matches)-len(ambiguous)}")

# verify with degree sequence (isomorphism invariant, cheap and strong)
print("\n--- degree-sequence verification ---")
ok = 0
for g, name in sorted(matches.items()):
    a = sorted(d for _, d in sage[g].degree())
    b = sorted(d for _, d in ohm[name][1].degree())
    good = a == b
    ok += good
    print(f"  g{g:<3d} <-> {name:28s} degree-seq {'MATCH' if good else 'MISMATCH'}")
print(f"degree-sequence confirmations: {ok}/{len(matches)}")

# resolve ambiguities by degree sequence
for g, cand in sorted(ambiguous.items()):
    a = sorted(d for _, d in sage[g].degree())
    hits = [c for c in cand if sorted(d for _, d in ohm[c][1].degree()) == a]
    print(f"  ambiguous g{g}: {cand} -> degree-seq compatible: {hits}")
    if len(hits) == 1:
        matches[g] = hits[0]

json.dump(matches, open("/home/claude/matches.json", "w"))
print(f"\nFINAL matched graphs: {len(matches)}/24")
for g in sorted(matches):
    print(f"  g{g:<3d} = {matches[g]}")

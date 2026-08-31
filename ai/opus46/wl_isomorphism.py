"""Recover node-index -> Entrez gene ID by colour-refinement (1-WL) canonical
labelling. Both sides are the same graph up to permutation, so WL colours
correspond. Nodes in singleton colour classes on BOTH sides map uniquely.
Verified afterwards by exact edge-set reconstruction."""
import json, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
OH = "/home/claude/ohm/bio-tissue-networks/bio-tissue-networks/"
matches = {int(k): v for k, v in json.load(open("/home/claude/matches.json")).items()}

gid = np.load(D + "train_graph_id.npy")
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")


def pruned(name):
    H = nx.read_edgelist(OH + name + ".edgelist", nodetype=int)
    H.remove_nodes_from([n for n in H if H.degree(n) - 2 * H.has_edge(n, n) == 0])
    return H


def wl_colors(H, rounds=6):
    col = {n: hash((H.degree(n), H.has_edge(n, n))) for n in H}
    for _ in range(rounds):
        col = {n: hash((col[n], tuple(sorted(col[m] for m in H.neighbors(n)))))
               for n in H}
    return col


mapping_all, stats = {}, []
for g in sorted(matches):
    name = matches[g]
    sub = G.subgraph(np.where(gid == g)[0].tolist())
    H = pruned(name)
    cs, ch = wl_colors(sub), wl_colors(H)
    gs, gh = collections.defaultdict(list), collections.defaultdict(list)
    for n, c in cs.items():
        gs[c].append(n)
    for n, c in ch.items():
        gh[c].append(n)
    mp = {gs[c][0]: gh[c][0] for c in gs
          if len(gs[c]) == 1 and len(gh.get(c, [])) == 1}
    # verify: every mapped edge must exist in H
    sube = [(u, v) for u, v in sub.edges() if u in mp and v in mp]
    good = sum(H.has_edge(mp[u], mp[v]) for u, v in sube)
    stats.append((g, name, sub.number_of_nodes(), len(mp), len(sube), good))
    mapping_all[g] = mp
    print(f"  g{g:<3d} {name:24s} n={sub.number_of_nodes():5d} "
          f"mapped={len(mp):5d} ({100*len(mp)/sub.number_of_nodes():5.1f}%)  "
          f"edge-check {good}/{len(sube)} "
          f"{'OK' if good == len(sube) else 'FAIL'}")

tot_n = sum(s[2] for s in stats); tot_m = sum(s[3] for s in stats)
tot_e = sum(s[4] for s in stats); tot_g = sum(s[5] for s in stats)
print(f"\nTOTAL mapped {tot_m}/{tot_n} = {100*tot_m/tot_n:.2f}% of nodes in "
      f"the {len(stats)} identified graphs")
print(f"edge verification: {tot_g}/{tot_e} = {100*tot_g/tot_e:.4f}% of mapped "
      f"edges reconstruct exactly")

json.dump({str(g): {str(k): v for k, v in m.items()} for g, m in mapping_all.items()},
          open("/home/claude/node2gene.json", "w"))

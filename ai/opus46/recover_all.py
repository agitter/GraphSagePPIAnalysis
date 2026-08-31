"""
Recover node-index -> Entrez gene ID for all 24 graphs.

Three tools, applied in order:
  1. WL colour refinement: maps nodes with unique structural fingerprints.
  2. Feature agreement: nodes sharing the same WL colour can be disambiguated
     if they have different 50-dim feature vectors (features are gene-level).
  3. Process of elimination: if only one candidate remains, assign it.
"""
import json, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
OH = "/home/claude/ohm/bio-tissue-networks/bio-tissue-networks/"
matches = {int(k): v for k, v in
           json.load(open("/home/claude/matches_full.json")).items()}

gid = np.load(D + "train_graph_id.npy")
feats = np.load(D + "ppi-feats.npy")
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")
N = feats.shape[0]


def get_lcc(H):
    lcc = max(nx.connected_components(H), key=len)
    return H.subgraph(lcc).copy()


def wl_colors(H, rounds=6):
    col = {n: hash((H.degree(n), H.has_edge(n, n))) for n in H}
    for _ in range(rounds):
        col = {n: hash((col[n], tuple(sorted(col[m] for m in H.neighbors(n)))))
               for n in H}
    return col


# Feature fingerprint for each node
fp = {i: feats[i].astype(np.uint8).tobytes() for i in range(N)}

total_nodes = 0
total_mapped_wl = 0
total_mapped_feat = 0
total_mapped_elim = 0
total_unmapped = 0
mapping_all = {}

for g in sorted(matches):
    name = matches[g]
    nodes = sorted(np.where(gid == g)[0])
    sub = G.subgraph(nodes)
    
    H = nx.read_edgelist(OH + name + ".edgelist", nodetype=int)
    
    # Work on LCC of both
    lcc_s = get_lcc(sub)
    lcc_o = get_lcc(H)
    
    sage_nodes = sorted(lcc_s.nodes())
    ohm_nodes = sorted(lcc_o.nodes())
    
    cs = wl_colors(lcc_s)
    ch = wl_colors(lcc_o)
    
    # Group by WL colour
    gs = collections.defaultdict(list)
    gh = collections.defaultdict(list)
    for n, c in cs.items():
        gs[c].append(n)
    for n, c in ch.items():
        gh[c].append(n)
    
    mp = {}
    n_wl = 0
    n_feat = 0
    n_elim = 0
    
    # Pass 1: WL singletons
    for c in gs:
        if len(gs[c]) == 1 and len(gh.get(c, [])) == 1:
            mp[gs[c][0]] = gh[c][0]
            n_wl += 1
    
    # Pass 2: within same WL colour class, disambiguate by feature vector
    for c in gs:
        if c not in gh:
            continue
        s_unmapped = [n for n in gs[c] if n not in mp]
        o_unmapped = [n for n in gh[c] if n not in mp.values()]
        if not s_unmapped or not o_unmapped:
            continue
        
        # Group by feature fingerprint within each side
        s_by_feat = collections.defaultdict(list)
        o_by_feat = collections.defaultdict(list)
        for n in s_unmapped:
            s_by_feat[fp[n]].append(n)
        for n in o_unmapped:
            # OhmNet nodes don't have features, but we can use the already-mapped
            # nodes to build a gene->feature table, then look up each gene
            pass  # Need different approach
        
        # Actually: features are on the GraphSAGE side only. OhmNet has gene IDs
        # but no features. So I need to match WITHIN the WL colour class using
        # the constraint that two sage nodes with different features must map
        # to different genes.
        #
        # Better approach: if |s_unmapped| == |o_unmapped| == 2 and the two sage
        # nodes have different features, I need another signal. But I can use
        # features across GRAPHS: the same gene has the same feature vector
        # everywhere. So if gene X is already mapped in another graph, I know
        # its feature vector and can match it here.
    
    # Pass 2 revised: use cross-graph feature knowledge
    # Build gene->feature from already mapped nodes across ALL graphs so far
    gene2feat = {}
    for gg, mm in mapping_all.items():
        for node, gene in mm.items():
            gene2feat[gene] = fp[node]
    # Also from this graph's already-mapped nodes
    for node, gene in mp.items():
        gene2feat[gene] = fp[node]
    
    for c in gs:
        if c not in gh:
            continue
        s_unmapped = [n for n in gs[c] if n not in mp]
        o_unmapped = [n for n in gh[c] if n not in mp.values()]
        if not s_unmapped or not o_unmapped:
            continue
        
        # For each unmapped OhmNet gene, do we know its feature from another graph?
        o_known = {n: gene2feat[n] for n in o_unmapped if n in gene2feat}
        
        for og, ofeat in o_known.items():
            candidates = [sn for sn in s_unmapped if fp[sn] == ofeat and sn not in mp]
            if len(candidates) == 1:
                mp[candidates[0]] = og
                s_unmapped.remove(candidates[0])
                n_feat += 1
    
    # Pass 3: elimination — if only one candidate remains in a colour class
    for c in gs:
        if c not in gh:
            continue
        s_unmapped = [n for n in gs[c] if n not in mp]
        o_unmapped = [n for n in gh[c] if n not in mp.values()]
        if len(s_unmapped) == 1 and len(o_unmapped) == 1:
            mp[s_unmapped[0]] = o_unmapped[0]
            n_elim += 1
    
    # Also map non-LCC nodes if they form matching small components
    # (pairs connected by an edge with matching degrees)
    non_lcc_s = set(nodes) - set(lcc_s.nodes())
    non_lcc_o = set(H.nodes()) - set(lcc_o.nodes())
    # Small components are typically isolated nodes or pairs
    # Skip for now — LCC covers the vast majority
    
    n_unmapped = len(sage_nodes) - len(mp)
    total_nodes += len(sage_nodes)
    total_mapped_wl += n_wl
    total_mapped_feat += n_feat
    total_mapped_elim += n_elim
    total_unmapped += n_unmapped
    mapping_all[g] = mp
    
    split = "TEST " if g >= 23 else ("VAL  " if g >= 21 else "train")
    print("g%-3d %-24s %s LCC=%d  WL=%d feat=%d elim=%d total=%d unmapped=%d" % (
        g, name, split, len(sage_nodes), n_wl, n_feat, n_elim,
        len(mp), n_unmapped))

print("\n" + "=" * 70)
print("TOTALS (LCC nodes only)")
print("  WL singleton:    %6d" % total_mapped_wl)
print("  Feature disamb:  %6d" % total_mapped_feat)
print("  Elimination:     %6d" % total_mapped_elim)
print("  TOTAL MAPPED:    %6d / %d = %.2f%%" % (
    total_mapped_wl + total_mapped_feat + total_mapped_elim,
    total_nodes,
    100 * (total_mapped_wl + total_mapped_feat + total_mapped_elim) / total_nodes))
print("  Unmapped:        %6d" % total_unmapped)

# Verify edge reconstruction
print("\n--- Edge verification ---")
total_e = 0; total_ok = 0
for g, mp in mapping_all.items():
    name = matches[g]
    H = nx.read_edgelist(OH + name + ".edgelist", nodetype=int)
    sub = G.subgraph(list(mp.keys()))
    edges = [(u, v) for u, v in sub.edges() if u in mp and v in mp]
    ok = sum(1 for u, v in edges if H.has_edge(mp[u], mp[v]))
    total_e += len(edges); total_ok += ok
    if ok != len(edges):
        print("  g%d: %d/%d edges OK" % (g, ok, len(edges)))

print("Total: %d/%d edges verified (%.4f%%)" % (total_ok, total_e, 100*total_ok/total_e))

# Save
out = {str(g): {str(k): v for k, v in mp.items()} for g, mp in mapping_all.items()}
json.dump(out, open("/home/claude/node2gene_full.json", "w"))

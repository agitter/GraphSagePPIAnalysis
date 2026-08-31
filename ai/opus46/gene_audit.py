"""With real Entrez gene IDs recovered for 16 of 24 graphs, run the audit that
the feature-fingerprint method could not do."""
import json, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
n2g = {int(g): {int(k): v for k, v in m.items()}
       for g, m in json.load(open("/home/claude/node2gene.json")).items()}
matches = {int(k): v for k, v in json.load(open("/home/claude/matches.json")).items()}

feats = np.load(D + "ppi-feats.npy")
gid = np.load(D + "train_graph_id.npy")
cm = json.load(open(D + "ppi-class_map.json"))
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")
N = feats.shape[0]
labels = np.zeros((N, 121), np.int8)
for k, v in cm.items():
    labels[int(k)] = v
is_val = np.zeros(N, bool); is_test = np.zeros(N, bool)
for n, d in G.nodes(data=True):
    if d.get("val"): is_val[n] = True
    if d.get("test"): is_test[n] = True
is_train = ~(is_val | is_test)

split_of = {g: ("test" if is_test[np.where(gid == g)[0][0]] else
                "val" if is_val[np.where(gid == g)[0][0]] else "train")
            for g in np.unique(gid)}

id_tr = [g for g in n2g if split_of[g] == "train"]
id_te = [g for g in n2g if split_of[g] == "test"]
print(f"identified graphs: {len(id_tr)} train ({sorted(matches[g] for g in id_tr)})")
print(f"                   {len(id_te)} test  ({[matches[g] for g in id_te]})")
print(f"NOTE: only {len(id_tr)}/20 train graphs identified -> all overlap figures "
      f"below are LOWER BOUNDS.\n")

train_genes = set()
for g in id_tr:
    train_genes |= set(n2g[g].values())
print(f"distinct genes across the {len(id_tr)} identified TRAIN graphs: {len(train_genes)}")

print("\n" + "=" * 70)
print("GENE-LEVEL TRAIN/TEST OVERLAP (the real audit)")
print("=" * 70)
for g in id_te:
    genes = set(n2g[g].values())
    ov = genes & train_genes
    print(f"  {matches[g]} (g{g}, test): {len(genes)} mapped genes; "
          f"{len(ov)} also in identified train graphs = {100*len(ov)/len(genes):.2f}%")

# ---- do the same gene's FEATURES agree across graphs? (validates gene-level claim)
print("\n" + "=" * 70)
print("IS THE FEATURE VECTOR A GENE-LEVEL PROPERTY?")
print("=" * 70)
gene2feat = collections.defaultdict(set)
gene2lab = collections.defaultdict(set)
gene2nodes = collections.defaultdict(list)
for g, m in n2g.items():
    for node, gene in m.items():
        gene2feat[gene].add(feats[node].astype(np.uint8).tobytes())
        gene2lab[gene].add(labels[node].tobytes())
        gene2nodes[gene].append((g, node))

multi = {k: v for k, v in gene2nodes.items() if len(v) > 1}
print(f"genes appearing in >1 identified graph: {len(multi)} "
      f"(of {len(gene2nodes)} total)")
print(f"  mean graphs per such gene: "
      f"{np.mean([len(v) for v in multi.values()]):.2f}")
same_f = sum(1 for k in multi if len(gene2feat[k]) == 1)
print(f"  identical FEATURE vector in every graph: {same_f}/{len(multi)} "
      f"= {100*same_f/len(multi):.2f}%")

print("\n" + "=" * 70)
print("ARE THE 121 LABELS TISSUE-SPECIFIC OR TISSUE-AGNOSTIC?")
print("=" * 70)
same_l = sum(1 for k in multi if len(gene2lab[k]) == 1)
print(f"  identical LABEL vector in every graph: {same_l}/{len(multi)} "
      f"= {100*same_l/len(multi):.2f}%")

# per-bit agreement for the same gene across graphs
ag = []
for k, v in multi.items():
    L = np.stack([labels[n] for _, n in v])
    ag.append((L == L[0]).mean())
print(f"  mean per-bit agreement across occurrences: {100*np.mean(ag):.2f}%")

# ---- the honest test-set question
print("\n" + "=" * 70)
print("TEST-GRAPH GENES: SEEN IN TRAINING WITH THE SAME ANSWER?")
print("=" * 70)
for g in id_te:
    same, diff, unseen = 0, 0, 0
    for node, gene in n2g[g].items():
        occ = [(gg, nn) for gg, nn in gene2nodes[gene] if gg in id_tr]
        if not occ:
            unseen += 1
        elif any(labels[nn].tobytes() == labels[node].tobytes() for _, nn in occ):
            same += 1
        else:
            diff += 1
    tot = same + diff + unseen
    print(f"  {matches[g]} (n={tot} mapped test genes)")
    print(f"    gene seen in train WITH IDENTICAL 121-label vector : "
          f"{same:5d} ({100*same/tot:.2f}%)")
    print(f"    gene seen in train, different label vector         : "
          f"{diff:5d} ({100*diff/tot:.2f}%)")
    print(f"    gene not in any identified train graph             : "
          f"{unseen:5d} ({100*unseen/tot:.2f}%)")

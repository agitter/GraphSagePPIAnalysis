"""Does an exact WL-1 structural match in TRAIN carry the test node's answer?
Pure lookup: no learning, no parameters, no gradient."""
import json, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
feats = np.load(D + "ppi-feats.npy")
gid = np.load(D + "train_graph_id.npy")
cm = json.load(open(D + "ppi-class_map.json"))
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")
N = feats.shape[0]; L = 121
labels = np.zeros((N, L), np.int8)
for k, v in cm.items():
    labels[int(k)] = v
is_val = np.zeros(N, bool); is_test = np.zeros(N, bool)
for n, d in G.nodes(data=True):
    if d.get("val"): is_val[n] = True
    if d.get("test"): is_test[n] = True
is_train = ~(is_val | is_test)

A = nx.to_scipy_sparse_array(G, nodelist=range(N), format="csr")
A.setdiag(0); A.eliminate_zeros()
fp = np.array([b.tobytes() for b in np.packbits(feats.astype(np.uint8), 1)], object)

ip, ind = A.indptr, A.indices
wl1 = np.array([hash((fp[i], tuple(sorted(fp[j] for j in ind[ip[i]:ip[i+1]]))))
                for i in range(N)], object)
# WL-2: refine using neighbours' WL-1 colours
wl2 = np.array([hash((wl1[i], tuple(sorted(wl1[j] for j in ind[ip[i]:ip[i+1]]))))
                for i in range(N)], object)


def micro_f1(y, p):
    tp = int((y & p).sum()); fp_ = int((~y & p).sum()); fn = int((y & ~p).sum())
    pr = tp / (tp + fp_) if tp + fp_ else 0.0
    rc = tp / (tp + fn) if tp + fn else 0.0
    return 2 * pr * rc / (pr + rc) if pr + rc else 0.0


yte = labels[is_test].astype(bool)
freq = labels[is_train].mean(0) > 0.5
test_idx = np.where(is_test)[0]

for name, col in [("WL-1", wl1), ("WL-2", wl2)]:
    # collision null
    dup = tot = 0
    for g in np.unique(gid):
        m = gid == g
        c = collections.Counter(col[m].tolist())
        dup += sum(v - 1 for v in c.values() if v > 1); tot += m.sum()
    lut = collections.defaultdict(list)
    for i in np.where(is_train)[0]:
        lut[col[i]].append(i)
    pred = np.zeros_like(yte); hits = 0
    for r, i in enumerate(test_idx):
        if col[i] in lut:
            pred[r] = labels[lut[col[i]]].mean(0) > 0.5; hits += 1
        else:
            pred[r] = freq
    res = np.array([col[i] in lut for i in test_idx])
    print(f"{name}: within-graph collision null = {100*dup/tot:5.2f}%  |  "
          f"test nodes matched in train = {hits}/{len(test_idx)} "
          f"({100*hits/len(test_idx):.2f}%)")
    print(f"      lookup micro-F1 (all test)      = {micro_f1(yte, pred):.4f}")
    print(f"      lookup micro-F1 (matched only)  = {micro_f1(yte[res], pred[res]):.4f}")
    print(f"      constant baseline (matched only)= "
          f"{micro_f1(yte[res], np.tile(freq,(res.sum(),1))):.4f}\n")

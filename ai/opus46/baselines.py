"""
Gene-identity baselines on the PPI benchmark.

We have Entrez gene IDs for 2,258 of 2,300 nodes in the midbrain test graph.
Every one of those genes also appears in the 15 identified training graphs,
carrying an identical 121-label vector.

Baselines:
  B1. Gene lookup: copy the label vector from training.
  B2. Gene-majority vote: if a gene appears in multiple training graphs,
      take the majority vote per label. (Moot here since labels are identical,
      but tests the mechanism.)
  B3. Nearest-gene by feature Hamming distance: for unmapped nodes, find
      the training node with the most similar 50-dim binary feature vector.
  B4. Nearest-gene by label frequency: predict the k most common labels.
  B5. Combined: gene lookup where available, nearest-feature fallback.

All results on the FULL test set (both test graphs, 5,524 nodes).
"""
import json, collections, numpy as np
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
n2g = {int(g): {int(k): v for k, v in m.items()}
       for g, m in json.load(open("/home/claude/node2gene.json")).items()}

feats = np.load(D + "ppi-feats.npy")
gid = np.load(D + "train_graph_id.npy")
cm = json.load(open(D + "ppi-class_map.json"))
G_data = json.load(open(D + "ppi-G.json"))
G = json_graph.node_link_graph(G_data, edges="links")

N = feats.shape[0]; L = 121
labels = np.zeros((N, L), np.int8)
for k, v in cm.items(): labels[int(k)] = v

is_val = np.zeros(N, bool); is_test = np.zeros(N, bool)
for n, d in G.nodes(data=True):
    if d.get("val"): is_val[n] = True
    if d.get("test"): is_test[n] = True
is_train = ~(is_val | is_test)

train_idx = np.where(is_train)[0]
test_idx = np.where(is_test)[0]
yte = labels[is_test].astype(bool)


def micro_f1(y, p):
    tp = (y & p).sum(); fp = (~y & p).sum(); fn = (y & ~p).sum()
    pr = tp / (tp + fp) if tp + fp else 0.0
    rc = tp / (tp + fn) if tp + fn else 0.0
    return 2 * pr * rc / (pr + rc) if pr + rc else 0.0


# Build gene -> label table from identified training graphs
gene2lab = {}
for g, m in n2g.items():
    if not is_train[np.where(gid == g)[0][0]]:
        continue
    for node, gene in m.items():
        gene2lab[gene] = labels[node]

# Build gene -> node index for test
test_gene = {}  # test row index -> gene
for g, m in n2g.items():
    if not is_test[np.where(gid == g)[0][0]]:
        continue
    for node, gene in m.items():
        r = np.searchsorted(test_idx, node)
        if r < len(test_idx) and test_idx[r] == node:
            test_gene[r] = gene

# Feature vectors of training nodes for nearest-feature fallback
train_feats = feats[is_train].astype(np.uint8)
train_labels = labels[is_train]
freq = labels[is_train].mean(0)

print("=" * 70)
print("BASELINES ON FULL TEST SET (n=%d)" % len(test_idx))
print("=" * 70)

# B0: constant baselines
const_med = (freq > 0.5)
const_all = np.ones(L, bool)
print(f"\nB0a. Predict all ones              : micro-F1 = {micro_f1(yte, np.tile(const_all, (len(test_idx),1))):.4f}")
print(f"B0b. Predict train freq > 0.5      : micro-F1 = {micro_f1(yte, np.tile(const_med, (len(test_idx),1))):.4f}")

# B1: Gene lookup (gene identity where known, constant fallback)
pred_b1 = np.tile(const_med, (len(test_idx), 1))
b1_hits = 0
for r in range(len(test_idx)):
    if r in test_gene and test_gene[r] in gene2lab:
        pred_b1[r] = gene2lab[test_gene[r]] > 0
        b1_hits += 1

print(f"\nB1.  Gene lookup + constant fallback")
print(f"     coverage: {b1_hits}/{len(test_idx)} = {100*b1_hits/len(test_idx):.1f}%")
print(f"     micro-F1 (full test set)      : {micro_f1(yte, pred_b1):.4f}")
covered = np.array([r in test_gene and test_gene[r] in gene2lab for r in range(len(test_idx))])
print(f"     micro-F1 (covered only)       : {micro_f1(yte[covered], pred_b1[covered]):.4f}")

# B3: Nearest training node by Hamming distance on 50-dim binary features
# For EVERY test node (not just unmapped ones)
test_feats = feats[is_test].astype(np.uint8)
print(f"\nB3.  Nearest-feature lookup (Hamming, all test nodes)")
# This is O(n_test * n_train * 50) -- tractable at these sizes
pred_b3 = np.zeros_like(yte)
dists = []
BATCH = 500
for start in range(0, len(test_idx), BATCH):
    end = min(start + BATCH, len(test_idx))
    chunk = test_feats[start:end]
    # Hamming via XOR
    d = np.unpackbits(
        np.bitwise_xor(
            np.packbits(chunk, axis=1)[:, np.newaxis, :],
            np.packbits(train_feats, axis=1)[np.newaxis, :, :]
        ), axis=2, count=50
    ).sum(axis=2)
    nn = d.argmin(axis=1)
    pred_b3[start:end] = train_labels[nn] > 0
    dists.extend(d.min(axis=1).tolist())

dists = np.array(dists)
print(f"     Hamming distances: mean={dists.mean():.1f} median={np.median(dists):.0f} "
      f"zero={int((dists==0).sum())} ({100*(dists==0).mean():.1f}%)")
print(f"     micro-F1 (full test set)      : {micro_f1(yte, pred_b3):.4f}")

# B4: combine gene lookup where available, nearest-feature elsewhere
pred_b4 = pred_b3.copy()
for r in range(len(test_idx)):
    if r in test_gene and test_gene[r] in gene2lab:
        pred_b4[r] = gene2lab[test_gene[r]] > 0

print(f"\nB4.  Gene lookup + nearest-feature fallback")
print(f"     gene coverage: {b1_hits}/{len(test_idx)} = {100*b1_hits/len(test_idx):.1f}%")
print(f"     micro-F1 (full test set)      : {micro_f1(yte, pred_b4):.4f}")

# B5: how good is nearest-feature on the midbrain graph alone?
# (where we know ground truth gene identity)
mid_mask = gid[is_test] == 24
print(f"\nB5.  Nearest-feature on midbrain only (n={mid_mask.sum()})")
print(f"     micro-F1                      : {micro_f1(yte[mid_mask], pred_b3[mid_mask]):.4f}")
print(f"     gene lookup micro-F1          : {micro_f1(yte[mid_mask], pred_b1[mid_mask]):.4f}")

# Summary table
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"{'Method':<45s} {'Full test':>10s} {'Midbrain':>10s}")
print("-" * 70)
for name, p in [
    ("B0a. Predict all ones", np.tile(const_all, (len(test_idx), 1))),
    ("B0b. Predict freq > 0.5", np.tile(const_med, (len(test_idx), 1))),
    ("B1.  Gene lookup + constant fallback", pred_b1),
    ("B3.  Nearest-feature lookup", pred_b3),
    ("B4.  Gene lookup + nearest-feature", pred_b4),
]:
    f_full = micro_f1(yte, p)
    f_mid = micro_f1(yte[mid_mask], p[mid_mask])
    print(f"{name:<45s} {f_full:>10.4f} {f_mid:>10.4f}")
print("-" * 70)
print("Published GAT (Velickovic et al. 2018)        0.9730")
print("Published Cluster-GCN (Chiang et al. 2019)    0.9936")

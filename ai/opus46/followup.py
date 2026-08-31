"""
Follow-ups after the fingerprint result came back negative.

(1) DECISIVE CONTROL: is cross-graph label agreement for same-fingerprint nodes
    any higher than WITHIN-graph agreement for same-fingerprint nodes? If not,
    the agreement is explained purely by feature-conditioning, not gene identity.
(2) Do structural signatures resolve nodes where features cannot?
(3) Where does the 0.99 actually come from -- measure the jump from raw features
    to neighbourhood-aggregated features.
(4) Label homophily.
"""
import json, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph
import scipy.sparse as sp

D = "/home/claude/dl/PPI-Inductive/ppi/"
rng = np.random.default_rng(0)

feats = np.load(D + "ppi-feats.npy")
gid = np.load(D + "train_graph_id.npy")
class_map = json.load(open(D + "ppi-class_map.json"))
G = json_graph.node_link_graph(json.load(open(D + "ppi-G.json")), edges="links")
N, F = feats.shape
L = 121
labels = np.zeros((N, L), np.int8)
for k, v in class_map.items():
    labels[int(k)] = v

is_val = np.zeros(N, bool); is_test = np.zeros(N, bool)
for n, d in G.nodes(data=True):
    if d.get("val"): is_val[n] = True
    if d.get("test"): is_test[n] = True
is_train = ~(is_val | is_test)

zero = (feats == 0).all(1)
nz = ~zero
fp = np.array([b.tobytes() for b in np.packbits(feats.astype(np.uint8), 1)], object)


def hdr(s): print("\n" + "=" * 74 + "\n" + s + "\n" + "=" * 74)


# ---------------------------------------------------------------------------
hdr("H. DECISIVE CONTROL: within-graph vs cross-graph label agreement")
print("If these match, same-fingerprint label agreement reflects feature")
print("conditioning, NOT recurrence of the same protein across tissues.\n")

within, cross = [], []
by_fp = collections.defaultdict(list)
for i in np.where(nz)[0]:
    by_fp[fp[i]].append(i)

for k, idx in by_fp.items():
    if len(idx) < 2:
        continue
    idx = np.array(idx)
    g = gid[idx]
    # sample pairs to keep it tractable
    for _ in range(min(200, len(idx))):
        a, b = rng.choice(len(idx), 2, replace=False)
        agree = (labels[idx[a]] == labels[idx[b]]).mean()
        (within if g[a] == g[b] else cross).append(agree)

print(f"same fingerprint, SAME graph  : n={len(within):6d}  "
      f"per-bit agreement = {100*np.mean(within):.2f}%")
print(f"same fingerprint, DIFF graph  : n={len(cross):6d}  "
      f"per-bit agreement = {100*np.mean(cross):.2f}%")
a = rng.choice(np.where(nz)[0], 50000); b = rng.choice(np.where(nz)[0], 50000)
print(f"random pairs (control)        : n={len(a):6d}  "
      f"per-bit agreement = {100*(labels[a]==labels[b]).mean():.2f}%")
print("\n=> difference cross-minus-within = "
      f"{100*(np.mean(cross)-np.mean(within)):+.2f} percentage points")

# ---------------------------------------------------------------------------
hdr("I. STRUCTURAL SIGNATURES: can topology identify nodes where features can't?")
A = nx.to_scipy_sparse_array(G, nodelist=range(N), format="csr")
A.setdiag(0); A.eliminate_zeros()
deg = np.asarray(A.sum(1)).ravel().astype(int)
print(f"degree: mean={deg.mean():.1f} median={np.median(deg):.0f} max={deg.max()}")
for g in [1, 23, 24]:
    m = gid == g
    print(f"  graph {g}: n={m.sum()} mean deg={deg[m].mean():.1f}")

# WL-1 colour = (own fingerprint, sorted multiset of neighbour fingerprints)
def wl1():
    out = np.empty(N, object)
    idxptr, ind = A.indptr, A.indices
    for i in range(N):
        nb = ind[idxptr[i]:idxptr[i + 1]]
        out[i] = hash((fp[i], tuple(sorted(fp[j] for j in nb))))
    return out

wl = wl1()
print(f"\ndistinct WL-1 colours: {len(set(wl.tolist()))} (vs {len(set(fp.tolist()))} "
      f"raw fingerprints, {N} nodes)")
dup = 0; tot = 0
for g in np.unique(gid):
    m = gid == g
    c = collections.Counter(wl[m].tolist())
    dup += sum(v - 1 for v in c.values() if v > 1); tot += m.sum()
print(f"WL-1 within-graph duplicate rate (collision null): "
      f"{dup}/{tot} = {100*dup/tot:.2f}%")
tr = set(wl[is_train].tolist())
hit = np.array([w in tr for w in wl[is_test]])
print(f"WL-1 colours of TEST nodes also seen in TRAIN: "
      f"{hit.sum()}/{is_test.sum()} = {100*hit.mean():.2f}%")

# ---------------------------------------------------------------------------
hdr("J. LABEL HOMOPHILY")
src, dst = A.nonzero()
same = (labels[src] == labels[dst]).mean()
print(f"per-bit label agreement across edges : {100*same:.2f}%")
a = rng.choice(N, 200000); b = rng.choice(N, 200000)
print(f"per-bit agreement, random node pairs : {100*(labels[a]==labels[b]).mean():.2f}%")
jac_e = ((labels[src] & labels[dst]).sum(1) /
         np.maximum((labels[src] | labels[dst]).sum(1), 1)).mean()
jac_r = ((labels[a] & labels[b]).sum(1) /
         np.maximum((labels[a] | labels[b]).sum(1), 1)).mean()
print(f"mean label Jaccard across edges      : {jac_e:.4f}")
print(f"mean label Jaccard, random pairs     : {jac_r:.4f}")

# ---------------------------------------------------------------------------
hdr("K. WHERE DOES THE SIGNAL LIVE? raw features vs neighbourhood aggregation")
from sklearn.linear_model import LogisticRegression
from sklearn.multioutput import MultiOutputClassifier


def micro_f1(y, p):
    tp = int((y & p).sum()); fp_ = int((~y & p).sum()); fn = int((y & ~p).sum())
    pr = tp / (tp + fp_) if tp + fp_ else 0.0
    rc = tp / (tp + fn) if tp + fn else 0.0
    return 2 * pr * rc / (pr + rc) if pr + rc else 0.0


Dg = sp.diags(1.0 / np.maximum(deg, 1))
P = Dg @ A                      # row-normalised propagation
X1 = P @ feats                  # 1-hop mean
X2 = P @ X1                     # 2-hop mean
reps = {"raw feats (50d)": feats,
        "feats + 1hop (100d)": np.hstack([feats, X1]),
        "feats + 1hop + 2hop (150d)": np.hstack([feats, X1, X2])}

ytr = labels[is_train].astype(bool); yte = labels[is_test].astype(bool)
for name, X in reps.items():
    Xtr, Xte = X[is_train], X[is_test]
    pred = np.zeros_like(yte)
    for c in range(L):
        col = ytr[:, c]
        if col.all() or not col.any():
            pred[:, c] = col.all()
            continue
        clf = LogisticRegression(max_iter=200, n_jobs=-1)
        clf.fit(Xtr, col)
        pred[:, c] = clf.predict(Xte)
    print(f"{name:30s} logistic-regression micro-F1 = {micro_f1(yte, pred):.4f}")

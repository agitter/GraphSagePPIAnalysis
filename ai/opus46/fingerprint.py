"""
PPI (GraphSAGE 24-tissue) fingerprint / leakage audit.

Core idea: the 50 features are MSigDB gene-set membership indicators, a GENE-level
property. If so, the same protein carries a bit-identical feature vector in every
tissue graph it appears in, so the feature vector acts as a protein fingerprint --
no gene-ID mapping needed.

Calibration: a gene appears at most once per tissue graph (each graph is an induced
subgraph of one global network). Therefore ANY repeated fingerprint *within* a
single graph is a false positive by construction. That gives an empirical null.

Assumptions are stated and tested rather than assumed.
"""
import json, collections, numpy as np, networkx as nx
from networkx.readwrite import json_graph

D = "/home/claude/dl/PPI-Inductive/ppi/"
rng = np.random.default_rng(0)


def hdr(s):
    print("\n" + "=" * 74 + "\n" + s + "\n" + "=" * 74)


# ---------------------------------------------------------------- load
feats = np.load(D + "ppi-feats.npy")
gid = np.load(D + "train_graph_id.npy")
class_map = json.load(open(D + "ppi-class_map.json"))
G_data = json.load(open(D + "ppi-G.json"))
G = json_graph.node_link_graph(G_data, edges="links")

N, F = feats.shape
# labels in node-index order (id_map verified to be identity)
L = len(next(iter(class_map.values())))
labels = np.zeros((N, L), dtype=np.int8)
for k, v in class_map.items():
    labels[int(k)] = v
assert set(np.unique(labels).tolist()) <= {0, 1}, "labels not binary"

is_val = np.zeros(N, bool)
is_test = np.zeros(N, bool)
for n, d in G.nodes(data=True):
    if d.get("val"):
        is_val[n] = True
    if d.get("test"):
        is_test[n] = True
is_train = ~(is_val | is_test)
split = np.where(is_test, "test", np.where(is_val, "val", "train"))

hdr("A. SPLIT / GRAPH STRUCTURE")
print(f"nodes={N} feats={F} labels={L} graphs={len(np.unique(gid))}")
print(f"train={is_train.sum()} val={is_val.sum()} test={is_test.sum()}")

# Is each tissue graph wholly within one split? (tests the 20/2/2 claim)
g2s = collections.defaultdict(set)
for g, s in zip(gid, split):
    g2s[int(g)].add(s)
mixed = {g: s for g, s in g2s.items() if len(s) > 1}
print("graphs spanning >1 split:", mixed if mixed else "none (clean 20/2/2)")
cnt = collections.Counter(next(iter(g2s[g])) for g in g2s)
print("graphs per split:", dict(cnt))

# ---------------------------------------------------------------- features
hdr("B. FEATURE MATRIX")
zero = (feats == 0).all(1)
print(f"all-zero feature rows: {zero.sum()} ({100*zero.mean():.2f}%)")
for s in ["train", "val", "test"]:
    m = split == s
    print(f"  {s:5s}: {zero[m].sum():6d}/{m.sum():6d} = {100*zero[m].mean():.2f}% zero")
print("\nper-graph zero-feature share (graph, split, n, %zero):")
for g in sorted(np.unique(gid)):
    m = gid == g
    print(f"  g{g:<3d} {next(iter(g2s[int(g)])):5s} n={m.sum():5d} "
          f"zero={100*zero[m].mean():6.2f}%")

# ---------------------------------------------------------------- fingerprints
bits = np.packbits(feats.astype(np.uint8), axis=1)
fp = np.array([b.tobytes() for b in bits], dtype=object)

hdr("C. FINGERPRINT RESOLVING POWER")
nz = ~zero
print(f"distinct fingerprints, all nodes      : {len(set(fp.tolist()))}")
print(f"distinct fingerprints, non-zero nodes : {len(set(fp[nz].tolist()))} "
      f"(of {nz.sum()} nodes)")

# Within-graph duplicates == collisions (null). Restricted to non-zero rows.
tot_nz, dup_nz = 0, 0
per_graph = []
for g in np.unique(gid):
    m = (gid == g) & nz
    c = collections.Counter(fp[m].tolist())
    d = sum(v - 1 for v in c.values() if v > 1)
    tot_nz += m.sum(); dup_nz += d
    per_graph.append(100 * d / max(m.sum(), 1))
print(f"\nWITHIN-GRAPH duplicate rate (collision null, non-zero rows): "
      f"{dup_nz}/{tot_nz} = {100*dup_nz/tot_nz:.3f}%")
print(f"  per-graph range: {min(per_graph):.3f}% - {max(per_graph):.3f}%")

# ---------------------------------------------------------------- overlap
hdr("D. CROSS-SPLIT FINGERPRINT OVERLAP")
train_fps = set(fp[is_train & nz].tolist())
for s in ["val", "test"]:
    m = (split == s)
    mnz = m & nz
    hit = np.array([f in train_fps for f in fp[mnz]])
    print(f"{s}: non-zero nodes {mnz.sum():5d}/{m.sum():5d}; "
          f"fingerprint seen in TRAIN: {hit.sum():5d} = {100*hit.mean():.2f}% of non-zero, "
          f"{100*hit.sum()/m.sum():.2f}% of all {s} nodes")

# How many distinct graphs does each fingerprint span?
fp2g = collections.defaultdict(set)
for f_, g in zip(fp[nz], gid[nz]):
    fp2g[f_].add(int(g))
spans = np.array([len(v) for v in fp2g.values()])
print(f"\ndistinct non-zero fingerprints: {len(spans)}")
print(f"  graphs spanned per fingerprint: mean={spans.mean():.2f} "
      f"median={np.median(spans):.0f} max={spans.max()}")
print(f"  fingerprints appearing in >1 graph: {(spans>1).sum()} "
      f"({100*(spans>1).mean():.1f}%)")

# ---------------------------------------------------------------- labels
hdr("E. LABEL STRUCTURE")
pos = labels.sum(1)
print(f"positive labels per node: mean={pos.mean():.2f} median={np.median(pos):.0f} "
      f"min={pos.min()} max={pos.max()} (of {L})")
print(f"overall label density: {100*labels.mean():.2f}%")
print(f"  zero-feature nodes : mean positives={pos[zero].mean():.2f}")
print(f"  non-zero-feat nodes: mean positives={pos[nz].mean():.2f}")
print(f"distinct label vectors: {len(set(map(bytes, np.packbits(labels,axis=1))))}")

hdr("F. LABEL AGREEMENT ACROSS GRAPHS FOR SAME FINGERPRINT")
print("(tests whether labels are tissue-agnostic -> direct answer leakage)")
groups = collections.defaultdict(list)
for i in np.where(nz)[0]:
    groups[fp[i]].append(i)
multi = {k: v for k, v in groups.items() if len(set(gid[v].tolist())) > 1}
print(f"fingerprints occurring in >1 graph: {len(multi)}")
ident, agrees, sizes = 0, [], []
for k, idx in multi.items():
    lv = labels[idx]
    sizes.append(len(idx))
    if (lv == lv[0]).all():
        ident += 1
    agrees.append((lv == lv[0]).mean())
print(f"  groups with byte-identical label vectors: {ident} "
      f"({100*ident/max(len(multi),1):.2f}%)")
print(f"  mean per-bit agreement with first member: {100*np.mean(agrees):.2f}%")
print(f"  mean occurrences per multi-graph fingerprint: {np.mean(sizes):.2f}")

# Random-pair control: agreement between unrelated non-zero nodes
a = rng.choice(np.where(nz)[0], 20000)
b = rng.choice(np.where(nz)[0], 20000)
print(f"  CONTROL random node pairs, per-bit agreement: "
      f"{100*(labels[a]==labels[b]).mean():.2f}%")

# ---------------------------------------------------------------- baselines
hdr("G. BASELINES (test micro-F1)")


def micro_f1(y, p):
    tp = int((y & p).sum()); fp_ = int((~y & p).sum()); fn = int((y & ~p).sum())
    pr = tp / (tp + fp_) if tp + fp_ else 0.0
    rc = tp / (tp + fn) if tp + fn else 0.0
    return 2 * pr * rc / (pr + rc) if pr + rc else 0.0


ytr = labels[is_train].astype(bool)
yte = labels[is_test].astype(bool)

# G1 constant: predict every label whose train frequency > 0.5
freq = ytr.mean(0)
const = np.tile(freq > 0.5, (yte.shape[0], 1))
print(f"G1 constant (train freq>0.5, {int((freq>0.5).sum())} labels on): "
      f"micro-F1 = {micro_f1(yte, const):.4f}")

# G1b predict ALL ones
print(f"G1b predict-all-ones                        : "
      f"micro-F1 = {micro_f1(yte, np.ones_like(yte)):.4f}")

# G2 fingerprint lookup (pure memorisation, no learning, no graph)
lut = {}
for k, idx in collections.defaultdict(list, {
        k: v for k, v in ((k, [i for i in v if is_train[i]])
                          for k, v in groups.items()) if v}).items():
    lut[k] = labels[idx].mean(0) > 0.5

test_idx = np.where(is_test)[0]
pred = np.zeros_like(yte)
hits = 0
for r, i in enumerate(test_idx):
    if not zero[i] and fp[i] in lut:
        pred[r] = lut[fp[i]]; hits += 1
    else:
        pred[r] = freq > 0.5
print(f"G2 fingerprint lookup (memorisation)        : "
      f"micro-F1 = {micro_f1(yte, pred):.4f}   [{hits}/{len(test_idx)} "
      f"= {100*hits/len(test_idx):.1f}% resolved by lookup]")

# G2b: restricted to the resolved subset only
res = np.array([(not zero[i]) and (fp[i] in lut) for i in test_idx])
if res.sum():
    print(f"    on resolved subset only                 : "
          f"micro-F1 = {micro_f1(yte[res], pred[res]):.4f}")
    print(f"    constant baseline on same subset        : "
          f"micro-F1 = {micro_f1(yte[res], const[res]):.4f}")

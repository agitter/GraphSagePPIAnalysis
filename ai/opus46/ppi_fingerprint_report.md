# PPI (GraphSAGE 24-tissue) — Initial Fingerprinting & Leakage Audit

Run on the real `ppi-G.json` / `ppi-feats.npy` / `ppi-class_map.json` / `ppi-id_map.json`
release. All numbers below are measured, not cited.

---

## 0. Provenance and integrity

`snap.stanford.edu` and the DGL S3 bucket are both blocked by this environment's network
policy, so the files were obtained from a GitHub mirror (`sufeidechabei/PPI-Inductive`,
`ppi/`). I could not checksum against the canonical host. Instead I verified every
published statistic, and all of them match exactly:

| Property | Measured | Published |
|---|---|---|
| Nodes | 56,944 | 56,944 |
| Edges | 818,716 | 818,716 |
| Features | 50 | 50 |
| Labels | 121 | 121 |
| Graphs | 24 | 24 |
| Train / val / test nodes | 44,906 / 6,514 / 5,524 | 44,906 / 6,514 / 5,524 |
| Graphs per split | 20 / 2 / 2 | 20 / 2 / 2 |

SHA-256 recorded in `verify.py` output for future comparison against the canonical release.
**Treat provenance as strong but not cryptographically established.**

### Three assumptions falsified at this stage

1. **Connected components ≠ tissue graphs.** The union graph has **295** connected
   components, not 24. Tissue membership comes from `train_graph_id.npy` — a misleadingly
   named file that actually contains a graph id for **all 56,944 nodes**, not just training
   ones. Anyone partitioning PPI by connected components is silently working with a
   different dataset.
2. **`ppi-id_map.json` contains no gene information.** It is the identity permutation over
   numeric strings (`"0"→0, "1"→1, …`). This confirms GraphSAGE Issue #188 directly: the
   release carries *zero* protein identifiers.
3. **The graph has 25,084 self-loops.** Removed before all degree/adjacency computations
   here; not all published pipelines do this.

No tissue graph spans more than one split — the 20/2/2 partition is clean.

---

## 1. The feature-fingerprint method **fails**. My earlier recommendation was wrong.

I proposed identifying proteins across tissue graphs by their 50-dim binary feature vector,
calibrating against within-graph duplicates (a gene appears at most once per tissue graph,
so within-graph duplicates are collisions by construction).

The calibration null kills the method:

```
distinct feature vectors, all 56,944 nodes : 592
distinct feature vectors, non-zero nodes   : 591 (of 32,733)
WITHIN-GRAPH duplicate rate (non-zero)     : 23,210 / 32,733 = 70.91%
  per-graph range                          : 54.57% – 73.94%
```

A 70.9% false-positive rate means the feature vector cannot identify a protein. The
nominally striking overlap number — **99.90% of non-zero-feature test nodes have a feature
vector present in training** — is therefore meaningless as evidence of leakage. It is what
you get when 32,733 nodes are drawn from 591 possible states.

The calibration step is what caught this. Without it, that 99.90% would have looked like a
headline finding.

### The failure is itself a finding

The features carry **≈9.2 bits of information across the entire dataset**. Combined with
**42.52% of nodes having an all-zero feature vector** (train 42.33%, val 43.38%,
test 43.03% — uniform across all 24 graphs, range 37.1%–44.5%), the input features are
close to uninformative.

### Label agreement: partially confounded, small real effect

| Comparison | Per-bit label agreement |
|---|---|
| Same fingerprint, **same** graph | 64.41% |
| Same fingerprint, **different** graph | 78.34% |
| Random node pairs | 62.48% |

The within-graph figure is the correct control (same feature state, definitionally different
proteins). Cross-graph exceeds it by **+13.9 points**, so there is *some* genuine
cross-tissue recurrence signal — but it is far smaller than the raw 93% agreement figure
suggested before controlling, and cannot support a strong leakage claim on its own.

---

## 2. Structural signatures **do** work — and this is where the leakage is

WL-1 colour = `(own feature vector, sorted multiset of neighbour feature vectors)`.

```
distinct WL-1 colours                      : 34,172   (vs 592 raw fingerprints)
WITHIN-GRAPH collision null                : 2.52%    (vs 70.91% for raw features)
TEST nodes whose WL-1 colour appears in TRAIN : 3,084 / 5,524 = 55.83%
```

A 2.52% collision null means WL-1 is a usable identity proxy. And **55.83% of the
"completely unseen" inductive test set carries a 1-hop structural signature that already
appears verbatim in the training graphs.**

### Does the structural match carry the answer? Yes.

Pure lookup table — no parameters, no learning, no gradients. For each test node with a
WL-1 match in training, copy the majority label vector of matching training nodes; otherwise
fall back to the constant baseline.

| Method | Test micro-F1 |
|---|---|
| Constant (train freq > 0.5, 15 labels on) | 0.3935 |
| Predict all ones | 0.4608 |
| **WL-1 lookup, all test nodes** | **0.7274** |
| **WL-1 lookup, matched subset only (55.8%)** | **0.9791** |
| Constant baseline, *same* matched subset | 0.4041 |
| WL-2 lookup, matched subset (2.23% match rate) | 0.9341 |

**0.9791 from a lookup table**, against 0.4041 for the constant baseline on the identical
subset. Over half the inductive test set is answerable by exact structural memorisation of
training data.

### Confound: matched nodes are systematically low-degree

This must be stated alongside the result.

| Test-node degree | n | Matched in train |
|---|---|---|
| 0–1 | 262 | 98.09% |
| 2–3 | 447 | 91.95% |
| 4–7 | 850 | 85.76% |
| 8–15 | 1,178 | 71.14% |
| 16–31 | 1,286 | 51.01% |
| 32–63 | 855 | 20.00% |
| 64+ | 646 | 3.41% |

Matched nodes: mean degree 12.2, mean 31.3 positive labels.
Unmatched: mean degree 51.0, mean 42.4 positive labels.

Exact WL-1 matching is far easier for small neighbourhoods, so the matched subset is
easier a priori. The 0.9791-vs-0.4041 comparison is *within* that subset, which controls for
subset difficulty — but the 55.83% headline rate should be read as "concentrated in the
low-degree half of the test set," not as uniform coverage.

---

## 3. Two claims from my earlier literature report are now contradicted

### Homophily (previous "Mechanism E") — **wrong on this dataset**

| Measure | Across edges | Random pairs |
|---|---|---|
| Per-bit label agreement | 62.05% | 63.30% |
| Mean label Jaccard | 0.3460 | 0.2315 |

By per-bit agreement, connected nodes are **no more label-similar than random pairs** —
marginally less. Jaccard shows a moderate real effect (0.35 vs 0.23), driven by hub nodes
carrying many labels. The "high functional homophily makes PPI trivially easy" mechanism I
cited from the general PPI literature **does not hold for this benchmark** and should be
withdrawn.

### The published feature-only baseline is below trivial

GraphSAGE reports 0.422 for raw-feature logistic regression; the MLP literature reports
~0.46. Measured here:

| Representation | Logistic regression micro-F1 |
|---|---|
| Raw features (50d) | 0.4316 |
| Features + 1-hop mean (100d) | 0.4607 |
| Features + 1-hop + 2-hop mean (150d) | 0.4771 |
| *Predict all ones (no model at all)* | *0.4608* |

**Predicting all ones beats the published raw-feature baseline and matches the 2-hop
linear model.** With 30.52% label density (mean 36.93 of 121 labels positive per node,
max 101), micro-F1 has a trivial floor near 0.46. Every "features-only baseline" number in
this literature sits at or below the value of having no model whatsoever.

---

## 4. Revised account of where the 99.5 comes from

Simple structure-aware linear models reach only 0.477. So the jump to 0.97–0.99 is **not**
explained by feature memorisation (features are near-empty) or by simple neighbourhood
smoothing (worth ~4 points) or by label homophily (absent).

What the evidence supports: high-capacity GNNs memorise **local structural signatures**, and
because all 24 tissue graphs are induced subgraphs of one global PPI network, a majority of
test-node neighbourhoods recur verbatim in training. A parameter-free lookup already
extracts 0.9791 on the matched 55.8%. This is consistent with the width-sensitivity
documented in the literature (Cluster-GCN needing 2,048 hidden units for 99.36) — capacity
buys memorisation of a training-set structure table, not generalisation.

The benchmark's "inductive" framing describes graph disjointness, not neighbourhood
disjointness. Those are not the same thing, and the difference is worth ~0.98 F1 on half the
test set.

---

## 5. Next steps

1. **Degree-stratified evaluation.** Report GNN micro-F1 separately for WL-matched vs
   unmatched test nodes, and by degree band. Predicted: near-ceiling on matched, sharply
   lower on the high-degree unmatched half.
2. **Better cross-graph collision null.** The 2.52% figure is within-graph; construct a
   cross-graph null by matching between two *training* graphs to test whether tissue context
   changes the base rate.
3. **Leakage-controlled re-split.** Partition graphs so that WL-1 signature overlap between
   train and test is minimised, then re-run a standard GAT. The delta is the leakage estimate.
4. **Recover gene identity properly.** Feature fingerprints cannot do it; WL-1 plausibly can.
   Match the 24 graphs against the OhmNet `bio-tissue-networks.tar.gz` release (which does
   carry Entrez IDs) by `(n_nodes, n_edges)` signature and degree sequence, then verify by
   subgraph isomorphism. Test the cheap hypothesis first: if the graphs were merged with
   `networkx.disjoint_union`, node ordering within each block may be preserved.
5. **Replace micro-F1** with macro-F1 and information-content-weighted Fmax, given the
   0.46 trivial floor.

## Caveats

- Mirror provenance verified by statistics, not by checksum against SNAP.
- WL-1 collision null measured within-graph; cross-graph base rate may differ.
- The matched/unmatched degree confound is real and quantified above; the headline 55.83%
  is concentrated in low-degree nodes.
- Logistic-regression baselines are one-vs-rest with default regularisation and no tuning;
  they establish a floor, not a ceiling.
- No GNN was trained here — the 0.97–0.99 figures are from the literature, not reproduced.

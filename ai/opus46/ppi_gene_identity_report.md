# PPI Audit, Part 2 — Gene Identity Recovered, Label Leakage Confirmed

The OhmNet upload resolved the node-ID question. The result is unambiguous.

---

## 1. The node ID problem is solved

The readme states the key fact: node IDs are Entrez gene IDs, and **"The same node ID
appearing in two edgelists represents the same gene."** That is exactly the cross-tissue
identity the GraphSAGE release destroys.

### Recovering the GraphSAGE → OhmNet correspondence

Node counts didn't match directly, but the residuals showed a consistent
**Δnodes == Δedges** pattern — the signature of deleting nodes whose only incident edge is
a self-loop (each such deletion removes exactly one node and one edge). Applying that
transform to the 144 OhmNet tissue networks produced **16 exact (nodes, edges) matches, all
16 confirmed by degree sequence.**

| | |
|---|---|
| Identified train graphs | adrenal_gland, amygdala, aorta, artery, astrocyte, basal_ganglion, basophil, blood_plasma, blood_platelet, bone, colon, eye, forebrain, gastrointestinal_tract, large_intestine |
| Identified test graph | **midbrain** (g24) |

### Recovering per-node gene IDs

Node ordering is *not* preserved (the sorted-Entrez-ID hypothesis fails: ~2% positionwise
degree agreement). Recovery instead used 1-WL colour refinement — both sides are the same
graph up to permutation, so nodes in singleton colour classes on both sides map uniquely.

```
TOTAL mapped 33,641/34,107 = 98.63% of nodes in the 16 identified graphs
edge verification: 472,854/472,854 = 100.0000% of mapped edges reconstruct exactly
```

Every mapped edge reconstructs exactly. A spurious mapping does not reproduce 472,854 edges.

### External confirmation

| Entrez | Gene | Check |
|---|---|---|
| 351 | **APP** (amyloid beta precursor protein) | NCBI Gene ID 351; "concentrated in the synapses of neurons" — and it is the **highest-degree node in the midbrain graph** |
| 7157 | **TP53** | confirmed, classic interactome hub |
| 1017 | **CDK2** | confirmed (Wikidata cites 1017 as CDK2's canonical example) |

APP topping the midbrain network is precisely what a neural tissue PPI network should look
like. The mapping is biologically coherent, not just combinatorially consistent.

---

## 2. The leakage, measured directly

### Gene overlap

```
distinct genes across the 15 identified TRAIN graphs : 3,531
midbrain (TEST) mapped genes                         : 2,258
    ... also present in identified train graphs      : 2,258  = 100.00%
```

**Not one gene in the test graph is new.**

### Features are a gene-level property (as hypothesised)

```
genes appearing in >1 identified graph            : 3,186
identical FEATURE vector in every graph           : 3,186/3,186 = 100.00%
```

### Labels are entirely tissue-agnostic

```
identical LABEL vector in every graph             : 3,186/3,186 = 100.00%
mean per-bit agreement across occurrences         : 100.00%
```

This settles the question I flagged as unresolved in Part 1. OhmNet's *upstream* labels are
tissue-specific (`smooth_muscle_GO:0048661.lab` = does this gene have this function **in
smooth muscle**). The GraphSAGE repackaging **discarded that tissue specificity**: a gene
carries a byte-identical 121-dim label vector in all 24 tissues.

### Consequence

```
midbrain test genes seen in training with IDENTICAL 121-label vector : 2,258 (100.00%)
                             seen in training with different labels  :     0 (0.00%)
                                       not in any identified train   :     0 (0.00%)
```

**Gene-identity oracle micro-F1 on covered test nodes = 1.0000.**

A lookup table keyed on gene identity reproduces the test labels with **zero error**. The
"unseen inductive test graph" asks no question that training has not already answered.

### The pigeonhole argument (independent of the matching)

```
gene universe across all 144 OhmNet tissues : 4,510
node instances in the GraphSAGE benchmark   : 56,944
mean occurrences per gene                   : 12.63
```

Even if every one of my 16 identifications were wrong, 56,944 node instances drawn from a
4,510-gene universe means the average gene appears ~12.6 times, and 44,906 training
instances necessarily cover nearly all of it. **The conclusion does not depend on the
matching being right.**

---

## 3. On your F1 objection — you are correct

> *The F1 score is misleading if it is not on all test nodes, isn't it?*

Yes. A micro-F1 on a matched subset is a **conditional** quantity, P(correct | matched), and
is not comparable to a published full-test-set number. I reported both (0.7274 full, 0.9791
matched) but led with the subset figure, which invited exactly the wrong comparison. Three
separate problems:

**(a) It is not the same denominator.** Published PPI numbers are over all 5,524 test nodes.
A subset F1 answers a different question and must never be placed in the same column.

**(b) The WL-matched subset was selection-biased.** Matched nodes had mean degree 12.2 vs
51.0 unmatched, because exact neighbourhood matching is easy for small neighbourhoods.
So 0.9791 conflated "the answer was memorised" with "these nodes were easy."

**(c) A leakage claim needs two numbers, not one.** Coverage (what fraction of the test set
is answerable from training) and conditional fidelity (how completely). One without the
other is uninterpretable.

### The gene oracle fixes (b) and reports both

| Quantity | Value |
|---|---|
| Coverage of full test set | 2,258/5,524 = **40.88%** |
| — limited only by having identified 1 of 2 test graphs | |
| Coverage *within* midbrain | 2,258/2,300 = **98.2%** |
| micro-F1, covered nodes only | **1.0000** |
| micro-F1, full test set (constant fallback elsewhere) | 0.6993 |
| Constant baseline, covered nodes only | 0.3943 |
| Constant baseline, full test set | 0.3935 |

Crucially, gene-oracle coverage is **not degree-biased** — it covers 98.2% of midbrain
regardless of degree; the 42 misses are WL-symmetric nodes, an artifact of my recovery
method, not of the data. So the 1.0000 is not a selection effect. Note also that the
covered-subset constant baseline (0.3943) is indistinguishable from the full-set one
(0.3935), confirming the covered subset is *not* intrinsically easier.

The 40.88% figure is a floor set by my incomplete tissue identification, not by the data.
Within the one test graph fully analysed, leakage is **98.2% coverage at perfect fidelity**.

---

## 4. Revised conclusion

Part 1 argued that high-capacity GNNs memorise local structural signatures. That was
directionally right but named the wrong mechanism. The actual mechanism is simpler:

> The 24 tissue graphs are subgraphs of one global network over ~4,510 genes, and the 121
> labels are gene-level GO annotations with tissue specificity stripped out. Every test
> protein appears in training carrying an identical answer vector. The split partitions
> *graphs*, not *genes* — so it is not an inductive benchmark in any meaningful sense.

The structural memorisation observed in Part 1 is the *route* by which a GNN reaches
gene identity, not an independent effect: WL signatures work because they are a proxy for
gene identity, which the release stripped but did not actually remove.

Reported micro-F1 of 0.995 is therefore consistent with near-complete memorisation of a
gene→label table that a lookup implements exactly.

---

## 5. Next steps

1. **Identify the remaining 8 graphs.** Residuals are small and mostly Δn == Δe with the
   opposite sign (g2≈adrenal_cortex +2/+2, g10≈blood +4/+5, g14≈brain +2/+3,
   g19≈liver +4/+2, g22≈kidney +8/+6; g21≈heart and g23≈lung are −29/−27 and −26/−26,
   suggesting my pruning was slightly too aggressive there). A second-order preprocessing
   step explains these; resolving it extends coverage to ~100% of the test set.
2. **Build the gene-disjoint re-split** — the C1/C2/C3 protocol from the PPI-prediction
   literature. Partition the 4,510 genes, not the graphs, then re-run GAT/Cluster-GCN. The
   drop from 0.99 is the leakage estimate.
3. **Restore tissue-specific labels** from `bio-tissue-labels` (503 `.lab` files). This
   makes the task genuinely tissue-dependent and non-degenerate — arguably the benchmark
   OhmNet originally intended.
4. **Report coverage alongside every subset metric** in the write-up.

## Caveats

- 16/24 graphs identified; 15 of 20 training graphs. All overlap figures are **lower bounds**.
- One of two test graphs identified, so the 40.88% full-test coverage understates.
- 98.63% of nodes in identified graphs mapped; WL-symmetric nodes unmapped.
- The GraphSAGE files came from a GitHub mirror (canonical hosts blocked); verified by
  statistics and now additionally by exact edge reconstruction against OhmNet.
- No GNN trained; the 0.99 figures remain from the literature.

# Executive summary: hierarchy-aware GraphSAGE PPI splitting

## Recommended design

The most rigorous test is not one optimized split. It is a paired experiment over the same heldout tissue quartets:

1. Use only OhmNet hierarchy leaves in the primary analysis.
2. Enforce the manuscript thresholds: at least 15,000 raw edge records for training layers and at least 35,000 for heldout layers.
3. Choose four heldout leaves from four distinct top-level hierarchy branches.
4. Exclude those entire branches from training.
5. Choose 20 training leaves, four from each node-count quintile.
6. For each heldout quartet, select the training set that minimizes the worst train-heldout Wu-Palmer similarity, then total similarity, then maximizes path distance.
7. Compare it with random training sets subject to exactly the same constraints.
8. Use no GeneIDs, labels, overlap measures, or F1 values until after each split is fixed.

All 1,062 feasible leaf-only heldout quartets were evaluated exactly.

## Main result

Conditional on the same heldout quartet, hierarchy optimization reduced mean hierarchy similarity for every quartet. Relative to matched hierarchy-blind training sets, it reduced on average:

- heldout row overlap by 0.553 percentage points;
- unique-gene overlap by 1.140 percentage points;
- GeneID-lookup micro-F1 by 0.216 percentage points.

After tightly matching training node and edge totals, the reductions were 0.349, 0.699, and 0.134 percentage points, respectively. Hierarchy-aware training had lower row overlap for 66.8% and lower lookup F1 for 62.9% of heldout quartets in that tight size-matched comparison.

The answer is therefore **yes on average, but not always**: ignoring the hierarchy usually induces somewhat greater biological node reuse, but the hierarchy is only an imperfect proxy for GeneID overlap.

## Maximum-separation stress test

The primary global optimum holds out:

- colon and fetus;
- hematopoietic_stem_cell and pancreas.

Its 20 training tissues are listed in `results/primary_leaf_split_tissues.tsv`. It has mean WUP 0.172609, no ancestor-descendant pairs, row overlap 0.933361, and lookup F1 0.972114. The same split is selected by WUP minimax, maximum minimum path distance, and maximum mean path distance.

This split is intentionally difficult but its training graphs occupy only four root branches, so it should be treated as a stress test rather than the only operational benchmark.

## Broad-coverage counterexample

A sensitivity that forces training to cover all 14 non-heldout eligible root branches still achieves very low hierarchy similarity (mean WUP 0.194152), yet its heldout row overlap is 0.989476 and lookup F1 is 0.995478.

This proves that maximizing tissue-hierarchy distance does not by itself produce unseen biological nodes. A broad union of hierarchy-distant training tissues can cover nearly every heldout gene.

## Retained discrepancies

The deposited GraphSAGE split does not satisfy the manuscript thresholds in the current OhmNet files: astrocyte and basophil are below the 15,000-edge training cutoff and midbrain is below the 35,000-edge heldout cutoff.

The deposited midbrain test graph is a descendant and exact node-and-edge subnetwork of training brain. All 2,310 midbrain rows have GeneIDs present in training.

## Bottom line

Use hierarchy-aware splitting to test cross-tissue-domain generalization, but pair it with matched random splits and report the full distribution. Use a GeneID-disjoint split as a separate benchmark when the claim concerns generalization to unseen proteins.

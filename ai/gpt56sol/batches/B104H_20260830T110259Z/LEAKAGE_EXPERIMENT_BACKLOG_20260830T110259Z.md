# Saved experiment: identity leakage under randomized biological labels

## Question

Can a trivial model attain essentially the same PPI test score after destroying all biological meaning in the labels, solely because the same Entrez Gene IDs occur in training and test tissue graphs with repeated labels?

## Primary model

A zero-parameter lookup model:

1. Collapse training rows by Entrez Gene ID.
2. Store that gene's 121-dimensional training label vector.
3. For validation or test rows whose GeneID appeared in training, return the stored vector.
4. For unseen GeneIDs, return an all-zero vector.

On the original labels, the current full mapping gives test micro-F1 approximately `0.9971784`, with 5,490 of 5,524 test rows corresponding to GeneIDs observed in training.

## Primary randomized-label control

Randomize at the **gene level**, not independently at the row level:

1. Construct one 121-vector per distinct Entrez Gene ID.
2. Randomly permute complete label vectors among GeneIDs, or generate random vectors with matched column prevalences and label-count distribution.
3. Copy each randomized gene vector to every tissue occurrence of that GeneID.
4. Keep the original graph split unchanged.
5. Fit the same lookup model.

This destroys the gene-to-function biology while preserving the data-leakage mechanism: a repeated identity has the same target in every tissue.

## Exact-performance control

For an exact algebraic control, permute complete label vectors separately within:

- genes appearing in training; and
- genes absent from training but appearing in validation/test.

This preserves the number of positive cells in the seen and unseen strata, so the lookup model's micro-F1 is exactly unchanged while gene-function correspondence is randomized.

## Essential controls

- Independent row-wise label randomization: should destroy the lookup advantage and demonstrates why gene-level repetition is the relevant unit.
- Gene-disjoint split: no Entrez Gene ID may occur in more than one split. Lookup performance should collapse.
- Tissue split with anonymous row IDs only: verifies that the leakage is biological identity reuse, not literal node-number reuse.
- Multiple random seeds and confidence intervals.
- Report micro-F1, macro-F1, Hamming loss, exact-match accuracy, seen-gene metrics, and unseen-gene metrics separately.
- Compare against a constant-label baseline, feature-only linear model, and a small GNN under both the original and gene-disjoint splits.

## Interpretation boundary

The randomized-label experiment is designed to show that high benchmark scores can arise from identity reuse even when the labels have no biological meaning. It does not imply that every published model literally memorized GeneIDs, but it establishes that the benchmark permits and rewards that solution.

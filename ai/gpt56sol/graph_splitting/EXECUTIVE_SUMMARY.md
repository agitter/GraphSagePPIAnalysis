# Executive summary

The GraphSAGE supplement says that 20 PPI networks with at least 15,000 edges were randomly selected for training and that four networks with at least 35,000 edges were selected for validation and testing. The deposited split is not compatible with those thresholds in the currently released OhmNet files: `astrocyte` and `basophil` fall below the training threshold and `midbrain` falls below the heldout threshold under raw edge counts; `basophil` remains too small even under directed-arc interpretations.

The released test set contains 5,524 tissue-instance rows, of which 5,490 (99.3845%) have a recovered Entrez GeneID already present in training. A forensic diagnostic that memorizes each training gene's 121-label vector and predicts zeros for unseen genes obtains micro-F1 0.9971784.

Alternative graph assignments show both robustness and sensitivity:

- Exhaustively reallocating the same 24 tissues across all 63,756 possible 20/2/2 splits gives median F1 0.998563; 92.95% exceed 0.99 and 22.09% are perfect. The released split is at only the 30.99th percentile.
- Sampling 1,000,000 threshold-valid splits from all 144 released layers gives median F1 0.995031 to 0.995452, depending on whether heldouts or training are selected first. The released value is at the 66.57th to 69.44th contextual percentile, although the released split itself is outside this universe.
- Restricting the universe to the 107 hierarchy leaves lowers median F1 to 0.990422 to 0.991944. The released value is at the 89.75th to 92.91st contextual percentile, but the released split contains 10 internal hierarchy layers and cannot occur in this universe.

The tissue hierarchy contains 107 leaf networks and 37 internal networks among the 144 files. All 167 network-file ancestor-descendant pairs show exact descendant-to-ancestor node and edge containment. GraphSAGE placed nine internal layers in training; the complete `midbrain` test graph is contained in the training `brain` graph, giving midbrain a perfect lookup score.

Hierarchy similarity is not a strong substitute for direct entity overlap. Across all 144 network pairs, Wu-Palmer hierarchy similarity and gene-set Jaccard similarity have Spearman rho 0.185. Across the all-144 sampled splits, F1 and row overlap have rho about 0.994, while F1 and mean hierarchy similarity have rho only 0.105 to 0.117.

The defensible conclusion is that the precise historical randomization remains unresolved, but the near-perfect lookup behavior is a robust structural property of these overlapping and hierarchically nested tissue networks rather than a one-off consequence of the deposited split.

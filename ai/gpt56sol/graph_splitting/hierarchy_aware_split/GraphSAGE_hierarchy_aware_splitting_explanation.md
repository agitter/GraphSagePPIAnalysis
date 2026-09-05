## The biological and benchmarking context

The GraphSAGE PPI benchmark is not one large network whose individual nodes were divided into training, validation, and test sets. It is a collection of **separate tissue-specific protein-interaction networks**. A node in one network represents a gene active in that tissue, and the same Entrez GeneID appearing in two tissue networks represents the same biological gene. fileciteturn5file1L61-L69

GraphSAGE trained on 20 tissue networks, used two networks for validation, and evaluated on two test networks. The manuscript describes the test graphs as entirely unseen and states that training networks had at least 15,000 edges, while validation and test networks had at least 35,000 edges. fileciteturn5file0L17-L25 fileciteturn5file4L260-L265

Holding out a tissue network, however, does not necessarily hold out the genes in that network. A gene can occur in many tissues. A model may therefore be evaluated on a nominally new graph while seeing mostly the same biological entities it encountered during training. The purpose of the hierarchy-aware analysis was to ask whether the original hierarchy-blind selection of tissues made this reuse worse than necessary.

There is an additional complication: the current OhmNet release does not contain only specific tissues. Its 144 network files comprise:

- 107 **leaf layers**, representing relatively specific tissues or cell types at the ends of the hierarchy;
- 37 **internal layers**, representing broader categories such as `brain`, `blood`, or `kidney`.

In the released OhmNet data, a descendant layer can be completely contained within its broader ancestor layer. Placing `brain` in training and `midbrain` in testing therefore creates two nominally separate graph files even though the test network is nested inside the training network. This is exactly what occurred in GraphSAGE: all 2,310 midbrain test rows correspond to genes already present in training, and the complete midbrain node and edge set is contained in the training brain layer. fileciteturn6file0L55-L68

## What “hierarchy-aware splitting” is intended to accomplish

The OhmNet hierarchy is a tree connecting broad and specific tissue concepts. Two tissues are considered close when they share a recent common ancestor in that tree. For example, two brain regions are more closely related than a brain region and a gastrointestinal tissue.

A hierarchy-aware split tries to place the validation and test tissues far from the training tissues in that tree. This is intended to create a stronger test of **generalization to biologically different tissue contexts**.

That goal is distinct from creating a **gene-disjoint** split:

- A hierarchy-aware split asks whether the model transfers to distant tissues.
- A gene-disjoint split asks whether the model transfers to proteins never encountered during training.

Distant tissues can still share many genes, so the first design does not guarantee the second.

## Why there is no single uniquely correct hierarchy-aware split

Several reasonable rules could be used, and they test slightly different notions of generalization.

A weak rule would simply forbid direct ancestor-descendant relationships. This would prevent `brain` from being paired with `midbrain`, but it could still place close sibling tissues on opposite sides of the split.

A stronger rule would place heldout tissues in different major branches of the hierarchy and exclude those whole branches from training. This creates a clearer domain shift, but it can narrow the biological diversity of the training set.

One could optimize the **average** distance between training and heldout tissues. However, an ordinary average can conceal one severe relationship. A split might look distant overall while still containing a pair such as training `brain` and test `midbrain`. For that reason, the primary procedure protects against the **closest** train-heldout pair first and considers average separation only afterward. fileciteturn6file0L77-L100

A further alternative is to require training to cover as many major hierarchy branches as possible. That produces a more biologically diverse training set, but the union of many diverse tissues can contain almost every gene found in the heldout tissues.

These are not interchangeable objectives. We therefore treated the maximum-separation and broad-coverage constructions as complementary analyses rather than declaring one universally optimal split. fileciteturn6file0L113-L129

## The primary split-construction procedure, in plain language

### 1. Apply the manuscript’s size rules

Only networks with at least 15,000 raw OhmNet edge records were eligible for training. Only networks with at least 35,000 records were eligible for validation or testing.

Under the current files, this produced 69 training-eligible leaf networks and 15 heldout-eligible leaf networks. The leaf-only restriction was introduced to avoid comparing a broad internal layer with one of its own descendant subnetworks. It is an analysis choice, not a claim that the GraphSAGE authors used only leaves. The analysis was also repeated with all 144 layers as a sensitivity check. fileciteturn6file0L59-L75

The deposited GraphSAGE split itself cannot be generated under this literal interpretation because `astrocyte` and `basophil` fall below the stated training threshold and `midbrain` falls below the stated heldout threshold in the present OhmNet release. Consequently, the released split was used only as an external comparison, not treated as one of the valid candidate splits. fileciteturn6file0L43-L53

### 2. Choose four biologically separated heldout tissues

Four heldout leaves were chosen, with each required to lie in a different top-level branch of the tissue hierarchy.

This prevents a heldout set composed entirely of closely related tissues, such as several neighboring brain regions. Of the 1,365 possible sets of four heldout leaves, 1,062 met this rule and left enough eligible training tissues to construct the required training set. All 1,062 were evaluated. fileciteturn6file0L131-L144

### 3. Exclude the heldout branches from training

After choosing the four heldout tissues, every candidate training tissue in those same four major branches was removed.

This is stronger than merely removing the four exact tissue names. It prevents the training set from including a parent, sibling, or close relative of a heldout tissue. Without this step, a nominally hierarchy-aware split could still place very similar tissue networks on opposite sides.

### 4. Prevent graph size from driving the result

Larger tissue networks contain more genes and are more likely to overlap other networks simply because of their size. A hierarchy-optimized algorithm could therefore appear to reduce gene overlap merely by selecting unusually small training networks.

To control this, eligible training tissues were sorted by node count and divided into five groups ranging from smaller to larger networks. Every training set had to contain four tissues from each group, for 20 tissues total.

A second, stricter comparison matched optimized and random training sets on their total numbers of nodes and edges. In the tightest comparison, the mean difference was only about 13 nodes and 138 edges across entire 20-network training sets. fileciteturn6file0L146-L165 fileciteturn6file0L200-L202

### 5. Maximize the separation of the closest pair

For each fixed group of four heldout tissues, the algorithm selected the 20 training tissues in this order of priority:

1. Make the **most similar** train-heldout pair as dissimilar as possible.
2. Among tied solutions, reduce overall train-heldout similarity.
3. Among remaining ties, increase ordinary tree-path distance.
4. Use a fixed alphabetical/index rule only for exact residual ties.

This is a “protect the weakest point first” strategy. It ensures that no heldout tissue is sacrificed merely to improve the average for the other three.

The optimization used only the tissue hierarchy and graph-size information. It did **not** inspect gene identities, labels, node overlap, or F1 scores. That separation is essential: otherwise the split would be selected using the outcome it was later supposed to evaluate. fileciteturn6file0L102-L111

### 6. Divide the four heldout tissues into validation and test pairs

Validation data are used to choose hyperparameters and models, so validation tissues should be protected from training similarity just as test tissues are.

For each quartet, the three possible ways to divide four tissues into two pairs were considered. The hierarchy-only rule selected the division that made the less-protected pair as distant from training as possible and kept validation and test difficulty reasonably balanced.

The labels “validation” and “test” should then be assigned using a predeclared rule or a recorded random draw—not by observing which pair gives more favorable model performance. fileciteturn6file0L167-L176

## How the hierarchy-aware rule was compared fairly with random splitting

A single highly optimized split is not enough to show that hierarchy awareness matters. It might simply contain unusually easy or difficult heldout tissues.

The stronger analysis was therefore paired:

- Fix one particular quartet of heldout tissues.
- Construct the best hierarchy-aware training set for that quartet.
- Generate many hierarchy-blind training sets for the **same quartet**.
- Require the random sets to obey the same leaf restriction, branch exclusions, edge thresholds, and size balance.
- Compare their biological overlap only after all splits have been selected.

This is analogous to giving two methods the same test cases and changing only how the training examples are chosen. Because the heldout tissues and all non-hierarchy constraints are identical, the comparison isolates the contribution of hierarchy-based selection much more cleanly than comparing unrelated splits. fileciteturn6file0L178-L198

## What the outcome measures mean

Three different biological outcomes were evaluated after each split was fixed.

**Heldout row overlap** is the percentage of tissue-specific gene occurrences in validation or test whose GeneID appears somewhere in training. A value of 99% means that 99 of every 100 heldout rows represent genes already encountered in another tissue during training.

**Unique-gene overlap** counts each GeneID once, regardless of how many heldout tissues contain it. This prevents frequently repeated genes from dominating the calculation.

**GeneID-lookup F1** is a diagnostic rather than a trained GraphSAGE result. For a test gene seen in training, it simply reuses that gene’s training label vector. A near-perfect score therefore indicates that biological identity alone nearly determines the heldout labels; it does not demonstrate that the lookup method learned tissue-network structure. The existing GraphSAGE split gives this diagnostic an F1 of approximately 0.997 because 5,490 of 5,524 test rows have GeneIDs already present in training. fileciteturn5file6L396-L416

For hierarchy similarity, we used a standard tree-based score in which larger values mean closer tissue relationships, together with ordinary path length as a sensitivity measure. We recorded both average similarity and the largest single train-heldout similarity. fileciteturn6file0L77-L98

## What the paired results mean

After very close matching of training-set node and edge totals, hierarchy-aware training sets had, on average:

- **0.349 percentage points less heldout-row overlap**;
- **0.699 percentage points less unique-gene overlap**;
- **0.134 percentage points lower lookup F1**.

The row-overlap reduction corresponds to roughly 35 fewer previously seen tissue-gene rows per 10,000 heldout rows.

Hierarchy-aware selection reduced row overlap for 66.8% of heldout quartets and reduced lookup F1 for 62.9%. Thus the direction was more often favorable than unfavorable, but it was not universal. fileciteturn6file0L244-L257

The correct interpretation is:

> Ignoring the tissue hierarchy usually caused somewhat more reuse of the same genes across training and heldout tissues, but the average effect was modest. Tissue hierarchy explains only part of gene sharing, and the outcome distributions for hierarchy-aware and hierarchy-blind splitting still overlap substantially.

Repeating the analysis while allowing all 144 layers, including internal layers, gave the same qualitative result, so the conclusion was not dependent on the leaf-only restriction. fileciteturn6file0L259-L269

## Why the two example splits are both important

### Maximum-separation stress test

The strongest hierarchy-separated split held out:

- `colon`
- `fetus`
- `hematopoietic_stem_cell`
- `pancreas`

Its closest training tissue was still far away in the hierarchy, and it contained no ancestor-descendant train-heldout pairs. Yet **93.3% of heldout rows still represented genes seen during training**, and the lookup diagnostic remained 0.972. fileciteturn6file0L291-L349

This is a useful difficult benchmark, but its training tissues occupy only four major branches. Some of its difficulty arises because training-domain coverage is deliberately narrow. It should therefore be described as a stress test, not as the one definitive replacement split.

### Broad-coverage hierarchy-aware split

A second construction required training to cover all 14 eligible major branches not used for heldout tissues. It remained exceptionally well separated in the hierarchy, but:

- 98.95% of heldout rows represented genes already seen in training;
- unique-gene overlap was 97.16%;
- lookup F1 was 0.9955.

The training set was biologically diverse, but the union of many different tissues collectively covered almost every heldout gene. fileciteturn6file0L355-L384

This counterexample is central to the interpretation:

> A tissue split can be rigorously separated in the OhmNet hierarchy and still provide almost no separation of biological node identity.

## Report-level conclusion

A hierarchy-aware splitting algorithm is appropriate when the intended claim is that a model generalizes to **different tissue contexts**. The rigorous evaluation is not one selected split, but a paired distribution comparing hierarchy-aware and hierarchy-blind training sets for the same heldout tissues and under the same graph-size constraints.

The analysis indicates that hierarchy-blind selection produces somewhat greater gene reuse on average. However, hierarchy distance is not a sufficient safeguard against identity leakage. A separate GeneID-disjoint benchmark is required when the intended claim is that the model generalizes to **previously unseen proteins**. The released brain-to-midbrain relationship demonstrates why both average hierarchy measures and the phrase “unseen graph” can be misleading when used without this additional biological context.

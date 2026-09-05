# GraphSAGE PPI graph-split, entity-overlap, and tissue-hierarchy investigation

**Investigation date:** 2026-09-04  
**Status:** completed side investigation; not part of the planned core reproduction package  
**Primary question:** How much could the GraphSAGE PPI train-test gene overlap and a GeneID-lookup diagnostic have changed under other graph splits consistent, or approximately consistent, with the published splitting description?

## Executive summary

The GraphSAGE supplement says that 20 PPI networks with at least 15,000 edges were randomly selected for training, while four networks with at least 35,000 edges were selected for validation and testing (two each). The deposited split cannot be generated literally from the currently released OhmNet edgelists under any ordinary edge-count convention we tested: the training tissues `astrocyte` and `basophil` have only 9,635 and 4,193 raw edge records, and the test tissue `midbrain` has 31,665. Even after counting undirected edges as directed arcs, `basophil` remains below 15,000. Therefore the historical candidate pool, source release, edge-count convention, or exception policy remains unresolved.

The central leakage result is nevertheless robust. The released split has 5,490 of 5,524 test tissue-instance rows (99.3845%) whose Entrez GeneID occurs in training. A diagnostic that memorizes each training gene's 121-label vector and predicts zeros for genes absent from training obtains pooled test micro-F1 = **0.9971784**. This is not the reported GraphSAGE model score and is not a fair baseline available from the anonymous release alone; it is a forensic diagnostic of entity reuse after recovering biological identities.

The effect of alternative graph choices depends strongly on what is treated as the candidate universe:

| Counterfactual universe and procedure | Analysis | Mean lookup F1 | 5th percentile | Median | Fraction F1 >= 0.99 | Fraction >= released F1 | Released F1 percentile |
|---|---:|---:|---:|---:|---:|---:|---:|
| The 24 deposited tissues; all 2-test/2-validation allocations | exact, 63,756 | 0.997242 | 0.988888 | 0.998563 | 92.95% | 69.01% | 30.99th |
| All 144 released layers; held-out first, uniform over valid raw-threshold splits | 1,000,000 samples | 0.993390 | 0.981837 | 0.995031 | 77.98% | 30.56% | 69.44th (external) |
| All 144 released layers; training first | 1,000,000 samples | 0.993914 | 0.983157 | 0.995452 | 80.80% | 33.43% | 66.57th (external) |
| Actual four held-outs fixed; random raw-threshold training set | 1,000,000 samples | 0.994591 | 0.987704 | 0.995562 | not tabulated | 21.43% | 78.57th (external) |
| The 107 hierarchy-leaf layers; held-out first | 1,000,000 samples | 0.987703 | 0.967355 | 0.990422 | 52.27% | 7.09% | 92.91st (external) |
| The 107 hierarchy-leaf layers; training first | 1,000,000 samples | 0.989842 | 0.974214 | 0.991944 | 61.73% | 10.25% | 89.75th (external) |

"External" means the released split is not a member of that counterfactual universe, so its percentile is contextual rather than the probability of selecting the released split.

Within the fixed set of 24 deposited tissues, the selected test pair explains **97.34%** of the variance in lookup F1; validation selection explains only 2.66% after conditioning on the test pair. The actual test pair, `lung` plus `midbrain`, is ordinary when averaged over all validation choices (49.82nd percentile among 276 test pairs), but the actual validation pair, `heart` plus `kidney`, is the eighth-lowest of 231 validation choices for that test pair. Withholding heart and kidney removes useful identity coverage: treating both validation graphs as additional training graphs would raise the diagnostic F1 from 0.997178 to 0.998854.

The supplied hierarchy contains 219 nodes and 218 edges. Exactly 107 of the 144 released network files correspond to hierarchy leaves and 37 to internal hierarchy nodes. GraphSAGE selected 14 leaf layers and 10 internal layers; nine internal layers are in training and `kidney` is an internal validation layer. Every one of the 167 ancestor-descendant pairs for which both nodes have network files satisfies exact node-set and edge-set containment from descendant to ancestor. Most strikingly, the test network `midbrain` is a descendant and exact subnetwork of the training network `brain`; all 2,310 midbrain genes are already in training, and its lookup F1 is 1.0.

Hierarchy proximity alone is not a dependable proxy for gene overlap. Across all 10,296 pairs of released networks, Wu-Palmer hierarchy similarity and gene-set Jaccard similarity have Spearman rho = 0.185; among the 36 raw-large networks, rho = 0.113. In the one-million-split all-144 simulations, lookup F1 correlates almost perfectly with row overlap (Spearman rho about 0.994) but only weakly with mean hierarchy similarity (rho 0.105 to 0.117) or the best hierarchy match in training (rho 0.223 to 0.233). The actual test pair is relatively far apart in the hierarchy, but the training set still provides broad gene coverage because it contains large and internal layers.

## 1. Scope and relation to the larger project

This is an exploratory side analysis. It does not alter the agreed scope of the clean GraphSAGE/DGL reproduction package. It builds on the earlier complete GraphSAGE node-to-Entrez mapping and exact 121-label reconstruction. Those prior results have their own attribution and evidence record in the project handoff; this report does not reassign credit for them.

The present work contributes:

1. exact identification of all 24 deposited GraphSAGE components against the supplied OhmNet network archive;
2. an audit of the published edge thresholds under several counting conventions;
3. exact enumeration of every 20/2/2 assignment among the deposited 24 tissues;
4. two one-million-draw simulations over all 144 released layers;
5. a fixed-heldout counterfactual and several edge-count sensitivity analyses;
6. two one-million-draw simulations restricted to the 107 hierarchy leaves;
7. a complete analysis of hierarchy distance, Wu-Palmer similarity, ancestor-descendant containment, and internal-layer selection;
8. deterministic seeds, hashes, independent checks, and compact result tables.

The DGL archive was not required and was not used in this analysis.

## 2. Source statements and evidence categories

### 2.1 What the GraphSAGE paper documents

Appendix B of the supplied GraphSAGE manuscript states:

- 20 PPI networks with at least 15,000 edges were randomly selected for training;
- four large networks with at least 35,000 edges were selected for validation and testing;
- two of those were used for validation and two for testing;
- the same validation networks and random training set were used in all experiments.

The main text describes training on 20 graphs and evaluation on two validation and two test graphs. It calls the held-out objects "entirely unseen graphs."

The wording leaves several unresolved choices: the universe of eligible networks, whether validation/test selection was random, whether training or held-out selection occurred first, what counted as an edge, and whether the released 144-file archive is the exact historical input.

### 2.2 What the OhmNet materials document

The supplied README says that each edgelist is one tissue-specific human PPI network, node IDs are Entrez GeneIDs, and the same ID in two layers denotes the same gene. It also supplies a tissue hierarchy. The official OhmNet materials describe 107 tissue-specific layers at the leaves of a 219-node hierarchy.

### 2.3 Evidence terminology used here

- **Documented:** directly stated by the GraphSAGE manuscript, OhmNet README, or official OhmNet materials.
- **Data-level exact:** calculated exactly from supplied bytes or established by exact comparison.
- **Exhaustive:** every member of a finite stated universe was evaluated.
- **Monte Carlo:** a fixed-seed sample from a stated probability mechanism; sampling error is reported or bounded.
- **Model-derived:** depends on the recovered global label-generation transformation outside the 4,301 genes observed in GraphSAGE.
- **Open:** the current evidence does not identify the historical mechanism uniquely.

## 3. Inputs, integrity, and topology matching

All supplied archives passed structural integrity checks. Source sizes and SHA-256 hashes are recorded in `results/source_manifest.csv` and `results/archive_integrity.csv`.

The 24 connected components in `graphsage_ppi.zip` were compared against all OhmNet edgelists as undirected edge multisets. Every component has one exact named match, with zero missing or extra records. The 24 matches total 818,716 edge records on both sides. The exact mapping and split assignment are in `results/graphsage_ohmnet_topology_validation.csv`.

The deposited split is:

- **Training (20):** adipose_tissue, adrenal_cortex, adrenal_gland, amygdala, aorta, astrocyte, artery, basal_ganglion, basophil, blood, blood_plasma, blood_platelet, bone, brain, colon, eye, forebrain, large_intestine, liver, gastrointestinal_tract.
- **Validation (2):** heart, kidney.
- **Test (2):** lung, midbrain.

## 4. Why the published threshold rule does not reproduce the deposited split

### 4.1 Literal raw-edgelist counts

Counting each line in a supplied undirected edgelist, including self-loops, produces 105 training-eligible networks with at least 15,000 records and 36 heldout-eligible networks with at least 35,000 records.

Three deposited assignments fail:

| Tissue | Deposited split | Raw records | Required | Result |
|---|---|---:|---:|---|
| astrocyte | training | 9,635 | >=15,000 | fails |
| basophil | training | 4,193 | >=15,000 | fails |
| midbrain | test | 31,665 | >=35,000 | fails |

### 4.2 Alternative ordinary edge-count conventions

Five count conventions were tested:

| Count convention | Training pool >=15k | Heldout pool >=35k | Deposited failures |
|---|---:|---:|---|
| undirected records including loops | 105 | 36 | training: astrocyte, basophil; heldout: midbrain |
| undirected non-self edges | 99 | 33 | training: adrenal_cortex, astrocyte, basophil; heldout: midbrain |
| directed non-self arcs | 131 | 89 | training: basophil |
| directed arcs, loop counted once | 133 | 90 | training: basophil |
| directed arcs, loop counted twice | 134 | 90 | training: basophil |

No ordinary interpretation makes the deposited assignment compliant. `basophil` is the decisive contradiction: even doubling every non-loop edge and counting loops twice gives only 8,386 arcs.

Nor can the same 24 tissues simply be reassigned to satisfy the thresholds. Under raw counts, `astrocyte` and `basophil` are too small for training and far too small for heldout assignment. Under directed conventions, `basophil` is still too small for either role.

### 4.3 Defensible interpretation

This analysis does not prove that the manuscript statement was false at the time. Plausible unresolved explanations include:

- a different historical OhmNet release;
- a count made before a filtering step not represented in the current archive;
- an undocumented exception or manual inclusion;
- a broader interaction count attached to a tissue before final network construction;
- imprecise reporting of the final selection rule.

Because the literal deposited split is outside the current threshold-defined universes, there is no unique historically correct null distribution. The multiple counterfactual universes below are therefore sensitivity analyses, not reconstructions of a proven randomization process.

## 5. Lookup diagnostic and label treatment

### 5.1 Diagnostic definition

For a proposed split:

1. collect the Entrez GeneIDs appearing in the 20 training tissues;
2. for every test occurrence of a seen GeneID, predict that gene's exact 121-dimensional label vector;
3. for every test occurrence of an unseen GeneID, predict the all-zero vector;
4. calculate pooled micro-F1 over all test rows and all 121 labels.

The recovered GraphSAGE labels are invariant across occurrences of the same GeneID. Consequently, this idealized diagnostic has no false positives: all error arises from positive labels on test genes absent from training. If `TP` is the number of covered positive cells and `FN` the number on unseen genes,

```
micro-F1 = 2*TP / (2*TP + FN)
```

This diagnostic should not be described as a model trained from GraphSAGE's anonymous public inputs. It uses the forensic node-to-GeneID recovery and measures how informative cross-split entity identity would be.

### 5.2 Extending labels to all 144 OhmNet layers

The 144 network files contain 4,510 distinct Entrez GeneIDs. Of these, 4,301 occur in the deposited GraphSAGE data and have directly validated released label vectors. For the remaining 209, labels were generated with the same recovered global GOA release-159 transformation that exactly reproduces every observed GraphSAGE label cell. One gene, 10159, is not covered by the historical mapping path and receives the all-zero vector; that is also consistent with its observed GraphSAGE behavior where applicable.

This extension is **model-derived**, not a comparison to unseen deposited targets. To bound its influence, a sensitivity run assigned zero positive-label weight to all 209 external genes while preserving every split and overlap calculation. Mean F1 rose by only 0.000710 in the heldout-first sample and 0.000702 in the training-first sample. The main conclusions are therefore not driven by extrapolated labels. Raw row-overlap and unique-gene-overlap statistics do not depend on label extension at all.

## 6. The released split

### 6.1 Direct overlap and F1

| Metric | Released value |
|---|---:|
| Test tissue-instance rows | 5,524 |
| Test rows whose GeneID occurs in training | 5,490 |
| Test row overlap | 0.993845 |
| Distinct test genes | 3,267 |
| Distinct test genes seen in training | 3,233 |
| Distinct-gene overlap | 0.989593 |
| Positive test label cells | 200,096 |
| Positive cells missed by unseen genes | 1,126 |
| Positive-cell coverage | 0.994373 |
| Lookup micro-F1 | **0.997178** |

The test tissues differ sharply:

| Test tissue | Rows | Rows seen in training | Unseen genes | Lookup F1 |
|---|---:|---:|---:|---:|
| lung | 3,214 | 3,180 | 34 | 0.995077 |
| midbrain | 2,310 | 2,310 | 0 | 1.000000 |

### 6.2 Effect of the validation choice

The validation tissues are not part of the training lookup. Twenty-one of the 34 test genes missing from the actual training set occur in heart or kidney. Counterfactually adding validation layers to training gives:

| Added to training | Unseen test rows | False-negative label cells | Lookup F1 |
|---|---:|---:|---:|
| none | 34 | 1,126 | 0.997178 |
| heart | 24 | 754 | 0.998112 |
| kidney | 17 | 575 | 0.998561 |
| heart and kidney | 13 | 458 | 0.998854 |

## 7. Exact enumeration within the deposited 24 tissues

This universe ignores the edge thresholds because no allocation of the 24 current tissues can satisfy them. It asks a clean finite question: if the same 24 named graphs had been allocated differently, how much would overlap and lookup F1 vary?

Every choice of two test graphs, two validation graphs from the remainder, and 20 training graphs was evaluated:

```
C(24,2) * C(22,2) = 63,756 splits
```

### 7.1 Distribution

| Statistic | Lookup micro-F1 | Test row overlap |
|---|---:|---:|
| minimum | 0.978316 | 0.946651 |
| 1st percentile | 0.985927 | 0.965098 |
| 5th percentile | 0.988888 | 0.971719 |
| 10th percentile | 0.991244 | 0.975972 |
| 25th percentile | 0.996676 | 0.991283 |
| median | 0.998563 | 0.996902 |
| 75th percentile | 0.999821 | 0.999416 |
| 90th percentile | 1.000000 | 1.000000 |
| maximum | 1.000000 | 1.000000 |

Additional exact facts:

- 92.954% of allocations have lookup F1 >= 0.99.
- 22.087% have F1 exactly 1.0.
- 69.012% have F1 at least as high as the released split.
- The released F1 is at the 30.989th midrank percentile.
- The released row overlap is at the 38.130th percentile.

The unique worst split tests on `blood` and `brain`, validates on `heart` and `lung`, and obtains F1 0.978316 with row overlap 0.946651. Even this worst allocation remains extremely high by ordinary classification standards.

### 7.2 Which random decision matters most

Grouping all 63,756 allocations by the identity of the two test tissues shows:

- 97.337% of F1 variance lies between the 276 possible test pairs;
- 2.663% remains among validation choices conditional on a test pair.

The `lung` plus `midbrain` test pair has mean F1 0.998497 across its 231 possible validation pairs, placing it at the 49.819th percentile among test-pair means. The test choice itself is therefore not unusually favorable within the selected 24.

The actual validation pair `heart` plus `kidney`, however, is the eighth-lowest of 231 choices for this fixed test pair and has a midrank position of 3.247%. It is a leakage-reducing validation choice because those two networks contain many of the lung genes that otherwise would be in training.

## 8. All 144 released network files under raw edge thresholds

Under literal raw line counts, 105 networks are eligible for training and 36 are eligible for heldout use. Two plausible randomization orders were simulated.

### 8.1 Procedure A: heldout first, uniform over valid labeled splits

1. select four distinct networks uniformly from the 36 large networks;
2. designate the first two as test and the next two as validation;
3. select 20 training networks uniformly from the remaining 101 training-eligible networks.

This is uniform over the stated valid labeled split space. The number of possible splits is:

```
C(36,2) * C(34,2) * C(101,20)
= 236,206,084,725,724,247,919,269,100
```

One million splits were sampled with seed 202609040001.

### 8.2 Procedure B: training first

1. select 20 networks uniformly from the 105 training-eligible networks;
2. select four remaining large networks;
3. designate two test and two validation.

This follows the sentence order in the supplement but is not uniform over final splits: training sets that contain different numbers of large networks leave different numbers of heldout choices. One million splits were sampled with seed 202609040002.

### 8.3 Results

| Statistic | Heldout first | Training first |
|---|---:|---:|
| mean F1 | 0.993390 | 0.993914 |
| SD | 0.006154 | 0.005714 |
| minimum sampled | 0.924701 | 0.922474 |
| 1st percentile | 0.971777 | 0.974124 |
| 5th percentile | 0.981837 | 0.983157 |
| 10th percentile | 0.985652 | 0.986661 |
| median | 0.995031 | 0.995452 |
| 90th percentile | 0.999258 | 0.999385 |
| 95th percentile | 1.000000 | 1.000000 |
| fraction F1 >= 0.99 | 77.977% | 80.803% |
| fraction F1 = 1 | 5.862% | 6.352% |
| mean test-row overlap | 0.983187 | 0.984469 |

The released score of 0.997178 is at the 69.437th and 66.566th contextual percentiles. Equivalently, 30.563% and 33.434% of sampled splits equal or exceed it. The released row overlap is at the 73.693rd and 71.001st contextual percentiles.

These distributions show that the particular graph choice can materially change leakage: sampled F1 extends down to about 0.92, while many splits are perfect. At the same time, the entire distribution remains heavily concentrated near one, so high entity-reuse performance is a general property of these tissue networks rather than a peculiarity of the one released split.

### 8.4 Fixed actual heldouts; only training randomized

Because the supplement explicitly uses "randomly" for the 20 training graphs but only says the four heldout graphs were "selected," a separate counterfactual fixed the released validation and test tissues and sampled one million 20-graph training sets from the raw >=15,000 pool.

The mean F1 was 0.994591, the median 0.995562, and the 5th to 95th percentile interval 0.987704 to 0.998371. The released training set's score lies at the 78.572nd contextual percentile; only 21.428% of sampled eligible training sets perform at least as well. Its row overlap is at the 83.841st percentile.

This comparison remains counterfactual because the released training set itself contains `astrocyte` and `basophil`, which are outside the current raw-threshold pool.

### 8.5 Edge-count sensitivity

Using undirected non-self edges gives nearly the same F1 distribution as raw counts. Treating each undirected edge as two directed arcs enlarges both pools and shifts the released score toward the middle of the sampled distribution:

| Count convention | Procedure | Mean F1 | Released contextual percentile |
|---|---|---:|---:|
| raw records | heldout first | 0.993390 | 69.44th |
| raw records | training first | 0.993914 | 66.57th |
| undirected non-self | heldout first | 0.993455 | 70.58th |
| undirected non-self | training first | 0.994053 | 67.23rd |
| directed non-self arcs | heldout first | 0.995014 | 51.93rd |
| directed non-self arcs | training first | 0.995117 | 51.39th |
| directed arcs, loop once | heldout first | 0.994936 | 52.30th |
| directed arcs, loop once | training first | 0.995020 | 51.78th |
| directed arcs, loop twice | heldout first | 0.994881 | 52.61st |
| directed arcs, loop twice | training first | 0.994974 | 52.15th |

None of these conventions admits the complete deposited split because `basophil` remains ineligible. The sensitivity demonstrates that conclusions about whether the released score was "high" relative to random splits depend on the historical edge definition, even though the broader conclusion of pervasive overlap does not.

## 9. The 107 leaf-layer sensitivity analysis

The supplied 144 network names map one-to-one to hierarchy nodes. Exactly 107 are leaves and 37 are internal nodes. This matches the official OhmNet description of 107 tissue-specific layers at hierarchy leaves and offers a natural candidate universe for a historical interpretation.

Under raw thresholds within the 107 leaves:

- 69 are training eligible;
- 15 are heldout eligible;
- there are 8,190 labeled choices of two test and two validation tissues;
- the uniform valid space contains 232,101,356,008,762,749,600 splits.

One million heldout-first and one million training-first splits were sampled with seeds 202609046001 and 202609046002.

| Statistic | Leaf heldout first | Leaf training first |
|---|---:|---:|
| mean F1 | 0.987703 | 0.989842 |
| SD | 0.009617 | 0.008017 |
| 5th percentile | 0.967355 | 0.974214 |
| median | 0.990422 | 0.991944 |
| 95th percentile | 0.997637 | 0.998092 |
| fraction F1 >= 0.99 | 52.271% | 61.733% |
| fraction F1 = 1 | 0.248% | 0.323% |
| mean row overlap | 0.970093 | 0.975052 |

The released F1 is at the 92.914th and 89.752nd external percentiles. The released row overlap is at the 95.602nd and 93.198th percentiles.

The released split is not a member of this universe: it uses 10 internal hierarchy layers. The leaf-only analysis therefore does not estimate the chance of obtaining the deposited split. Instead, it shows that including broad internal layers, especially in training, is associated with substantially greater identity coverage than would usually arise from 20 leaf-only training tissues.

## 10. Tissue hierarchy and network containment

### 10.1 Hierarchy representation

The supplied `tissue.edges` is a rooted tree with 219 nodes and 218 edges. Every released network name maps uniquely after case/underscore normalization. Two complementary metrics were calculated:

- **shortest-path distance:** number of tree edges between two tissues;
- **Wu-Palmer similarity:** `2*depth(LCA)/(depth(a)+depth(b))`, using root depth 1.

Wu-Palmer similarity accounts for unequal depths, but both metrics are reported.

### 10.2 Leaf and internal network files

Among the 144 edgelists:

- 107 correspond to leaves;
- 37 correspond to internal hierarchy nodes;
- GraphSAGE selected 14 leaves and 10 internal layers.

The selected internal layers are:

```
adrenal_gland, basal_ganglion, blood, bone, brain, eye,
gastrointestinal_tract, kidney, large_intestine, liver
```

Nine are training layers; kidney is validation.

For all 167 ancestor-descendant pairs where both hierarchy nodes have edgelists:

- the descendant node set is a subset of the ancestor node set;
- the descendant edge set is a subset of the ancestor edge set.

For every internal network, the union of its descendant leaf edges is contained in the internal network. The internal network also contains additional edges, so it is not merely the exact union of released descendant leaves.

### 10.3 Actual split hierarchy relationships

Relevant paths are:

```
heart:    Root > CardiovascularSystem > Heart
kidney:   Root > UrogenitalSystem > UrinarySystem > Kidney
lung:     Root > RespiratorySystem > Lung
midbrain: Root > NervousSystem > CentralNervousSystem > Brain > BrainStem > Midbrain
brain:    Root > NervousSystem > CentralNervousSystem > Brain
```

The two test tissues are separated by seven hierarchy edges and have Wu-Palmer similarity 0.2222. Relative to all fixed-24 allocations, that test-pair distance is at the 71.92nd percentile and its similarity at the 27.72nd percentile. Relative to all-144 raw-threshold simulations, the distance is at about the 66.8th percentile and the similarity at about the 23.8th percentile. Thus the two test tissues are themselves relatively dissimilar.

Their relationship to training is mixed:

| Test tissue | Best hierarchy match in training | Best Wu-Palmer | Key gene-set fact |
|---|---|---:|---|
| lung | gastrointestinal_tract by WUP (0.4); brain covers most genes | 0.4 | 92.875% of lung genes occur in training brain alone |
| midbrain | brain | 0.8 | 100% of midbrain genes and edges are contained in training brain |

`midbrain` is a descendant network of `brain`, and its complete 2,310-node set is also covered by `forebrain` despite forebrain being on a different branch below Brain. This explains its perfect lookup score without requiring close similarity for every test tissue.

Across all four heldout tissues, the actual mean heldout-to-training hierarchy distance is relatively small in the all-144 context (about the 15th percentile), while mean Wu-Palmer similarity is mildly above median (about the 59th to 60th percentile). In the leaf-only context, the same actual configuration is extreme: heldout-to-training Wu-Palmer similarity is around the 98.5th percentile. Again, this is because the released training set includes internal layers that do not exist in the leaf-only universe.

### 10.4 Hierarchy similarity versus gene similarity

Across all 10,296 pairs of the 144 network files:

- Spearman rho(Wu-Palmer, gene Jaccard) = 0.1849;
- Spearman rho(Wu-Palmer, smaller-set gene coverage) = 0.3009.

Among the 36 raw-large networks:

- rho(Wu-Palmer, gene Jaccard) = 0.1134;
- rho(Wu-Palmer, smaller-set coverage) = 0.3512.

The hierarchy captures some biological organization, but gene sets overlap broadly even across distant top-level branches. Median smaller-set coverage for different top-level branches is 0.9227. That broad reuse limits the hierarchy's ability to predict leakage by itself.

In the all-144 split simulations:

| Relationship | Heldout first Spearman rho | Training first Spearman rho |
|---|---:|---:|
| F1 vs test-row overlap | 0.9942 | 0.9942 |
| F1 vs unique-gene overlap | 0.9915 | 0.9917 |
| F1 vs mean test-training WUP | 0.1053 | 0.1171 |
| F1 vs best test-training WUP | 0.2228 | 0.2330 |

Hierarchy similarity is therefore a weak, secondary predictor; direct entity overlap is the operative quantity.

### 10.5 Ancestor-descendant and containment effects

In the all-144 heldout-first simulation:

- 39.72% of splits have neither test tissue in an ancestor-descendant relation with training;
- 46.91% have one related test tissue;
- 13.37% have both related.

But exact gene-set containment is more decisive:

- 62.13% have neither test gene set fully contained in the training union;
- 32.01% have one contained;
- 5.86% have both contained.

Whenever both test gene sets are contained in the training union, lookup F1 is exactly 1.0 by construction. Hierarchy relation is neither necessary nor sufficient for full containment because broad gene reuse crosses hierarchy branches.

### 10.6 Internal training layers

With the actual four heldouts fixed and raw-threshold training sets sampled, the expected number of internal hierarchy layers among 20 training networks is 6.927. The released training set contains nine, at the 85.35th midrank percentile of this counterfactual distribution.

Mean F1 increases descriptively with internal-layer count:

| Internal training layers | Mean F1 |
|---:|---:|
| 4 | 0.992952 |
| 6 | 0.994225 |
| 7 | 0.994751 |
| 8 | 0.995209 |
| 9 | 0.995625 |
| 10 | 0.995952 |
| 12 | 0.996568 |

This is an association, not an isolated causal effect: internal layers are also larger and differ in identity. It nevertheless supports the mechanistic interpretation that broad aggregate layers in training increase coverage of heldout genes.

## 11. Which tissue choices tend to lower or raise lookup F1

In the all-144 raw-threshold simulations, test-tissue marginal means were lowest for broad networks such as `fetus`, `blood`, `nervous_system`, `central_nervous_system`, and `brain`. These networks contain large and comparatively difficult-to-cover gene sets. Marginal means were highest for more specific or heavily nested tissues such as `corpus_striatum`, `basal_ganglion`, `temporal_lobe`, `hippocampus`, and `cerebral_cortex`.

This pattern explains why test identity dominates variance. Choosing a broad aggregate network for testing exposes many genes that may not be represented by 20 training layers; choosing a specific nested network often makes the test set nearly or completely contained in training.

The statement is descriptive under the sampled procedures. It should not be interpreted as a claim about biological difficulty or model generalization independent of the lookup diagnostic.

## 12. Monte Carlo and computational quality control

### 12.1 Sampling implementation

The primary one-million-draw simulations used a compact C++ sampler with:

- xoshiro256** pseudorandom generation;
- fixed 64-bit seeds;
- rejection-based unbiased bounded selection;
- direct bitset union and membership calculations;
- fixed packed 76-byte output records;
- no dependence on platform hash iteration order.

The raw sample binaries are not included in the compact bundle because they total hundreds of megabytes. Their filenames, sizes, seeds, formats, and SHA-256 hashes are retained in `sampling/RAW_SAMPLE_MANIFEST.csv`; source code and generated data headers are included so they can be regenerated.

### 12.2 Deterministic reruns

The two primary all-144 sample files and their extrema files were regenerated from the same seeds and matched byte-for-byte. The two leaf-107 samples were independently rerun and also matched byte-for-byte. Hashes are in:

- `results/sampling_deterministic_rerun_sha256.txt`
- `results/leaf107_deterministic_rerun_sha256.txt`

### 12.3 Independent record verification

A separate Python implementation reconstructed the first 1,000 records from each primary sample and verified all integer, tissue-identity, overlap, F1, and hierarchy fields; floating-point hierarchy summaries agreed within the expected float32 representation. Both the 144-layer and 107-leaf checks passed.

### 12.4 Uniformity and convergence checks

For the heldout-first samples, all heldouts were distinct and in the required pool. Test-pair and test-tissue marginal count chi-square tests showed no evidence against the expected uniformity. The same checks were performed for the training-first draws conditional on their mechanism.

Each one-million sample was divided into 20 consecutive blocks of 50,000. Means and key tail probabilities were stable across blocks. The standard error of mean F1 was about 0.000006 in the all-144 samples and at most about 0.000010 in the leaf-only samples. For a probability near 0.3, one million draws gives a binomial standard error near 0.00046, or about 0.09 percentage points for an approximate 95% interval.

Sample minima and maxima are observed sample extrema, not proven global extrema. Only the 63,756 fixed-24 analysis is exhaustive.

## 13. Conclusions

1. **The historical split-generation rule is still open.** The literal published thresholds do not admit the deposited split from the currently released OhmNet files under any ordinary edge count.

2. **The leakage phenomenon is not fragile.** Across all plausible universes, most test rows have a GeneID already present in training and the identity-lookup diagnostic is usually near perfect.

3. **Graph choice still matters.** Under the all-144 raw-threshold models, the sampled F1 range extends from about 0.92 to 1.0, and broad test layers produce the lowest coverage. The released value is moderately high in those contexts, but not exceptional.

4. **Within the deposited 24, the released split is actually less leakage-favorable than most alternatives.** Its F1 is only at the 31st percentile, largely because heart and kidney are withheld for validation.

5. **The candidate universe changes the interpretation.** Relative to leaf-only splits, the released F1 is near the 90th to 93rd percentile, but the released split uses 10 internal hierarchy layers and therefore is not a member of that universe.

6. **The test designation "unseen graphs" is formally true at the component-name level but incomplete at the biological-entity and topology levels.** The underlying genes are overwhelmingly reused, and the entire midbrain test network is an ancestor-contained subnetwork of the training brain layer.

7. **Hierarchy distance is not the main leakage mechanism.** It weakly predicts overlap at best. The dominant causes are repeated Entrez entities, broad cross-tissue gene activity, nested internal layers, and the choice of broad versus specific test networks.

## 14. Recommended report wording

A concise, defensible summary for the later project report is:

> The GraphSAGE supplement describes a random selection of 20 training tissue networks above a 15,000-edge threshold and four larger heldout networks above 35,000 edges, but the deposited split is incompatible with those thresholds in the currently released OhmNet files. Across exact and sampled alternative splits, train-test GeneID reuse remains pervasive: a forensic GeneID-lookup diagnostic is usually near perfect, although its value varies with the test tissues and is especially elevated when broad or ancestor networks occur in training; in the released split, the entire midbrain test network is contained in the training brain layer.

## 15. File guide

Start with:

- `EXECUTIVE_SUMMARY.md` - compact findings.
- `results/headline_results.csv` - primary distribution summary.
- `results/criteria_interpretations.csv` - exact definition and status of each universe.
- `results/graphsage_ohmnet_topology_validation.csv` - exact 24-graph mapping.
- `results/fixed24_all_63756_splits.csv.gz` - complete finite enumeration.
- `results/full144_sampling_summary.json` - all-144 sample distributions.
- `results/leaf107_sampling_summary.json` - leaf-only sample distributions.
- `results/hierarchy_layer_and_containment_validation.json` - leaf/internal and containment facts.
- `results/full144_metric_correlations.csv` - overlap/hierarchy/F1 relationships.
- `CLAIMS.csv` - claim-status and caveat register.
- `SHA256SUMS.txt` - integrity inventory for the compact bundle.

# Hierarchy-aware train-validation-test splitting for the GraphSAGE PPI benchmark

## Purpose

This investigation asks whether the GraphSAGE PPI graphs could have been split more rigorously by using the OhmNet tissue hierarchy, and whether ignoring that hierarchy increases reuse of the same biological nodes across training and heldout graphs.

The GraphSAGE supplement states that 20 PPI networks with at least 15,000 edges were selected for training, while four networks with at least 35,000 edges were selected for validation and testing, two for each role. The supplied OhmNet documentation states that network nodes are Entrez GeneIDs and that the same node ID in two tissue layers is the same gene. These two facts make it possible to define a hierarchy-aware split without using gene identities during split selection, and then measure biological identity reuse after the split is fixed.

The investigation has two distinct goals:

1. Construct defensible train-validation-test splits that minimize tissue-hierarchy similarity while respecting the manuscript's stated graph-size thresholds.
2. Compare hierarchy-aware and hierarchy-blind procedures under matched constraints to determine whether hierarchy awareness reduces biological node reuse and the GeneID-lookup diagnostic F1.

The second goal is the stronger inferential question. A single globally optimized split is useful as a stress test, but a paired ensemble is more rigorous because it separates the effect of the training-selection rule from the identity of the heldout tissues.

## Source-derived constraints and retained discrepancies

### Manuscript split criteria

The GraphSAGE supplement states:

- 20 training PPI networks, each with at least 15,000 edges;
- four large heldout networks, each with at least 35,000 edges;
- two validation and two test networks;
- one fixed random training set and one fixed validation pair across model development.

This investigation interprets an edge as one record in the released undirected OhmNet edgelist. That is the most literal interpretation of the supplied files and produces a nonempty valid universe.

### OhmNet node semantics

The OhmNet README states that each tissue edgelist is a tissue-specific human PPI layer, nodes are Entrez GeneIDs, and the same node ID in different edgelists represents the same gene. Therefore, overlap of Entrez IDs across split roles is biological entity reuse, not merely an accidental collision of anonymous identifiers.

### Threshold discrepancy retained for the project report

The deposited GraphSAGE split is not a member of the literal threshold-defined universe in the current OhmNet release:

| Tissue | Deposited role | Raw edge records | Stated cutoff |
|---|---:|---:|---:|
| astrocyte | training | 9,635 | at least 15,000 |
| basophil | training | 4,193 | at least 15,000 |
| midbrain | test | 31,665 | at least 35,000 |

This is a confirmed manuscript/data discrepancy for the current released files. It does not prove what happened historically: an earlier network release, a different counting stage, an undocumented exception, or imprecise reporting remain possible.

### Midbrain containment retained for the project report

The deposited midbrain test layer is a descendant of the brain training layer in the OhmNet hierarchy. Its complete 2,310-node and edge set is contained in brain, and all 2,310 midbrain tissue-instance rows have GeneIDs present in training. This is more informative than saying only that the graphs are separate components: the test graph is nominally unseen but is not biologically or topologically independent of training.

## Why the primary analysis uses hierarchy leaves

The current OhmNet release contains 144 edgelists mapped to the supplied hierarchy:

- 107 leaf layers;
- 37 internal hierarchy layers.

Internal layers are broad ancestors such as brain, blood, or kidney. Prior exact comparison found that every released descendant layer is a node-and-edge subnetwork of its released ancestor layer. If an internal layer is placed in training and one of its descendants is held out, the hierarchy-aware objective is defeated by construction and topology can be nested across the split.

The primary clean universe therefore contains only hierarchy leaves. This is a methodological restriction, not a claim about the undocumented historical selection process. A separate all-144 sensitivity analysis allows internal layers to mirror the released collection more closely.

Under literal raw-edgelist thresholds:

| Universe | Training-eligible layers (>=15,000) | Heldout-eligible layers (>=35,000) |
|---|---:|---:|
| All 144 layers | 105 | 36 |
| 107 hierarchy leaves only | 69 | 15 |

## Hierarchy metrics

The supplied hierarchy is treated as a rooted tree.

For tissues x and y, Wu-Palmer similarity is calculated from their depths and lowest common ancestor:

    WUP(x,y) = 2 * (depth(LCA(x,y)) + 1)
               / ((depth(x) + 1) + (depth(y) + 1))

The added one gives the root a nonzero contribution and matches the implementation used throughout this analysis. Higher values indicate greater hierarchy similarity.

The sensitivity metric is ordinary tree-path distance in edges:

    distance(x,y) = depth(x) + depth(y) - 2 * depth(LCA(x,y))

The analysis records:

- mean train-heldout WUP similarity over all 20 x 4 tissue pairs;
- maximum train-heldout WUP similarity, which protects against one severe near-neighbor relationship;
- mean and minimum path distance;
- ancestor-descendant pair count;
- the closest training tissue for each heldout tissue.

A minimax objective is preferable to mean similarity alone. The released split illustrates why: an average can look ordinary while a single test tissue, midbrain, is directly nested inside training brain.

## Split-construction principles

A rigorous split procedure should satisfy six principles.

1. **No outcome leakage during split selection.** Tissue hierarchy, leaf/internal status, graph sizes, and the manuscript thresholds may be used. Gene identities, labels, node overlap, and lookup F1 may not be used until after a split is fixed.
2. **Validation and test receive equal protection.** Validation influences model and hyperparameter selection, so both heldout roles must be distant from training.
3. **Worst-case protection.** Minimize the closest train-heldout hierarchy relationship before optimizing averages.
4. **Block broad hierarchy branches.** A heldout tissue's entire root branch is excluded from training, not merely the exact tissue name.
5. **Control graph size.** Training graph sizes are stratified so hierarchy optimization cannot win primarily by choosing unusually small graphs.
6. **Use paired nulls.** Compare optimized and random training sets for the same heldout quartet under the same non-hierarchy constraints.

## Alternative valid formulations considered

Several procedures satisfy the manuscript thresholds but answer somewhat different questions:

| Procedure | Additional rule | Strength | Main weakness |
|---|---|---|---|
| Threshold-only random | No hierarchy rule | Closest current-release analogue to the manuscript prose | Can place parents, descendants, or close siblings across roles |
| Leaf-only random | Exclude internal hierarchy layers | Removes direct ancestor-layer supersets | Still ignores distances among leaves |
| Ancestor-blocked | Forbid direct ancestor-descendant train-heldout pairs | Prevents the most obvious containment | Close siblings can remain; internal layers still broaden gene coverage |
| Distinct-root-branch | Put heldouts in four root branches and exclude those branches from training | Strong domain blocking and simple interpretation | May concentrate training into a few remaining branches |
| Mean-WUP optimization | Minimize average train-heldout similarity | Maximizes overall separation | Can hide one very close train-heldout pair |
| Minimax-WUP optimization | Minimize the single closest train-heldout relationship first | Protects every heldout graph | May sacrifice average separation or domain coverage |
| Path-distance optimization | Maximize minimum or average tree distance | Metric sensitivity independent of WUP scaling | Shares the hierarchy's biological limitations |
| Size-stratified hierarchy optimization | Balance training node counts while optimizing hierarchy | Reduces graph-size confounding | Does not force broad training-domain coverage |
| Broad-coverage hierarchy optimization | Represent every remaining root branch in training | Maintains training diversity | A broad union of tissues can recreate high GeneID overlap |

The primary method combines leaves, root-branch blocking, size stratification, and minimax WUP. The all-144, mean-distance, maximin-distance, and broad-coverage methods are reported as sensitivities rather than silently blended into one score.

## Primary hierarchy-aware algorithm

The confirmatory construction is called `branch_distinct_node_stratified_minimax_wup`.

### Step 1: choose the four heldout tissues

Choose four threshold-eligible hierarchy leaves, each from a different top-level branch immediately below the hierarchy root.

Of the 15 heldout-eligible leaves, there are:

- C(15,4) = 1,365 possible heldout quartets;
- 1,062 quartets satisfying the four-distinct-root-branches rule and leaving enough training candidates in every size stratum.

All 1,062 feasible quartets were evaluated exactly.

### Step 2: block heldout branches from training

Exclude every training candidate in any of the four heldout root branches. This prevents selection of a parent, sibling, or close cousin from the same broad branch.

### Step 3: balance training graph sizes

Sort the 69 training-eligible leaves by node count and divide them deterministically into five near-equal strata. Select exactly four training graphs from each stratum, for 20 total.

This controls node-count composition while retaining the manuscript edge thresholds. Edge totals are recorded and additionally controlled in a post hoc nearest-size sensitivity analysis.

### Step 4: optimize the training set lexicographically

For each fixed heldout quartet, select the training set that:

1. minimizes the largest WUP similarity among all 80 train-heldout pairs;
2. among ties, minimizes the total WUP similarity;
3. among remaining ties, maximizes total tree-path distance;
4. uses stable tissue-index order only as a final deterministic tie-break.

Because the size-stratum constraints are separable, this optimum was enumerated exactly. An independent implementation recomputed all 1,062 primary rows with zero failures.

### Step 5: partition heldouts into validation and test pairs

There are three unordered 2+2 partitions of each heldout quartet. Select the partition that, using hierarchy information only:

1. minimizes the worse of the two train-to-role mean WUP similarities;
2. minimizes the difference between validation and test similarity;
3. minimizes validation-test cross-similarity;
4. uses stable IDs only as the final tie-break.

The two resulting pairs are symmetric. For a released benchmark, pair A can be called validation and pair B test by a pre-registered lexical rule or a recorded fair coin flip. Their names must not be assigned after seeing model performance.

## Hierarchy-blind comparison procedures

### Manuscript-like uniform null

For each universe, select a heldout quartet uniformly from all threshold-eligible quartets, then choose 20 training layers uniformly from the remaining training-eligible pool. This is the closest current-release analogue to an edge-threshold-only random split.

### Matched hierarchy-blind null

Use exactly the same leaf restriction, four-distinct-heldout-branch rule, heldout-branch blocking, and four-per-node-size-stratum constraints as the hierarchy-aware algorithm, but select the 20 training graphs randomly within strata.

This is the primary comparator. Any difference is attributable to hierarchy-based training selection rather than to different eligibility, branch blocking, or node-size composition.

### Conditional paired comparison

For every one of the 1,062 feasible heldout quartets:

- retain the exact hierarchy-optimized training set;
- average the matched-random outcomes observed for that same quartet;
- subtract matched random from optimized.

The one-million matched-null stream yields a median of 942 random training sets per heldout quartet, with 842 to 1,039 per quartet. This conditional analysis controls heldout-tissue identity.

### Tight size-matched sensitivity

For each heldout quartet, compare the hierarchy-optimal training set with the 25, 50, 100, or 200 matched-random training sets nearest in standardized total node and edge count. At K=25, mean training-size differences are only -12.5 nodes and +137.5 edges, making residual size confounding very small.

## Exact and sampled search scope

Exact enumeration covered:

| Universe | All threshold-eligible heldout quartets | Feasible branch-distinct quartets |
|---|---:|---:|
| All 144 layers | 58,905 | 26,830 |
| Leaf-only | 1,365 | 1,062 |

Four fixed-seed Monte Carlo streams each contained 1,000,000 splits:

| Universe | Null | Seed |
|---|---|---:|
| all 144 | uniform | 2026090401 |
| all 144 | matched stratified | 2026090501 |
| leaf-only | uniform | 2026090402 |
| leaf-only | matched stratified | 2026090502 |

The sampler used xoshiro256** with unbiased bounded selection. Each stream was regenerated byte-for-byte from its recorded seed. An independent implementation recomputed the first 1,000 records of every stream within a maximum absolute difference of 5.3e-7. Twenty consecutive 50,000-split blocks were used to inspect convergence.

## Primary result: hierarchy awareness modestly lowers biological node reuse on average

The strongest comparison is conditional on the same heldout quartet.

### Leaf-only paired effects

The following values are hierarchy-optimized minus matched-random means across all 1,062 heldout quartets:

| Metric | Mean paired difference | Fraction with optimized value lower |
|---|---:|---:|
| Mean WUP similarity | -0.03058 | 100.0% |
| Maximum WUP similarity | -0.03535 | 96.4% |
| Heldout row overlap | -0.00553 | 69.6% |
| Unique-gene overlap | -0.01140 | 69.7% |
| GeneID-lookup micro-F1 | -0.00216 | 65.2% |

Thus hierarchy optimization always improved its intended average hierarchy objective, and usually reduced biological identity reuse, but did not do so for every heldout quartet.

The modest Spearman association between the change in hierarchy similarity and the change in lookup F1 is 0.279 (p = 1.84e-20). Since both improvements are negative differences, larger hierarchy reductions tend to accompany larger F1 reductions, but hierarchy distance explains only part of the identity-overlap behavior.

### Tight size-matched result

At K=25 nearest random controls per heldout quartet:

| Metric | Mean paired difference | Fraction with optimized value lower |
|---|---:|---:|
| Mean WUP similarity | -0.02983 | 100.0% |
| Heldout row overlap | -0.00349 | 66.8% |
| Unique-gene overlap | -0.00699 | 66.7% |
| GeneID-lookup micro-F1 | -0.00134 | 62.9% |
| Training node total | -12.5 | 74.6% |
| Training edge total | +137.5 | 31.4% |

The effect attenuates but remains when node and edge totals are nearly identical. This supports a real algorithmic association between hierarchy-aware selection and reduced entity reuse, while also showing that part of the unadjusted difference was related to graph composition.

### All-144 sensitivity

Allowing internal hierarchy layers gives the same qualitative result across 26,830 feasible heldout quartets:

- mean WUP difference: -0.03700;
- row-overlap difference: -0.00294;
- unique-overlap difference: -0.00690;
- lookup-F1 difference: -0.00118;
- optimized F1 lower for 67.9% of heldout quartets.

The conclusion therefore does not depend on the leaf-only restriction, although the all-144 universe permits nested ancestor layers and is less clean biologically.

## Distribution-level comparison

For the primary leaf-only matched null, one million splits produced:

| Metric | Mean | 5th percentile | Median | 95th percentile |
|---|---:|---:|---:|---:|
| Mean WUP | 0.23369 | 0.21308 | 0.23380 | 0.25426 |
| Heldout row overlap | 0.97563 | 0.95323 | 0.97820 | 0.98994 |
| Unique-gene overlap | 0.93976 | 0.89544 | 0.94393 | 0.97194 |
| Lookup micro-F1 | 0.99018 | 0.98088 | 0.99121 | 0.99606 |

Across the 1,062 exact hierarchy-optimized training sets, the means were:

- mean WUP: 0.20309;
- heldout row overlap: 0.97010;
- unique-gene overlap: 0.92834;
- lookup micro-F1: 0.98802.

These ensemble means agree with the paired conclusion: hierarchy awareness shifts the distributions toward less similarity and less identity reuse, but the distributions overlap substantially.

## A single maximum-separation stress-test split

The global primary optimum is useful when one particularly difficult hierarchy-separated split is needed. It was selected without using genes or labels.

### Heldout quartet

- colon
- fetus
- hematopoietic_stem_cell
- pancreas

Balanced hierarchy-only role partition:

- pair A: colon, fetus
- pair B: hematopoietic_stem_cell, pancreas

### Training tissues

- aorta
- b_lymphocyte
- caudate_nucleus
- corpus_callosum
- forebrain
- frontal_lobe
- hypophysis
- medulla_oblongata
- midbrain
- myometrium
- prostate_gland
- renal_glomerulus
- spermatid
- temporal_lobe
- testis
- thalamus
- thymocyte
- umbilical_cord
- uterine_cervix
- uterine_endometrium

### Metrics

| Metric | Value |
|---|---:|
| Mean train-heldout WUP | 0.172609 |
| Maximum train-heldout WUP | 0.200000 |
| Mean path distance | 9.70 |
| Minimum path distance | 8 |
| Ancestor-descendant pairs | 0 |
| Training nodes | 40,713 |
| Training edge records | 553,047 |
| Heldout row overlap | 0.933361 |
| Unique-gene overlap | 0.854394 |
| Combined lookup micro-F1 | 0.972114 |
| Pair A row overlap / F1 | 0.929504 / 0.971108 |
| Pair B row overlap / F1 | 0.937632 / 0.973205 |

No one of the one million matched-null draws had mean WUP as low; the standard plus-one Monte Carlo estimate is at most approximately 1.0e-6. Its row overlap and F1 are each near the 0.86th percentile of the matched null.

The same exact split is selected when the primary objective is changed from WUP minimax to either maximum minimum path distance or maximum mean path distance. This metric agreement is strong evidence that the separation is not an artifact of one hierarchy formula.

### Limitation of the global optimum

The 20 training tissues occupy only four root branches. This is a legitimate hierarchy optimum but a potentially unrealistic benchmark: it obtains separation partly by concentrating the training domain. It should be treated as a maximum-separation stress test, not as the only recommended operational split.

## Broad-coverage sensitivity

A more domain-diverse construction requires at least one training leaf from every non-heldout eligible root branch while preserving the four-per-size-stratum constraint. The objective then minimizes total/mean WUP, followed by worst WUP and path distance.

One of two tied hierarchy optima has heldouts:

- colon
- pancreas
- prostate_gland
- temporal_lobe

Its training set spans all 14 remaining eligible root branches. It has:

| Metric | Value |
|---|---:|
| Mean train-heldout WUP | 0.194152 |
| Maximum train-heldout WUP | 0.250000 |
| Mean path distance | 8.60 |
| Minimum path distance | 6 |
| Training nodes | 42,350 |
| Training edge records | 575,499 |
| Heldout row overlap | 0.989476 |
| Unique-gene overlap | 0.971563 |
| Lookup micro-F1 | 0.995478 |

Hierarchy separation remains extreme: only 143 of one million matched-null splits had mean WUP at or below this value. Yet row overlap is at the 93.8th percentile and F1 at the 91.0th percentile of that null.

This is a decisive counterexample to the idea that hierarchy distance alone guarantees biological independence. A training set can span many distant branches and collectively cover almost every heldout gene.

The second tied hierarchy optimum replaces prostate_gland with testis in the heldout quartet. It has the same hierarchy and training-size metrics but row overlap 0.984102 and F1 0.993282. Even tied hierarchy optima can have materially different biological overlap, reinforcing the requirement not to select among ties using gene or label outcomes.

## Context for the deposited GraphSAGE split

Because the deposited split violates the current-release thresholds, it is an external comparison rather than a possible draw from the literal valid universe.

For its four heldout tissues (heart, kidney, lung, midbrain):

| Metric | Value |
|---|---:|
| Mean train-heldout WUP | 0.281429 |
| Maximum train-heldout WUP | 0.800000 |
| Minimum path distance | 2 |
| Ancestor-descendant pairs | 1 |
| Heldout row overlap | 0.988786 |
| Unique-gene overlap | 0.972018 |
| Combined lookup micro-F1 | 0.995307 |

For the actual two test tissues alone (lung and midbrain):

- row overlap: 0.993845;
- unique-gene overlap: 0.989593;
- lookup micro-F1: 0.997178;
- all 2,310 midbrain rows are seen in training.

The released split is not globally exceptional by average WUP or average overlap when compared externally with all-144 random draws. Its main concern is local and severe: the maximum WUP of 0.8 and the exact brain-to-midbrain containment relationship. This is why worst-case metrics and explicit ancestor blocking are necessary.

## Interpretation

### Answer to the main question

Under matched thresholds, leaf status, branch blocking, and graph-size strata, hierarchy-aware training selection causes less biological node reuse on average than hierarchy-blind random selection. The mean reduction is modest after tight size matching: about 0.35 percentage points in heldout row overlap and 0.13 percentage points in lookup micro-F1.

The effect is not universal. Roughly one third of heldout quartets show equal or greater identity reuse after hierarchy optimization, and the broad-coverage counterexample shows that very distant tissues can still share nearly all genes.

### Why the relationship is limited

The OhmNet hierarchy encodes anatomical or cell-type relatedness, while gene reuse is driven by several additional processes:

- widely active genes occur in many tissues;
- broad training coverage accumulates a large union of genes;
- some internal layers are exact supersets of descendants;
- graph size affects how many distinct genes can be covered;
- hierarchy-distant tissues can still have high gene-set overlap.

Hierarchy is therefore a useful pre-outcome design variable but an imperfect proxy for entity disjointness.

### What a hierarchy-aware benchmark can and cannot claim

A hierarchy-aware split can claim stronger separation of tissue domains and avoidance of direct ancestor-descendant leakage. It cannot claim that heldout nodes are unseen biological entities unless GeneID disjointness is enforced directly.

If the scientific target is generalization to new tissues containing many familiar proteins, hierarchy-aware graph splitting is appropriate. If the target is generalization to new proteins, a gene-disjoint or grouped-by-GeneID split is required. These are different tasks and should be reported separately.

## Recommended experimental design

The most rigorous assessment is a paired multi-split experiment, not a comparison of one favorable split with the deposited split.

1. Freeze the literal edge-count definition, the leaf-only primary universe, hierarchy metric, root-branch definition, graph-size strata, and tie-breaking rules before inspecting model outcomes.
2. Enumerate all 1,062 feasible heldout quartets.
3. For every quartet, generate the exact hierarchy-optimal training set and multiple matched hierarchy-blind training sets.
4. Assign the two heldout role pairs using hierarchy only, then fix validation/test orientation by a recorded rule.
5. Train the same model and hyperparameter procedure on paired splits with identical random seeds and compute budgets.
6. Report paired distributions of model F1, node overlap, graph sizes, and hierarchy metrics, including confidence intervals and the fraction of quartets favoring each procedure.
7. Include the maximum-separation split as a stress test and the broad-coverage split as a domain-diversity sensitivity analysis.
8. Keep the all-144 analysis as a release-mirroring sensitivity, while clearly labeling its ancestor-layer contamination risk.
9. Never choose a split or resolve a hierarchy tie based on gene overlap, lookup F1, or downstream model performance.
10. Report a gene-disjoint evaluation separately if claims concern unseen proteins rather than unseen tissue graphs.

For a manageable model-training study, a pre-registered stratified sample of heldout quartets can replace all 1,062, but the hierarchy-aware and matched-random sets must remain paired on the same quartet. The exact enumeration tables in this bundle permit such a sample to be drawn reproducibly without rerunning the optimization.

## Lookup diagnostic definition and caveat

The GeneID-lookup predictor memorizes each training GeneID's multilabel vector and applies it to the same GeneID in heldout tissues; an unseen GeneID receives the zero vector. Micro-F1 is computed over all heldout label cells. It is a diagnostic of entity reuse, not a trained GraphSAGE result.

The 144 OhmNet layers contain 4,510 genes, 209 more than the 4,301-gene deposited GraphSAGE universe. Label weights for those 209 genes were generated by the same recovered global GO transformation that exactly reproduces the deposited labels for the observed universe. Their hypothetical labels are model-derived rather than independently validated against a GraphSAGE target. Node-count, row-overlap, and unique-gene-overlap results do not depend on this extension. A prior sensitivity that assigned the extra genes zero positive-label weight changed mean F1 by only about 0.0007.

## Quality control

The following checks passed:

- independent recomputation of every one of the 1,062 primary leaf rows: zero failures;
- independent recomputation of all reported exact global optima: maximum numeric discrepancy below 5e-13;
- deterministic regeneration of all four one-million-split streams: byte-identical SHA-256 values;
- independent recomputation of the first 1,000 records in every stream: zero failures at tolerance 2e-6;
- conditional matched-null coverage: 842 to 1,039 draws per heldout quartet;
- blockwise convergence inspection across twenty 50,000-split blocks;
- fixed source-file hashes and an append-only claims/discrepancy register.

The optimizer used no gene identity or label information. Biological overlap was evaluated only after each hierarchy-defined split was fixed.

## Scope boundaries

This investigation does not recover the authors' historical random split or prove which edge-count interpretation they used. It does not estimate GraphSAGE neural-model performance under the new splits. It quantifies hierarchy separation, biological entity reuse, and the lookup diagnostic under explicit current-release procedures.

## Principal conclusion

A hierarchy-aware split is more rigorous than an edge-threshold-only random split because it can prevent direct ancestor-descendant relationships and reduce biological node reuse on average. However, the reduction is modest under tight controls and is not guaranteed: hierarchy distance and GeneID overlap are related but distinct. The defensible benchmark design is therefore a pre-registered paired ensemble of hierarchy-aware and matched hierarchy-blind splits, supplemented by direct GeneID-disjoint evaluation when the scientific claim concerns unseen proteins.

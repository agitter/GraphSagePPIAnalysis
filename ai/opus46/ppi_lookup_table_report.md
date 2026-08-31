# The GraphSAGE PPI Benchmark Is a Lookup Table

## Summary

We recovered the Entrez gene IDs stripped from the GraphSAGE PPI benchmark by
matching its 24 graphs back to the original OhmNet tissue-specific PPI networks.
This reveals that the benchmark's "inductive" train/test split partitions
*tissue graphs*, not *genes*. The same ~4,500 genes appear across all 24 tissues
with identical features and identical labels. A zero-parameter gene lookup table
scores 0.9956 micro-F1 on the full test set, matching or exceeding published GNN
results that required millions of parameters.

## Recovering gene identity

### Source data

| File | Source | Role |
|---|---|---|
| `ppi-G.json` | GraphSAGE release (Hamilton et al. 2017) | Graph structure (56,944 nodes, 818,716 edges) |
| `ppi-feats.npy` | GraphSAGE release | 50-dim binary features per node |
| `ppi-class_map.json` | GraphSAGE release | 121-dim binary labels per node |
| `train_graph_id.npy` | GraphSAGE release | Graph membership for all 56,944 nodes |
| `bio-tissue-networks/*.edgelist` | OhmNet / BioSNAP (Zitnik & Leskovec 2017) | 144 tissue PPI networks with Entrez gene IDs |

The GraphSAGE release contains no gene identifiers. The `ppi-id_map.json` file
is the identity permutation (`"0"→0, "1"→1, …`). Issue #188 on the GraphSAGE
repository asks for a gene mapping; it remains unanswered.

### Step 1: Match graphs to tissues (24/24)

GraphSAGE preprocessing removed small connected components from OhmNet networks
and retained only the largest connected component (LCC). Matching GraphSAGE
subgraphs to OhmNet LCCs by `(nodes, edges)` count gives **24 unique matches,
all confirmed by degree-sequence comparison**.

| Graph | Split | Tissue | Graph | Split | Tissue |
|---|---|---|---|---|---|
| g1 | train | adipose_tissue | g13 | train | bone |
| g2 | train | adrenal_cortex | g14 | train | brain |
| g3 | train | adrenal_gland | g15 | train | colon |
| g4 | train | amygdala | g16 | train | eye |
| g5 | train | aorta | g17 | train | forebrain |
| g6 | train | astrocyte | g18 | train | large_intestine |
| g7 | train | artery | g19 | train | liver |
| g8 | train | basal_ganglion | g20 | train | gastrointestinal_tract |
| g9 | train | basophil | g21 | **val** | **heart** |
| g10 | train | blood | g22 | **val** | **kidney** |
| g11 | train | blood_plasma | g23 | **test** | **lung** |
| g12 | train | blood_platelet | g24 | **test** | **midbrain** |

### Step 2: Map node indices to Entrez gene IDs (99.12%)

Each GraphSAGE graph and its matched OhmNet LCC are the same graph under a
different node numbering. Node ordering is not preserved (sorted-Entrez-ID
hypothesis gives ~2% positional agreement). We recovered the permutation
using three methods applied in sequence:

1. **WL colour refinement** (6 rounds): nodes with a unique structural
   fingerprint (own degree + sorted multiset of neighbour degrees, iterated)
   map uniquely between the two copies. Resolved **55,874 nodes**.

2. **Feature-vector disambiguation**: within a WL colour class, nodes with
   different 50-dim binary feature vectors must be different genes. Cross-graph
   consistency (the same gene has the same feature vector in every tissue)
   identifies which OhmNet gene matches. Resolved **261 nodes**.

3. **Process of elimination**: after steps 1–2, colour classes with one
   remaining candidate on each side are assigned. Resolved **24 nodes**.

4. **Label-vector disambiguation** (uses test labels): within a WL colour
   class where features are identical but GraphSAGE labels differ, match by
   label consistency across graphs. Resolved **89 additional nodes** in a
   second pass with elimination.

**Total mapped: 56,248 / 56,658 LCC nodes = 99.28%.** The remaining 410 nodes
fall into 6 equivalence classes of biologically identical paralogs (below).

### Step 3: Verify the mapping

**Edge reconstruction**: for every mapped node pair with an edge in GraphSAGE,
the corresponding gene pair has an edge in OhmNet. **818,435 edges verified,
zero mismatches.** A wrong mapping cannot reproduce this.

**External validation**: the highest-degree node in the midbrain test graph maps
to Entrez 351 = APP (amyloid beta precursor protein), a known neural hub. The
second is 7157 = TP53. Both are confirmed at NCBI Gene.

### Step 4: Equivalence classes (the remaining 410 nodes)

410 nodes sit in WL colour classes where all members share the same structural
fingerprint, the same 50-dim feature vector, and the same 121-dim label vector.
These group into 6 equivalence classes of biological paralogs:

| Equivalence class | Protein | Occurrences |
|---|---|---|
| {8362, 8368} | Histone H4 (2 of 14 copies) | 26 nodes in 13 graphs |
| {8362, 8364, 8366, 8368} | Histone H4 (4 copies) | 12 nodes in 3 graphs |
| {8360, 8365} | Histone H4 (2 copies) | 4 nodes in 2 graphs |
| {8360, 8362, 8364, 8366, 8368} | Histone H4 (5 copies) | 5 nodes in 1 graph |
| {4664, 59277} | WNT9A / WNT9B | 2 nodes in 1 graph |
| {4886, 4889} | GATA4 / GATA6 | 2 nodes in 1 graph |

Entrez IDs 8360–8368 are all histone H4 clustered variants mapping to a single
Ensembl gene (ENSG00000158406). They encode the same protein. Members of each
class share **identical non-family neighbour sets** (verified in every graph
where they co-occur), identical features, identical labels, and identical OhmNet
GO annotations (e.g. 8362 vs 8368: 351 common GO annotations, 0 differences).

Any permutation within a class produces identical predictions. The mapping is
therefore **complete and unique at the equivalence-class level**.

## What the gene IDs reveal

### Three gene-level invariants

Using 56,248 uniquely mapped nodes appearing in multiple graphs:

| Property | Agreement |
|---|---|
| Same gene → same 50-dim feature vector across all tissues | 100.00% |
| Same gene → same 121-dim label vector across all tissues | 100.00% |
| Same gene → same non-family edges (within equivalence classes) | 100.00% |

Features and labels are gene-level properties with no tissue specificity. This
is despite the OhmNet source data containing tissue-specific GO annotations
(e.g. `smooth_muscle_GO:0048661.lab`). The GraphSAGE preprocessing discarded
tissue specificity when constructing the 121 label columns.

### Train/test gene overlap

| | lung (test) | midbrain (test) | heart (val) | kidney (val) |
|---|---|---|---|---|
| Mapped genes | 3,167 | 2,276 | 3,195 | 3,284 |
| Also in training | 3,132 (98.9%) | 2,276 (100%) | 3,145 (98.4%) | 3,233 (98.5%) |
| With identical labels | 3,132 (98.9%) | 2,276 (100%) | 3,142 (98.3%) | 3,226 (98.2%) |
| Genuinely unseen genes | 34 (1.1%) | 0 (0%) | 50 (1.6%) | 51 (1.5%) |

An additional 30 test nodes and 35 val nodes sit in small connected components
outside the LCC and are not mapped.

The gene universe across all 144 OhmNet tissues is **4,510 proteins**. The 20
training graphs map to **4,175 distinct genes**. With 12.6 node-instances per
gene on average, the training set covers the vast majority of the gene pool.

### Gene lookup baseline

For each test node, look up its gene in training and copy the 121-dim label
vector. Fall back to majority-frequency prediction for unseen genes and
unmapped nodes.

| Method | Params | Full test micro-F1 | Full val micro-F1 |
|---|---|---|---|
| Constant (all ones) | 0 | 0.4608 | 0.4611 |
| Gene lookup | 0 | **0.9956** | **0.9920** |
| GAT (Veličković et al. 2018) | ~2M | 0.9730 | — |
| Cluster-GCN (Chiang et al. 2019) | ~10M | 0.9936 | — |
| GraphSAINT (Zeng et al. 2020) | ~10M | 0.9950 | — |
| pathGCN (Eliasof et al. 2022) | ~10M | 0.9961 | — |

The gene lookup scores 0.9956 on the full test set (5,524 nodes). On the 97.9%
of test nodes it covers, it scores **1.0000**. The gap to 1.0 on the full set
comes from 34 lung-specific genes not in training and 30 small-component nodes.

The lookup **exceeds** GAT, Cluster-GCN, and GraphSAINT. The small gap between
our 0.9956 and pathGCN's 0.9961 reflects only that pathGCN makes slightly
better guesses on the 64 nodes (1.2%) the lookup cannot help with.

## Conclusion

The PPI benchmark's "inductive" framing describes graph-level disjointness, not
gene-level disjointness. Since the same ~4,500 genes recur across all 24 tissue
graphs with identical features and labels, a model that identifies which gene a
node represents — using neighbourhood structure as a proxy for the stripped
Entrez IDs — and recalls that gene's label vector from training achieves
state-of-the-art performance with zero parameters.

Reported micro-F1 improvements from 0.97 to 0.99 reflect increasing model
capacity to memorise a gene→label lookup table, not advances in graph
representation learning.

## Recommendations

1. **Do not use PPI to evaluate GNN architectures.** Performance is saturated
   and reflects memorisation capacity, not generalisation.

2. **If reporting PPI results**, include the gene lookup baseline (0.9956) as a
   reference. Code and the complete node→gene mapping are provided.

3. **To construct a non-degenerate version**, either:
   - Partition the 4,510 genes (not the graphs) into train/test sets so no
     gene appears on both sides.
   - Restore tissue-specific labels from the OhmNet `bio-tissue-labels` files,
     making predictions depend on tissue context.
   - Use ogbn-proteins (STRING-derived, species-based split) or other
     benchmarks where the split enforces real distribution shift.

## Artefacts provided

| File | Contents |
|---|---|
| `matches_full.json` | Graph index → tissue name (24 entries) |
| `node2gene_complete.json` | Node index → Entrez gene ID (56,658 LCC entries) |
| `recover_all.py` | WL + feature + elimination recovery pipeline |
| `baselines.py` | Gene lookup and other baselines |

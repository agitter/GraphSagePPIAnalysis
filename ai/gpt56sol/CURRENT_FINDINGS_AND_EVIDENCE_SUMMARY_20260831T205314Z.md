# GraphSAGE PPI Provenance Investigation: Current Findings and Evidence Summary

**Snapshot date:** 2026-08-31 15:53 CDT  
**UTC datestamp:** 2026-08-31T20:53:14Z  
**Status:** Internal results summary as of this date; the project is continuing  

## Purpose of this report

This document summarizes what is currently known about the GraphSAGE PPI benchmark, the evidence supporting each major conclusion, the places where the released data differ from the literature or repository descriptions, and the questions that remain open. It is intended to inform a collaborator, not yet to serve as the final fully audited publication.

## 1. Executive summary

The investigation has reconstructed the complete released dataset at the level of observed outputs.

Under a deterministic row-to-Entrez mapping inferred from the original OhmNet edgelist order and legacy CPython 2 dictionary behavior:

- all 56,944 GraphSAGE node rows receive an Entrez GeneID;
- all 24 graph node and edge sets match named OhmNet tissue networks;
- all 2,847,200 feature cells match an MSigDB-derived reconstruction;
- all 6,890,224 label cells match a single GOA-derived reconstruction;
- the DGL archive is reproduced exactly, including graph IDs, standardized features, labels, and directed edge sets.

The released train/test split is graph-wise but not gene-wise. Of 5,524 test rows, 5,490 have a GeneID already present in training. A trivial gene-identity lookup reaches test micro-F1 **0.9971784** without learning from graph structure or features.

Several important distinctions remain:

- the observed outputs are reproduced exactly;
- some historical mechanisms are strongly inferred rather than documented by original source code;
- six GO column identities remain provisionally oriented because three term pairs have identical membership vectors on the entire 4,301-gene graph universe;
- the exact procedure for selecting the 24 tissues and assigning the 20/2/2 split remains unknown;
- the literal upstream label product could have been GOA release 159 or an equivalent Entrez-native/processed source.

## 2. Evidence categories used in this project

| Category | Meaning |
|---|---|
| **Byte-level exact** | Direct file equality, edge-set equality, or complete matrix reconstruction with zero differences. |
| **Strongly inferred** | A historically plausible mechanism explains all available data and survives independent or negative controls, but original source code or an ordered intermediate artifact is absent. |
| **Documented** | Explicitly stated in a manuscript, official archive, source page, or historical software source. |
| **Open** | Multiple histories remain compatible with the observed output, or a decisive intermediate artifact is missing. |

Many conclusions have two parts. For example, the node mapping is an exact table with complete downstream validation, while the claim that the original code literally used the reconstructed Python 2/NetworkX path is strongly inferred.

## 3. Released dataset anatomy

The GraphSAGE PPI archive contains:

| Item | Value |
|---|---:|
| Node rows | 56,944 |
| Stored undirected links | 818,716 |
| Input features | 50 |
| Output labels | 121 |
| Graph blocks | 24 |
| Distinct Entrez GeneIDs across those blocks | 4,301 |
| Training rows | 44,906 |
| Validation rows | 6,514 |
| Test rows | 5,524 |

Each row is a gene occurrence in one tissue graph. The same biological GeneID can therefore occur in many rows and in more than one split.

## 4. Graph topology and external PPI provenance

### 4.1 Immediate graph source

**Current status: byte-level exact.**

All 24 anonymous GraphSAGE graph blocks have exact node-set and edge-set matches to named OhmNet tissue edgelists. Under the complete row mapping, all 818,716 stored GraphSAGE links are accounted for and no link joins different tissue blocks.

This establishes OhmNet as the immediate source of the released GraphSAGE topology.

### 4.2 External source chain

A separate BioSNAP audit showed:

- the union of all 144 supplied OhmNet tissue edgelists contains 70,338 unique undirected gene pairs;
- that union exactly equals the BioSNAP combined OhmNet network;
- all 70,338 pairs occur in the 342,353-edge BioSNAP global human interactome;
- the global BioSNAP resource contains 21,557 nodes.

The supported chain is therefore:

```text
BioSNAP global human interactome
21,557 nodes / 342,353 unique edges
        ↓ contains every tissue edge
BioSNAP combined OhmNet tissue network
70,338 unique undirected pairs
        ↓ exactly equals
union of the 144 supplied OhmNet tissue edgelists
        ↓ GraphSAGE selected 24 tissues
GraphSAGE PPI benchmark topology
```

### 4.3 What should be said about BioGRID and Menche et al.

The topology should not be described as “a BioGRID network.” OhmNet's upstream global network was a composite experimentally supported interactome assembled from several resources; BioGRID was one contributor.

Menche et al. 2015 reports 13,460 proteins and 141,296 interactions, which differs substantially from the OhmNet/BioSNAP 21,557-node, 342,353-edge network. The Menche network is best treated as part of the source or methodological lineage, not as the exact unchanged file used by OhmNet.

Per-edge attribution to individual constituent databases is not available in the flattened released network.

### 4.4 Remaining topology questions

The exact identities of the 24 tissues are resolved, but the rule that selected them from 144 and assigned 20 train, 2 validation, and 2 test graphs is still open. Simple edge-count and “largest 24” rules do not reproduce the selection.

## 5. Complete node-to-Entrez mapping

### 5.1 Reconstructed mechanism

**Output status: complete and exact. Historical mechanism: strongly inferred.**

For each selected OhmNet tissue edgelist, the successful model:

1. reads edges in original line order;
2. leaves node tokens as strings;
3. inserts them into a simulated 64-bit, unrandomized CPython 2.7 dictionary as the graph is built;
4. emits nodes by scanning occupied dictionary slots from low to high.

This deterministic sequence is used as the tissue-local GraphSAGE row order.

### 5.2 Independent evidence

Before using the ordering model, 56,411 row identities had been established independently:

- 55,878 rows from graph topology alone;
- 533 additional rows from topology plus independently identified MSigDB feature membership.

The string-key Python 2 order agreed with **all 56,411** of those identities and disagreed with none. Competing orders performed extremely poorly:

| Candidate node order | Agreements among 56,411 independently anchored rows |
|---|---:|
| First appearance in edgelist | 24 |
| Sorted numeric GeneID | 30 |
| Reverse numeric order | 17 |
| Python 2 dictionary with integer keys | 23 |
| **Python 2 dictionary with string keys** | **56,411** |

The same CPython 2 dictionary simulator also reproduced the full irregular order of all 56,944 keys in `ppi-class_map.json`, providing an independent environment fingerprint.

A second implementation reproduced all 56,944 row assignments exactly.

### 5.3 Mapping product

The canonical table records, for each row:

- GraphSAGE node ID and feature/label row index;
- graph index and tissue name;
- local node index;
- scalar split value (`train`, `validation`, or `test`);
- Entrez GeneID.

A companion evidence table includes neighbor-set hashes, dictionary slots, source-edgelist hashes, observed and reconstructed feature hashes, observed and reconstructed label hashes, and evidence tier.

### 5.4 Remaining caveat

The original preprocessing script has not been found. The mapping is therefore best described as a **forensic reconstruction with exact out-of-sample confirmation**, not direct documentary proof that the authors' code used exactly this sequence of NetworkX calls.

## 6. The 50 input features

### 6.1 Observed feature identity

**Matrix reproduction: byte-level exact. Selection procedure: strongly supported behavioral reconstruction.**

The deposited feature matrix contains:

- 30 MSigDB C1 positional gene-set columns;
- 20 MSigDB C3 motif gene-set columns;
- no observed C7 immunological-signature column.

All 56,944 × 50 = 2,847,200 feature cells are reproduced exactly, with no conflicts among repeated appearances of the same GeneID.

Observed and reconstructed row-major matrix SHA-256:

```text
274ecfee66596b7e9dfb19b71a7fb39a2611a3c140bc199c25062c6ca75bfca1
```

### 6.2 Strongly supported candidate selection procedure

A simple procedure yields the exact 50 columns and order:

1. preserve source-file row order;
2. process C1, then C3;
3. retain sets above a source-level size boundary near 200 unique Entrez IDs;
4. stop once 50 columns have been selected.

This selects 30 C1 sets followed by 20 C3 sets. The next qualifying C3 set would come after the cap, so no C7 set is reached if C7 was nominally included later in the collection sequence.

The released data do not distinguish:

- `size >= 200` from `size > 200`, because no relevant set lies on the discriminating boundary;
- an explicit cap of 50 from a predefined list containing the same 50 sets;
- a C1→C3-only loop from a C1→C3→C7 loop that reaches the cap during C3.

### 6.3 The all-zero `chryq11` column

Feature column 10 is identified as `chryq11` under the candidate procedure.

- `chryq11` contains 204 Entrez GeneIDs in the full MSigDB source collection;
- none is among the 4,301 genes in the GraphSAGE benchmark;
- the released column is therefore entirely zero.

This is strong evidence that feature sets were filtered by their full source size before intersection with graph genes and that empty columns were not removed afterward. It is not, by itself, proof of the literal source code.

### 6.4 MSigDB version ambiguity

MSigDB versions 5.0, 5.1, 5.2, and 6.0 all produce the same selected membership vectors under this rule. The feature matrix therefore does not identify which of those versions was used.

The earlier claim that v5.2 was uniquely determined has been withdrawn.

## 7. The 121 GO labels

### 7.1 Exact output reconstruction

**Matrix reproduction: byte-level exact under the reconstructed row mapping. Literal upstream product: strongly supported but not proven.**

One fixed transformation reproduces every released label cell:

```text
GO annotations:       GOA human release 159
GeneID crosswalk:     May-2016 gp2protein.geneid
GPI metadata:         release-159 UniProt accession and primary symbol
Evidence retained:    EXP, IDA, IEP, IGI, IMP, ISS
Negation:             exclude NOT
Default relations:    involved_in, part_of, enables
Excluded qualifiers:  colocalizes_with, contributes_to
GO normalization:     canonicalize alternate IDs
Propagation:          direct term plus transitive is_a ancestors
Ontology part_of:     do not propagate
Per-gene tuning:      none
Per-column tuning:    none
```

Complete comparison:

| Quantity | Result |
|---|---:|
| Rows | 56,944 |
| Label columns | 121 |
| Binary cells | 6,890,224 |
| Mismatched cells | **0** |
| Exact columns | **121/121** |
| Repeated-gene label conflicts | **0** |

Observed and reconstructed row-major matrix SHA-256:

```text
677cc50459190ba22afc0762356e573ca56e422c80fd36204669a81694afa78d
```

A second independently implemented set-based pipeline produced the same complete matrix.

### 7.2 Why qualifier handling was decisive

An earlier model treated every positive GAF relation as ordinary membership and produced 901 excess positive cells.

Those differences decomposed as:

- 501 dependent on `colocalizes_with`;
- 387 dependent on `contributes_to`;
- 13 caused by one identifier component.

Keeping only the default aspect relations removed the first 888 differences without losing any observed positive label.

### 7.3 Ambiguous Entrez–UniProt mappings

Mappings were not forced one-to-one. Historical many-to-many components were preserved and resolved only with explicit semantic evidence.

The key exceptional component was:

```text
Q9Y620 / RAD54B → GeneID 25788 / RAD54B
O95073 / FSBP   → GeneID 100861412 / FSBP
```

Official reviewed UniProt releases from April, May, and June 2016 also carried an O95073→25788 cross-reference. That historical edge was real, probably reflecting the nested relationship of FSBP within the RAD54B locus.

The GraphSAGE labels nevertheless behave as if O95073/FSBP annotations were not transferred to the node representing GeneID 25788/RAD54B. This behavior can arise from a symbol-aware component resolution or from an Entrez-native annotation source. The exact original mechanism remains open.

Other many-to-many mappings, including hemoglobin, calmodulin, HSPA1A/HSPA1B, and histone H4 components, were retained rather than discarded.

### 7.4 ATP6AP2

GeneID 10159, ATP6AP2, is absent from the release-159 GPI projection used here. Its GraphSAGE label vector is all zero, and the fixed predictor also yields all zero, so this provenance gap does not hide a label mismatch.

### 7.5 GOA source date

GOA release 159 is uniquely exact among tested releases 158–169 when all other transformation choices are held fixed.

| Release | Exact columns | Total gene-label differences |
|---:|---:|---:|
| 158 | 8 | 846 |
| **159** | **121** | **0** |
| 160 | 15 | 652 |
| 161 | 9 | 863 |
| 162 | 1 | 1,182 |
| 163–169 | 0–1 | 1,526–9,782 |

This is strong evidence for the release-159 annotation state. It does not prove that the original code literally opened `goa_human.gaf.159.gz`; an Entrez-native or processed product could encode the same effective associations.

### 7.6 Label aspects and selection

The candidate terms comprise:

- 85 Biological Process terms;
- 26 Cellular Component terms;
- 10 Molecular Function terms.

The common description of all 121 labels as Biological Process terms is incorrect.

Applying the exact transformation to the full historical human annotation universe shows that the candidate set is exactly:

- the 121 most prevalent propagated terms; and
- all propagated terms with approximately 1,000 or more annotated human genes or proteins in release 159.

The release-159 boundary is sharp:

- rank 121: 1,007 mapped human GeneIDs;
- rank 122: 997.

The data do not distinguish “top 121” from “at least 1,000,” or whether the original implementation described the counted entities as genes or proteins.

### 7.7 Direct MSigDB and OhmNet label sources were rejected

Direct memberships from supplied MSigDB releases 5.0 through 6.0 do not reproduce a complete GraphSAGE label column. The supplied OhmNet tissue-specific label files also fail under individual, same-tissue, and union transformations.

Thus, a literal interpretation that the released label vectors were copied directly from one tested MSigDB collection is contradicted by the data. MSigDB may still have supplied term names or a candidate list.

## 8. GO column identities and order

### 8.1 What is uniquely identifiable

The 121 columns correspond to 121 exact GO candidate terms but only 118 distinct membership vectors. Three column pairs contain two GO terms with identical membership on all 4,301 graph genes:

| Columns | Indistinguishable candidate terms |
|---|---|
| 24 and 71 | `GO:0043228` and `GO:0043232` |
| 39 and 63 | `GO:0006464` and `GO:0036211` |
| 48 and 70 | `GO:0043230` and `GO:1903561` |

Membership alone uniquely determines 115 column identities. Six identities require ordering evidence.

### 8.2 Python 2 dictionary-order evidence

A large dictionary keyed by GO ID strings gives a strong match to the observed column order. Across hundreds of plausible construction variants, the best result has approximately:

- longest common subsequence: 94/121;
- Kendall tau: 0.766;
- pairwise concordance: 88.3%;
- exact absolute positions: 13/121;
- exact initial prefix: first five columns.

Ordinary biological orders and small rebuilt dictionaries perform much worse.

Every high-scoring model and every tested GOA release 158–169 gives the same orientation:

- column 24 → `GO:0043228`; column 71 → `GO:0043232`;
- column 39 → `GO:0006464`; column 63 → `GO:0036211`;
- column 48 → `GO:1903561`; column 70 → `GO:0043230`.

This orientation is strongly supported and useful as a provisional complete map, but it is not proven because the exact original 121-column sequence has not been reproduced.

## 9. Exact GraphSAGE-to-DGL conversion

**Current status: byte-level exact.**

The DGL archive is reproduced by:

1. assigning each tissue's largest connected component to that tissue's own graph ID;
2. leaving non-largest components in the first graph ID of the corresponding split;
3. concatenating rows by graph ID while preserving source row order within each ID;
4. fitting `sklearn.preprocessing.StandardScaler` on all GraphSAGE training rows in float64;
5. transforming validation and test with the same scaler;
6. replacing every non-loop undirected edge with both directions;
7. retaining exactly one self-loop per node.

The original GraphSAGE JSON does not contain the DGL graph-ID field, so the conversion is more specific than a simple stable sort.

The supplied DGL features are float64, not float32. Labels, graph IDs, transformed features, and complete directed edge multisets all match exactly.

## 10. Split semantics and leakage

### 10.1 The split is graph-wise

Each node row has one scalar split inherited from its tissue graph. A GeneID can occur in many rows and multiple splits.

Among 4,301 distinct genes:

| Split membership | Genes |
|---|---:|
| Train only | 620 |
| Validation only | 69 |
| Test only | 13 |
| Train + validation | 345 |
| Train + test | 166 |
| Validation + test | 21 |
| Train + validation + test | 3,067 |

A total of 3,599 genes occur in more than one split.

### 10.2 Trivial identity lookup

For the test split:

- 5,524 node rows total;
- 5,490 rows have a GeneID that occurs in training;
- overlap fraction: 99.3845%;
- every training-seen test GeneID has the same 121-label vector as in training.

A lookup model that memorizes the label vector for each training GeneID and predicts all zero for unseen GeneIDs obtains:

```text
test micro-F1 = 0.9971784
true-positive label cells = 198,970
false-positive label cells = 0
false-negative label cells = 1,126
```

The false negatives come from the 34 test GeneIDs absent from training.

This does not yet replace a full leakage study. The planned next controls are:

- randomized label vectors assigned at the GeneID level and copied to every tissue occurrence;
- a permutation that preserves overlap strata so the lookup score remains algebraically unchanged while biological meaning is destroyed;
- a genuinely gene-disjoint split;
- ordinary machine-learning baselines for comparison.

## 11. Main differences from manuscripts and repository descriptions

| Literature or common description | Current evidence |
|---|---|
| The topology is a BioGRID PPI network. | The immediate source is OhmNet; the upstream network is a composite interactome containing BioGRID among several resources. |
| Features include C1, C3, and C7. | The deposited matrix contains exactly 30 C1 and 20 C3 columns and no observed C7 column. A 50-column cap is a plausible explanation, not documentary proof. |
| The 121 GO labels were collected from MSigDB. | No tested MSigDB release reproduces a complete label vector directly; the fixed GOA release-159 transformation reproduces every cell. MSigDB may have supplied term identities rather than memberships. |
| The labels are GO Biological Process terms. | The reconstructed set is 85 BP, 26 CC, and 10 MF. |
| Test graphs are unseen. | Correct graph-wise, but most underlying GeneIDs recur in training; the task is not entity-disjoint. |
| OhmNet experiments use 107 tissue layers. | The released OhmNet archive contains 144 tissue edgelists; GraphSAGE selected 24. The relation among these counts and the exact selection rule remains open. |
| The Menche 2015 interactome is the OhmNet network. | The dimensions differ substantially; Menche is part of the source/method lineage, not the unchanged OhmNet/BioSNAP file. |

## 12. Current open questions

The major output matrices are no longer open. The remaining uncertainty concerns historical provenance and selection choices:

1. **Original preprocessing code:** not found for either GraphSAGE data preparation or the precise label/feature construction.
2. **Literal label source:** GOA release 159 is uniquely exact under the fixed transformation, but an equivalent Entrez-native or processed product may have been used.
3. **Six GO column names:** the three duplicate-vector orientations are strongly supported but not documentary proof.
4. **Exact label-column order:** a strong Python 2 dictionary fingerprint exists, but no plausible model reproduces every position.
5. **Feature selection syntax:** the observed rule is exact behaviorally, but `>=200` versus `>200`, explicit cap versus fixed list, C7 traversal, and MSigDB release remain indistinguishable.
6. **Network selection and split assignment:** the choice of 24 tissues and their 20/2/2 allocation is unresolved.
7. **Historical nested-gene handling:** the data require FSBP annotations not to be assigned to RAD54B, but the original implementation may have used Entrez-native annotations rather than our symbol-aware component resolution.
8. **ATP6AP2 provenance:** the all-zero row is reproduced, but its absence from GPI159 remains a small mapping-history question.
9. **Downstream impact:** usage of this exact PPI benchmark across papers, packages, tutorials, and derivative archives has not yet been systematically cataloged.
10. **Leakage consequences:** the lookup result is measured, but randomized-label and gene-disjoint experiments remain to be implemented.

## 13. Recommended claim language at this stage

A balanced summary is:

> We reconstructed a deterministic mapping from every one of the 56,944 GraphSAGE PPI node rows to an Entrez GeneID using the original OhmNet edgelist order and legacy CPython 2 string-dictionary behavior. The mapping agrees with all 56,411 identities established independently from topology and features and exactly validates all 24 graph topologies, all 50 feature columns, and all 121 label columns. A single fixed GOA release-159 transformation reproduces the complete 56,944 × 121 label matrix with zero differences. The original preprocessing source code, the exact feature-selection syntax, the network-selection rule, and three duplicate-vector GO-column orientations remain unproven.

Avoid claiming that the original source code was recovered or that every historically equivalent implementation choice has been uniquely identified.

## 14. Highest-value next work

The project is ready to shift from open-ended source hunting toward reproducible documentation and consequence testing:

1. create a clean external-data-to-GraphSAGE-to-DGL reproduction workflow;
2. define `sources.yaml` with source URLs, dates, checksums, licenses, and archival mirrors;
3. define `claims.csv` with evidence categories and literature comparisons;
4. archive external data where licensing permits;
5. produce the workflow diagram and current-results report;
6. catalog downstream use of the benchmark;
7. implement randomized-label and gene-disjoint leakage controls;
8. freeze a versioned release before responding to GraphSAGE GitHub issues.

## 15. Principal internal artifacts

- Complete row mapping: `gpt56sol/batches/B104H_20260830T110259Z/graphsage_ppi_node_to_entrez_20260830T110259Z.tsv.gz`
- Evidence-rich row mapping: `gpt56sol/batches/B104H_20260830T110259Z/graphsage_ppi_node_to_entrez_evidence_20260830T110259Z.tsv.gz`
- Full row/feature/label report: `gpt56sol/batches/B104G_20260829T150633Z/B104G_REPORT_20260829T151452Z.md`
- GO-label transformation breakthrough: `gpt56sol/batches/B104A_20260828T145842Z/B104A_REPORT_20260828T145842Z.md`
- GOA date and feature analysis: `gpt56sol/batches/B104E_20260829T121535Z/B104E_GOA_DATE_SCREEN_UNIPROT_AND_FEATURE_RULE_REPORT_20260829T121535Z.md`
- Column-order analysis: `gpt56sol/batches/B104D_20260829T010311Z/B104D_COLUMN_ORDER_DATE_RANGE_AND_UNIPROT_REPORT_20260829T010311Z.md`
- PPI and split provenance: `gpt56sol/batches/B104I_20260830T114918Z/B104I_PPI_PROVENANCE_MENCHE_AND_SPLIT_REPORT_20260830T114918Z.md`
- Literature comparison register: `gpt56sol/batches/B104I_20260830T114918Z/B104I_EVIDENCE_VS_LITERATURE_REGISTER_20260830T114918Z.csv`
- Current source ledger: `gpt56sol/results/source_ledger_through_B104I_FINAL_20260830T114918Z.csv`
- Current provenance events: `gpt56sol/results/provenance_events_through_B104I_20260830T114918Z.csv`

## 16. Bottom line

The benchmark's released data are now understood at a much deeper level than the accompanying documentation provides. The external graph source, node identities, feature values, GO label values, and DGL conversion can all be reproduced. The main remaining challenge is no longer discovering what the files contain; it is documenting which conclusions are exact, which historical mechanisms are inferred, and how the graph-wise split affects the validity of machine-learning evaluation.

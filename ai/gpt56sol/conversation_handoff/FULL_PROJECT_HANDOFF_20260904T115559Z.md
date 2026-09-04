# GraphSAGE PPI provenance/reproduction project — full handoff snapshot

**Snapshot time:** 2026-09-04 06:55 America/Chicago / 2026-09-04T11:55:59Z  
**Purpose:** Seed a new ChatGPT conversation with enough scientific, technical, provenance, and project-management context to continue the work without losing the substance of the current long conversation.  
**Important:** This is a functional reconstruction of the conversation and project state, not a verbatim export of every chat message. It preserves the important questions, decisions, corrections, results, evidence, artifacts, and next steps. The nested `gpt56sol_complete_artifact_snapshot_20260830T124820Z.zip` preserves the large body of assistant-generated artifacts through batch B104I.

---

## 1. Project objective

The project investigates the widely used GraphSAGE PPI benchmark and its DGL derivative. The goals evolved from reverse-engineering the benchmark into a rigorous data-provenance and reproducibility study.

The current principal goals are:

1. Identify the biological identities of the anonymous GraphSAGE PPI nodes.
2. Determine the external origins and construction of the graph topology, 50 input features, and 121 output labels.
3. Reconstruct the deterministic GraphSAGE PPI data from upstream sources.
4. Reconstruct the DGL PPI derivative from the reconstructed GraphSAGE data.
5. Build a clean, readable, independently runnable reproduction package from scratch.
6. Programmatically verify reconstruction at the strongest meaningful level: exact node identities, graph structure, features, labels, split metadata, and DGL transforms.
7. Carefully distinguish exact reconstruction from historically inferred mechanisms and from unresolved details.
8. Document where the released data disagree with or are more specific than the GraphSAGE/OhmNet literature and project webpages.
9. Later, demonstrate the severity of entity leakage in the benchmark with trivial lookup and randomized-label controls.
10. Eventually publish a clean report, reproducible repository, source archive strategy, and helpful responses to old GraphSAGE GitHub issues.

The user strongly values scientific caution, provenance, exact hashes, reproducibility, and explicit correction of earlier overclaims.

---

## 2. People/agents and attribution context

There are three principal intellectual contributors in the investigation:

- **User / human investigator**: research framing, source acquisition, insistence on rigorous biological identifier treatment, skepticism, historical checking, storage-aware execution, and repeated pressure to distinguish proof from hypothesis.
- **Claude / Opus 4.6**: earlier agent whose work predated much of the current ChatGPT investigation. The user keeps that work under `opus46/` in the repository. Claude did not observe the current conversation and was not synchronized with later results.
- **GPT-5.6 Sol**: current assistant, which independently re-derived and substantially extended the earlier work.

Important attribution notes:

- The initial state supplied to GPT-5.6 already contained major groundwork: matching GraphSAGE tissue graphs to OhmNet, Weisfeiler-Lehman structural node matching, partial biological ID recovery, C1/C3 feature-family identification, an initial DGL transformation, and the key observation that train/test graphs reuse the same underlying genes, enabling a lookup baseline.
- GPT-5.6's major original conceptual contributions include the GAF relation/qualifier policy, exact evidence-code/propagation policy, full-human GO prevalence term selection, CPython 2 dictionary fingerprint for label-column order, exact MSigDB feature-generation hypothesis, and especially the CPython 2 string-key/legacy NetworkX node-order reconstruction that resolves all 56,944 rows.
- The user repeatedly caught or prevented overclaims. A particularly important correction was insisting that the early "121/121 exact" result did not yet cover every anonymous row individually.
- A later attribution discussion estimated conceptual shares approximately GPT-5.6 61%, Claude/Opus 22%, user 17%, but those numbers are subjective and not themselves a scientific result.

---

## 3. Initial state provided to GPT-5.6

The initial uploaded investigation summary reported approximately the following state:

- GraphSAGE PPI had 56,944 tissue-instance nodes, 818,716 link records, 50 features, 121 labels, and 24 graph components/tissues.
- 24 GraphSAGE graphs had been matched to OhmNet tissue-specific PPI networks.
- A Weisfeiler-Lehman/color-refinement approach matched most anonymous GraphSAGE nodes to Entrez Gene IDs from OhmNet.
- MSigDB C1 positional and C3 motif sets matched the feature columns after some ambiguous swaps.
- An initial DGL reverse engineering had been done.
- A zero-parameter gene lookup across splits already achieved extremely high micro-F1 because the same genes recur in different tissue graphs and retain identical labels.
- GOA release 159 was a promising source for the labels, but earlier analysis only achieved near matches and was blocked by identifier mapping/evidence-code complications.
- Several earlier interpretations were later corrected, including the claim that the PPI topology was simply "BioGRID", the claim that features uniquely identified MSigDB v5.2, and the claim that labels were only Biological Process terms.

GPT-5.6 was instructed to independently verify, not merely trust the prior agent's conclusions.

---

## 4. Core released dataset facts

### GraphSAGE PPI

Current reconstructed dimensions:

- **56,944 node rows**. A row is a gene occurrence in one tissue graph, not a globally unique biological entity.
- **24 graphs/tissues**.
- **818,716 raw link records** in the released GraphSAGE representation.
- **50 binary input feature columns**.
- **121 binary output label columns**.
- **4,301 distinct Entrez Gene IDs** across the 56,944 tissue-instance rows after the complete node-order reconstruction.

### Split structure

Graph split counts:

- train: 20 tissue graphs, 44,906 rows
- validation: 2 tissue graphs, 6,514 rows
- test: 2 tissue graphs, 5,524 rows

Distinct-gene split membership from the full 4,301-gene mapping:

- train only: 620 genes
- validation only: 69
- test only: 13
- train + validation: 345
- train + test: 166
- validation + test: 21
- train + validation + test: 3,067

Thus 3,599 / 4,301 distinct genes occur in multiple splits, and 3,067 occur in all three.

This is crucial for the later leakage analysis: the benchmark is graph/tissue-disjoint but emphatically not entity/gene-disjoint.

---

## 5. Graph topology provenance

### Immediate source

The 24 GraphSAGE graphs are exact matches to named tissue-specific OhmNet edgelists.

The released OhmNet data README states that each edgelist is a tissue-specific human PPI layer and uses Entrez Gene IDs. The same Gene ID across different edgelists represents the same gene.

### Upstream global interactome

A later BioSNAP audit established:

- union of all 144 supplied OhmNet tissue edgelists = **70,338 unique undirected gene pairs**;
- that union exactly equals the BioSNAP combined OhmNet network;
- every one of those pairs is contained in the BioSNAP global physical interactome with **342,353 unique edges**;
- no missing/extra edges in the comparisons.

The global interactome has 21,557 proteins and 342,353 interactions, matching the OhmNet manuscript's stated global network dimensions.

### Important literature correction

Do **not** describe the GraphSAGE PPI graphs simply as BioGRID networks.

The GraphSAGE dataset webpage links only to BioGRID when describing the PPI data, which is incomplete/misleading. The OhmNet manuscript describes a composite experimentally supported physical interactome sourced from multiple resources/antecedents including MIntAct/IntAct, Rolland et al., BioGRID, HPRD, CORUM, and Menche et al. The tissue graphs are derived from that composite network.

The defensible wording is approximately:

> GraphSAGE's 24 PPI graphs are exact tissue-specific OhmNet subnetworks. OhmNet constructed tissue-specific layers by filtering a larger composite experimentally supported human interactome; BioGRID was one contributor, not the sole source.

### Menche clarification

Menche et al. 2015 used a smaller interactome (13,460 proteins, 141,296 interactions), so the exact OhmNet/BioSNAP 21,557 / 342,353 network is not simply the unchanged Menche network. Menche is part of the lineage/methodological ancestry.

---

## 6. Complete node-to-Entrez mapping — major breakthrough

### Earlier partial mapping

The earlier agent used Weisfeiler-Lehman graph-color refinement to align GraphSAGE nodes with OhmNet Entrez nodes. This recovered the vast majority, and MSigDB feature vectors resolved many structural equivalence classes.

Immediately before the final breakthrough, 56,411 / 56,944 rows had independent biological identities and 533 rows remained ambiguous within 183 topology/feature-equivalence classes.

### CPython 2 string-key / legacy NetworkX ordering reconstruction

GPT-5.6 discovered that the remaining row order can be exactly reconstructed by modeling the likely legacy Python/NetworkX graph construction:

1. Read each original OhmNet edgelist in line order.
2. Keep Entrez node identifiers as **strings**, rather than converting them to integers.
3. Insert those strings into a graph/dictionary as edges are processed.
4. Simulate **64-bit, unrandomized CPython 2.7 dictionary** behavior: historical string hash, probing, resize, and table iteration.
5. Iterate the occupied dictionary table in legacy order to obtain local GraphSAGE node order.

This global mechanism:

- agrees with **all 56,411 previously independently anchored rows**;
- resolves the remaining 533 rows without using label-specific optimization;
- yields a complete mapping of all 56,944 rows to 4,301 distinct Entrez Gene IDs;
- reproduces graph neighbor sets, feature vectors, and label vectors for all rows.

Obvious alternatives were dramatically worse (first appearance, sorted numeric GeneID, reverse sorted, Python 2 integer keys, etc.).

### Independent positive control

A separate CPython 2 dictionary simulator inserting the string keys `"0"` through `"56943"` reproduced the complete serialized key order of `ppi-class_map.json` exactly under a 64-bit unrandomized Python 2 model. A 32-bit model did not.

This is strong independent evidence that legacy CPython 2 dictionary iteration materially shaped GraphSAGE artifact ordering.

### Confidence wording

The data-level mapping is exact under the recovered mechanism. The historical mechanism remains **strongly inferred**, because the original preprocessing script has not been found.

Use wording such as:

> The released mapping deterministically reconstructs all 56,944 GraphSAGE PPI node identities and exactly reproduces every deposited edge, feature value, and label value. The historical explanation — legacy NetworkX graph construction with string node IDs and 64-bit unrandomized CPython 2 dictionary iteration — is strongly supported by exact out-of-sample agreement and a class-map serialization fingerprint, but is not directly documented by the original preprocessing code.

Key current files are included in this handoff bundle under `critical_artifacts/`:

- `graphsage_ppi_node_to_entrez_20260830T110259Z.tsv.gz`
- `graphsage_ppi_node_to_entrez_evidence_20260830T110259Z.tsv.gz`
- `B104G_full_graphsage_row_to_entrez_mapping.csv.gz`
- `B104G_full_4301_gene_universe.csv.gz`
- `B104H_MAPPING_VALIDATION_20260830T110259Z.json`
- independent validation JSON files from B104G

---

## 7. Feature reconstruction

### Exact observed feature identities

The 50 binary GraphSAGE feature columns are exactly reproducible as:

- **30 MSigDB C1 positional gene sets**
- **20 MSigDB C3 motif/regulatory-target gene sets**
- **0 C7 immunological-signature sets** in the released matrix

This differs from a literal reading of GraphSAGE documentation that mentions C1, C3, and C7 feature families.

### Strong reconstruction hypothesis for feature selection

A simple source-level selection rule exactly reproduces the full feature matrix:

1. Process C1 gene sets in source GMT/XML order.
2. Then process C3 sets in source order.
3. Retain gene sets meeting an approximately 200-member **source-level** size threshold.
4. Append qualifying sets until a global maximum of 50 columns is reached.
5. Project those sets onto the GraphSAGE gene universe.

This produces 30 C1 + 20 C3 columns exactly.

The exact historical operator cannot be distinguished between `>= 200` and `> 200` because no relevant set lies at the discriminating boundary. The procedure itself is therefore a **strongly supported exact reconstruction hypothesis**, not source-code proof.

### `chryq11` as strong supporting evidence

Feature column 10 is the C1 set `chryq11`.

- It has 204 Entrez Gene IDs in the full MSigDB source set, so it passes the approximately 200-member source threshold.
- None of those genes occurs in the resolved GraphSAGE universe.
- The deposited GraphSAGE feature column is therefore entirely zero.

This is strong evidence that filtering happened on full-source gene-set size **before** projection to GraphSAGE genes, and that no later cleanup removed empty columns. It also suggests little manual inspection after matrix construction.

### MSigDB version is not identifiable from the feature matrix

The same procedure produces the same 50 observed vectors for supplied MSigDB versions 5.0, 5.1, 5.2, and 6.0.

Therefore, earlier claims that the features uniquely identified v5.2 are withdrawn.

For the new clean reproduction, the plan is to use **MSigDB v6.0 C1/C3**, because the observed feature vectors are unchanged and v6.0 has a more permissive license regime. The user wants to verify licensing and, if permitted, archive only the minimal original-order C1/C3 inputs needed by the workflow.

A helper script can transform a privately downloaded original MSigDB package into a minimal reproducibility input. Public CI should not depend on automating the historical MSigDB login wall.

Current critical feature artifacts included here:

- `B104E_exact_MSigDB52_feature_generation_rule_20260829T121535Z.csv`
- `B104E_exact_MSigDB52_feature_generation_validation_20260829T121535Z.json`
- `B104E_independent_MSigDB52_feature_validation_20260829T121535Z.json`
- `B104E_MSigDB_feature_version_screen_PROVENANCE_SAFE_20260829T121535Z.csv`

Full feature reconstruction over all rows:

- 56,944 x 50 = 2,847,200 binary cells
- zero differences under the final row map

---

## 8. GO label reconstruction

### Current exact transformation

The 121 GraphSAGE labels are reproduced using a single global policy based on GOA human release 159 and historical mapping/ontology data.

Core sources used in the successful reconstruction:

- `goa_human.gaf.159.gz` — GOA human release 159, generated 2016-07-04
- `goa_human.gpi.159.gz`
- historical `2016-06-01-gp2protein.geneid.gz` (internally generated May 2016)
- June 1, 2016 `go.obo`

Exact annotation policy:

**Evidence codes retained:**

- EXP
- IDA
- IEP
- IGI
- IMP
- ISS

**Negation:**

- exclude `NOT`

**Ordinary GAF relations retained:**

- Biological Process: `involved_in`
- Cellular Component: `part_of`
- Molecular Function: `enables`

**Qualified relations not treated as ordinary binary membership:**

- `colocalizes_with`
- `contributes_to`

**GO ID processing:**

- canonicalize alternate GO IDs

**Ontology propagation:**

- direct annotated term + transitive `is_a` ancestors only
- do **not** propagate ontology `part_of`

### Crucial qualifier breakthrough

A near-exact reconstruction initially had 901 excess cells. Those decomposed into:

- 501 from `colocalizes_with`
- 387 from `contributes_to`
- 13 from the O95073/RAD54B mapping issue

Removing the qualified relations as ordinary membership and resolving the O95073 component eliminated the discrepancies.

### Labels span all three GO namespaces

The reconstructed 121 candidate GO terms are:

- 85 Biological Process
- 26 Cellular Component
- 10 Molecular Function

Do not describe all labels as BP.

This result was discovered in the GPT-5.6 analysis and later communicated to Claude; Claude did not originate it.

### Term selection

On the full historical human annotation universe, the exact candidate set is recovered as the **top 121 GO terms by prevalence**.

Equivalent observed rule under release 159:

- all terms with at least approximately **1,000 annotated distinct historical human Gene IDs/proteins**.

Boundary values previously observed:

- rank 120 GO:0051254 — 1,016
- rank 121 GO:0031399 — 1,007
- rank 122 GO:0043167 — 997
- rank 123 GO:0051049 — 991

Thus top-121 and >=1000 are observationally equivalent for the successful snapshot. The exact original source-code formulation remains unknown.

### GOA release date screen

A controlled screen of GOA human releases 158 through 169 held mapping, ontology, and transformation policy fixed.

Release 159 was uniquely exact under the tested model.

Selected results:

- release 158: 8 exact columns; 21 FP, 825 FN, 846 total differences
- release **159**: 121 exact columns; 0 FP, 0 FN
- release 160: 15 exact columns; 438 FP, 214 FN, 652 total differences
- release 161: 9 exact; 863 total differences
- release 162: 1 exact; 1,182 total differences
- release 168: 0 exact; 9,225 total differences
- release 169: 0 exact; 9,782 total differences

Release 159 also uniquely yields the exact 121 candidate terms together with the natural approximately 1,000-member threshold behavior among the tested nearby releases.

Critical file included:

- `B104E_GOA_release158_169_validated_summary_20260829T121535Z.csv`

### Exactness scope — important correction

An early claim said "121/121 exact" when the comparison covered only 4,268 individually resolved Entrez genes / 56,411 rows.

The user correctly objected that this was not yet a complete row-level reconstruction.

The analysis then progressed through:

1. class-multiset equality for all 183 unresolved equivalence classes; and later
2. the CPython 2 string-key row-order solution.

Current result:

> Under the full 56,944-row node map, all 121 label columns and all 6,890,224 binary label cells are reproduced exactly by one fixed global transformation with no per-gene or per-column tuning.

Independent validation JSON is included in `critical_artifacts/`.

---

## 9. Identifier mapping policy and O95073 / FSBP / RAD54B

### General policy

The correct historical identifier treatment preserves many-to-many mappings rather than forcing one-to-one matches.

Key principles:

1. Preserve historical GeneID–UniProt edges as a bipartite graph.
2. Analyze complete connected mapping components before restricting to the GraphSAGE gene universe.
3. Use primary symbols only for unique semantic resolution/bijection where justified.
4. Do not use broad synonyms as an aggressive fallback.
5. Keep unresolved mappings ambiguous rather than selecting whatever improves the label match.
6. Never tune mapping decisions per gene against the observed label matrix.

### O95073 historical audit

Date-matched reviewed Swiss-Prot records for releases 2016_04, 2016_05, and 2016_06 showed:

- O95073 / FSBP carried GeneIDs **100861412 and 25788** in all three releases.
- Q9Y620 / RAD54B carried GeneID **25788**.

Therefore, the O95073→25788 edge was a real historical UniProt cross-reference, not a malformed mapping row.

The correct interpretation is narrower:

- FSBP annotations associated with O95073 should not be projected onto the GraphSAGE node representing RAD54B / GeneID 25788.
- The semantic mapping used for reconstruction treats Q9Y620 as RAD54B→25788 and O95073 as FSBP→100861412 for annotation projection.

This eliminates 13 false-positive label cells.

The original preprocessing mechanism remains unknown: it may have used Entrez-native annotations, HGNC/symbol-aware mapping, or another intermediate that naturally separated the nested/host-gene relationship.

Critical file included:

- `B104E_UniProt_O95073_Q9Y620_flatfile_record_comparison_20260829T121535Z.csv`

---

## 10. Label-column identities and order

### 121 terms versus 118 distinct vectors

The 121 reconstructed candidate GO terms correspond to only 118 distinct binary membership vectors over the GraphSAGE gene universe.

Three pairs are membership-indistinguishable:

- columns 24 and 71: GO:0043228 / GO:0043232
- columns 39 and 63: GO:0006464 / GO:0036211
- columns 48 and 70: GO:0043230 / GO:1903561

Therefore, membership alone cannot prove which member of each pair belongs to which column.

### Python 2 label-order fingerprint

Conventional orderings performed poorly:

- GO ID order
- GO term-name alphabetic order
- GraphSAGE positive counts
- full-human prevalence
- OBO stanza order
- first accepted GAF occurrence
- MSigDB XML order

The strongest consistent signal comes from simulating a large 64-bit unrandomized CPython 2 dictionary keyed by GO ID strings.

A strong model based on a large GO dictionary built during GAF processing reached approximately:

- Kendall tau about 0.77
- pairwise concordance about 88.6%
- LCS 92-94 of 121 depending on model variant
- exact prefix first 5

A fresh dictionary containing only the selected 121 terms destroys the signal, which suggests the selected labels were filtered while iterating a larger existing dictionary rather than copied to a fresh 121-key dict and iterated again.

### Provisional duplicate-vector orientation

Across all high-scoring plausible Python 2 order models, orientation was stable:

- col 24 -> GO:0043228, col 71 -> GO:0043232
- col 39 -> GO:0006464, col 63 -> GO:0036211
- col 48 -> GO:1903561, col 70 -> GO:0043230

These assignments are **strongly supported but not proven**.

Critical file included:

- `B104C_inferred_unique_121_GO_column_order_20260828T194921Z.csv`

Further combinatorial search for exact historical label order is now considered low priority unless an ordered label-name artifact is found.

---

## 11. DGL reconstruction

The DGL PPI dataset was reverse engineered as a transformation of GraphSAGE PPI.

Current broad understanding includes:

- node grouping/reordering by graph ID and connected-component handling;
- StandardScaler fit on training rows only;
- feature transform in float64 then cast/represented appropriately for DGL outputs;
- undirected GraphSAGE edges converted to directed edges;
- self-loops added;
- split packaging into DGL graph objects/files.

Earlier Claude work established much of the broad transformation; GPT-5.6 independently corrected/verified details.

The new reproduction package should reconstruct DGL from the **reconstructed GraphSAGE artifacts**, never from the downloaded DGL target.

Validation should compare:

- graph membership/counts
- row order
- labels
- feature values and dtype/tolerance
- directed edges
- self-loops
- split assignments

Exact binary archive serialization is not necessarily scientifically meaningful if the logical content is exact.

---

## 12. Leakage finding and later planned experiment

The benchmark splits tissues/graphs rather than genes. Because the same genes appear in train, validation, and test and have tissue-invariant label vectors, a trivial GeneID lookup can predict held-out graph labels almost perfectly.

Under the full node map, prior measured test behavior was approximately:

- test rows: 5,524
- rows whose GeneID appears in training: 5,490
- overlap: 99.3845%
- GeneID lookup micro-F1: approximately **0.9971784**

At an earlier partial-mapping stage, values around 0.9956 were reported. Use the full-row value for current work.

The user wants a later controlled ML demonstration that:

1. fits a trivial model that exploits only repeated GeneID identity;
2. achieves essentially the same performance on **gene-level randomized labels**, provided each gene keeps the same randomized vector across tissue occurrences;
3. contrasts this with a true gene-disjoint split;
4. potentially uses a non-lookup baseline as an additional demonstration.

A particularly strong randomization control is to permute complete label vectors among genes while preserving overlap strata so the lookup score remains unchanged despite destroying biological meaning.

This experiment is intentionally deferred until after the clean provenance/reproduction workflow is built.

Critical planning file included:

- `LEAKAGE_EXPERIMENT_BACKLOG_20260830T110259Z.md`

---

## 13. Evidence categories agreed with user

The final project should classify claims into four main evidence levels:

1. **Byte-level exact / data-level exact** — directly reproduced or compared exactly.
2. **Strongly inferred** — one mechanism explains the observed data and survives strong controls, but original source code/documentation is absent.
3. **Documented** — stated by a manuscript, repository, or official source.
4. **Open** — multiple histories remain observationally equivalent or evidence is insufficient.

Generated comparisons may use finer technical labels such as:

- byte_exact
- array_exact
- structurally_exact
- semantically_exact
- numeric_tolerance
- strongly_inferred
- documented_only
- unresolved

The user wants later `claims.csv` to record claim, evidence, supporting artifacts, literature statement, discrepancy, and remaining uncertainty.

---

## 14. Major corrections and superseded claims

Preserve these explicitly so a new agent does not regress to outdated conclusions.

### Superseded: incomplete node map

Earlier state: about 4,268-4,278 independently resolved genes with several hundred ambiguous rows.

Current state: complete 56,944-row map to 4,301 distinct Entrez Gene IDs under the strong legacy Python ordering reconstruction.

### Superseded: "121/121 exact" implied all rows

Earlier exact label result covered only the independently resolved subset.

Current state: all 56,944 rows / 6,890,224 label cells are exact under the final full row map.

### Superseded: labels are all BP

Current: 85 BP, 26 CC, 10 MF candidate terms.

### Superseded: features uniquely identify MSigDB v5.2

Current: supplied MSigDB 5.0, 5.1, 5.2, and 6.0 produce the same 50 observed feature vectors under the reconstructed rule.

### Superseded: PPI source is simply BioGRID

Current: exact immediate source is OhmNet tissue networks; OhmNet uses a larger composite experimentally supported global interactome in which BioGRID is one contributor.

### Superseded: O95073→25788 was a corrupt mapping

Current: historical Swiss-Prot proves that cross-reference was real in 2016. The important behavior is that FSBP annotations were not projected onto the RAD54B node.

### Corrected MSigDB label-presence parser bug

An earlier B104B parser counted GO IDs anywhere in XML rather than strictly within C5.

Correct C5 candidate GO-ID presence among the 121 reconstructed terms:

- v5.0: 57 / 121
- v5.1: 57 / 121
- v5.2: 6 / 121
- v6.0: 6 / 121

The conclusion that direct MSigDB membership does not generate the 121 label columns remains unchanged.

---

## 15. Current discrepancies with literature/project documentation

These should later be summarized in a dedicated report.

### GraphSAGE PPI source webpage

The GraphSAGE dataset webpage links only to BioGRID for PPI source. Our evidence shows the released topology is exact OhmNet tissue networks from a composite upstream interactome. Calling the webpage wholly false is too strong; call it incomplete/misleading provenance.

### GraphSAGE feature description

Documentation mentions positional gene sets, motif gene sets, and immunological signatures. The released 50 columns are exactly 30 C1 + 20 C3 and no C7.

A plausible explanation is a global 50-column cap reached during C3, but the literal original feature-selection code is not known.

### GraphSAGE label description

The paper says 121 GO gene sets were collected from MSigDB. Direct memberships from tested MSigDB versions do not reproduce the labels. GOA release 159 with the recovered policy does reproduce them exactly.

MSigDB may have supplied a term list/metadata rather than final memberships; exact historical provenance remains open.

### "Unseen graphs"

The test graphs are indeed unseen tissue graphs, but the underlying genes are overwhelmingly seen during training. The distinction is critical for understanding the later leakage analysis.

### OhmNet network count

OhmNet manuscript experiments discuss 107 tissue layers, while a released archive contains 144 tissue edgelists. GraphSAGE uses 24. The exact rule selecting those 24 and assigning 20/2/2 splits remains unresolved.

Critical register included:

- `B104I_EVIDENCE_VS_LITERATURE_REGISTER_20260830T114918Z.csv`

---

## 16. GitHub issues to revisit after publication

The user plans to respond helpfully to old GraphSAGE GitHub issues once the clean repository/report is stable.

High-value issues already recorded:

- #188 — PPI ID to protein-name correspondence: directly addressed by the complete node-to-Entrez map.
- #190 — how C1/C3/C7 became 50 dimensions: addressed by exact feature reconstruction and the 30 C1 + 20 C3 result.
- #78 — exact PPI source, features, preprocessing.
- #86 — dataset description.
- #16 — PPI F1 evaluation.
- #32 — reproducing reported results/hyperparameters.

Tracker files included:

- `B104I_GRAPHSAGE_GITHUB_ISSUES_TRACKER_20260830T114918Z.md`
- `B104I_GRAPHSAGE_GITHUB_ISSUES_TRACKER_20260830T114918Z.csv`

Do not post until there is a stable public reproduction repository/tag and durable links.

---

## 17. Repository organization already agreed locally

The user's current repository root is approximately:

```text
GraphSagePPIAnalysis/
├── README.md
├── data/          # ignored; external/source data
├── gpt56sol/      # GPT-5.6 artifacts, scripts, reports, user-run outputs
├── opus46/        # Claude/Opus 4.6 artifacts
└── papers/        # ignored; local manuscript PDFs/supplements
```

The user moved the working directory contents into the repository and organized files manually.

The root `.gitignore` should broadly ignore:

- `/data/`
- `/papers/`
- redundant snapshot ZIP `gpt56sol_complete_artifact_snapshot_*.zip`
- Python caches/environments
- temp/download fragments

But should **not** broadly ignore `.gz`, `.json`, `.csv`, `.npy`, `.zip`, etc., because compact derived artifacts may be important evidence.

The user prefers erring on the side of including compact scientific/provenance outputs in Git.

---

## 18. Current TODO list

User's project TODOs:

- Export Claude artifacts.
- Create a summary of manuscript/repository claims that are incorrect or incomplete and the evidence.
- Create a workflow diagram showing source data and transformations into GraphSAGE and DGL files.
- Create a GitHub Actions workflow that starts with external data, reconstructs GraphSAGE/DGL, and runs validation tests.
- Prepare a report describing the data workflow, evidence, and uncertainty.
- Catalog usage of the PPI dataset to estimate impact.
- Create machine-learning controls demonstrating leakage, including randomized-label predictive performance.
- Respond to GraphSAGE GitHub issues with stable reproducibility resources.
- Archive external data to protect against future unavailability, if licenses permit.

The immediate priority is the clean core reproduction package.

---

## 19. Clean reproduction package — agreed scope

The exact agreed scope statement is:

> **Reconstruct every deterministic artifact needed by supervised GraphSAGE and DGL, but do not attempt to regenerate the original stochastic `ppi-walks.txt` byte-for-byte.**

The clean reproduction does **not** need to reproduce every exploratory analysis performed during this investigation.

It should:

1. Freshly acquire/validate external source files and released reference targets.
2. Reconstruct the selected OhmNet topology and complete node order.
3. Reconstruct the 50 MSigDB feature columns.
4. Reconstruct the 121 GO labels.
5. Assemble deterministic GraphSAGE artifacts needed by supervised GraphSAGE.
6. Compare reconstructed GraphSAGE outputs to independently downloaded released targets.
7. Transform reconstructed GraphSAGE data into the DGL representation.
8. Compare reconstructed DGL outputs to released DGL targets.
9. Generate a compact run/provenance/validation summary.

Optional historical audits are deliberately omitted from the production package.

---

## 20. Agreed reproduction technology choices

### Workflow manager

Use **Snakemake**.

Reason:

- explicit dependency DAG;
- readable rule-based workflow;
- strong fit for data acquisition -> transformation -> validation;
- familiar to the user;
- suitable for later GitHub Actions.

Start with **one Snakefile**. Split into included `.smk` files only if it becomes genuinely difficult to read.

### Environment

Use **Pixi** as the canonical environment manager.

Do **not** add `environment.yml`; the user explicitly wants to avoid competing canonical environment definitions.

Use a committed `pixi.lock`.

Linux x86-64 is the preferred canonical exact-regression platform; WSL2 is appropriate on the user's Windows machine.

Do not install actual Python 2.7 merely for historical ordering. Reimplement the relevant CPython 2 behavior explicitly in readable modern Python and unit test it.

### Testing

Use **pytest** for unit and integration/regression tests.

### CLI

Avoid a large central CLI with many modes.

User-facing commands should be simple Pixi tasks such as:

```bash
pixi run reproduce
pixi run test
pixi run clean
```

Snakemake is the workflow interface. Individual Python modules may expose small internal `argparse` entry points for Snakemake rules.

---

## 21. Agreed simplified reproduction tree

The latest design is intentionally compact:

```text
reproduction/
├── README.md
├── pixi.toml
├── pixi.lock
├── pyproject.toml
├── Snakefile
│
├── config/
│   ├── sources.tsv
│   └── reconstruction.yaml
│
├── resources/
│   ├── selected_graphs.tsv
│   ├── label_columns.tsv
│   ├── identifier_decisions.tsv
│   └── expected.yaml
│
├── src/
│   └── graphsage_ppi_repro/
│       ├── __init__.py
│       ├── sources.py
│       ├── manifest.py
│       ├── legacy_dict.py
│       ├── topology.py
│       ├── features.py
│       ├── labels.py
│       ├── graphsage.py
│       ├── dgl.py
│       └── validate.py
│
├── tests/
│   ├── data/
│   ├── test_legacy_dict.py
│   ├── test_identifier_mapping.py
│   ├── test_go_labels.py
│   ├── test_features.py
│   └── test_reproduction.py
│
├── docs/
│   ├── architecture.md
│   ├── data-sources.md
│   └── limitations.md
│
├── build/              # generated, gitignored
└── results/            # generated, gitignored
```

The standalone agreed design is included as:

- `REPRODUCTION_REPOSITORY_OUTLINE_CORE_20260901T164642Z.md`

### Directory meaning

- `config/`: small human-maintained operational/scientific configuration.
- `resources/`: small evidence-backed recovered specifications and expected invariants. These are not ordinary knobs users should casually edit.
- `src/`: clean modern Python implementation of core transformations.
- `tests/`: synthetic/unit and full-reproduction validation.
- `docs/`: narrowly technical software/data documentation, not the eventual scientific narrative.
- `build/`: regenerated intermediate data, ignored.
- `results/`: regenerated run outputs/validation reports, ignored and published only as release artifacts when appropriate.

Scientific reports and `claims.csv` belong at repository level, outside the reproduction software package.

---

## 22. Proposed source manifest simplification

The user rejected a deeply nested `sources.yaml` as unnecessarily complex.

Use a flat `config/sources.tsv`, with columns approximately:

```text
source_id
role
filename
url
archive_url
sha256
size_bytes
access
license_note
description
```

Possible `access` values:

- `public`
- `project_archive`
- `manual`

Long licensing discussion belongs in `docs/data-sources.md`, not the TSV.

A generated run manifest can record actual retrieval time, redirects, observed hashes, etc.

---

## 23. `reconstruction.yaml`

Keep this small and readable. It records the canonical reconstruction policy, not machine-specific execution settings.

Conceptual content:

```yaml
node_order:
  method: cpython27_string_dict
  word_size_bits: 64
  hash_randomization: false

features:
  collections: [C1, C3]
  minimum_source_members: 200
  maximum_columns: 50

labels:
  evidence_codes: [EXP, IDA, IEP, IGI, IMP, ISS]
  allowed_relations:
    biological_process: involved_in
    cellular_component: part_of
    molecular_function: enables
  exclude_not: true
  excluded_relations: [colocalizes_with, contributes_to]
  propagation: [is_a]

resources:
  selected_graphs: resources/selected_graphs.tsv
  label_columns: resources/label_columns.tsv
  identifier_decisions: resources/identifier_decisions.tsv
```

Note the unresolved historical ambiguity around the exact feature threshold operator; the production code may use one explicit equivalent policy while documentation records the non-identifiability.

---

## 24. Core `resources/`

### `selected_graphs.tsv`

Frozen, evidence-backed selection/order/split of the 24 known GraphSAGE tissue graphs.

We know the 24 identities and split assignments exactly from released data, but not the original algorithm that selected them from all OhmNet networks.

### `label_columns.tsv`

Stores the 121 recovered GO column identities.

- 115 have unique membership identities.
- 6 rows belong to the three duplicate-vector pairs and should carry a provisional/strongly-supported status.

### `identifier_decisions.tsv`

Stores the few semantic mapping decisions that cannot be expressed by a simple generic historical crosswalk rule, including O95073/FSBP/RAD54B and other historically ambiguous components.

Avoid unexplained hard-coded special cases in Python.

### `expected.yaml`

Small committed invariants, for example:

- 24 graphs
- 56,944 rows
- 818,716 GraphSAGE link records
- feature shape 56,944 x 50
- label shape 56,944 x 121
- 4,301 distinct Gene IDs
- split graph counts 20/2/2
- 30 C1, 20 C3
- zero feature column 10
- 121 label columns / 118 distinct membership vectors

Source file SHA-256 values belong in `sources.tsv`, not duplicated here.

---

## 25. Minimal core Python modules

### `sources.py`

Download/validate external files according to `sources.tsv`.

Requirements:

- download to `.part`/temporary path;
- verify expected hash and size;
- reject HTML error/login pages masquerading as data;
- check magic bytes/archive integrity;
- atomically rename only after validation;
- never silently update a changed upstream hash.

### `manifest.py`

Write run-level provenance:

- source IDs and resolved URLs;
- observed hashes;
- Git commit;
- Pixi lock hash;
- software versions;
- output hashes;
- validation summary.

### `legacy_dict.py`

Readable modern implementation of the necessary 64-bit unrandomized CPython 2.7 string hash/dict behavior.

Must have focused unit tests for collisions, resizes, iteration order, and the known GraphSAGE class-map positive control.

### `topology.py`

- parse selected OhmNet tissue edgelists;
- reconstruct local row ordering using `legacy_dict.py`;
- build node-to-Entrez mapping;
- build deterministic graph topology/split representation.

### `features.py`

- parse minimal MSigDB C1/C3 source in original order;
- apply documented source-size rule and 50-column cap;
- project memberships to ordered GraphSAGE rows;
- emit deterministic 56,944 x 50 matrix.

### `labels.py`

- parse GOA GAF/GPI, historical gp2protein mapping, and GO ontology;
- build ambiguity-preserving identifier components;
- apply documented mapping decisions;
- filter evidence/relations/NOT;
- canonicalize alt IDs;
- propagate `is_a` only;
- select/reconstruct label membership;
- project to ordered GraphSAGE rows and recovered column order.

### `graphsage.py`

Assemble deterministic released GraphSAGE components used by supervised GraphSAGE.

### `dgl.py`

Transform reconstructed GraphSAGE logical data into DGL logical data.

### `validate.py`

Only module that compares reconstructed outputs to independently downloaded reference GraphSAGE/DGL targets.

This separation protects against circular reconstruction.

---

## 26. Circularity boundary

A central design principle is:

> Reconstruction must not read the released GraphSAGE/DGL reference targets. Validation may read both reconstructed and reference outputs.

The recovered `resources/` are allowed because they are explicitly documented results of the forensic discovery phase.

A useful integration test is to make the GraphSAGE/DGL reference-target directory temporarily unavailable and verify that reconstruction still completes; only comparison/validation should fail or remain pending.

---

## 27. MSigDB licensing and reproducibility plan

Historical MSigDB releases had login/licensing constraints that complicate public unattended GitHub Actions.

Important empirical result:

- v5.0, v5.1, v5.2, and v6.0 produce the same observed 50 GraphSAGE feature vectors under the recovered procedure.

Current plan:

1. Use **MSigDB v6.0 C1/C3** as the canonical reproducibility source.
2. Document that earlier tested releases yield the same result and that the actual historical release used by GraphSAGE is not identifiable from the matrix.
3. Verify the v6.0 licensing terms for redistribution of the needed C1/C3 content.
4. If permitted, archive only the minimal original-order C1/C3 files needed by the pipeline in a durable public archive.
5. Also provide a small script that takes a user's privately downloaded original MSigDB package and extracts/transforms the minimal reproduction input.
6. Avoid private CI if a legally redistributable v6.0 minimal source can be archived publicly.

Do not commit historical login-protected MSigDB packages to the public repository.

---

## 28. `ppi-walks.txt` — agreed scope and open side investigation

### Scope decision

`ppi-walks.txt` is **not required** for the supervised GraphSAGE PPI GO prediction task.

GraphSAGE repository documentation describes `<train_prefix>-walks.txt` as optional and used only for the unsupervised version.

The supplied `graphsage.utils.py` source confirms:

- `load_data(..., load_walks=False)` by default;
- walks are only read when `load_walks=True`;
- `run_random_walks` uses `random.choice` from Python's `random` module;
- constants `WALK_LEN=5`, `N_WALKS=50`;
- main-mode node list selects nodes that are neither validation nor test;
- graph is restricted to the training-node subgraph;
- 50 walks are run from each source node;
- in each length-5 walk, co-occurrence pairs `(source_node, curr_node)` are emitted after each non-initial step.

The GraphSAGE manuscript says 50 random walks of length 5 from each node were used for the unsupervised objective.

The original Python utility does not explicitly call `random.seed`, making byte-for-byte regeneration dependent on undocumented Python random state plus graph/neighbor iteration order.

Therefore the core reproduction package will omit it as a reconstruction target.

### Open exploratory side task

Immediately before this handoff, the user asked GPT-5.6 to explore whether the original `ppi-walks.txt` could nonetheless be reproduced by exhaustive/pruned search using:

- the original GraphSAGE code (uploaded as `GraphSAGE-a0fdef9.zip`);
- the supplied `graphsage.utils.py`;
- likely Python/NetworkX versions;
- the observed `ppi-walks.txt` target;
- arbitrary/common random seeds;
- alternative plausible graph/node/neighbor orderings;
- early-abort matching against the output prefix so many candidate seeds/configurations can be tested quickly.

That analysis had **not yet been performed** when the conversation hit its length limit.

This is a side investigation only. It must not delay the agreed deterministic reproduction package.

The new conversation can either perform this exploratory walk search first or proceed directly to implementing the reproduction skeleton.

External input `GraphSAGE-a0fdef9.zip` is intentionally not included in this handoff bundle because the user said external/input data can be resupplied.

---

## 29. Why identity features are a later leakage concern

The GraphSAGE repository recommends using `--identity_dim` when the task does not require generalization to unseen data.

The user correctly noticed this is potentially severe for PPI because test tissue graphs reuse nearly all the same genes as training.

This is not part of the deterministic data reproduction. It should be explored later in the leakage/evaluation phase, especially because identity-sensitive models can exploit repeated biological entities across graph splits.

---

## 30. Current confidence / unresolved issues

### Essentially exact at the data level

- 24 GraphSAGE graph identities/topology versus OhmNet.
- complete 56,944 row-to-Entrez mapping under the strong legacy Python order mechanism.
- 50 feature matrix over all rows.
- 121 label matrix over all rows.
- external OhmNet/BioSNAP edge-set chain.
- split membership and gene reuse.

### Strongly inferred historical mechanisms

- exact legacy NetworkX/CPython 2 route that produced node order.
- feature-selection implementation (source order + approximately 200 threshold + cap 50).
- label-column order mechanism via a larger Python 2 GO dictionary.
- semantic mapping treatment used by the original authors for nested/host-gene cross-references.

### Open / observationally equivalent

- exact original MSigDB release among tested versions producing identical feature vectors.
- `>=200` versus `>200` feature threshold operator.
- exact source-code wording of top-121 versus >=1000 GO term selection.
- exact historical source file versus an Entrez-native intermediate equivalent to GOA release 159.
- six duplicate-vector term-to-column identities are strongly supported but not proven.
- why exactly those 24 OhmNet tissue graphs were selected and the precise 20/2/2 split-selection algorithm.
- original random state for `ppi-walks.txt`.

None of these open items prevents the core deterministic reproduction strategy.

---

## 31. Progress reports already prepared

Included at the handoff root:

- `INVESTIGATION_HISTORY_AND_EXPLORATION_LOG_20260831T205314Z.md`
- `CURRENT_FINDINGS_AND_EVIDENCE_SUMMARY_20260831T205314Z.md`
- `REPRODUCTION_REPOSITORY_OUTLINE_CORE_20260901T164642Z.md`

These are useful human-readable snapshots immediately preceding implementation.

---

## 32. Master assistant-artifact archive

Included under `archive/`:

- `gpt56sol_complete_artifact_snapshot_20260830T124820Z.zip`
- its SHA-256 text file

This is the earlier consolidated assistant-generated artifact snapshot containing hundreds of reports, scripts, results, provenance tables, and historical batch materials through B104I. It intentionally excludes large user-provided raw inputs.

The snapshot was previously reported as:

- 597 unique assistant-generated artifacts
- 900 recorded original occurrences/aliases
- 15 historical delivery bundles expanded/deduplicated

It is the archival evidence base if a later agent needs a detailed intermediate result that is not loose in `critical_artifacts/`.

---

## 33. Critical loose artifacts included in this handoff

The `critical_artifacts/` directory intentionally duplicates a small set of the most useful current outputs so a new conversation does not have to unpack the historical archive first.

### Identity/topology

- complete row-to-Entrez maps
- evidence-rich row map
- 4,301-gene universe
- mapping validation JSON
- independent row-order validation

### Features

- exact feature-generation rule table
- primary and independent validations
- cross-version screen

### Labels

- GOA release 158-169 screen
- inferred 121-column GO map
- independent full-label validation
- UniProt O95073/Q9Y620 historical comparison

### Claims/provenance

- B104G claim-status table
- B104I evidence-versus-literature table
- latest source ledger, input manifest, and provenance events

### Leakage/future work

- gene-level split membership
- leakage experiment backlog

### Human reports

- B104G full-row resolution report
- B104I PPI provenance/Menche/split report
- GraphSAGE GitHub issue tracker

---

## 34. External inputs intentionally omitted from this handoff

The user explicitly said raw/input files can be resupplied as needed. Therefore this handoff does not attempt to bundle large external data.

Likely future inputs include some subset of:

- `graphsage_ppi.zip`
- `dgl_ppi.zip`
- `bio-tissue-networks.tar.gz`
- `bio-tissue-readme.txt`
- GOA release 159 GAF/GPI files
- June 2016 `gp2protein.geneid.gz`
- June 2016 `go.obo`
- minimal MSigDB v6.0 C1/C3 source or original private package for extraction
- GraphSAGE repository/source snapshot `GraphSAGE-a0fdef9.zip`
- `graphsage.utils.py`
- BioSNAP combined/global PPI files if external topology provenance checks are rerun

The new conversation should ask for a missing raw input only when actually needed.

---

## 35. Provenance discipline

Throughout the project the user required:

- SHA-256 tracking of inputs and outputs;
- distinction between user-provided bytes and downloaded bytes;
- source URLs and release metadata;
- append-only provenance/source ledgers;
- explicit deletion clearances when large temporary uploads were no longer needed;
- no guessing about deleted local files;
- no silent mutation of source files;
- careful handling of historical mapping ambiguity;
- exact distinction between direct evidence and inference.

The clean new reproduction should preserve these values but with a much smaller, maintainable implementation rather than carrying every exploratory ledger.

---

## 36. Storage constraints that influenced prior workflow

The user had limited local/cloud upload capacity and frequently deleted large historical inputs after compact derivatives or extraction results were validated.

A screenshot showed 72.3 MB / 512 MB visible Library use even while upload-limit warnings occurred; the visible Library meter was not necessarily the limiting quota.

Large historical UniProt reviewed archives were processed sequentially and deleted after target records were extracted and verified.

The new reproduction should still avoid unnecessary duplication of large files and should support durable local caching of validated sources.

---

## 37. User's desired final products

Eventually the project should have:

1. Clean external-source manifest.
2. Controlled Pixi environment and lock file.
3. Readable Snakemake workflow.
4. Clean Python reconstruction package written from scratch.
5. Pytest suite with synthetic and full integration tests.
6. Exact source -> GraphSAGE -> DGL reproduction.
7. Workflow diagram.
8. Scientific report with confidence categories and literature discrepancies.
9. `claims.csv` machine-readable claim register.
10. Archived source copies where licensing permits.
11. Impact/citation catalog for downstream uses of the PPI benchmark.
12. Leakage/randomized-label experiments.
13. Stable GitHub release and later responses to old GraphSAGE issues.

---

## 38. Recommended immediate next steps in a new conversation

### Main path

Begin implementation of `reproduction/` from scratch, starting with the skeleton and source acquisition layer, then topology/node-order reconstruction.

The user asked explicitly for readable, modular, heavily commented code suitable for human reviewers and GitHub Actions.

Do not reuse the old exploratory code verbatim; use it only as evidence/reference while writing clean implementations.

### Optional side path

Before implementation, the user is curious whether `ppi-walks.txt` can be reproduced or its random seed/configuration inferred.

If doing that side task:

- inspect full uploaded `GraphSAGE-a0fdef9.zip` and target walk file from `graphsage_ppi.zip`;
- replicate likely Python 2 / NetworkX <=1.11 semantics;
- derive exact node list and neighbor iteration order;
- use target-prefix matching to reject candidates very early;
- search plausible seed families and random-state assumptions efficiently;
- consider Python version/state differences and whether code was invoked in a larger process after prior random calls;
- document evidence even if exact seed remains unrecoverable;
- stop before this becomes a distraction from the core reproduction.

---

## 39. Ready-to-paste seed for a new ChatGPT conversation

Use the separate `NEW_CONVERSATION_SEED_20260904T115559Z.md` for a concise first message. The full handoff file is the detailed reference.

---

## 40. Final current-state statement

The most important current scientific result is:

> A complete, globally consistent biological reconstruction now exists for the deterministic GraphSAGE PPI benchmark: the 24 graph topologies map exactly to OhmNet tissue networks; all 56,944 anonymous row identities are resolved to 4,301 Entrez Gene IDs under a strongly supported legacy NetworkX/CPython 2 string-dictionary ordering mechanism; all 2,847,200 feature cells are exactly reproduced from a simple MSigDB C1/C3 source-order selection hypothesis; and all 6,890,224 GO-label cells are exactly reproduced by a fixed GOA release-159 mapping/filtering/`is_a`-propagation policy. The DGL dataset is a deterministic downstream transformation to be reimplemented cleanly. The remaining uncertainties concern historical implementation details and observationally equivalent choices rather than unexplained data mismatches.

The next major engineering/scientific task is to turn this forensic reconstruction into a clean, minimal, independently runnable reproduction package.

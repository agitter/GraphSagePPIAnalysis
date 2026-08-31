# GraphSAGE PPI Provenance Investigation: History and Exploration Log

**Snapshot date:** 2026-08-31 15:53 CDT  
**UTC datestamp:** 2026-08-31T20:53:14Z  
**Status:** Rough internal project snapshot; not the final external report  

## Purpose of this report

This document gives a human colleague the history of the investigation: what information was available at the beginning, what questions we were trying to answer, which analyses were attempted in sequence, what failed, what succeeded, and how the working conclusions changed. It deliberately preserves dead ends and superseded interpretations because those explain why the current reconstruction is credible and where caution is still required.

A separate companion report, `CURRENT_FINDINGS_AND_EVIDENCE_SUMMARY_20260831T205314Z.md`, summarizes only the present state of knowledge.

## 1. Starting point

The project began with a provenance and reproducibility problem involving the protein–protein interaction benchmark released with GraphSAGE and later repackaged by DGL. The benchmark contains:

- 24 graph blocks interpreted as tissue-specific human interaction networks;
- 56,944 node rows;
- 818,716 stored undirected links;
- 50 binary input features;
- 121 binary output labels;
- train, validation, and test partitions defined at the graph level.

The main motivation was that the benchmark appeared to allow very high prediction accuracy through repeated biological entities across graph splits. Understanding that issue required recovering the hidden biological identity of the anonymous node rows and reconstructing the provenance of the graph edges, features, labels, and DGL conversion.

### 1.1 Information initially supplied

The initial investigation summary from another agent proposed several apparently resolved claims:

- the anonymous GraphSAGE nodes could mostly be mapped to Entrez Gene IDs using graph topology and MSigDB features;
- the 24 graph topologies came from OhmNet tissue networks;
- the 50 features came from MSigDB C1 and C3 rather than C7;
- DGL transformed the GraphSAGE archive by reordering nodes, standardizing features, directing edges, and adding self-loops;
- the graph-level split caused extensive repeated-gene leakage;
- the 121 labels were probably GO annotations, but the exact source and transformation remained unresolved.

Those claims were treated as hypotheses, not accepted inputs.

### 1.2 Files and source families available

Over the investigation, the supplied and locally acquired materials included:

- `graphsage_ppi.zip` and `dgl_ppi.zip`;
- OhmNet tissue-network, tissue-label, hierarchy, and README archives;
- the GraphSAGE, OhmNet, Greene et al. 2015, and related manuscripts and supplements;
- Greene supplementary Tables 6 and 9;
- MSigDB releases 5.0, 5.1, 5.2, the 5.2 chip files, and 6.0;
- GOA human annotation releases around 2016, including the GAF, GPAD/GPA, and GPI formats;
- the June 2016 GO ontology in OBO and OWL forms;
- `gp2protein.geneid` and `gp2protein.human` mappings;
- historical and current NCBI `gene2go` snapshots or derivatives;
- Bioconductor `org.Hs.eg.db` releases 3.0.0, 3.1.2, 3.3.0, and 3.4.0;
- HumanBase/GIANT files and sample data;
- date-matched UniProt Swiss-Prot releases for targeted identifier audits;
- a historical commit-pinned human GO annotation derivative from the dhimmel `gene-ontology` repository;
- source listings, web-archive leads, GitHub issues, and prior investigator notes.

### 1.3 Initial goals

The requested work had five central requirements:

1. Independently rederive every transformation claimed to be resolved.
2. Continue searching for the exact Entrez-to-GO label-generation procedure without assuming that a missing `gene2go` snapshot was the answer.
3. Handle ordinary bioinformatics ambiguity, including many-to-many Entrez–UniProt mappings, without arbitrarily forcing one-to-one mappings.
4. Track every input, URL, release date, checksum, generated artifact, failed acquisition, and deletion event so a later reproduction package could be built.
5. Determine whether the train/test design leaked gene identity strongly enough to explain the benchmark's unusually high predictive performance.

## 2. Investigation principles that emerged

Several operating principles became important as the work progressed.

### 2.1 Do not import conclusions from the earlier summary

The prior node mapping, feature mapping, tissue mapping, and label mapping were not used as ground truth. They were independently reconstructed from raw files and only compared afterward.

### 2.2 Separate exact output reproduction from historical proof

A transformation can reproduce every observed byte or matrix cell without proving that the original authors used exactly the same code or source file. The project now distinguishes:

- **Byte-level exact:** direct equality or complete reconstruction from retained files.
- **Strongly inferred:** one historically plausible mechanism explains all observations and survives independent controls.
- **Documented:** explicitly stated in a manuscript, source archive, official page, or historical source code.
- **Open:** multiple histories remain observationally equivalent or an intermediate artifact is missing.

### 2.3 Preserve mapping ambiguity

Entrez Gene IDs and protein accessions are not one-to-one. Mapping was therefore represented as a bipartite graph of components. Ambiguous components were retained, examined before filtering to GraphSAGE genes, and resolved only when primary symbols or a unique component structure justified doing so.

### 2.4 Use batch uploads with deletion receipts

Because large historical data files could not all be uploaded at once, the work adopted batch IDs, SHA-256 verification, compact retained derivatives, provenance events, and explicit “safe to delete” notices. This was also used to correct earlier incomplete manifests and to distinguish local holdings from files actually materialized in the analysis runtime.

## 3. Sequential exploration

## 3.1 Initial independent reproduction and correction of the project infrastructure

The first pass inspected the GraphSAGE, DGL, OhmNet, MSigDB, and manuscript files. It reproduced substantial portions of the earlier investigation but also exposed problems in the first generated deliverables:

- the original master report referred to diagnostics that were not included;
- the manifest did not record URLs for downloaded data;
- actual inputs, candidate inputs, generated outputs, and web references had been mixed together;
- several conclusions were stated more strongly than the supporting files justified.

The reporting and provenance structure was rebuilt. Separate artifacts were created for actual inputs, candidate sources, generated outputs, execution diagnostics, and append-only provenance events.

Scientific findings at this stage included:

- every GraphSAGE graph block matched a supplied OhmNet tissue network;
- topology alone resolved most node identities;
- 49 nonzero feature columns matched MSigDB C1/C3 sets, while one column was all zero and could not yet be named from its restricted vector alone;
- the DGL archive could be reproduced, but the previously stated description of its graph IDs and feature dtype was wrong;
- a conservative gene-identity lookup already achieved approximately 0.994 micro-F1 on test data.

## 3.2 B101: GOA release 159 formats and direct annotation screening

The first focused GO batch contained:

- `goa_human.gaf.159.gz`;
- `goa_human.gpa.159.gz`;
- `goa_human.gpi.159.gz`.

The analysis first reconciled the formats rather than assuming they were interchangeable:

- GAF and GPAD represented the same projected assertion set after translating ECO evidence terms and relations;
- the extra GPAD rows were explained by multiple ECO terms collapsing to the same traditional GAF evidence code;
- every annotated object appeared in GPI;
- GPI's `DB_Xrefs` and `Parent_Object_ID` fields were empty, so GPI did not itself provide the Entrez crosswalk.

A large direct-label grid varied evidence filters, mappings, term restrictions, and comparison universes. Direct, non-propagated GOA annotations did not reproduce any label column at 95% agreement. This ruled out a simple direct-membership construction but did not yet identify the missing transformation.

## 3.3 B102: historical `gp2protein` mapping

The next batch added the June 2016 GO `gp2protein` files and README.

Key results:

- `gp2protein.human` was a historical human UniProt accession set, essentially a UniProt self-map rather than an Entrez mapping;
- `gp2protein.geneid` was an all-species GeneID–UniProt mapping;
- filtering through the human accession set yielded a historical human many-to-many crosswalk;
- most GraphSAGE genes were covered, but several accessions and GeneIDs were ambiguous or missing;
- a direct GOA label screen still failed badly.

This established that identifier conversion mattered but was not sufficient. It also motivated explicit handling of many-to-many components instead of dropping ambiguous accessions.

## 3.4 B103: ontology propagation and evidence-code search

The June 2016 GO ontology was then parsed for:

- primary and alternate GO IDs;
- namespaces;
- obsolete terms and replacements;
- `is_a` relations;
- `part_of` and other relation variants.

A broad search varied:

- evidence-code subsets;
- direct annotations versus propagated annotations;
- `is_a` versus `is_a + part_of` propagation;
- alternative mapping policies.

The first major breakthrough was that the following policy gave a near-exact fit:

- evidence codes `EXP, IDA, IEP, IGI, IMP, ISS`;
- `is_a` propagation only;
- an ambiguity-aware historical GeneID–UniProt projection.

It produced 89 exact columns and 901 excess predicted gene-label cells, with no missing positive labels. This also showed that the labels were not all Biological Process terms: the candidate set included Biological Process, Cellular Component, and Molecular Function.

At this point the leading hypothesis was that ontology drift or an earlier annotation snapshot explained the remaining 901 differences.

## 3.5 B104: release 158 comparison and residual analysis

GOA release 158 was compared with release 159 under the same mapping and ontology logic.

Release 159 was much closer:

- release 158 had both false positives and hundreds of false negatives;
- release 159 removed all false negatives but retained 901 false positives.

The remaining differences were concentrated in broad, shallow Cellular Component and Molecular Function terms. This initially strengthened the ontology-drift hypothesis because many errors arose through ancestor propagation. Alternative tests involving annotation dates, assigned-by sources, NOT annotations, and `part_of` propagation did not resolve them.

## 3.6 B104A: qualifier semantics and the O95073/Q9Y620 mapping component

This was the decisive label-reconstruction breakthrough.

The 901 excess predictions were decomposed by GAF relation semantics:

- 501 depended on `colocalizes_with`;
- 387 depended on `contributes_to`;
- 13 remained under default relations.

Treating only the default aspect relations as ordinary binary membership—`involved_in`, `part_of`, and `enables`—removed 888 differences without losing any observed positive labels.

The final 13 all came from one historical mapping component:

- `Q9Y620`, symbol RAD54B, mapped to GeneID 25788;
- `O95073`, symbol FSBP, mapped historically to both GeneID 100861412 and GeneID 25788.

The initial component builder had filtered to GraphSAGE GeneIDs too early, hiding GeneID 100861412 and incorrectly transferring FSBP annotations to the GraphSAGE RAD54B node. Resolving the full component by concordant primary symbols removed the remaining 13 cells.

On the then-resolved 4,268-gene subset, all 121 label columns and all 516,428 cells matched exactly. A separate set-based implementation independently reproduced the result.

This phase also corrected an important interpretation: broad ontology terms accumulated more errors because qualified descendant annotations propagated to them, but ontology depth was an amplifier rather than the root cause.

## 3.7 B104B: recovering the GO-term selection rule

The next question was why these 121 terms had been selected.

Counting annotations only within the GraphSAGE gene universe did not explain the selection. The key was to apply the exact transformation to the full historical human GOA universe.

Under the accepted release-159 transformation:

- the candidate columns were exactly the 121 most prevalent full-human propagated GO terms;
- they were also exactly the terms with at least approximately 1,000 annotated human genes or proteins;
- the boundary was sharp: rank 121 had 1,007 mapped genes and rank 122 had 997.

A grid of alternative evidence, qualifier, and propagation policies showed that only the same policy that matched the label matrix also recovered the exact 121-term prevalence set.

The matrix contained 121 columns but only 118 distinct membership vectors. Three pairs of GO terms happened to have identical membership on the GraphSAGE gene universe, leaving six column identities unresolved by values alone.

Direct MSigDB memberships from versions 5.1, 5.2, and 6.0 still failed to reproduce the labels. This suggested that MSigDB might have supplied names or a candidate list, while memberships were regenerated from GOA or an equivalent annotation product.

## 3.8 B104C: MSigDB 5.0 and the label-column ordering fingerprint

MSigDB 5.0 was added. It also failed to reproduce any complete label column directly, extending the negative result across versions 5.0 through 6.0.

Attention then shifted to the irregular order of the 121 label columns.

A simulator implemented the historical 64-bit CPython 2.7 string hash, dictionary probing, resizing, and table-order iteration. As an independent positive control, inserting string keys `"0"` through `"56943"` reproduced the full irregular order of all 56,944 keys in `ppi-class_map.json` exactly.

Ordinary biological orderings—GO number, name, prevalence, ontology order, first GAF appearance, and MSigDB order—performed poorly. Large Python 2 dictionaries keyed by GO ID strings reproduced much of the label order, with a best longest common subsequence of 94/121 after later model expansion.

All high-scoring models chose the same orientation for the three duplicate-vector pairs. This yielded a fully named provisional 121-column map, but the orientation remained strongly supported rather than proven because no model reproduced all absolute positions.

## 3.9 B104D: expanded ordering search and low-storage source-date scripts

Hundreds of additional order models were tested, varying:

- key representation;
- direct versus propagated dictionaries;
- row ordering;
- ancestor insertion order;
- mapping scope;
- nested dictionary and set constructions;
- accession-, symbol-, GeneID-, and tuple-keyed alternatives.

No plausible model perfectly matched the 121-column order. The best remained LCS 94/121, with strong pairwise concordance but only 13 exact absolute positions.

This strengthened the claim that the column order has a legacy Python dictionary fingerprint, while leaving the exact dictionary construction open.

Low-storage scripts were created to:

- screen GOA releases 158–169 one pair at a time;
- download large historical UniProt Swiss-Prot releases sequentially, extract only O95073 and Q9Y620, record provenance, and delete each parent archive after verification.

## 3.10 B104E: GOA date screen, date-matched UniProt audit, and feature selection

The GOA date screen showed that release 159 was uniquely exact among releases 158–169 under the fixed transformation. The closest alternative, release 160, still had 652 gene-label differences.

Release 159 was also the only tested release that simultaneously produced:

- the exact label matrix;
- the exact top-121 candidate set;
- exactly 121 full-human terms at the natural approximately 1,000-gene threshold.

The date-matched UniProt audit corrected the interpretation of O95073:

- official reviewed UniProt releases 2016_04, 2016_05, and 2016_06 all mapped O95073/FSBP to both GeneID 100861412 and GeneID 25788;
- Q9Y620/RAD54B mapped to GeneID 25788;
- therefore O95073→25788 was a genuine historical cross-reference, not a parser error;
- the exact GraphSAGE labels nevertheless behaved as if FSBP annotations were not assigned to the RAD54B node.

This left the original mechanism open: the authors may have used an Entrez-native source, a symbol-aware mapping, or another curated intermediate.

The same phase found a simple candidate procedure that reproduced the 50 feature matrix exactly on all then-resolved rows:

- process large C1 sets in source order, then C3 sets;
- use a source-level size cutoff near 200 genes;
- stop after 50 sets.

This selected 30 C1 and 20 C3 sets. The all-zero column was identified as `chryq11`, which has 204 source genes but none in the GraphSAGE gene universe. That zero column strongly supports source-level filtering and little or no post-construction removal of empty columns.

However, the same 50 feature vectors were produced by MSigDB 5.0, 5.1, 5.2, and 6.0. The earlier claim that v5.2 was uniquely identifiable was withdrawn.

## 3.11 B104F: clarification for independent verification and the remaining 533 rows

An independent agent reproduced only part of the label result and reported 12 apparent error genes. Investigation showed that the 121/121 statement had initially been misunderstood:

- 56,411 rows had individually resolved Entrez identities;
- 533 rows remained in 183 topology/feature-equivalence classes;
- within each class, observed and predicted label-vector multisets matched exactly, but individual row identities were not yet assigned.

The fixed GOA transformation was documented in sufficient detail for independent reproduction, including exact file hashes, evidence codes, qualifier handling, and mapping rules.

Four earlier “TSV” downloads were found to be saved GitHub HTML pages and were rejected. A validated raw downloader was supplied, and genuine commit-pinned TSVs were later analyzed.

## 3.12 B104G: complete node-row reconstruction

The remaining row ambiguity was resolved by reconstructing the legacy node order.

For each selected OhmNet tissue edgelist:

1. read edges in original line order;
2. retain node tokens as strings;
3. insert them into a simulated 64-bit, unrandomized CPython 2 dictionary as the graph is built;
4. iterate occupied dictionary slots in table order.

This ordering agreed with all 56,411 previously independent row identities and disagreed with none. In contrast, first-appearance order, numeric sorting, reverse sorting, and integer-key dictionary order matched only a few dozen rows.

The model deterministically assigned Entrez GeneIDs to all 56,944 rows and resolved the 533 formerly ambiguous rows without consulting labels.

Under the full mapping:

- all 24 node and edge sets matched their OhmNet tissue networks exactly;
- all 2,847,200 feature cells matched exactly;
- all 6,890,224 label cells matched exactly;
- repeated occurrences of the same GeneID had no conflicting feature or label vectors.

Independent implementations validated row ordering, features, and labels.

The genuine dhimmel Entrez-native annotation files were also tested. None reproduced any full label column exactly. They were useful controls, especially for showing that an Entrez-native source naturally separates RAD54B and FSBP, but they were not the complete label source.

## 3.13 B104H: exportable node mapping and remaining-source prioritization

A compact node-to-Entrez table and an evidence-rich companion table were produced. The mapping records:

- GraphSAGE node/row ID;
- tissue graph;
- local row index;
- split;
- Entrez GeneID;
- topology, feature, label, and dictionary-order evidence.

The evidence was divided into tiers:

- 55,878 rows resolved by topology alone;
- 533 more resolved by topology plus independently identified MSigDB features;
- the remaining 533 resolved by the inferred Python 2 string-dictionary order and then verified against all downstream data.

A broader search for alternate public copies containing label names did not find an ordered 121-name artifact. Several old GraphSAGE GitHub issues were cataloged for later response.

A leakage-experiment plan was saved, including gene-level label randomization and a truly gene-disjoint split.

## 3.14 B104I: PPI provenance, Menche comparison, and split semantics

A small BioSNAP audit was run against external source files.

It established that:

- the union of all 144 supplied OhmNet tissue edgelists contains 70,338 unique undirected pairs;
- that union exactly equals the BioSNAP combined OhmNet network;
- all 70,338 pairs occur in the 342,353-edge BioSNAP global human interactome.

This externally verifies the chain from the released OhmNet tissue files to the broader BioSNAP interactome.

The Menche et al. 2015 interactome has different dimensions—13,460 proteins and 141,296 interactions—so it was not used unchanged as the OhmNet/BioSNAP 21,557-node, 342,353-edge network. It is best understood as part of the methodological or source lineage of a later enlarged composite interactome.

This also corrected the shorthand “BioGRID network.” The upstream interactome was assembled from multiple resources; BioGRID was one contributor.

Split analysis showed that `split` is one scalar value per tissue-specific node row, but the same Entrez GeneID often appears in multiple splits. Of 4,301 distinct genes, 3,599 occur in more than one split and 3,067 occur in all three.

## 3.15 Repository and provenance organization

A complete `gpt56sol/` artifact snapshot was assembled with:

- reports;
- scripts;
- retained normalized inputs;
- result tables;
- validations;
- source ledgers;
- provenance-event histories;
- deletion receipts;
- superseded artifacts retained for auditability;
- an inventory and file-placement plan.

The project repository was then organized into:

```text
README.md
data/
gpt56sol/
opus46/
papers/
```

Bulk external data and papers are intended to be ignored by Git, while scripts, compact results, provenance tables, reports, and validation artifacts remain versioned. External data will be archived separately where licensing permits.

## 4. Important negative results and dead ends

The following routes were tested and did not explain the released data by themselves:

- direct MSigDB label memberships in versions 5.0, 5.1, 5.2, or 6.0;
- direct use of supplied OhmNet tissue-specific label files;
- unions of OhmNet label files across selected or all tissues;
- Greene Tables 6 or 9 as direct gene-to-GO membership matrices;
- direct, non-propagated GOA membership;
- GO ontology `part_of` propagation;
- treating `colocalizes_with` and `contributes_to` as ordinary binary membership;
- adding IPI to the successful six-code evidence set;
- arbitrary one-to-one reduction of ambiguous Entrez–UniProt mappings;
- dropping all ambiguous mappings;
- broad synonym-based identifier mapping;
- annotation-date or assigned-by-source filters as the main explanation;
- ontology drift as the primary explanation for the former 901 differences;
- ordinary node orders such as sorted GeneID or first appearance;
- ordinary label orders such as GO number, GO name, prevalence, OBO order, or MSigDB order;
- the dhimmel direct/inferred human GO tables as the exact label source;
- simple edge-count thresholds as a complete explanation of why 24 of 144 tissues were selected;
- general web and repository searches for a public copy containing the ordered 121 GO label names.

These negative results are valuable because they constrain what the original preprocessing could have done.

## 5. Major corrections made during the investigation

| Earlier interpretation | Current interpretation |
|---|---|
| The PPI topology was a BioGRID network. | The immediate source is OhmNet; the upstream global interactome is a composite resource that includes BioGRID among several sources. |
| Approximately 4,510 genes and 4,278 identified genes. | The 24 selected tissues contain 4,301 distinct Entrez GeneIDs; all 56,944 rows now have a deterministic reconstructed mapping. |
| All 50 features were uniquely identifiable from v5.2. | The full matrix is exact, but the apparent procedure is strongly inferred; `chryq11` is identified through the source-level rule, and versions 5.0–6.0 produce the same selected vectors. |
| The 121 labels were all Biological Process terms. | The candidate set contains 85 BP, 26 CC, and 10 MF terms. |
| A 121/121 label match covered the entire released matrix. | Initially it covered 56,411 individually resolved rows; after the legacy node-order reconstruction it covers all 56,944 rows. |
| O95073→25788 was probably a bad historical mapping. | It was an official UniProt cross-reference in multiple 2016 releases; the label generator nevertheless separated FSBP annotations from the RAD54B node. |
| DGL used float32 features and simply stable-sorted an existing GraphSAGE graph ID. | DGL stores float64 standardized features and constructs graph IDs using largest-connected-component assignments plus aggregation of non-largest components. |
| GOA release 159 was merely one close match. | Under the fixed transformation, release 159 is uniquely exact among tested releases 158–169. |
| High performance was approximately 0.995 from a partial mapping. | With the full row map, a trivial GeneID lookup reaches test micro-F1 0.9971784. |

## 6. Current transition to the next phase

The exploratory phase has largely achieved the central provenance reconstruction. The project is now moving toward:

1. a clean end-to-end reproduction script;
2. `sources.yaml` with URLs, versions, hashes, licenses, and archival status;
3. `claims.csv` with claim text, evidence category, supporting artifacts, literature comparison, and uncertainty;
4. a workflow diagram from external sources through GraphSAGE and DGL;
5. a report distinguishing exact output reproduction from inferred historical mechanisms;
6. external-data archiving where licenses permit;
7. a catalog of downstream PPI benchmark usage;
8. leakage controls using original labels, gene-level randomized labels, and gene-disjoint splits;
9. eventual responses to long-standing GraphSAGE GitHub provenance questions.

## 7. Recommended companion artifacts

The most useful current internal references are:

- `gpt56sol/batches/B104G_20260829T150633Z/B104G_REPORT_20260829T151452Z.md`
- `gpt56sol/batches/B104H_20260830T110259Z/GRAPHSAGE_NODE_MAPPING_README_20260830T110259Z.md`
- `gpt56sol/batches/B104I_20260830T114918Z/B104I_PPI_PROVENANCE_MENCHE_AND_SPLIT_REPORT_20260830T114918Z.md`
- `gpt56sol/batches/B104I_20260830T114918Z/B104I_EVIDENCE_VS_LITERATURE_REGISTER_20260830T114918Z.csv`
- `gpt56sol/batches/B104E_20260829T121535Z/B104E_GOA_DATE_SCREEN_UNIPROT_AND_FEATURE_RULE_REPORT_20260829T121535Z.md`
- `gpt56sol/batches/B104D_20260829T010311Z/B104D_COLUMN_ORDER_DATE_RANGE_AND_UNIPROT_REPORT_20260829T010311Z.md`
- `gpt56sol/batches/B104A_20260828T145842Z/B104A_REPORT_20260828T145842Z.md`
- `gpt56sol/results/source_ledger_through_B104I_FINAL_20260830T114918Z.csv`
- `gpt56sol/results/provenance_events_through_B104I_20260830T114918Z.csv`

## 8. Bottom line of the exploration history

The investigation moved from a partially mapped anonymous benchmark with an unresolved GO-label source to a complete forensic reconstruction of every GraphSAGE node row and exact reproduction of all graph, feature, label, and DGL outputs. The strongest discoveries came not from finding the original preprocessing script, but from combining historical biological data with deterministic behavior of the legacy Python/NetworkX environment.

The remaining uncertainty is mostly historical rather than numerical: several exact outputs can be reproduced, but some original source-code choices and ordered intermediate artifacts are no longer available.

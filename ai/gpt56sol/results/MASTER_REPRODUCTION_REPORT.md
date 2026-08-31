# GraphSAGE PPI benchmark provenance — corrected independent reproduction report

Generated: 2026-08-24T16:10:11.352732+00:00

## Status of the previous deliverables

The previous `MASTER_REPRODUCTION_REPORT.md`, `RUN_STATUS.md`, and manifest should not be used. The report contained dangling references to error sections that were not embedded, and its manifest omitted source URLs. This corrected report is self-contained, links every table it relies on, and includes a complete execution-diagnostics section. The corrected `actual_input_file_manifest.csv` contains only actual inputs; the broader `source_ledger.csv` separately records web references, historical candidates, and generated outputs.

## Scope and independence rules

The supplied investigation summary was used only as a list of hypotheses to test. It was never imported as a node map, feature map, tissue map, or label map. GraphSAGE labels were not used to resolve gene identities. Gene identities were inferred from graph topology and then, only for residual topology-equivalent classes, from independently identified MSigDB v5.2 feature memberships.

The corrected run directly analyzes the supplied GraphSAGE archive, DGL archive, OhmNet network and label archives, three MSigDB releases, and the Greene supplementary tables. Historical GOA/NCBI/Bioconductor candidates are not described as tested unless their bytes were available.

## Executive verdicts

| Claim from prior investigation | Corrected verdict | Independent result |
| --- | --- | --- |
| Dataset dimensions | Verified | 56,944 rows, 818,716 stored undirected links, 50 features, 121 labels, 24 blocks. |
| Immediate graph source | Verified | Each of the 24 anonymous GraphSAGE blocks has a unique exact node/edge-statistic match to one supplied OhmNet tissue network; mapped edges then verify with zero mismatches. |
| Upstream PPI source = BioGRID alone | Not verified; prior wording too narrow | The OhmNet paper describes a composite physical-interaction network assembled from several resources. The supplied edge files verify OhmNet as the immediate source, not BioGRID as the sole upstream database. |
| Gene identity recovery | Partially verified | Topology alone identifies 55,878 rows and 4,232 genes. MSigDB features resolve another 533 rows, yielding 4,268 of the 4,301-gene 24-tissue universe. The earlier 4,278/~4,510 claim is not reproduced. |
| Feature provenance | Verified with one qualification | All 49 nonzero columns uniquely match MSigDB v5.2 sets: 29 C1 and 20 C3. Column 10 is all-zero and therefore cannot be assigned a unique gene-set name from the observations; '50/50 uniquely matched' was overstated. |
| DGL transformation | Verified, earlier description corrected | All labels, graph IDs, transformed features, and complete directed edge sets reproduce exactly. DGL stores float64 features, and its graph-ID construction assigns each tissue LCC separately while aggregating non-LCC components into the first graph ID of each split. |
| Leakage mechanism | Verified conservatively | Using only identities resolved without consulting labels, 5,430/5,524 test rows (98.2983%) occur in training and all have identical label vectors. Zero-filled unresolved/unseen predictions give micro-F1 0.9940678477. |
| Labels are direct MSigDB gene sets | Rejected for supplied releases | No exact match across MSigDB 5.1, 5.2, or 6.0. v5.2/v6.0 have only one column at ≥99% agreement and two at ≥95% when all collections are searched. |
| Labels are supplied OhmNet tissue labels | Rejected for tested transformations | No exact, ≥99%, or ≥95% matches for individual label files, same-tissue comparisons, all-tissue GO unions, or selected-24-tissue GO unions. |
| Exact public GO/gene2go/Bioconductor source | Open | The historical files are recorded with candidate URLs but were not materialized in this runtime. The corrected report does not repeat the earlier claim that a missing mid-2016 gene2go file is necessarily the cause. |

## 1. Inputs and traceability

Every actual input was hashed. Public inputs have a direct or canonical URL where one exists; MSigDB uses the official authenticated download page. The complete URL and retrieval-status fields are in `actual_input_file_manifest.csv` and `source_ledger.csv`.

| Input | Role | Direct/canonical source | SHA-256 |
| --- | --- | --- | --- |
| `graphsage_ppi.zip` | primary dataset under investigation | https://snap.stanford.edu/graphsage/ppi.zip | `53aeb76e54fd41b645e7edb48b62929240b89839495396b048086fd212503fbd` |
| `dgl_ppi.zip` | downstream DGL dataset to reproduce | https://data.dgl.ai/dataset/ppi.zip | `1f5b2b09ac0f897fa6aa1338c64ab75a5473674cbba89380120bede8cddb2a6a` |
| `bio-tissue-networks.tar.gz` | OhmNet tissue PPI source candidate | https://snap.stanford.edu/ohmnet/bio-tissue-networks.tar.gz | `2c79e17f4a7c8680a7cbf8b20cef4acf356a7523c9a75fce586494153c0603d1` |
| `bio-tissue-labels.tar.gz` | OhmNet tissue-specific GO labels | https://snap.stanford.edu/ohmnet/bio-tissue-labels.tar.gz | `6abf272940d2407849bd779e5f85c0377a2fb07c2351d1ebc82e3d06a46bc11d` |
| `bio-tissue-hierarchy.tar.gz` | OhmNet tissue hierarchy; retained for later work | https://snap.stanford.edu/ohmnet/bio-tissue-hierarchy.tar.gz | `c4568a68bb83319bff854eecf73a93f698fe2c41ed6e95639af974dd024ffef7` |
| `bio-tissue-readme.txt` | OhmNet data documentation | https://snap.stanford.edu/ohmnet/bio-tissue-readme.txt | `7f3372f8ae3a90852951c73b18980386ceb4ad2f5d32d81366adf22fd75e2b20` |
| `msigdb_v5.1_files_to_download_locally.zip` | historical MSigDB comparison | https://www.gsea-msigdb.org/gsea/downloads.jsp | `5a8b3f10ea92f8e71eaaa0705ab9d3a5229d838864eec9699544b75884bc9e29` |
| `msigdb_v5.2_files_to_download_locally.zip` | feature recovery and historical label comparison | https://www.gsea-msigdb.org/gsea/downloads.jsp | `a618c1c60b11570036034e6357e73e80ee43065ec7a57c1dbd238f205405fbdb` |
| `msigdb_v5.2_chip_files_to_download_locally.zip` | historical chip mappings; retained for later GO work | https://www.gsea-msigdb.org/gsea/downloads.jsp | `252befc853b5e01cfe99439ec50a7ab8de747cd0abe6224a93077f2d9b0b20fc` |
| `msigdb_v6.0_files_to_download_locally.zip` | later MSigDB comparison | https://www.gsea-msigdb.org/gsea/downloads.jsp | `39fa82c4cedc9183c532afb1c1431683536b5945e50fb8be4f5bcce3ac136edf` |
| `Greene2015.pdf` | Greene et al. manuscript | https://doi.org/10.1038/ng.3259 | `15c734d37bf63dc586d9bfb95673612209a2f2d298a0a1dc84fa63a1d7a17ce2` |
| `Greene2015_sup.pdf` | Greene et al. supplementary information | https://www.nature.com/articles/ng.3259#Sec23 | `89e84c545590a3d34890f24cf6543a336b59cace8926a83701f503afcd979ed9` |
| `Greene2015_Table6.xlsx` | expert-curated GO term list | https://www.nature.com/articles/ng.3259#Sec23 | `691b9d895ac6d0f6ed7abedb96d9b206965fe221e3ccdae940b4daa5db50533e` |
| `Greene2015_Table9.xlsx` | GO-to-tissue mapping | https://www.nature.com/articles/ng.3259#Sec23 | `18ae68f28d9b84f4b1cb7f7c7c1cc8eb76716414de2089a329f070a5aeca6cd5` |
| `OhmNet.pdf` | OhmNet manuscript | https://doi.org/10.1093/bioinformatics/btx252 | `e60daf8341d0e322ce58e7c6ad194f7e4b573df7c8aba1716ad78c98992b02fe` |
| `investigation_summary_2026_08_23.md` | prior hypotheses to independently test | Supplied local artifact; no external URL provided | `fe23f5d35c1c3a21bba13c2241e0c8c783c31a0a31026c8e5f1c0d7a8c320d16` |
| `Pasted markdown(1).md` | copied EBI GOA archive listing | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | `f856e6f29f02a17c83650724f2b8434383dc7952c03e5200d4bb9c47c1c8782e` |
| `historical_go_mapping_inventory.md` | prior source-discovery inventory | Supplied local artifact; no external URL provided | `ad10c6ed2919409a3dda0f2e89f4a5f7709cb739f9e3129297b1fa9f497f4b0c` |

Input verification result: all 18 actual inputs were present and their recorded checksums verified. See `input_verification_log.csv`.

## 2. GraphSAGE archive structure

The supplied archive contains:

- 56,944 node rows.
- 818,716 stored undirected links, including 25,084 self-loops.
- Feature array shape (56944, 50), dtype `float64`.
- Label array shape (56944, 121), dtype `uint8`.
- Split counts: 44,906 training, 6,514 validation, and 5,524 test rows.
- 295 connected components.
- The node JSON has only `id, test, val` fields; it has no `graph_id` field.

The edge-crossing scan finds exactly 24 contiguous blocks separated by cuts crossed by no non-loop edge. Matching each block's node and edge counts against all supplied OhmNet tissues gives one global one-to-one assignment.

## 3. Exact 24-block tissue partition

| Block | Split | GraphSAGE rows | Nodes | Undirected links incl. loops | OhmNet tissue | Counts exact |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | train | 0–1553 | 1554 | 16912 | adipose_tissue | True |
| 2 | train | 1554–2941 | 1388 | 15517 | adrenal_cortex | True |
| 3 | train | 2942–5214 | 2273 | 30839 | adrenal_gland | True |
| 4 | train | 5215–7563 | 2349 | 33731 | amygdala | True |
| 5 | train | 7564–9150 | 1587 | 18871 | aorta | True |
| 6 | train | 9151–10179 | 1029 | 9635 | astrocyte | True |
| 7 | train | 10180–12014 | 1835 | 23041 | artery | True |
| 8 | train | 12015–14514 | 2500 | 36276 | basal_ganglion | True |
| 9 | train | 14515–15116 | 602 | 4193 | basophil | True |
| 10 | train | 15117–18442 | 3326 | 54554 | blood | True |
| 11 | train | 18443–20849 | 2407 | 33201 | blood_plasma | True |
| 12 | train | 20850–22733 | 1884 | 24016 | blood_platelet | True |
| 13 | train | 22734–24558 | 1825 | 23741 | bone | True |
| 14 | train | 24559–28055 | 3497 | 54824 | brain | True |
| 15 | train | 28056–30862 | 2807 | 43854 | colon | True |
| 16 | train | 30863–33195 | 2333 | 30993 | eye | True |
| 17 | train | 33196–35859 | 2664 | 39671 | forebrain | True |
| 18 | train | 35860–38688 | 2829 | 43964 | large_intestine | True |
| 19 | train | 38689–41868 | 3180 | 48409 | liver | True |
| 20 | train | 41869–44905 | 3037 | 46956 | gastrointestinal_tract | True |
| 21 | valid | 44906–48115 | 3210 | 50061 | heart | True |
| 22 | valid | 48116–51419 | 3304 | 52126 | kidney | True |
| 23 | test | 51420–54633 | 3214 | 51666 | lung | True |
| 24 | test | 54634–56943 | 2310 | 31665 | midbrain | True |

The statistic assignment is only a candidate identification step. It is then tested structurally through joint Weisfeiler–Lehman refinement and, after mapping, through exact edge membership. Every one of the 816,950 GraphSAGE edges whose endpoints are independently mapped is present in the assigned OhmNet network; mismatches: 0.

This verifies **OhmNet as the immediate graph source**. It does not verify that BioGRID alone was the upstream source. The OhmNet manuscript describes its global physical interactome as a combination of several interaction resources; the tissue networks are induced from that global physical network using tissue-activity information.

## 4. Gene identity reconstruction

### 4.1 Topology-only stage

For each assigned tissue, the anonymous GraphSAGE block and the Entrez-keyed OhmNet graph were refined jointly. Initial signatures use non-loop degree and self-loop presence; subsequent signatures use the current color and the multiset of neighbor colors until stable.

- Topology-unique rows: 55,878/56,944.
- Unique Entrez IDs identified at this stage: 4,232.
- Residual topology-equivalent classes: 335 classes containing 1,066 rows before feature disambiguation.

### 4.2 Independent feature disambiguation

The topology-unique genes were used to identify the observed feature columns against all Entrez GMT files in supplied MSigDB v5.2. The 49 nonzero columns each have a unique exact observed gene-set membership. Those 49 independently named sets were then used to subdivide residual WL equivalence classes.

- Additional rows resolved: 533.
- Total mapped rows: 56,411/56,944.
- Distinct mapped Entrez IDs: 4,268.
- Distinct Entrez IDs in the assigned 24 OhmNet tissues: 4,301.
- Residual rows: 533 in 183 equivalence classes.
- Candidate genes participating in residual classes: 117; because candidate sets overlap across tissues, this is not the same quantity as the 33-gene gap between mapped and complete tissue-universe counts.

The complete row-level map, candidate sets, method, and confidence are in `graphsage_row_to_entrez_topology_features.csv`. No arbitrary choice is made inside unresolved equivalence classes.

### Corrected verdict

The independent run supports 4,268 mapped Entrez IDs, not the earlier 4,278 figure. It also establishes a 4,301-gene universe across the selected OhmNet tissues, not the approximate 4,510 figure. The earlier counts may have depended on an undisclosed mapping or manual swaps, but they are not accepted here without the corresponding artifact and derivation.

## 5. Feature provenance

| Column | Status | Collection | Gene-set name | Observed positives in topology-unique universe | Best mismatches |
| --- | --- | --- | --- | --- | --- |
| 0 | exact_unique | C1 | chr8q24 | 30 | 0 |
| 1 | exact_unique | C1 | chr14q11 | 30 | 0 |
| 2 | exact_unique | C1 | chr12q24 | 49 | 0 |
| 3 | exact_unique | C1 | chr19p13 | 115 | 0 |
| 4 | exact_unique | C1 | chr11q12 | 24 | 0 |
| 5 | exact_unique | C1 | chr22q13 | 43 | 0 |
| 6 | exact_unique | C1 | chr22q11 | 23 | 0 |
| 7 | exact_unique | C1 | chr1q21 | 50 | 0 |
| 8 | exact_unique | C1 | chr6p21 | 68 | 0 |
| 9 | exact_unique | C1 | chr5q31 | 41 | 0 |
| 10 | all_zero_unidentifiable | — | — | 0 | 0 |
| 11 | exact_unique | C1 | chr19q13 | 113 | 0 |
| 12 | exact_unique | C1 | chr8p23 | 10 | 0 |
| 13 | exact_unique | C1 | chr17p13 | 46 | 0 |
| 14 | exact_unique | C1 | chr14q32 | 21 | 0 |
| 15 | exact_unique | C1 | chr3p21 | 52 | 0 |
| 16 | exact_unique | C1 | chr21q22 | 37 | 0 |
| 17 | exact_unique | C1 | chr9q34 | 48 | 0 |
| 18 | exact_unique | C1 | chr11q13 | 74 | 0 |
| 19 | exact_unique | C1 | chr17q25 | 37 | 0 |
| 20 | exact_unique | C1 | chrxp11 | 32 | 0 |
| 21 | exact_unique | C1 | chr1p36 | 82 | 0 |
| 22 | exact_unique | C1 | chr17q21 | 76 | 0 |
| 23 | exact_unique | C1 | chr11p15 | 63 | 0 |
| 24 | exact_unique | C1 | chr20q13 | 39 | 0 |
| 25 | exact_unique | C1 | chr12q13 | 57 | 0 |
| 26 | exact_unique | C1 | chr16p13 | 51 | 0 |
| 27 | exact_unique | C1 | chr7q11 | 20 | 0 |
| 28 | exact_unique | C1 | chr16p11 | 25 | 0 |
| 29 | exact_unique | C1 | chr12p13 | 31 | 0 |
| 30 | exact_unique | C3 | AAAYRNCTG_UNKNOWN | 128 | 0 |
| 31 | exact_unique | C3 | V$MYOD_01 | 87 | 0 |
| 32 | exact_unique | C3 | V$E47_01 | 94 | 0 |
| 33 | exact_unique | C3 | V$CMYB_01 | 91 | 0 |
| 34 | exact_unique | C3 | V$AP4_01 | 97 | 0 |
| 35 | exact_unique | C3 | AACTTT_UNKNOWN | 666 | 0 |
| 36 | exact_unique | C3 | V$ELK1_01 | 100 | 0 |
| 37 | exact_unique | C3 | V$SP1_01 | 87 | 0 |
| 38 | exact_unique | C3 | V$ATF_01 | 86 | 0 |
| 39 | exact_unique | C3 | V$ELK1_02 | 67 | 0 |
| 40 | exact_unique | C3 | V$RSRFC4_01 | 82 | 0 |
| 41 | exact_unique | C3 | V$CETS1P54_01 | 86 | 0 |
| 42 | exact_unique | C3 | V$P300_01 | 98 | 0 |
| 43 | exact_unique | C3 | V$NFE2_01 | 96 | 0 |
| 44 | exact_unique | C3 | V$CREB_01 | 84 | 0 |
| 45 | exact_unique | C3 | V$CREBP1CJUN_01 | 84 | 0 |
| 46 | exact_unique | C3 | V$SOX5_01 | 93 | 0 |
| 47 | exact_unique | C3 | V$E4BP4_01 | 68 | 0 |
| 48 | exact_unique | C3 | V$E2F_02 | 83 | 0 |
| 49 | exact_unique | C3 | V$NFKAPPAB65_01 | 95 | 0 |

The result is:

- 29 uniquely identified nonzero C1 columns.
- 20 uniquely identified nonzero C3 columns.
- Column 10 is zero for every GraphSAGE row. Many external gene sets can induce an all-zero vector after restriction to this graph universe, so the column cannot be uniquely named from the dataset itself.

Thus the observable feature provenance is C1+C3. The stronger earlier statement that all 50 columns were uniquely matched is not justified; 49 are uniquely identified and one is observationally unidentifiable.

## 6. Exact DGL transformation

The DGL archive is reproduced by this algorithm:

1. Within each split, initially assign all rows to that split's first graph ID.
2. For every tissue block, assign its largest connected component to that tissue's own graph ID.
3. Leave every non-largest component in the first graph ID of its split.
4. Concatenate by graph ID while preserving original row order inside each ID.
5. Fit `sklearn.preprocessing.StandardScaler` on all GraphSAGE training rows in float64 and transform every split.
6. Convert every non-loop undirected edge into both directions and retain exactly one self-loop per node.

This is more specific than “stable sort by graph_id”: the original GraphSAGE JSON has no graph-ID field, and the DGL graph IDs encode a largest-connected-component rule. The supplied DGL feature arrays are float64, not float32.

| Split | Rows | Archive feature dtype | Max \|difference\| vs float64 StandardScaler | Labels exact | Graph IDs exact | Directed edges | Missing edges | Extra edges |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| train | 44906 | float64 | 1.14e-13 | True | True | 1271274 | 0 | 0 |
| valid | 6514 | float64 | 1.14e-13 | True | True | 205434 | 0 | 0 |
| test | 5524 | float64 | 1.14e-13 | True | True | 167500 | 0 | 0 |

All DGL checks pass. Complete per-component assignments are in `dgl_component_assignment.csv`; complete split checks are in `dgl_split_verification.csv` and `dgl_transformation_verification.json`.

## 7. Leakage measurement

This measurement uses only gene identities obtained from topology and independently identified features. Labels are not used to choose a gene ID.

- Test rows: 5,524.
- Test rows with resolved genes: 5,465.
- Test rows whose Entrez ID occurs in training: 5,430 (98.2983% of all test rows).
- Seen test rows with byte-identical label vectors: 5,430/5,430.
- Gene lookup prediction, with unresolved or unseen rows predicted as all-zero: micro-F1 = 0.9940678477; TP=197,736, FP=0, FN=2,360.

This verifies the leakage mechanism. The earlier 98.8% and 0.9956 values are plausible under a larger node map, but they were not reproduced independently from the current supplied files and therefore are not reported as verified.

## 8. Gene-to-label investigation using locally available sources

Comparison universe: 4,268 independently mapped Entrez IDs and 121 GraphSAGE label columns.

### 8.1 MSigDB screens

| Version | Search scope | Exact | ≥99% | ≥95% | Median best agreement | Minimum best mismatches |
| --- | --- | --- | --- | --- | --- | --- |
| 5.1 | all_collections | 0 | 0 | 0 | 0.807170 | 307 |
| 5.1 | c5_bp | 0 | 0 | 0 | 0.804358 | 312 |
| 5.2 | all_collections | 0 | 1 | 2 | 0.855904 | 10 |
| 5.2 | c5_bp | 0 | 0 | 0 | 0.846532 | 265 |
| 6.0 | all_collections | 0 | 1 | 2 | 0.855904 | 10 |
| 6.0 | c5_bp | 0 | 0 | 0 | 0.846532 | 265 |

No label column exactly equals a gene set in the supplied MSigDB releases. The near match in v5.2/v6.0 all-collection search is not sufficient to identify the 121-column matrix as direct MSigDB membership, and the GO Biological Process collection itself has no ≥95% matches.

### 8.2 OhmNet supplied label files

| Transformation | Exact | ≥99% | ≥95% | Median best agreement | Minimum best mismatches |
| --- | --- | --- | --- | --- | --- |
| OhmNet individual tissue label files, absent genes treated as 0 | 0 | 0 | 0 | 0.781396 | 327 |
| OhmNet union by GO term across all supplied tissues | 0 | 0 | 0 | 0.781396 | 327 |
| OhmNet union by GO term across selected 24 tissues | 0 | 0 | 0 | 0.780928 | 352 |

A separate same-tissue comparison tested 111 OhmNet label files belonging to the selected 24 tissues; no GraphSAGE column reached 95% agreement. Therefore the GraphSAGE labels are not direct copies of the supplied OhmNet tissue-label files under the tested absent-gene and union transformations.

### 8.3 Greene supplementary restrictions

The workbook parser finds 973 unique GO IDs in Table 6 and 6,172 in Table 9, with 276 in common. These tables define candidate term restrictions and tissue associations; they do not themselves supply the Entrez-to-GO membership matrix needed to reconstruct the 121 labels.

### 8.4 What remains open

The corrected run has **not** tested the historical GOA, GPAD, GPI, `gene2go`, GO ontology, `gp2protein`, UniProt mapping, or Bioconductor package combinations numerically because those candidate bytes were not available in the runtime. Their exact candidate URLs and materialization status are recorded in `source_ledger.csv` and `SOURCE_ACQUISITION.md`.

Consequently, this report does not endorse the earlier explanation that the problem is necessarily an unavailable July–August 2016 `gene2go` snapshot. That is one hypothesis, not a verified conclusion. The next source grid must test:

- GAF/GPAD/GPI releases 155–160 with release-matched object metadata.
- Direct GPI xrefs, release-matched `gp2protein.geneid`, historical UniProt mappings, and Bioconductor Entrez–UniProt maps, separately and in controlled unions.
- No propagation, `is_a` only, and a declared safe-relation closure from release-matched GO ontology files.
- Explicit evidence-code subsets, NOT qualifiers, obsolete/alternate GO IDs, annotation-extension handling, canonical-protein versus isoform collapse, and many-to-many identifier mappings.
- Whether the 121 terms were selected before or after restriction to the GraphSAGE/OhmNet gene universe.

## 9. Reproduction commands

From the corrected bundle root:

```bash
python -u scripts/run_core_verification.py \
  --input-dir /mnt/data \
  --work-dir work \
  --output-dir results

python -u scripts/run_local_label_source_screen.py \
  --input-dir /mnt/data \
  --work-dir work \
  --output-dir results

python scripts/build_source_manifest.py
python scripts/build_corrected_report.py
python scripts/build_full_source_ledger.py
python scripts/refresh_run_status.py
python scripts/validate_corrected_bundle.py

python scripts/download_or_verify_sources.py \
  --manifest results/actual_input_file_manifest.csv \
  --dest inputs \
  --verify-only \
  --log results/input_verification_log.csv
```

## 10. Execution diagnostics

### Final accepted runs

- Core verification: PASS, exit 0, 68.656 s, empty stderr.
- Local label-source screen: PASS, exit 0; one nonfatal openpyxl workbook-extension warning is preserved in stderr.
- Source-ledger construction: PASS, exit 0, empty stderr; the final ledger separately records actual inputs, web references, historical candidates, analysis scripts, and stable generated outputs.
- Actual-input verification: PASS, exit 0, 18/18 files verified.
- Corrected report generation: PASS, exit 0, empty stderr.
- Bundle validation: PASS; actual-input hashes, public-input URLs, report references, workbook readability, and absence of missing tracked outputs were checked programmatically.

### Superseded failures

None of the following partial outputs is used:

1. An initial core-verification invocation exceeded a 600-second wrapper timeout after archive extraction. It produced no accepted result. The unchanged analysis was rerun directly and completed in 68.656 seconds.
2. An early label-screen implementation raised `IndexError: list index out of range` while selecting a best C5-BP candidate after duplicate removal. Candidate construction was corrected before any result was accepted.
3. The next attempt raised `ValueError: No candidate records for MSigDB 5.1 C5 BP`. Historical archive-name recognition and specialized-collection retention were corrected.
4. A non-optimized exhaustive label-screen attempt timed out after the MSigDB comparisons but before the OhmNet same-tissue loop finished. Its partial output was discarded. The accepted bitset implementation completed successfully.

The corresponding commands, exit statuses, warnings, stdout/stderr files, and scope limitations are also collected in `EXECUTION_DIAGNOSTICS.md`.

## 11. Output guide

- `actual_input_file_manifest.csv`: **actual inputs only**, with local path, URL, source page, status, size, and SHA-256.
- `source_ledger.csv`: actual inputs plus web references, historical candidates, and generated outputs.
- `SOURCE_ACQUISITION.md`: download/verify commands and readable source tables.
- `EXECUTION_DIAGNOSTICS.md`: complete success, warning, failure, and scope-limit section.
- `core_verification_summary.json`: machine-readable central results.
- `tissue_partition.csv`: complete 24-block assignment.
- `graphsage_row_to_entrez_topology_features.csv`: row-level mapping and unresolved candidates.
- `feature_column_mapping.csv`: 50-column feature provenance evidence.
- `dgl_transformation_verification.json`: exact DGL reconstruction checks.
- `local_label_source_screen_summary.json`: machine-readable local label-source exclusions.
- `msigdb_label_source_screen.csv`, `ohmnet_global_label_source_screen.csv`, and `ohmnet_same_tissue_label_screen_best.csv`: per-column label comparisons.

## 12. Requested files for the next historical GO grid

The highest-priority missing bytes are:

1. `goa_human.gaf.159.gz`, `goa_human.gpa.159.gz`, and `goa_human.gpi.159.gz`.
2. The matching release 158 and 160 triplets.
3. `gene_association.goa_human.157.gz`, `gp_association.goa_human.157.gz`, and `gp_information.goa_human.157.gz`.
4. `2016-06-01-go.obo` or `go-basic.obo`, plus `gp2protein.geneid.gz`.
5. `org.Hs.eg.db_3.3.0.tar.gz` and `org.Hs.eg.db_3.4.0.tar.gz`.
6. Any independently archived human `gene2go` snapshot dated June–September 2016.

The candidate URLs are already in `source_ledger.csv`; upload is needed only for files the included downloader cannot materialize in the execution environment.

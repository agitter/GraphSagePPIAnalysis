# B104E — GOA source-date screen, date-matched UniProt mapping audit, and exact feature-generation rule

## Executive summary

This batch validates two user-side audit packages and adds a new exact reconstruction of the 50 GraphSAGE input features.

The principal findings are:

1. **GOA release 159 is the only tested GOA human release from 158 through 169 that reproduces the complete GraphSAGE label matrix exactly under the fixed, previously established transformation.** Release 159 gives 121/121 exact columns and zero differences. The closest alternative, release 160, has 652 gene-label differences.
2. **Changing the GOA release date does not solve the remaining column-order problem.** The best tested order model has LCS 94/121 for releases 158–166 and LCS 93/121 for releases 167–169. No tested release produces an exact 121-column order.
3. **The provisional orientation of the three duplicated-vector column pairs is stable across every tested release 158–169.** This increases confidence in the orientation but does not prove the original ordered source list.
4. **The official reviewed UniProt records for April, May, and June 2016 explicitly map O95073/FSBP to both GeneID 100861412 and GeneID 25788.** Q9Y620/RAD54B maps to GeneID 25788. Therefore, the historical O95073→25788 edge is real in UniProt and must not be dismissed as a parser artifact or a transient bad download.
5. **The exact label fit nevertheless requires FSBP annotations not to be assigned to the GraphSAGE node identified as GeneID 25788/RAD54B.** This is a semantic disambiguation of a genuine many-to-many historical component, not proof that the UniProt cross-reference itself was erroneous. An Entrez-native annotation source such as historical `gene2go`, or a symbol-aware mapping policy, could naturally produce this behavior.
6. **A simple MSigDB threshold-and-source-order rule reproduces all 50 GraphSAGE features exactly on all 56,411 resolved node rows: 2,820,550 binary cells with zero differences.** The rule selects 30 C1 sets and 20 C3 sets; no C7 set is reached before the 50-feature cap.
7. **The selected 50 feature memberships are identical in supplied MSigDB versions 5.0, 5.1, 5.2, and 6.0.** Consequently, the feature matrix does not uniquely identify MSigDB v5.2 as the source version.

## 1. Inputs and integrity

### Uploaded B104E packages

| Uploaded package | Bytes | SHA-256 | Integrity |
|---|---:|---|---|
| `uniprot_2016_mapping_audit_ledger.zip` | 30,449 | `4ac8cb1a900215ded9dc35a4fc44a4abaaff3e5e774cf62769592f1ab153a7b0` | ZIP test passed; 16 members |
| `goa_date_screen_results.zip` | 54,246 | `502c8ffdb7b809c1665e82d31db52b72e7e855e7e9fcc06a8a2f46a64bc30de9` | ZIP test passed; 15 members |

The UniProt package now contains the requested extracted `.dat`, summary `.tsv`, per-release provenance JSON, official metalink, and audit ledger for each of releases 2016_04, 2016_05, and 2016_06.

The GOA package contains one detailed JSON result for each release 158–169, a release summary, run metadata, and append-only source-integrity/analysis/deletion events.

A programmatic audit executed 81 checks across the two packages. All passed.

## 2. GOA release-date screen

### 2.1 What was held fixed

The screen deliberately varied only the GOA GAF/GPI release. It held the following transformation fixed:

- May-2016 historical `gp2protein` mapping;
- preservation of legitimate many-to-many mappings;
- the previously described semantic separation of O95073/FSBP from GeneID 25788/RAD54B for label projection;
- unique GPI primary-symbol fallback for otherwise unresolved graph genes;
- evidence codes `EXP`, `IDA`, `IEP`, `IGI`, `IMP`, and `ISS`;
- default GAF relations `involved_in`, `part_of`, and `enables` only;
- exclusion of `NOT`, `colocalizes_with`, and `contributes_to` from ordinary binary membership;
- propagation through `is_a` only;
- the archived 2016-06-01 ontology graph;
- the previously recovered candidate GO term for each GraphSAGE column.

This design isolates annotation-release effects. It is not a simultaneous search over GOA release, ontology release, and identifier mapping.

### 2.2 Membership results

| GOA release | Listed date | Exact columns | False positives | False negatives | Total differences | Terms ≥1,000 | Candidate overlap among ≥1,000 | Candidate overlap among top 121 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 158 | 2016-06-07 | 8 | 21 | 825 | 846 | 120 | 120 | 121 |
| **159** | **2016-07-04** | **121** | **0** | **0** | **0** | **121** | **121** | **121** |
| 160 | 2016-09-14 | 15 | 438 | 214 | 652 | 122 | 120 | 120 |
| 161 | 2016-10-03 | 9 | 625 | 238 | 863 | 123 | 120 | 120 |
| 162 | 2016-10-31 | 1 | 894 | 288 | 1,182 | 124 | 120 | 119 |
| 163 | 2016-11-28 | 1 | 1,231 | 295 | 1,526 | 125 | 120 | 119 |
| 164 | 2017-01-16 | 1 | 3,890 | 1,757 | 5,647 | 127 | 120 | 119 |
| 165 | 2017-02-13 | 0 | 4,458 | 2,059 | 6,517 | 130 | 119 | 118 |
| 166 | 2017-03-13 | 0 | 6,385 | 2,099 | 8,484 | 131 | 119 | 118 |
| 167 | 2017-04-10 | 0 | 6,634 | 2,143 | 8,777 | 131 | 119 | 118 |
| 168 | 2017-05-08 | 0 | 6,903 | 2,322 | 9,225 | 132 | 119 | 118 |
| 169 | 2017-06-05 | 0 | 7,456 | 2,326 | 9,782 | 133 | 119 | 118 |

Release 159 is uniquely exact in the tested date range. Release 160 is the next closest, but it still changes 652 of the 516,428 resolved gene-label cells.

This result materially strengthens the release-159 attribution. It should still be phrased as:

> Under the fixed historical mapping, qualifier, evidence, and ontology policy, GOA release 159 is the unique exact GOA source among releases 158–169.

It is not proof that the original authors opened these exact GAF/GPI files. An Entrez-native or preprocessed source could encode the same associations.

### 2.3 Term-selection result

The candidate set is exactly the top 121 full-human terms for releases 158 and 159. However, the natural round threshold of at least 1,000 mapped human genes selects exactly 121 candidate terms only in release 159:

- release 158: rank 121 has 995 genes; only 120 terms meet the threshold;
- release 159: rank 121 has 1,007 genes and rank 122 has 997; exactly 121 terms meet the threshold;
- release 160: 122 terms meet the threshold, and one candidate term drops from the top 121.

This makes the joint hypothesis “release 159 plus a 1,000-gene/protein cutoff” more specific than a generic “take the top 121” rule.

### 2.4 Column order across dates

| Release range | Best LCS | Exact absolute positions | Exact prefix | Duplicate orientation |
|---|---:|---:|---:|---|
| 158–166 | 94 / 121 | 13 | first 5 | `001` |
| 167–169 | 93 / 121 | 7 | first 5 | `001` |

No date in the full practical window produces a perfect order under the current Python-2 dictionary construction model. The order is almost invariant across 158–166, despite substantial changes in annotation membership.

Therefore, GOA source date was not the missing variable needed to recover all 121 positions.

The stable `001` duplicate orientation across all 12 releases is additional supporting evidence for the provisional assignments:

| Earlier column | Provisional term | Later column | Provisional term |
|---:|---|---:|---|
| 24 | `GO:0043228` | 71 | `GO:0043232` |
| 39 | `GO:0006464` | 63 | `GO:0036211` |
| 48 | `GO:1903561` | 70 | `GO:0043230` |

This remains strongly supported rather than proven because every order model still shares substantial assumptions about the original preprocessing environment.

### 2.5 Validation limitation

The user-side screening script downloaded, integrity-checked, analyzed, and deleted each release sequentially. The package preserves source URLs, sizes, SHA-256 values, headers, and deletion events. Releases 158 and 159 reproduce earlier analyses performed from uploaded raw files.

Raw releases 160–169 were not independently downloaded again in this runtime. An attempted independent download of releases 160 and 168 failed because the execution environment had temporary DNS resolution failure. This is recorded as an execution limitation, not as a problem with the supplied results.

## 3. Date-matched UniProt audit of O95073 and Q9Y620

### 3.1 Parent archive provenance

| Release | Date | Parent archive bytes | Official MD5 | Observed MD5 | Parent SHA-256 | Swiss-Prot records scanned | Parent deleted after success |
|---|---|---:|---|---|---|---:|---|
| 2016_04 | 2016-04-13 | 1,516,525,310 | `e607b83de1ac87e6f63b13715c049a3f` | match | `a75cd81a9114141aeb43d31f6a2742d41329d02cbada070015985191dc8ac754` | 550,960 | yes |
| 2016_05 | 2016-05-11 | 1,504,161,063 | `fe9525832026b03ab34f0971b43c0c81` | match | `f68bc63f9705e97ddfe2cb8c9dfad7a3e717cf5e49b050ec35c8b10a0bf45421` | 551,193 | yes |
| 2016_06 | 2016-06-08 | 1,504,963,399 | `e3a5ac5a166efc95e9ad06465d5bd2c4` | match | `ad6c96bc6d433ed99ef8722d7bf003003ada43c3ab57c005d0343910a515d484` | 551,385 | yes |

The small retained records match every hash declared in the per-release provenance JSON and ledger.

### 3.2 What the records say

The two reviewed records are stable across all three releases:

| UniProt accession | UniProt gene name | HGNC symbol | GeneID cross-references in 2016_04, 2016_05, and 2016_06 |
|---|---|---|---|
| `O95073` | FSBP | FSBP | `100861412`, `25788` |
| `Q9Y620` | RAD54B | RAD54B | `25788` |

The O95073 record also contains the historical comment, quoted exactly from the file:

```text
Intragenic, in the second intron of RAB54B gene.
```

The spelling above is preserved from the historical UniProt record. It should not be silently corrected in the provenance report.

### 3.3 Revised interpretation

The date-matched files require a correction to our earlier wording.

It is no longer defensible to describe `O95073 → 25788` as merely an erroneous historical edge. The official reviewed UniProt entries explicitly carried it in April, May, and June 2016, together with the FSBP-specific GeneID 100861412.

What the GraphSAGE labels establish is narrower:

> The label-generating transformation behaved as if the FSBP annotations attached to O95073 were not assigned to the GraphSAGE node representing GeneID 25788/RAD54B.

That behavior can be reconstructed by resolving the component according to the concordant primary symbols:

```text
Q9Y620 / RAD54B -> GeneID 25788 / RAD54B
O95073 / FSBP   -> GeneID 100861412 / FSBP
```

But the exact fit alone does not reveal which upstream mechanism the original authors used. Plausible mechanisms include:

1. an Entrez-native annotation source such as `gene2go`, where RAD54B and FSBP are separate genes;
2. a symbol-aware UniProt-to-GeneID disambiguation;
3. another curated mapping table that distinguishes the nested gene product from the host-gene cross-reference.

This increases, rather than decreases, the value of testing the May-2016 `gene2go` derivative and `org.Hs.eg.db` 3.3.0.

## 4. Exact reconstruction of the 50 input features

### 4.1 Exact rule found

The following simple procedure reproduces every feature value for every resolved GraphSAGE row:

1. Read MSigDB Entrez GMT collections in the order C1, C3, C7.
2. Within each GMT, preserve source-file row order.
3. Retain sets whose full MSigDB membership contains at least 200 unique Entrez GeneIDs.
4. Append retained sets until 50 features have been selected.
5. Encode each graph gene as binary membership in the selected sets.

The observed result is:

```text
30 C1 sets + 20 C3 sets + 0 C7 sets = 50 features
```

The 50th selected set is `V$NFKAPPAB65_01`; the next qualifying C3 set would be `V$CREL_01`. The cap is reached before any C7 set is considered.

### 4.2 Exact validation

| Quantity | Result |
|---|---:|
| GraphSAGE feature matrix | 56,944 × 50, `float64` binary values |
| Resolved node rows compared | 56,411 |
| Distinct resolved Entrez genes | 4,268 |
| Resolved feature cells compared | 2,820,550 |
| Differences | **0** |
| Exact columns | **50 / 50** |
| Conflicting feature vectors among repeated appearances of the same gene | **0** |

A second implementation using independent nested Python loops also checks all 2,820,550 resolved cells and finds zero differences.

### 4.3 The formerly unidentified all-zero feature

Feature column 10 is:

```text
chryq11
```

It contains 204 unique Entrez GeneIDs in the full MSigDB collection, which is why it passes the global size filter. None of those GeneIDs occurs among the 4,268 resolved GraphSAGE genes, so the deposited feature column is entirely zero.

### 4.4 C7 interpretation

The feature matrix itself establishes only that the selected sets are 30 C1 sets followed by 20 C3 sets. It cannot prove that the original code explicitly included C7 in its collection loop.

However, the exact threshold/cap behavior offers a coherent explanation for the project documentation stating that C1, C3, and C7 supplied features while the deposited matrix contains no C7 column: a global 50-feature cap is reached during C3.

This is a plausible reconciliation, not source-code proof.

### 4.5 Threshold identifiability

The selected sequence is unchanged for either of these conditions:

```text
set size >= 200
set size >= 201   # equivalently, set size > 200
```

The smallest selected set has 201 genes, and the largest excluded C1 set has 199. Therefore, the data identify a cutoff boundary around 200 but do not determine the exact comparison operator.

### 4.6 MSigDB version is not identifiable from the feature matrix

The same threshold/source-order rule gives the same 50 membership vectors in every supplied version:

| MSigDB version | Input used in this validation | Exact columns | Differences | Selected membership-sequence SHA-256 |
|---|---|---:|---:|---|
| 5.0 | hash-verified retained normalized derivative | 50 | 0 | `41c5d821e1b706ec4c8dceb47ab25c5dbad689998483e34e9e41d08094448101` |
| 5.1 | prior user-supplied raw archive | 50 | 0 | same |
| 5.2 | prior user-supplied raw archive | 50 | 0 | same |
| 6.0 | prior user-supplied raw archive | 50 | 0 | same |

MSigDB v6.0 changes some displayed set names, such as removing the `V$` prefix, but not the selected membership vectors.

The feature provenance conclusion must therefore be revised from “uniquely MSigDB v5.2” to:

> The 50 features are exactly recoverable from the supplied MSigDB versions 5.0 through 6.0 using the same simple source-order and size rule. The feature matrix alone does not distinguish which of these versions was used.

## 5. Column-order status after B104E

The date screen closes one major line of inquiry:

- releases 158–166 all yield the same best LCS of 94;
- releases 167–169 are slightly worse at 93;
- no tested date produces a perfect order;
- the duplicate orientation remains stable.

The strongest supported statement remains:

> The label order has a strong 64-bit, unrandomized CPython-2 dictionary-iteration fingerprint, but the exact dictionary construction has not been recovered.

A perfect fit can always be manufactured by choosing an artificial hash-table insertion history. Such an overfit would not be evidence. A scientifically useful perfect match must arise from a plausible source file and preprocessing sequence.

The highest-value remaining order hypotheses are:

1. an Entrez-native `gene2go` scan and its original row order;
2. an `org.Hs.eg.db` 3.3.0 map/table iteration order;
3. the exact late-June ontology key universe from B105;
4. separate direct-annotation and propagated-annotation dictionaries;
5. original preprocessing code or an ordered intermediate term list.

## 6. Claim-strength table

| Claim | Strength after B104E | Basis |
|---|---|---|
| Release 159 is the unique exact GOA release among 158–169 under the fixed transformation | **Very strong within the tested model** | 159 has zero differences; every other release has at least 652. Results contain source hashes and passed integrity/cleanup events. |
| The original authors necessarily read `goa_human.gaf.159.gz` | **Not established** | An Entrez-native or processed source could encode the same associations. |
| The ≥1,000 full-human prevalence rule was used | **Strongly supported** | Exactly 121 terms only in release 159, with a 1,007/997 boundary; not source-code proven. |
| O95073→25788 existed in official 2016 UniProt | **Established** | Present in reviewed records from 2016_04, 2016_05, and 2016_06. |
| O95073→25788 was simply a bad historical cross-reference | **Withdrawn / unsupported** | Date-matched UniProt explicitly retained it across releases. |
| FSBP annotations were not projected to RAD54B/GeneID 25788 in the label generator | **Established by matrix behavior** | Projecting them creates 13 false-positive cells; excluding them gives exact membership. |
| The original disambiguation method was symbol-aware | **Plausible, not proven** | It exactly reproduces the behavior, but Entrez-native annotation is an alternative. |
| The exact 50-feature selection rule was recovered | **Exact behavioral reconstruction; original code not proven** | Two implementations reproduce 2,820,550 resolved cells with zero differences. |
| MSigDB v5.2 was uniquely the feature source | **Withdrawn** | Versions 5.0, 5.1, 5.2, and 6.0 produce identical selected memberships. |
| The label columns use Python-2 dictionary iteration | **Strongly supported** | Independent 56,944-key positive control and high GO-order concordance; no perfect GO construction yet. |
| Duplicate-pair orientation `001` is correct | **Strongly supported, provisional** | Stable across all high-scoring model families and releases 158–169; no original ordered list. |

## 7. Recommended next analyses

The most informative next test is no longer another GOA release. The practical priorities are:

1. **Entrez-native annotation test.** Use the May-2016 `gene2go` derivative or a close snapshot to determine whether the exact label matrix and O95073/RAD54B separation arise without any UniProt mapping.
2. **Bioconductor 3.3.0 mapping audit.** Extract the Entrez–UniProt and Entrez–GO tables from `org.Hs.eg.db_3.3.0.tar.gz`, focusing on GeneIDs 25788 and 100861412 and on the full label reconstruction.
3. **Column-order screen from those Entrez-native table orders.** Test whether dictionary insertion while scanning `gene2go` or Bioconductor tables improves LCS 94.
4. **B105 ontology robustness.** Test the exact 2016-06-29 or nearest July-1 ontology while holding release 159 annotations fixed.
5. **Preprocessing-code search.** Search for constants or patterns corresponding to feature threshold 200/201, feature cap 50, label threshold 1,000, and direct iteration over a GO dictionary.

## 8. Provenance and deletion state

The user-side scripts deleted the multi-gigabyte UniProt parent archives and the downloaded GOA GAF/GPI pairs only after successful validation and retained-output hashing. Those deletion events are preserved.

The two B104E uploaded ZIPs are small wrappers around complete retained audit outputs. After this report, checksum inventory, and frozen bundle are validated, the conversation copies can be deleted without losing the extracted records or result JSONs.

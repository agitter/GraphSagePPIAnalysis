# B104C — MSigDB v5.0 and GraphSAGE GO-label column-order investigation

Generated: 2026-08-28T20:09:48Z  
Batch: B104C  
B105 status: deferred by user; no late-June ontology was used in this batch.

## Executive conclusions

1. **MSigDB v5.0 does not directly reproduce the GraphSAGE label memberships.** Across 1,454 C5 ontology sets and 10,348 total v5.0 gene sets, no set exactly matches any of the 121 GraphSAGE label columns on the 4,268 independently resolved Entrez genes. The closest v5.0 set still differs at 307 genes.
2. **The GraphSAGE column order has a strong legacy-Python fingerprint.** The complete JSON key order of `ppi-class_map.json` is reproduced exactly by a 64-bit, unrandomized CPython 2.7-style dictionary containing string keys `"0"` through `"56943"` inserted in ascending order. This gives an independent positive control that legacy Python 2 dictionary iteration was used when serializing at least one core dataset artifact.
3. **Conventional intentional GO orderings do not explain the 121 columns.** GO ID order, name order, prevalence, ontology-file order, first GAF occurrence, and MSigDB XML order all have weak or negligible agreement with the deposited order.
4. **A large CPython 2 GO-term dictionary explains most of the column order.** Simulating a dictionary built while scanning the accepted release-159 GOA rows and inserting direct terms plus `is_a` ancestors produces Kendall tau 0.7722, pairwise concordance 88.61%, and a longest common subsequence of 92 of 121 terms. A reduced grid of 48 plausible GAF-derived simulations reaches a longest common subsequence of 94 and is unanimous about the three otherwise ambiguous duplicate-vector pairs.
5. **The six duplicate-vector columns can now be assigned provisionally.** The order evidence strongly supports:

   - column 24 = `GO:0043228`; column 71 = `GO:0043232`
   - column 39 = `GO:0006464`; column 63 = `GO:0036211`
   - column 48 = `GO:1903561`; column 70 = `GO:0043230`

   These assignments are strongly supported by independent hash-table simulations, but remain provisional until original preprocessing code, an exact ordered source list, or an exact recreation of the complete dictionary construction is found.
6. **The likely ordering operation was direct filtering of a large dictionary.** The data are consistent with a Python 2 operation structurally similar to:

   ```python
   labels = [go_id for go_id in go_to_genes if len(go_to_genes[go_id]) >= 1000]
   ```

   Copying the selected terms into a new 121-key dictionary and then iterating it destroys the observed ordering signal.

## Inputs and integrity

### New raw input

| File | Bytes | SHA-256 | Verification |
|---|---:|---|---|
| `msigdb_v5.0_files_to_download_locally.zip` | 115,484,475 | `d372fc23f229cbb79656d824e0519587db6110963d22d1f4c95e5154963a32d2` | Matches the user's full local inventory; ZIP and XML parsed successfully |

The archive contains `msigdb_v5.0.xml`, whose internal build date is **Apr 27, 2015**. The public MSigDB release-notes page labels v5.0 as **Mar 2015**. These are recorded as different metadata fields rather than forced into one date.

### Reused verified inputs and retained derivatives

- GraphSAGE `ppi.zip`
- exact B104A GO-label reconstruction and 121-column candidate map
- normalized GOA human release-159 GAF
- normalized June 1, 2016 GO term and `is_a` tables
- full-human release-159 propagated GO prevalence counts
- previously supplied MSigDB v5.1, v5.2, and v6.0 archives for cross-version comparison

The user reported deleting prior MSigDB archives from local storage to free space. Because the exact deleted filenames were not enumerated, the provenance ledger records a collection-level user-local deletion report rather than silently changing each file's local status. Previously supplied conversation/runtime copies were still available for this analysis, and compact cross-version derivatives have now been retained.

## MSigDB v5.0 direct-membership test

The direct comparison used every v5.0 Entrez membership vector, restricted only to the same 4,268 resolved GraphSAGE genes used in the exact GOA reconstruction.

| Scope | Sets tested | Exact GraphSAGE columns | Closest mismatch | Median best mismatch | Columns at least 99% | Columns at least 95% |
|---|---:|---:|---:|---:|---:|---:|
| C5 ontology sets | 1,454 | 0 | 307 | 823 | 0 | 0 |
| All v5.0 collections | 10,348 | 0 | 307 | 823 | 0 | 0 |

The closest v5.0 C5 result is `RNA_BINDING` compared with GraphSAGE column 86, but it still disagrees for 307 of 4,268 genes.

Only 57 of the 121 exact GO candidate IDs are represented in v5.0 C5, and none of those 57 direct MSigDB memberships equals the exact GOA-derived membership.

### Cross-version C5 identity coverage

A corrected parser now counts GO IDs only when the corresponding XML row is actually in C5:

| MSigDB version | Internal build date | Total sets | C5 sets | Exact candidate GO IDs represented in C5 |
|---|---|---:|---:|---:|
| 5.0 | Apr 27, 2015 | 10,348 | 1,454 | 57 / 121 |
| 5.1 | Jan 27, 2016 | 13,311 | 1,454 | 57 / 121 |
| 5.2 | Sep 22, 2016 | 18,890 | 6,166 | 6 / 121 |
| 6.0 | Feb 22, 2017 | 18,643 | 5,917 | 6 / 121 |

The earlier B104B count of 23 target IDs in v5.2/v6.0 was too broad: it counted GO IDs found elsewhere in the XML, not only C5. This batch supersedes that count with the C5-scoped value of 6.

MSigDB v5.0 and v5.1 have 1,452 C5 standard names in common, and all 1,452 common sets have identical full Entrez membership. Each release has two C5 names absent from the other. Their complete C5 sequence is not byte-for-byte identical.

### Interpretation

The v5.0 result strengthens the previous conclusion:

- no tested MSigDB release from 5.0 through 6.0 directly supplies the deposited label memberships;
- v5.0/v5.1 cannot even supply 64 of the 121 GO identities through C5;
- v5.2/v6.0 C5 include only six of the 121 exact identities under their revised C5 organization;
- the exact memberships instead come from the historical GOA transformation already reconstructed.

MSigDB may still have influenced terminology, documentation, or another preprocessing stage, but the tested archives are not the direct label matrix source.

## Column-order investigation

### Positive control from `ppi-class_map.json`

The GraphSAGE class-map JSON contains 56,944 object keys. The observed first keys are not numeric order; they begin with values such as `50088`, `44884`, `11542`, and `11543`.

A simulator implementing the CPython 2.7 string hash, collision probing, resize rule, and low-to-high table scan was tested as follows:

1. Create an empty legacy dictionary.
2. Insert string keys `"0"`, `"1"`, ..., `"56943"` in ascending order.
3. Iterate the final internal table from low to high slot.

On a 64-bit model with hash randomization disabled, the simulated sequence matches **all 56,944 JSON keys exactly**. A 32-bit model matches only 5,410 positions and first diverges at position 24.

This does not by itself prove how the GO-term dictionary was built, but it establishes that the dataset contains a direct legacy-Python dictionary-order artifact and identifies the relevant hash width and default zero hash secret.

### Conventional order models

The deposited GO-column sequence is poorly explained by intentional biological or lexical orderings:

| Candidate order | Comparable terms | Kendall tau | Longest common subsequence |
|---|---:|---:|---:|
| GO ID ascending | 121 | -0.1008 | 18 |
| GO name alphabetical | 121 | 0.0121 | 20 |
| GraphSAGE positive count descending | 121 | -0.0802 | 14 |
| Full-human GOA prevalence descending | 121 | -0.0347 | 18 |
| June 1 OBO stanza order | 121 | -0.1008 | 18 |
| First accepted GAF occurrence | 121 | -0.0540 | 19 |
| MSigDB v5.0 XML order | 57 | -0.0501 | 12 |
| MSigDB v5.1 XML order | 57 | -0.0501 | 12 |

The MSigDB v5.2/v6.0 comparison has only six comparable C5 terms and therefore provides little ordering information.

### Legacy-hash models

| Model | Terms | Final table size | Kendall tau | Pairwise concordance | LCS | Exact positions | Exact prefix |
|---|---:|---:|---:|---:|---:|---:|---:|
| Raw 64-bit Python 2 GO-ID hash, low 15 bits | 121 | 32,768 conceptual slots | 0.7138 | 85.69% | 86 | 16 | 5 |
| Accepted GAF rows; insert direct terms and `is_a` ancestors in a large dictionary | 16,357 dictionary keys | 32,768 | 0.7722 | 88.61% | 92 | 6 | 5 |
| Best of 48 plausible GAF-derived variants | 13,302 dictionary keys | 32,768 | 0.7661 | 88.31% | 94 | 13 | 5 |
| Copy selected terms into a new 121-key dictionary | 121 | 512 | 0.0477 | 52.38% | 23 | 0 | 0 |
| Dictionary containing all 44,797 ontology terms | 44,797 | 131,072 | 0.2361 | 61.80% | 38 | 1 | 0 |

The full accepted-GAF model has two-sided Kendall p = `3.54e-36`. Even an intentionally overconservative Bonferroni adjustment for 10,000 hypothetical order models leaves p = `3.54e-32`.

The **32,768-slot table-size fingerprint** is important. It is the size expected for a growing CPython 2 dictionary with roughly 13,000–16,000 GO keys under the documented two-thirds fill and quadrupling resize policy. The strongest order signal occurs at 15 low hash bits, exactly corresponding to 32,768 slots.

### Best current construction hypothesis

The evidence is most consistent with this family of operations:

1. Build one large dictionary mapping GO IDs to human genes/proteins while processing annotations and ancestors.
2. Let CPython 2 determine dictionary slot placement.
3. Iterate that large dictionary directly.
4. Retain terms meeting the approximately 1,000-gene/protein criterion.
5. Append membership vectors in that iteration order.

This explains why:

- the selected **set** is exactly the full-human top 121 / at least 1,000 terms;
- the selected **order** strongly follows a 32,768-slot Python 2 dictionary;
- recreating a smaller dictionary from the selected 121 terms does not reproduce the order;
- the first five columns are reproduced exactly by the large-dictionary simulations.

The exact sequence is not yet completely reconstructed. Likely remaining variables include the exact late-June ontology key universe, the exact order in which direct terms and ancestors were first inserted, annotation preprocessing order, and the original data structure used by the authors.

## Disambiguating the three duplicate-vector pairs

The exact membership matrix alone cannot distinguish the two GO IDs in each pair because both terms have identical membership on the 4,268 resolved genes. Order evidence supplies an additional independent signal.

All 48 plausible GAF-derived Python 2 dictionary simulations selected the same orientation. In each pair, the GO term assigned to the earlier GraphSAGE column occupies the earlier dictionary slot.

| Earlier column | Earlier GO term | Later column | Later GO term | Order status |
|---:|---|---:|---|---|
| 24 | `GO:0043228` — non-membrane-bounded organelle | 71 | `GO:0043232` — intracellular non-membrane-bounded organelle | Strongly supported |
| 39 | `GO:0006464` — cellular protein modification process | 63 | `GO:0036211` — protein modification process | Strongly supported |
| 48 | `GO:1903561` — extracellular vesicle | 70 | `GO:0043230` — extracellular organelle | Strongly supported |

The relevant low-15-bit / final-table slots are:

- `GO:0043228`: 5,267; `GO:0043232`: 21,094
- `GO:0006464`: 10,702; `GO:0036211`: 20,179
- `GO:1903561`: 14,299; `GO:0043230`: 21,092

These assignments are incorporated into `B104C_inferred_unique_121_GO_column_order_20260828T194921Z.csv`. Rows with unique membership-based identities remain classified as exact; only these six rows are marked as order-inferred.

## What would resolve the remaining ordering uncertainty?

Highest-information tests are:

1. **B105, when storage permits:** repeat only the dictionary-construction/order simulation with the exact or nearest late-June 2016 ontology. Membership should remain fixed; changes in key universe and ancestor insertion order may refine the exact column sequence.
2. **Original preprocessing code or an intermediate GO-to-gene dictionary:** search for a Python 2 script that filters by `1000`, `>= 1000`, `121`, or iterates a GO dictionary without sorting.
3. **Original ordered term-list artifact:** a text, pickle, JSON, NumPy, or MATLAB file containing 121 GO terms would settle the six ambiguous assignments and the complete sequence.
4. **Exact annotation-row ordering and mapping universe:** reproduce whether the dictionary was populated before or after taxon filtering, identifier projection, qualifier/evidence filtering, and ontology propagation.
5. **Python runtime confirmation:** package metadata, logs, or environment files that establish Python 2.7, 64-bit architecture, and disabled hash randomization would independently confirm the serialization model already inferred from `ppi-class_map.json`.

## Sequential UniProt mapping audit

A single script is provided:

```text
download_extract_uniprot_2016_mapping_audit.py
```

It handles releases 2016_04, 2016_05, and 2016_06 one at a time. For each release it:

- downloads the official `RELEASE.metalink`;
- downloads/resumes the reviewed-only Swiss-Prot archive;
- checks the release identifier, exact byte size, and official MD5;
- records a local SHA-256;
- streams the complete Swiss-Prot flat file without unpacking it to disk;
- retains full records containing `O95073` or `Q9Y620`;
- writes a compact TSV, complete `.dat` records, JSON provenance, and an append-only CSV ledger;
- hashes all retained outputs; and
- deletes the approximately 1.5 GB archive only after every validation succeeds.

Any failure leaves the archive or partial download in place. A synthetic self-test scanned 100,004 records, extracted both targets, and exercised the post-validation deletion gate successfully. A dry run confirmed all three official URLs, sizes, and MD5 values.

Recommended command:

```bash
python download_extract_uniprot_2016_mapping_audit.py \
  --work-dir uniprot_audit_work \
  --output-dir uniprot_audit_results
```

Only the small `uniprot_audit_results` directory should later be uploaded. The large archives do not need to be uploaded or retained.

Official releases pinned by the script:

| Release | Release date | Archive bytes | MD5 |
|---|---:|---:|---|
| 2016_04 | 2016-04-13 | 1,516,525,310 | `e607b83de1ac87e6f63b13715c049a3f` |
| 2016_05 | 2016-05-11 | 1,504,161,063 | `fe9525832026b03ab34f0971b43c0c81` |
| 2016_06 | 2016-06-08 | 1,504,963,399 | `e3a5ac5a166efc95e9ad06465d5bd2c4` |

## Retained MSigDB v5.0 derivative

The raw 115 MB archive has been normalized into a compact, complete Entrez-membership table:

```text
B104C_msigdb_v5.0_normalized_entrez_gene_sets_20260828T194921Z.tsv.gz
```

Reconciliation:

- 10,348 gene-set rows
- 1,454 C5 rows
- 1,315,074 total unique membership entries summed across rows
- SHA-256 `6775fb96be44080c768bd5789a0dbb0c802a1a0faa45927aa2a07d70af9f7c1f`

The normalized file preserves all gene-set names, systematic IDs, collection/subcollection, GO ID when available, chip namespace, external details URL, source XML order, C5 order, and sorted unique Entrez memberships.

## Current scientific interpretation

- The GOA-based label **membership** reconstruction remains exact.
- The full-human **term-selection set** remains exactly top 121 / at least 1,000 genes or proteins.
- Direct MSigDB memberships from versions 5.0, 5.1, 5.2, and 6.0 do not reproduce the labels.
- The **column order** is very likely an unintentional CPython 2 dictionary iteration artifact.
- The order evidence strongly disambiguates all six duplicate-vector columns.
- Exact ordering reproduction remains open and should be treated separately from the already exact label membership reconstruction.

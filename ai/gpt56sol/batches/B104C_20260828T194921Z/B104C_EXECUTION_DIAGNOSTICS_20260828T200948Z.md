# B104C execution diagnostics

Generated: 2026-08-28T20:09:48Z

## Accepted executions

### Main MSigDB v5.0 and column-order analysis

Script:

```text
scripts/analyze_B104C_msigdb50_column_order.py
```

Status: completed successfully.  
Standard error: empty.  
Primary outputs:

- v5.0 direct-membership grid and summary
- normalized v5.0 gene-set metadata
- GraphSAGE class-map CPython 2 dictionary validation
- conventional and hash-based column-order comparisons
- duplicate-vector disambiguation
- cross-version MSigDB C5 identity metadata

Input checks performed by the script:

- v5.0 archive byte size matched 115,484,475
- v5.0 archive SHA-256 matched `d372fc23f229cbb79656d824e0519587db6110963d22d1f4c95e5154963a32d2`
- GraphSAGE archive and all retained B104A/B104B derivatives were hashed
- v5.0 ZIP and XML parsed successfully

### Extended GAF-derived Python 2 dictionary simulations

Script:

```text
scripts/explore_B104C_column_order_variants.py
```

Status: completed successfully.  
Standard error: empty.  
Grid size: 48 accepted simulations.

The reduced grid varied:

- annotation-row inclusion policy;
- source-row ordering;
- ancestor insertion order; and
- dictionary-key population details.

All 48 simulations ended with a 32,768-slot dictionary and selected the same orientation for all three duplicate-vector pairs.

### MSigDB v5.0 normalization

Script:

```text
scripts/normalize_msigdb50.py
```

Status: completed successfully.  
Standard error: empty.  
Checks:

- normalized rows: 10,348
- C5 rows: 1,454
- summed unique memberships: 1,315,074
- gzip integrity: passed
- normalized SHA-256: `6775fb96be44080c768bd5789a0dbb0c802a1a0faa45927aa2a07d70af9f7c1f`

### UniProt sequential audit script tests

Script:

```text
scripts/download_extract_uniprot_2016_mapping_audit.py
```

Compilation: passed.

Synthetic self-test:

- 100,004 Swiss-Prot-format records scanned;
- O95073 and Q9Y620 records extracted;
- complete `.dat`, compact TSV, and provenance outputs hashed;
- large synthetic archive deleted only after successful validation;
- standard error empty.

Dry run:

- all three official archive URLs printed;
- official sizes and MD5 values matched the pinned values;
- no network download initiated;
- standard error empty.

The actual approximately 1.5 GB UniProt archives were not downloaded in this runtime.

## Non-accepted or incomplete attempts

### Initial exhaustive order grid

Exploratory script:

```text
/mnt/data/ppi_repro_corrected/work_B104C/order_grid.py
```

The initial Cartesian search attempted substantially more combinations than were needed for the scientific question. It exceeded the practical execution window after approximately three minutes and produced no accepted output. Its stdout and stderr files were empty. Results from that attempt were not used.

It was replaced by the bounded 48-configuration simulation grid described above. The reduced grid preserved the major plausible axes of variation and completed successfully.

### Direct Git clone of the GraphSAGE repository

A container-side `git clone` attempt failed because the execution environment could not resolve `github.com`. No cloned source or result from that attempt was used. Public repository pages and official CPython source files were inspected with the web-access tool instead.

No publicly indexed original PPI preprocessing script was located during the web search. Therefore, the legacy-dictionary explanation remains an inference supported by the serialized data and CPython behavior rather than a quotation from the original preprocessing source.

## Corrected prior result

B104B reported 23 exact candidate GO IDs represented in C5 for MSigDB v5.2 and v6.0. That count was based on GO IDs detected anywhere in the XML rows. B104C reparsed all versions with an explicit `CATEGORY_CODE == C5` requirement and obtained:

```text
v5.0: 57 / 121
v5.1: 57 / 121
v5.2:  6 / 121
v6.0:  6 / 121
```

The B104C C5-scoped result supersedes the earlier count.

## Storage and availability notes

The user reported deleting previously held MSigDB archives locally to free disk space. Because that statement did not enumerate exact filenames, it is recorded as a collection-level local deletion event. It does not imply that every corresponding conversation attachment was deleted.

Previously supplied v5.1, v5.2, and v6.0 raw archives were still present in the analysis runtime and were used only for the documented cross-version comparisons. Compact metadata and summaries have been retained so that those conclusions do not require immediate re-upload.

B105 was explicitly deferred by the user and was not executed.

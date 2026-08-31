# B104D — Column-order limits, source-date plan, and UniProt audit review

## Executive conclusions

1. **No tested model reproduces the complete 121-column order.** B104D added 492 model configurations to the 48 previously tested in B104C. Across the combined 540 scored configurations, no model reached a longest common subsequence (LCS) of 100, let alone 121.
2. **The previous best model remains the best tested model.** It has LCS 94/121, Kendall tau 0.7661, 88.31% pairwise concordance, 13 columns in their exact absolute positions, and an exact five-column prefix.
3. **LCS 94 does not mean 94 absolute column positions are correct.** It means 94 terms can be retained in the same relative order after deleting 27 terms from both sequences. Only 13 absolute positions match.
4. **The duplicate-vector disambiguation remains strongly supported but not proven.** Every one of the 161 models with LCS at least 80 chose the same orientation (`001`) for the three duplicated vector pairs. This includes all 48 original GAF-derived models, every high-scoring targeted model, and all 36 outer/pair-inversion models.
5. **The most important untested axis is source date.** The GraphSAGE ZIP records an extended UTC timestamp of 2017-05-10 19:25:23 for `ppi-class_map.json`. GOA releases 160 through 168 predate that artifact, but only releases 158 and 159 are currently materialized here. A sequential low-storage screening script was produced for releases 158–169.
6. **The uploaded UniProt package is incomplete but useful.** Its ledger says all three official 2016 Swiss-Prot archives passed expected size and MD5 checks, both O95073 and Q9Y620 were found, and the large archives were deleted. The ZIP contains only the ledger, however, so the extracted `GN` and `DR   GeneID;` lines cannot yet be inspected independently.

## 1. Claim-strength framework

| Claim | Current strength | Reason |
|---|---|---|
| At least one core GraphSAGE artifact was serialized using 64-bit, unrandomized CPython-2 dictionary order | Very strong | The complete sequence of 56,944 keys in `ppi-class_map.json` is reproduced exactly by the CPython-2 dictionary simulator. |
| The GO-label column order arose from direct iteration over a large Python-2 hash table keyed by GO IDs | Strongly supported | Large GO-ID dictionaries give LCS 94 and 88.31% pairwise concordance; names, numeric IDs, MSigDB names, small rebuilt dictionaries, and ordinary biological orders perform poorly. |
| The exact original dictionary construction has been recovered | Not established | No tested model reproduces all 121 positions, and source releases 160–168 remain untested. |
| The six duplicated columns have the `001` orientation | Strongly supported, provisional | Every high-scoring model chooses it, but no original ordered term list or preprocessing code has been found. |
| GOA release 159 was necessarily the source release | Not yet established | Release 159 exactly reproduces membership under the fixed transformation, but later releases have not yet been screened. |

## 2. Expanded column-order model search

### Model families

B104D tested three additional model families on top of the original B104C grid:

- **256 targeted key/construction variants**
  - GO ID strings;
  - GO names;
  - MSigDB standard names;
  - numeric GO identifiers and alternate string encodings;
  - OBO-primary and OBO-primary-plus-alternate key universes;
  - direct-term dictionaries, propagated dictionaries, and different row/ancestor orders.
- **200 GraphSAGE-mapped variants** that varied whether the key universe was limited to annotations reaching resolved graph genes.
- **36 outer-dictionary and annotation-pair variants**
  - accession-to-term, symbol-to-term, and GeneID-to-term dictionaries followed by inversion;
  - Python-2 dictionary or set inner containers;
  - tuple-keyed annotation-pair sets in both element orders.

These are 492 new configurations. Together with the 48 prior B104C models, 540 configurations were scored. They are not statistically independent; many share data and key universes. The count is a coverage description, not a multiple-testing claim.

### Best-scoring model family

Several tied configurations attain the same top score. One representative best-scoring configuration is:

```text
key representation:        GO ID string, e.g. "GO:0050789"
annotation filter:         exact accepted evidence/relation policy
mapping scope:             accessions with a historical GeneID edge after component-aware correction
GAF row order:             annotation date, accession, GO ID
key insertion:             first direct GO occurrence, then sorted is_a ancestors
runtime model:             64-bit unrandomized CPython-2 dictionary
final table size:          32,768 slots
unique GO keys:            13,302
```

Its scores are:

```text
LCS:                       94 / 121
Kendall tau:               0.7661157025
pairwise concordance:      0.8830578512
exact absolute positions:  13 / 121
exact prefix:              5 columns
```

### Negative results that sharpen the interpretation

The strong signal is specific to GO-ID string keys in a large table.

- GO names and MSigDB standard names reach only approximately LCS 22.
- Numeric GO-ID encodings reach approximately LCS 18–23.
- Static ontology dictionaries reach at most LCS 38.
- Copying the 121 selected terms into a fresh small dictionary destroys the signal.
- Outer accession/symbol/GeneID dictionaries followed by inversion reach at most LCS 90.
- Annotation-pair sets do not improve on the large GO-ID dictionary.

These failures make a direct GO-ID-keyed hash table more plausible than an intentionally ordered biological list, a name-keyed table, or an MSigDB list order.

## 3. Duplicate-vector columns

The label matrix contains three pairs of columns whose values are identical across the resolved GraphSAGE gene universe. Membership alone cannot identify which exact GO term belongs to the earlier versus later column.

All models with LCS at least 80 choose the same orientation:

| Earlier column | Provisional GO term | Later column | Provisional GO term |
|---:|---|---:|---|
| 24 | GO:0043228 — non-membrane-bounded organelle | 71 | GO:0043232 — intracellular non-membrane-bounded organelle |
| 39 | GO:0006464 — cellular protein modification process | 63 | GO:0036211 — protein modification process |
| 48 | GO:1903561 — extracellular vesicle | 70 | GO:0043230 — extracellular organelle |

This orientation should be described as **strongly supported** rather than proven. A perfect recreation of the source dictionary, an ordered intermediate label file, or the preprocessing code would upgrade it to proven.

## 4. Source-date range

### GraphSAGE artifact timestamps

The ZIP extended timestamps are:

| Member | Extended UTC timestamp |
|---|---|
| `ppi-walks.txt` | 2017-05-09 03:20:50 |
| `ppi-class_map.json` | 2017-05-10 19:25:23 |
| `ppi-feats.npy` | 2017-05-10 19:26:10 |
| `ppi-id_map.json` | 2017-05-10 19:35:25 |
| `ppi-G.json` | 2017-05-29 16:12:53 |

The `ppi-class_map.json` timestamp is a strong practical upper bound on label creation, but ZIP timestamps can be preserved or rewritten and are not cryptographic proof.

The GraphSAGE v1 preprint was submitted on 7 June 2017. Taken together, the most relevant GOA source window is:

```text
release 159: 2016-07-04
release 160: 2016-09-14
release 161: 2016-10-03
release 162: 2016-10-31
release 163: 2016-11-28
release 164: 2017-01-16
release 165: 2017-02-13
release 166: 2017-03-13
release 167: 2017-04-10
release 168: 2017-05-08
```

Release 158 is a useful earlier control. Release 169, dated 5 June 2017, is a post-`class_map` negative control.

### What has actually been tested

Only releases 158 and 159 are materialized in this runtime.

Under the final exact transformation and fixed June-2016 ontology:

| Release | Exact columns | FP | FN | Total mismatches | Top-121 term-set overlap | Order LCS | Duplicate orientation |
|---:|---:|---:|---:|---:|---:|---:|---|
| 158 | 8 | 21 | 825 | 846 | 121/121 | 94 | `001` |
| 159 | 121 | 0 | 0 | 0 | 121/121 | 94 | `001` |

Thus membership sharply favors 159 over 158, while these two releases do not distinguish the ordering model.

### Low-storage date screen

`screen_goa_release_date_range.py` was produced to test 158–169 sequentially. It downloads one GAF/GPI pair, verifies gzip and SHA-256, analyzes it, writes compact results, and deletes script-downloaded inputs before continuing. It never deletes files supplied through a local directory.

The script reports for each release:

- exact GraphSAGE label columns and total FP/FN;
- the ≥1,000 and top-121 term-selection sets;
- source hashes and headers;
- mapping fallback counts;
- the fixed Python-2 order-model LCS, tau, exact positions, and duplicate orientation.

This date screen holds the May-2016 mapping and 1 June 2016 ontology fixed. B105 remains a separate ontology-version robustness test.

## 5. UniProt 2016 audit package

The uploaded ZIP has SHA-256:

```text
494845b09874e26ef95c724f2a904f3d8cfeadb9c93c662975dc0028fcd0c61a
```

It contains only:

```text
uniprot_2016_mapping_audit_ledger.csv
```

The ledger internally reports:

| Release | Archive bytes | Expected/observed MD5 | Records scanned | Targets found | Archive deleted |
|---|---:|---|---:|---|---|
| 2016_04 | 1,516,525,310 | match | 550,960 | O95073, Q9Y620 | yes |
| 2016_05 | 1,504,161,063 | match | 551,193 | O95073, Q9Y620 | yes |
| 2016_06 | 1,504,963,399 | match | 551,385 | O95073, Q9Y620 | yes |

This supports that the user-side script correctly downloaded and verified the official parent archives. It does **not** yet show what the extracted records say, because the following small outputs are absent from the ZIP:

- the three extracted `.dat` files;
- the three summary `.tsv` files;
- the three per-release provenance JSON files.

No 1.5-GB archive needs to be reacquired. `package_uniprot_audit_outputs.py` verifies the hashes already recorded in the ledger and creates a compact complete package.

Until those record files are reviewed, the O95073/Q9Y620 correction remains supported by the full historical mapping component and concordant symbols, but the date-matched UniProt flat-file confirmation is still pending.

## 6. ChatGPT storage message

The Library meter and chat attachment/upload limits are separate mechanisms. OpenAI's Library documentation explicitly says Library storage is separate from daily attachment/chat limits. The upload FAQ also documents rolling upload-rate limits and notes that ChatGPT does not show how much of that rolling quota remains.

Therefore, a Library display of 72.3 MB used does not rule out:

- a rolling upload-rate limit;
- a per-chat or Project file-count limit;
- a failed-upload attempt counting toward the rate cap;
- a transient service issue.

The screenshot does not establish which limit was encountered. Any earlier statement from this analysis that equated the warning with the visible Library storage meter was too strong.

## 7. Current conclusion

The Python-2 dictionary explanation is the best-supported ordering mechanism, but the exact original construction has not been recovered. The provisional duplicate orientation is substantially more reliable than any single model score because it is invariant across all high-scoring model families tested so far.

The highest-value next evidence is now:

1. sequential GOA release screening from 160 through 168;
2. complete small UniProt extracted-record package;
3. exact or nearest late-June 2016 ontology for B105;
4. original preprocessing code or an ordered intermediate GO-term list.

# B104D execution diagnostics

## Accepted executions

| Analysis | Exit status | Key output |
|---|---:|---|
| Targeted key/construction order grid | 0 | 256 models; no perfect order; best LCS 94 |
| Graph-mapped order grid | 0 | 200 models; no perfect order; best LCS 94 |
| Outer-dictionary/pair inversion grid | 0 | 36 models; no perfect order; best LCS 90 |
| GOA date-screen offline self-test | 0 | Release 158 and 159 results reproduced prior calculations |
| UniProt repackager synthetic test | 0 | Hash/path/integrity/deletion-safe packaging checks passed |
| B104D provenance build | 0 | Actual-input manifest, source ledger, and event log generated |

## Resource observations

The targeted order grid used approximately 689 MB maximum resident memory and 26 seconds wall time. The outer/pair grid used approximately 648 MB and 23 seconds. No accepted execution reported a nonzero exit status.

## Warnings and limitations

1. The 540 scored model configurations are not statistically independent; many share source rows, key universes, and scoring functions.
2. A longest common subsequence of 94 means relative-order agreement for a 94-term subsequence, not 94 correct absolute positions. The best model has 13 correct absolute positions.
3. Releases 160–168 were not materialized in this runtime. The supplied low-storage script is required to screen the full likely date range.
4. `uniprot_2016_mapping.zip` contains only the ledger. The extracted flat-file records cannot be audited until the small output directory is repackaged.
5. The GraphSAGE ZIP timestamps are useful practical bounds, not cryptographic proof of when the labels were computed.
6. No source-code or ordered intermediate label list has been located; duplicate-vector disambiguation remains strongly supported and provisional.

## Empty stderr files

The graph-mapped analysis, date-screen self-test, and UniProt packager self-test produced empty stderr. The nonempty `fast2.stderr` and `outerpair_reduced.stderr` contain GNU `time` resource summaries, not errors; both report exit status 0.

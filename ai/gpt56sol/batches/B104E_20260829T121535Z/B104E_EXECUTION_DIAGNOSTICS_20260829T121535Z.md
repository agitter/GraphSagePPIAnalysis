# B104E execution diagnostics

## Accepted executions

1. Verified SHA-256 and ZIP integrity for both uploaded packages.
2. Extracted 16 UniProt audit members and 15 GOA screen members.
3. Ran `validate_uploaded_audit_archives.py`: 81/81 checks passed.
4. Parsed the six UniProt target records directly from the retained flat files and independently confirmed stable `GN`, `HGNC`, and `DR   GeneID;` fields across 2016_04, 2016_05, and 2016_06.
5. Rebuilt the GOA release summary and a 121-column × 12-release mismatch trajectory from the per-release JSON files.
6. Reconstructed the GraphSAGE features from MSigDB source GMT order and compared 2,820,550 resolved binary cells; zero mismatches.
7. Ran a second nested-loop feature implementation; zero mismatches.
8. Screened supplied MSigDB versions 5.0, 5.1, 5.2, and 6.0 using a provenance-safe implementation. Version 5.0 was read from the retained normalized derivative, not from the logically deleted raw archive. All four versions produce the same selected membership sequence and all 50 columns exactly.
9. Evaluated feature-threshold identifiability. Thresholds `>=200` and `>=201` produce the same selected sequence; no other integer threshold from 1 through 1,000 does.

## Superseded execution

An initial cross-version feature screen read a residual runtime mount of the raw MSigDB v5.0 archive after the user had confirmed deletion of the B104C conversation copy. Its numerical result was correct, but that execution is not treated as provenance-authoritative. It was superseded by `screen_msigdb_feature_versions_provenance_safe.py`, which uses the hash-verified complete normalized v5.0 derivative retained before deletion. Superseded files were moved under `logs/superseded` and are excluded from the final analysis bundle's primary result set.

## Download limitation

An independent rerun of GOA releases 160 and 168 was attempted with the previously validated date-screen script. The execution environment failed DNS resolution for the EBI FTP host after three attempts. No partial source file was accepted. The failure is recorded in:

```text
logs/independent_date_screen_160_168.stderr
```

This does not invalidate the user-side results, whose source hashes, integrity events, and cleanup events were preserved. It means raw releases 160–169 were not independently re-downloaded in this runtime.

## Interpretation corrections

- The date-matched UniProt records show that O95073→25788 was an official reviewed cross-reference in all three audited releases. Earlier wording that treated it as simply an erroneous edge is superseded.
- The exact feature rule is compatible with MSigDB 5.0, 5.1, 5.2, and 6.0. Earlier wording that uniquely attributed features to v5.2 is superseded.
- The GOA date screen does not produce a perfect column order. The Python-2 dictionary mechanism remains strongly supported, but the exact source construction is unresolved.

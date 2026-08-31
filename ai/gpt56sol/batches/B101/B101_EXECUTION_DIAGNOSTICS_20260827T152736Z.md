# B101 execution diagnostics

Final analysis stamp: `20260827T152736Z`  
Diagnostics finalized: `2026-08-27T15:33:19.298657+00:00`

## Superseded attempts

| Attempt stamp | Outcome | Retained location | Interpretation |
|---|---|---|---|
| `20260827T152148Z` | execution timeout during unoptimized direct-label grid | `failed_attempt_20260827T152148Z_timeout/` | Partial outputs only; superseded. |
| `20260827T152423Z` | execution timeout after first optimization | `failed_attempt_20260827T152423Z_timeout/` | Partial outputs only; superseded. |
| `20260827T152614Z` | foreground execution-tool timeout before completion | `failed_attempt_20260827T152614Z_tool_timeout/` | Partial outputs only; superseded. |
| `20260827T152736Z` | completed successfully using packed-bit vectorization | final B101 artifacts | Authoritative B101 result. |

The failed-attempt directories are deliberately retained and marked as superseded. They are not included as scientific outputs.

## Completed-run checks

- The three uploaded files matched B000 sizes and SHA-256 hashes.
- All three passed `gzip -t`.
- Every raw data row was preserved in normalized derivatives.
- Re-reading each derivative reproduced its raw data-row SHA-256 exactly.
- GAF and GPAD projected assertion sets were identical.
- The complete direct grid was produced for four mapping strategies, nine evidence filters, five term scopes, and two comparison masks.
- A separate brute-force implementation recomputed the selected `fallback_synonym_union / all_except_IEA_ND_NAS / all_bp / full-universe` configuration. It agreed for all 121 label columns.
- The final ZIP passed `unzip -t` when rebuilt after report finalization.

## Warnings and limitations

- `openpyxl` emitted: `Unknown extension is not supported and will be removed` while reading a Greene workbook. The workbook remained readable and GO identifiers were extracted.
- Exact official remote bytes were not independently downloaded in this runtime. Official EBI URLs are recorded, while upload identity is tied to B000 by local size and SHA-256.
- B101 contains no ontology, so true-path propagation was not tested.
- B101 GPI has no populated `DB_Xrefs`; the Entrez mapping used for this direct screen is explicitly provisional.

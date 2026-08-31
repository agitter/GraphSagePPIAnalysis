# B104A execution diagnostics

Generated: 2026-08-28T15:20:26Z

## Accepted executions

### `finish_B104A_graph_summary.py`

- Exit status: 0
- Purpose: construct per-gene tissue-occurrence, graph-degree, feature-count, and label-count summaries after correcting the row-indexing bug in the continuation script.
- Accepted outputs:
  - `derived/B104A_all_resolved_gene_topology_feature_label_properties_20260828T145842Z.csv`
  - `analysis/B104A_special_gene_topology_feature_label_properties_20260828T145842Z.csv`
  - `analysis/B104A_special_gene_group_comparisons_20260828T145842Z.csv`
  - `B104A_gene_property_summary_20260828T145842Z.json`

### `finalize_B104A_exact_reconstruction.py`

- Exit status: 0
- Purpose: reconstruct the full historical O95073/Q9Y620 mapping component, apply the qualifier and mapping policies factorially, compare GOA releases 158 and 159, emit the reconstructed matrix, and perform cell-level/hash validation.
- Result:
  - release 159: 121 exact columns, 0 FP, 0 FN;
  - release 158 under the same final policy: 8 exact columns, 21 FP, 825 FN;
  - 516,428 binary cells compared, 0 differences;
  - deposited and reconstructed row-major uint8 matrix hashes identical.
- Standard error: empty.

### `validate_B104A_exact_reconstruction_set_based.py`

- Exit status: 0
- Purpose: independent validation using ordinary Python sets and a separate parser/predictor implementation; no pandas, numpy, or integer-bitset predictor.
- Result: 121 exact columns, 0 FP, 0 FN, identical matrix hash.
- Standard error: empty.

## Superseded or incomplete executions

### `analyze_B104A_residuals.py`

- Status: timed out at the execution limit after writing the completed relation-policy, ontology-depth, evidence-mask, and several residual-analysis outputs.
- Stdout/stderr files are empty because the process was terminated by the execution timeout before buffered output was flushed.
- Acceptance policy: individual outputs were accepted only after the continuation/finalization scripts successfully re-read them and their scientific baselines were independently reproduced. The timed-out process itself is not treated as a successful complete run.

### `continue_B104A.py`

- Status: partially completed, then exited with an exception.
- Completed before failure:
  - remaining-13 pair/witness analysis;
  - date/source/distance hypothesis screens;
  - immediate-edge impact and subset analyses.
- Failure:

```text
ValueError: Length of values (56944) does not match length of index (56411)
```

- Cause: the row-to-Entrez table contains only 56,411 resolved rows, whereas graph arrays contain 56,944 rows. The script attempted positional assignment rather than indexing by the explicit `graphsage_row` field.
- Correction: `finish_B104A_graph_summary.py` indexed graph arrays by `graphsage_row`, then completed and validated the gene-property outputs.
- No output written after the failing statement was accepted from this script.

## Input discipline

B104A used no new user-uploaded raw files. It consumed:

- the frozen B104 analysis bundle;
- the frozen corrected core-analysis bundle;
- retained normalized B101/B102/B103/B104 derivatives;
- the existing MSigDB v5.2 and GraphSAGE archives, whose hashes were already recorded in the input manifest.

The raw B104 conversation attachments remained logically deleted after the user's `Deleted B104` confirmation and were not used as B104A inputs.

# B104G execution diagnostics

## Accepted executions

### Primary analysis

Script:

```text
scripts/analyze_B104G_full_row_mapping_and_dhimmel.py
SHA-256: 39e1b09856d8d8fd43265093b6a292b1daa4aeda03680b70a1458f991f38bc51
```

Result:

- completed with exit status 0;
- `logs/analyze_B104G.stderr` is empty;
- summary written to `B104G_ANALYSIS_SUMMARY_20260829T150633Z.json`;
- all expected detailed CSV/JSON and retained-input files were created.

### Independent row-order implementation

```text
scripts/independent_validate_B104G_row_order.py
SHA-256: 4443a10e6cff1bd573526d92b607950b8d0320a62223959b8043208386a4adf1
```

Result:

- completed with exit status 0;
- stderr empty;
- 56,944/56,944 independently reconstructed assignments agree with the primary map;
- 24/24 node sets and 24/24 edge sets agree;
- no repeated-gene feature or label conflicts.

### Independent full-feature implementation

```text
scripts/independent_validate_B104G_full_features.py
SHA-256: a7d25ffb9b15627551056fdc09f42f3fff10bc36c9b9089f2631e81c8ba63c88
```

Result:

- completed with exit status 0;
- stderr empty;
- all 2,847,200 feature cells agree;
- observed and expected matrix SHA-256 both equal `274ecfee66596b7e9dfb19b71a7fb39a2611a3c140bc199c25062c6ca75bfca1`.

### Independent full-label implementation

```text
scripts/independent_validate_B104G_full_labels.py
SHA-256: a73f491f6f759d6f85c5361cd108e487205be3c7ab010dde6fdf478512b8dda5
```

Result:

- completed with exit status 0;
- stderr empty;
- all 6,890,224 label cells agree;
- observed and expected matrix SHA-256 both equal `677cc50459190ba22afc0762356e573ca56e422c80fd36204669a81694afa78d`.

### Prototype row-order screen

```text
scripts/prototype_test_node_order.py
SHA-256: dd3362faa2c3e341b8d1f4977a58d2b7b0ecb4392e019d47a252c607e6a6ad33
```

This was the exploratory implementation that identified the exact Python-2 string-dictionary model. The accepted result was reimplemented in both the primary analysis and the independent row-order validator.

## Input validation

The four uploaded annotation files were validated as actual TSV data rather than HTML pages. Each has the expected header and was retained as a deterministic gzip copy. Decompressing each retained file reproduces the raw upload byte-for-byte.

The reconciliation is recorded in:

```text
analysis/B104G_dhimmel_raw_to_retained_reconciliation.json
```

## Source-code compilation

All five B104G Python scripts passed `python -m py_compile` after analysis.

## Warnings and non-scientific issues

A display-only shell attempt used the optional `column` utility, which was unavailable in the runtime. No analysis depended on that command; the underlying CSV files were inspected directly instead. This did not change any result or output.

No accepted scientific execution reported stderr output, timeout, parser warning, malformed annotation row, failed integrity check, or incomplete output.

## Scope limitations

- The inferred row-order mechanism is not backed by the original preprocessing source code.
- The three duplicate GO-vector orientations remain provisional.
- The exact label-column insertion history remains unrecovered.
- The exact 24-network selection and 20/2/2 split procedure remain open.
- The four dhimmel controls are derived summary products; their nonmatch does not exclude every possible Entrez-native annotation pipeline.

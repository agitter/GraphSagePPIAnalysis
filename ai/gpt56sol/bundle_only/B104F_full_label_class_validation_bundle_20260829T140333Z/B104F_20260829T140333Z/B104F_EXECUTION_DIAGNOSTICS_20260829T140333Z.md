# B104F execution diagnostics

Generated: 2026-08-29T14:03:33Z

## Accepted analyses

### Full class-multiset validation

Command:

```bash
python scripts/analyze_full_label_class_multisets.py \
  --graphsage-zip /mnt/data/graphsage_ppi.zip \
  --inputs retained_inputs \
  --analysis-dir analysis \
  --summary-json B104F_FULL_LABEL_CLASS_VALIDATION_20260829T140333Z.json
```

Status: passed.

Key assertions:

- 56,411 individually resolved rows; zero row mismatches.
- 183 unresolved topology/feature classes; all 183 observed/predicted label-vector multisets identical.
- 533 unresolved rows; 438 unique by fixed GOA vector within class; 95 remain ambiguous.
- 4,301 graph candidate GeneIDs; only GeneID 10159 lacks a GPI159 mapping, and its fixed prediction is all zero.

### Independent matrix/signature validation

Command:

```bash
python scripts/independent_validate_full_label_classes_numpy.py ...
```

Final status: passed.

Independent implementation characteristics:

- pandas/NumPy matrix construction rather than per-gene Python sets;
- 121-dimensional vectors compared as packed-byte signatures;
- class equality checked by sorted signature sequences and independent counters.

Final output agrees exactly with the primary implementation:

- zero mismatches among 56,411 resolved rows;
- 183/183 unresolved classes exact;
- 438 uniquely assignable rows and 95 still ambiguous.

## Superseded attempt

The first independent NumPy implementation incorrectly represented `GO_ID -> column` as a one-to-one dictionary. The selected 121-column list contains three duplicated GO-vector/term situations, so repeated GO IDs in the provisional map were collapsed to one column. This produced 29,083 row mismatches and 68/183 exact classes.

The bug was corrected by representing `GO_ID -> list of columns`. The corrected implementation passed all validations. The failed result is not used scientifically.

## Uploaded dhimmel-file audit

Command:

```bash
python scripts/audit_uploaded_dhimmel_files.py --input-dir /mnt/data \
  --output analysis/B104F_uploaded_dhimmel_file_audit.csv
```

Status: completed; all four uploads rejected as scientific data.

Each file:

- begins with `<!DOCTYPE html>`;
- contains a GitHub HTML page title;
- does not contain a `go_id` TSV header;
- is approximately 147–151 KB.

No annotation rows were parsed from these files.

## Input-retention discipline

B104F consumed compact retained derivatives from frozen B104A/core artifacts. It did not rely on deleted B101–B104E conversation attachments. `graphsage_ppi.zip` remains an active original input. The MSigDB v5.2 raw archive was used only once to create a compact symbol-to-Entrez derivative included in B104F; the scientific class validation consumes that retained derivative.

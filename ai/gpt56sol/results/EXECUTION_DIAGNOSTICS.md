# Execution diagnostics — corrected run

Generated: 2026-08-24T16:10:11.352311+00:00

This is the error/diagnostic section that the previous report referred to but did not contain. Final results use only the successful runs below. Superseded failed attempts are documented separately so that no error is hidden.

## Final successful runs

### 1. Core verification

Command:

```bash
python -u scripts/run_core_verification.py \
  --input-dir /mnt/data \
  --work-dir work \
  --output-dir results
```

- Exit status: 0
- Reported status: PASS
- Runtime: 68.656 s
- Standard error: empty
- Complete standard output: `core_verification.stdout`
- Complete structured result: `core_verification_summary.json`

### 2. Local label-source screen

Command:

```bash
python -u scripts/run_local_label_source_screen.py \
  --input-dir /mnt/data \
  --work-dir work \
  --output-dir results
```

- Exit status: 0
- Reported status: PASS
- Structured runtime: 5.473 s; `/usr/bin/time` wall time: 8.22 s
- Warning: openpyxl reported an unknown workbook extension while reading one Greene supplementary XLSX; it removed that unsupported extension. Numeric parsing completed.
- Complete standard output: `local_label_source_screen.stdout`
- Complete standard error and timing: `local_label_source_screen.stderr`
- Complete structured result: `local_label_source_screen_summary.json`

### 3. Source-ledger construction

Command:

```bash
python scripts/build_source_manifest.py
```

- Exit status: 0
- Standard error: empty
- Output: `build_source_manifest.stdout`

### 4. Input checksum verification

Command:

```bash
python scripts/download_or_verify_sources.py \
  --manifest results/actual_input_file_manifest.csv \
  --dest inputs \
  --verify-only \
  --log results/input_verification_log.csv
```

- Exit status: 0
- 18/18 actual inputs: `verified_present`
- Standard error: empty
- Output: `input_verification.stdout`

### 5. Corrected report and manifest packaging

Command:

```bash
python scripts/build_corrected_report.py
```

- Exit status: 0
- Reported status: PASS
- Standard error: empty
- Output: `build_corrected_report.stdout`

### 6. Comprehensive source-ledger construction

Command:

```bash
python scripts/build_full_source_ledger.py
```

- Exit status: 0
- Reported status: PASS
- Records are separated by type: actual inputs, web references, historical candidates, analysis scripts, and stable generated outputs.
- Standard error: empty
- Output: `build_full_source_ledger.stdout`

### 7. Bundle validation

Command:

```bash
python scripts/validate_corrected_bundle.py
```

- Exit status: 0
- Reported status: PASS
- Checks include input existence and SHA-256, URL population for public inputs, historical candidate URLs, XLSX readability, presence of the embedded execution-diagnostics section, exact superseded-error text, absence of dangling report artifact references, and absence of missing tracked outputs.
- Structured result: `bundle_validation.json`
- Standard output: `bundle_validation.stdout`
- Standard error: empty

## Superseded attempts and errors

These attempts produced no accepted scientific result and were replaced by the successful runs above.

1. An initial core-verification invocation exceeded a 600-second command timeout after archive extraction. It produced no accepted output. The unchanged script was rerun under a direct unbuffered invocation and completed in 68.656 seconds; only that final run is used.
2. An early label-screen implementation raised `IndexError: list index out of range` while selecting a best C5-BP candidate after duplicate removal. The candidate-construction code was corrected.
3. The next label-screen attempt raised `ValueError: No candidate records for MSigDB 5.1 C5 BP`. Parsing was corrected to recognize the historical archive naming and to keep the specialized collection.
4. A non-optimized exhaustive label-screen attempt exceeded its time limit after completing MSigDB comparisons but before finishing the OhmNet same-tissue loop. Those partial files were discarded. Bitset-based comparisons and output reuse were implemented; the final accepted run completed successfully.

## Known scope limitations, not execution errors

- Historical GOA, GPAD, GPI, GO ontology, `gene2go`, `gp2protein`, and Bioconductor package bytes were not materialized in the corrected runtime. They are therefore not silently treated as tested.
- The manifest records candidate URLs but labels them `not_materialized` rather than implying that a download occurred.
- MSigDB direct downloads require authentication; the supplied archives are verified by SHA-256 and linked to the official download page.

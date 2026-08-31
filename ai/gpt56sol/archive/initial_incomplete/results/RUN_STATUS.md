# Run status

Generated: 2026-08-24T04:56:28.524240+00:00

- `core_summary.json`: MISSING
- `core_reconstruction_report.md`: MISSING
- `label_source_search_summary.json`: MISSING
- `label_source_search_report.md`: MISSING
- `MASTER_REPRODUCTION_REPORT.md`: present, 1,895 bytes
- `source_manifest.csv`: MISSING
- `priority_source_manifest.csv`: MISSING
- `repository_history_search.md`: present, 32 bytes
- `public_code_search.md`: present, 1,078 bytes
- `humanbase_discovery.md`: present, 695 bytes

## Nonempty stderr/error files

### `label_source_search.stderr`
```text
Traceback (most recent call last):
  File "/mnt/data/ppi_repro/scripts/label_source_search.py", line 44, in <module>
    if not coll.exists():raise FileNotFoundError('Run core_reproduction.py first: collapsed_gene_labels.csv absent')
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FileNotFoundError: Run core_reproduction.py first: collapsed_gene_labels.csv absent

```
# GraphSAGE PPI / DGL / gene-to-GO reconstruction

## Scope and independence

All numerical checks in this report were recomputed from the supplied GraphSAGE, DGL, OhmNet, MSigDB, Greene, and paper files. The prior investigation summary was treated only as a list of hypotheses. It was not used as a source of row mappings, tissue names, feature identities, or GO labels.

> **Core reconstruction did not complete.** See the error section and rerun `scripts/core_reproduction.py`.

## Gene-to-GO investigation

> **The historical GO search did not complete.** See the error section and rerun `scripts/label_source_search.py`.

## Source tracking and acquisition

`source_manifest.csv` and `source_manifest.json` are the authoritative input ledger. They include every supplied file and every attempted historical download, with URL, status, byte size, and SHA-256. Failed URLs are retained rather than omitted.

## Reproduction commands

Run from an environment with Python, NumPy, pandas, SciPy, scikit-learn, NetworkX, and openpyxl:

```bash
python scripts/inspect_inputs.py
python scripts/core_reproduction.py
python scripts/download_sources.py
python scripts/repo_history_search.py
python scripts/label_source_search.py
python scripts/build_master_report.py
```

The scripts never use the prior summary as a mapping input. Downloaded files are cached and hash-checked in `downloads/`.

## Execution diagnostics

### `label_source_search.stderr`
```text
Traceback (most recent call last):
  File "/mnt/data/ppi_repro/scripts/label_source_search.py", line 44, in <module>
    if not coll.exists():raise FileNotFoundError('Run core_reproduction.py first: collapsed_gene_labels.csv absent')
                         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
FileNotFoundError: Run core_reproduction.py first: collapsed_gene_labels.csv absent
```

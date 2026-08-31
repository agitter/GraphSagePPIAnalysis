# Recommended repository layout

```text
repo/
├── README.md
├── .gitignore
├── data/                         # local raw/reference data; normally not committed
│   ├── raw/
│   ├── derived/
│   ├── reference_papers/
│   ├── cache/
│   └── manifests/
├── gpt56sol/                     # this snapshot: scripts, reports, results, provenance
├── other_agent/                  # other investigator's reports/scripts/results
└── notebooks_or_analysis/        # your later clean external reproduction
```

Use `metadata/LOCAL_DIRECTORY_CLASSIFICATION_AND_MOVE_PLAN.csv` as the proposed path-by-path move plan for the directory listing supplied on 20260830T124820Z.

## Git recommendation

- Commit `gpt56sol/`, `other_agent/`, and small provenance manifests.
- Keep `data/raw/`, `data/cache/`, and most binary archives out of ordinary Git.
- Use download scripts plus hashes for reproducibility; use Git LFS only when a binary must be versioned.
- Do not commit credentials, authenticated MSigDB cookies, or transient download caches.
- Quarantine or delete the four `dhimmel-gene-ontology-...tsv` files that are actually saved GitHub HTML pages.

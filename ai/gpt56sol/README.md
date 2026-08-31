# gpt56sol complete investigation snapshot

Created: `20260830T124820Z`

This repo-ready snapshot contains every unique assistant-generated artifact currently recoverable from the working directories and historical delivery bundles. Historical bundles are expanded and deduplicated by SHA-256; the bundle ZIPs themselves are not nested.

## Contents

- **597 unique assistant artifacts** (78,078,137 bytes before snapshot metadata).
- **900 original occurrences/aliases** across direct files and prior bundles.
- **15 prior delivery bundles** with complete member-by-member SHA-256 coverage.
- **117 local paths** classified into `data/`, `gpt56sol/`, `other_agent/`, cache, or quarantine destinations.

## Start here

- `CURRENT_CANONICAL_INDEX.md`
- `RECOMMENDED_REPO_LAYOUT.md`
- `LOCAL_DIRECTORY_MOVE_PLAN.md`
- `metadata/GPT56SOL_COMPLETE_INVENTORY.xlsx`
- `metadata/ASSISTANT_ARTIFACT_INVENTORY.csv`
- `metadata/LOCAL_DIRECTORY_CLASSIFICATION_AND_MOVE_PLAN.csv`

## Snapshot policy

- All unique artifact bytes are retained once.
- Duplicate filenames and historical aliases are preserved in `ARTIFACT_OCCURRENCES_AND_ALIASES.csv`.
- Prior delivery bundles are omitted as redundant containers; `OMITTED_REDUNDANT_BUNDLE_COVERAGE.csv` proves that every bundle member is represented by SHA-256.
- Early incomplete or superseded deliveries are retained under `archive/` or marked in the inventory.
- Raw external datasets are not added to this snapshot; the local move plan places those under `data/`.

## Artifact status counts

| Status | Files |
|---|---:|
| `current_core` | 86 |
| `historical_foundation` | 80 |
| `historical_or_supporting` | 138 |
| `superseded_initial_delivery` | 5 |
| `superseded_or_early_summary` | 1 |
| `supporting_analysis` | 287 |

## Artifact category counts

| Category | Files |
|---|---:|
| `analysis_result_or_derived_data` | 186 |
| `batch_administration` | 16 |
| `documentation` | 37 |
| `other_artifact` | 46 |
| `provenance_manifest` | 87 |
| `provenance_workbook` | 1 |
| `reference_pack_or_nested_bundle` | 1 |
| `report_or_documentation` | 27 |
| `retained_normalized_input` | 52 |
| `script` | 84 |
| `validation_or_checksum` | 60 |

## Verification

Run:

```bash
python gpt56sol/scripts/verify_snapshot.py
```

`metadata/SNAPSHOT_FILE_MANIFEST.csv` intentionally excludes itself and the ZIP container to avoid a self-referential checksum.

# Local directory classification and move plan

This is a filename/path-based classification of the local `find .` listing supplied by the user. It does not silently assert hashes for files that were not uploaded or previously inventoried.

## Summary by origin class

| Origin class | Entries |
|---|---:|
| `external_current_mapping_output` | 1 |
| `external_raw_data` | 35 |
| `external_raw_data_downloaded_by_gpt56sol_script` | 4 |
| `external_raw_data_or_code` | 5 |
| `external_raw_or_derived_data` | 4 |
| `external_reference_document` | 5 |
| `gpt56sol_artifact_or_script_output` | 9 |
| `gpt56sol_script` | 6 |
| `gpt56sol_script_output` | 36 |
| `invalid_download_html_disguised_as_tsv` | 4 |
| `other_agent_or_prior_investigation` | 7 |
| `working_cache` | 1 |

## Important actions

- Move large external source files to `data/` and normally keep them out of ordinary Git.
- Move assistant scripts, reports, script outputs, and provenance records to `gpt56sol/`.
- Keep the other agent’s artifacts separate under `other_agent/` so claims remain attributable.
- Delete or quarantine the four `dhimmel-gene-ontology-...tsv` files; prior inspection established they are HTML pages, not annotation TSVs.
- Treat `goa_date_screen_cache/` as disposable working cache, not a canonical result.
- Preserve the misspelled `HuamnBase-kidney.dat` filename until content/hash evidence supports a safe rename.

The complete path-level plan is in `metadata/LOCAL_DIRECTORY_CLASSIFICATION_AND_MOVE_PLAN.csv` and the workbook.

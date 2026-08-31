# Corrected run status

Generated: 2026-08-24T16:10:16.142843+00:00

## Scientific stages

| Stage | Status | Basis |
|---|---|---|
| Actual-input checksum verification | PASS | 18/18 supplied inputs verified |
| GraphSAGE structure and 24-block derivation | PASS | Unique safe-cut partition and OhmNet statistic assignment |
| Immediate OhmNet topology provenance | PASS | 816,950/816,950 independently mapped edges match |
| Topology-only gene mapping | PASS with unresolved equivalence classes | 55,878 rows unique |
| Feature-assisted gene disambiguation | PASS with residual ambiguity | 4,268/4,301 unique genes mapped |
| MSigDB feature provenance | PASS with one unidentifiable zero column | 49 nonzero exact unique columns; column 10 all-zero |
| DGL transformation | PASS | All graph IDs, labels, float64 standardized features, and edge sets exact |
| Conservative leakage measurement | PASS | 5,430/5,524 test rows seen; lookup F1 0.9940678477 |
| Local label-source exclusion screen | PASS | MSigDB, OhmNet, and Greene restrictions tested |
| Historical GOA/gene2go/Bioconductor source grid | NOT RUN | Candidate bytes not materialized; no numerical claim made |
| Comprehensive source ledger | PASS | 99 records: {'actual_input': 18, 'analysis_script': 8, 'generated_output': 34, 'historical_candidate_not_materialized': 31, 'web_reference': 8} |
| Programmatic bundle validation | PASS | See `bundle_validation.json`; validation outputs are not self-hashed in this status table |

## Output integrity

| File | Status | Bytes | SHA-256 |
| --- | --- | --- | --- |
| `MASTER_REPRODUCTION_REPORT.md` | present | 25,609 | `0a1df00dbc4d6c5c3d3f15cfd0980632688445f69f9290d3559489b5508617b5` |
| `MASTER_REPRODUCTION_REPORT_v2.md` | present | 25,609 | `0a1df00dbc4d6c5c3d3f15cfd0980632688445f69f9290d3559489b5508617b5` |
| `EXECUTION_DIAGNOSTICS.md` | present | 4,657 | `0b651d11e14fe6f826d51ed447a01a712210ca0bfbaa5566d9dfbc0135747b02` |
| `SOURCE_ACQUISITION.md` | present | 12,529 | `f2a003e6391699b88f3c2ffb5cd018da29dcc16a8ebeeded3ada621be546bcd3` |
| `actual_input_file_manifest.csv` | present | 8,658 | `37d802b1283e0a8c68e626d64ac39ea0d77b3f6fad0323e4edb2b8734ec94a77` |
| `actual_input_file_manifest.json` | present | 15,537 | `59b7f1c6d189e3490b76eafaf4b6a13772bc87936fbc639a0de46d7181e4ddee` |
| `actual_input_file_manifest.xlsx` | present | 8,973 | `642b6c6ac9ef1fdb3a54ebd83184718ef81f0d5fea8417731f7def4765c1b1b1` |
| `source_ledger.csv` | present | 40,525 | `99d2e3cf69886206c449b7b3967e88e756dd5227c97b9b31e34b7e685b13c003` |
| `source_ledger.json` | present | 79,389 | `15b8c5c589ac090ff5cb9567d7cb4ad09ea546a678a4464e5c7b16809f0664ea` |
| `source_ledger.md` | present | 22,543 | `ed5b5517ca04d9d46aa1a40862aacbfbe03d7a6c9637a365c97ec6bff314736d` |
| `source_ledger.xlsx` | present | 18,780 | `f650aece2fc07d78d8a9a4cc04f0bb4bdefa934cdca6d745c4804b89bbcef972` |
| `input_verification_log.csv` | present | 4,107 | `c62c6d36b52177efea7e8ecd32ad0dc089ad5a8e7a5361a57632978893b1ce95` |
| `core_verification_summary.json` | present | 6,123 | `3d4e3624a4d00d6c8c85dcdaab01f951c55a0b6cdf338c97a492ce3930983ec0` |
| `tissue_partition.csv` | present | 1,550 | `064c396d394b90428ae2f73768c44e470a776259c1d69bea26d028fdd50fb059` |
| `wl_tissue_summary.csv` | present | 991 | `ca0faa1d0d7eaa8ef63a2f531dba557a3796992f45c0b68ec3acc34260326cb2` |
| `wl_ambiguous_topology_only.json` | present | 69,252 | `1f9c89ee87b6c5c8a5deae2340fba7237f5f5dbfec226a1907a22a9670c8d584` |
| `wl_residual_after_features.json` | present | 36,578 | `664ab9a55cb4087316a632be8eaa92bd5da72dff4a80fe4a94976f1a28dced93` |
| `feature_column_mapping.csv` | present | 7,670 | `be22af448afae0aebccafc16e94e3e916caf2cbf7d46e1686b7d475c29e4939a` |
| `graphsage_row_to_entrez_topology_features.csv` | present | 1,201,265 | `4aaada1751c958db2653e19e05ca5a8e2bdd2da86c1a4441bad26e0468ae7207` |
| `dgl_component_assignment.csv` | present | 861 | `b76360dbf149133385821fe981517c844b191406800e2f67d74659eb21eff249` |
| `dgl_split_verification.csv` | present | 701 | `5a0823ae7eb4f6a087b79ed3e2e7232e8d13d02c82d30178594d8e30d431d551` |
| `dgl_transformation_verification.json` | present | 2,582 | `b4d60f38f9ca2111ccb3b5eabe0022e5c1e067ebc4e3890864e4d07c1a707cda` |
| `collapsed_gene_labels_topology_features.csv` | present | 1,061,035 | `4fd88002b2600c0b7c5dcb076390f5b74a39edcb34184f1b0768d9926dcb1907` |
| `local_label_source_screen_summary.json` | present | 3,795 | `7e3cdf99fea507db9571f28847581492811645ffaaf378d7c618f223a0aca726` |
| `msigdb_label_source_screen.csv` | present | 213,568 | `3d11550adc4a24853719b44643bff95ac34370f36106ae18d76de286befebddd` |
| `ohmnet_global_label_source_screen.csv` | present | 99,457 | `257de22d3d83a4a8b0c3de4f738afe62471964f420290815712680042225a90b` |
| `ohmnet_same_tissue_label_screen_best.csv` | present | 12,371 | `cae495550f7700b009c1b0ddf4f12f24e747b26c624e53dc4d0ff79a90031101` |
| `core_verification.stdout` | present | 1,652 | `df2b49c9a05f7e960b7fd61d40607959d8de1bc77ee1b81e6a6d3a13084c08e1` |
| `core_verification.stderr` | present | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `local_label_source_screen.stdout` | present | 3,846 | `3ae6d4f049e5ae3aac84826e935d3f545053cee7842c41dce16f1121fcfaa9a9` |
| `local_label_source_screen.stderr` | present | 1,079 | `bfbaef32806945b073b4ec4009a435fe9a4565e64c3364a17e4519e35d1561f7` |
| `build_source_manifest.stdout` | present | 270 | `a5454693ab548e0703ca0272be8d1d502bda780c6d2a4127c81bcec9b509e6fb` |
| `build_source_manifest.stderr` | present | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `build_corrected_report.stdout` | present | 180 | `53691bcd56c7b89c6512aa369695c8ad7ec001194343ba6fc0a8dece24dd23ca` |
| `build_corrected_report.stderr` | present | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `build_full_source_ledger.stdout` | present | 205 | `c9b9cf9b3db4117158acb163bea1178315c0a487aaeba10904153bc114fab2d3` |
| `build_full_source_ledger.stderr` | present | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `input_verification.stdout` | present | 131 | `cd06fd1cfa0b00c7aeb62c0ad7631af16a1ded134d9f3b11f1ff454e24ee2f5b` |
| `input_verification.stderr` | present | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` |
| `REQUESTED_ADDITIONAL_INPUTS.md` | present | 832 | `9187d31017c77817bcc22337b231649a4ffdb5b9fdc53e4945431edf37b14553` |

## Diagnostics

`MASTER_REPRODUCTION_REPORT.md` embeds the exact superseded timeout and exception descriptions. `EXECUTION_DIAGNOSTICS.md` additionally records all accepted commands, exit statuses, warnings, logs, and scope limits.

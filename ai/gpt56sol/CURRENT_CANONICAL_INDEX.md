# Current canonical / start-here index

This index points to the latest core reconstruction and provenance artifacts. Earlier batches remain in the snapshot for auditability.

## Recommended reading order

1. `batches/B104G_20260829T150633Z/B104G_REPORT_20260829T151452Z.md` — complete row identity, features, labels, and ambiguity review.
2. `batches/B104H_20260830T110259Z/GRAPHSAGE_NODE_MAPPING_README_20260830T110259Z.md` — node-to-Entrez mapping schema and interpretation.
3. `batches/B104H_20260830T110259Z/B104H_SOURCE_PROVENANCE_AND_NEXT_PRIORITIES_20260830T110259Z.md` — source-provenance assessment.
4. `batches/B104I_20260830T114918Z/B104I_PPI_PROVENANCE_MENCHE_AND_SPLIT_REPORT_20260830T114918Z.md` — PPI provenance and split semantics.
5. `batches/B104I_20260830T114918Z/B104I_EVIDENCE_VS_LITERATURE_REGISTER_20260830T114918Z.csv` — exact/inferred/documented/open claim register.
6. `batches/B104A_20260828T145842Z/B104A_REPORT_20260828T145842Z.md` — exact GO-label transformation breakthrough.
7. `batches/B104E_20260829T121535Z/B104E_GOA_DATE_SCREEN_UNIPROT_AND_FEATURE_RULE_REPORT_20260829T121535Z.md` — GOA date and MSigDB feature results.
8. `batches/B104C_20260828T194921Z/B104C_MSIGDB50_COLUMN_ORDER_REPORT_20260828T200948Z.md` — label-column order fingerprint.

## Latest provenance snapshots

- `results/actual_input_file_manifest_through_B104I_20260830T114918Z.csv`
- `results/source_ledger_through_B104I_FINAL_20260830T114918Z.csv`
- `results/provenance_events_through_B104I_20260830T114918Z.csv`

## Current core files discovered automatically

- `batches/B104I_20260830T114918Z/B104I_DELETION_CLEARANCE_20260830T114918Z.md` — Batch-specific record authorizing deletion of raw conversation attachments after retention checks.
- `batches/B104I_20260830T114918Z/B104I_EVIDENCE_VS_LITERATURE_REGISTER_20260830T114918Z.csv` — Claim-by-claim register comparing file-derived evidence with OhmNet, GraphSAGE, and related literature.
- `batches/B104I_20260830T114918Z/B104I_GRAPHSAGE_GITHUB_ISSUES_TRACKER_20260830T114918Z.csv` — Tracker of GraphSAGE GitHub issues relevant to dataset provenance, mappings, features, and evaluation.
- `batches/B104I_20260830T114918Z/B104I_GRAPHSAGE_GITHUB_ISSUES_TRACKER_20260830T114918Z.md` — Tracker of GraphSAGE GitHub issues relevant to dataset provenance, mappings, features, and evaluation.
- `batches/B104I_20260830T114918Z/B104I_PPI_PROVENANCE_MENCHE_AND_SPLIT_REPORT_20260830T114918Z.md` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `batches/B104I_20260830T114918Z/derived/graphsage_ppi_entrez_split_membership_20260830T114918Z.tsv.gz` — One-row-per-Entrez summary of train/validation/test tissue-graph membership.
- `bundle_only/B104I_PPI_provenance_and_literature_bundle_20260830T114918Z/B104I_VALIDATION_20260830T114918Z.json` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `bundle_only/B104I_PPI_provenance_and_literature_bundle_20260830T114918Z/results/provenance_events_through_B104I_20260830T114918Z.csv` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `bundle_only/B104I_PPI_provenance_and_literature_bundle_20260830T114918Z/results/source_ledger_through_B104I_20260830T114918Z.csv` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `bundle_only/B104I_PPI_provenance_and_literature_bundle_20260830T114918Z/retained_inputs/README_20260830T113352Z.md` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `bundle_only/B104I_PPI_provenance_and_literature_bundle_20260830T114918Z/retained_inputs/biosnap_ppi_provenance_audit_20260830T113352Z.json` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `bundle_only/B104I_PPI_provenance_and_literature_bundle_20260830T114918Z/retained_inputs/biosnap_ppi_provenance_audit_20260830T113352Z.zip` — PPI topology provenance audit or report linking GraphSAGE, OhmNet, and BioSNAP edge sets.
- `results/B104I_FINAL_DELIVERY_VALIDATION_20260830T114918Z.json` — Machine-readable validation checks for the associated analysis or delivery.
- `results/actual_input_file_manifest_through_B104I_20260830T114918Z.csv` — Manifest of actual materialized inputs used by the analysis at that point in time.
- `results/provenance_events_through_B104I_20260830T114918Z.csv` — Append-only provenance-event history, including receipt, validation, and deletion states.
- `results/source_ledger_through_B104I_FINAL_20260830T114918Z.csv` — Append-only source and artifact provenance ledger.
- `batches/B104H_20260830T110259Z/B104H_MAPPING_VALIDATION_20260830T110259Z.json` — Machine-readable validation checks for the associated analysis or delivery.
- `batches/B104H_20260830T110259Z/B104H_SOURCE_PROVENANCE_AND_NEXT_PRIORITIES_20260830T110259Z.md` — assistant-generated investigation artifact: B104H_SOURCE_PROVENANCE_AND_NEXT_PRIORITIES_20260830T110259Z.md.
- `batches/B104H_20260830T110259Z/GRAPHSAGE_NODE_MAPPING_README_20260830T110259Z.md` — report or documentation artifact: GRAPHSAGE_NODE_MAPPING_README_20260830T110259Z.md.
- `batches/B104H_20260830T110259Z/LEAKAGE_EXPERIMENT_BACKLOG_20260830T110259Z.md` — Saved design for the gene-identity lookup, randomized-label, and gene-disjoint leakage experiments.
- `batches/B104H_20260830T110259Z/build_B104H_node_mapping_artifacts.py` — analysis, download, validation, or packaging script: build_B104H_node_mapping_artifacts.py.
- `batches/B104H_20260830T110259Z/download_and_verify_biosnap_ppi_sources.py` — Downloader and exact edge-set audit for BioSNAP OhmNet and global PPI sources.
- `batches/B104H_20260830T110259Z/graphsage_ppi_node_to_entrez_20260830T110259Z.tsv.gz` — Compact GraphSAGE row/node-to-Entrez mapping for all 56,944 node rows.
- `batches/B104H_20260830T110259Z/graphsage_ppi_node_to_entrez_evidence_20260830T110259Z.tsv.gz` — Evidence-rich GraphSAGE row/node-to-Entrez mapping with topology, feature, label, and legacy-dictionary checks.
- `bundle_only/B104H_node_mapping_provenance_bundle_20260830T110259Z/B104H/B104H_BUNDLE_README_20260830T110259Z.md` — B104H report or documentation artifact: B104H_BUNDLE_README_20260830T110259Z.md.
- `bundle_only/B104H_node_mapping_provenance_bundle_20260830T110259Z/B104H/B104H_FILE_CHECKSUMS_20260830T110259Z.csv` — B104H assistant-generated investigation artifact: B104H_FILE_CHECKSUMS_20260830T110259Z.csv.
- `bundle_only/B104H_node_mapping_provenance_bundle_20260830T110259Z/B104H/manifests/provenance_events_through_B104H_PREBUNDLE_20260830T110259Z.csv` — Append-only provenance-event history, including receipt, validation, and deletion states.
- `bundle_only/B104H_node_mapping_provenance_bundle_20260830T110259Z/B104H/manifests/source_ledger_through_B104H_PREBUNDLE_20260830T110259Z.csv` — Append-only source and artifact provenance ledger.
- `results/B104H_FINAL_DELIVERY_VALIDATION_20260830T110259Z.json` — Machine-readable validation checks for the associated analysis or delivery.
- `results/provenance_events_through_B104H_FINAL_20260830T110259Z.csv` — Append-only provenance-event history, including receipt, validation, and deletion states.
- `results/source_ledger_through_B104H_FINAL_20260830T110259Z.csv` — Append-only source and artifact provenance ledger.
- `batches/B104G_20260829T150633Z/B104G_CLAIM_STATUS_20260829T151452Z.csv` — analysis result or compact derived dataset: B104G_CLAIM_STATUS_20260829T151452Z.csv.
- `batches/B104G_20260829T150633Z/B104G_DELETION_CLEARANCE_20260829T151452Z.md` — Batch-specific record authorizing deletion of raw conversation attachments after retention checks.
- `batches/B104G_20260829T150633Z/B104G_INDEPENDENT_FULL_FEATURE_VALIDATION_20260829T150633Z.json` — Independent full 56,944-row feature reconstruction validation.
- `batches/B104G_20260829T150633Z/B104G_INDEPENDENT_FULL_LABEL_VALIDATION_20260829T150633Z.json` — Independent full 56,944-row GO-label reconstruction validation.
- `batches/B104G_20260829T150633Z/B104G_INDEPENDENT_ROW_ORDER_VALIDATION_20260829T150633Z.json` — Independent validation of the recovered CPython 2/NetworkX row order.
- `batches/B104G_20260829T150633Z/B104G_REPORT_20260829T151452Z.md` — Narrative scientific and provenance report for the associated analysis batch.
- `batches/B104G_20260829T150633Z/B104G_REPRODUCTION_INSTRUCTIONS_20260829T151452Z.md` — Instructions for reproducing the associated analysis from retained inputs.
- `batches/B104G_20260829T150633Z/B104G_UPDATED_AGENT_HANDOFF_20260829T151452Z.md` — Investigator handoff summarizing methods, findings, caveats, and next tests.
- `batches/B104G_20260829T150633Z/derived/B104G_full_4301_gene_universe.csv.gz` — Deduplicated 4,301-gene universe derived from the complete row mapping.
- `batches/B104G_20260829T150633Z/derived/B104G_full_graphsage_row_to_entrez_mapping.csv.gz` — Complete GraphSAGE row-to-Entrez mapping produced by the recovered legacy Python/NetworkX node order.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/B104G_BUNDLE_README_20260829T151452Z.md` — B104G report or documentation artifact: B104G_BUNDLE_README_20260829T151452Z.md.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/B104G_EXECUTION_DIAGNOSTICS_20260829T151452Z.md` — Execution diagnostics, including accepted runs, failures, retries, and limitations.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/B104G_FILE_CHECKSUMS_20260829T151452Z.csv` — B104G assistant-generated investigation artifact: B104G_FILE_CHECKSUMS_20260829T151452Z.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_all_OhmNet_network_sizes_and_GraphSAGE_selection.csv` — B104G analysis result or compact derived dataset: B104G_all_OhmNet_network_sizes_and_GraphSAGE_selection.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_dhimmel_annotation_file_summary.csv` — B104G analysis result or compact derived dataset: B104G_dhimmel_annotation_file_summary.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_dhimmel_column_order_controls.csv` — B104G analysis result or compact derived dataset: B104G_dhimmel_column_order_controls.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_dhimmel_per_column_comparison.csv.gz` — B104G analysis result or compact derived dataset: B104G_dhimmel_per_column_comparison.csv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_dhimmel_raw_to_retained_reconciliation.json` — B104G analysis result or compact derived dataset: B104G_dhimmel_raw_to_retained_reconciliation.json.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_dhimmel_selected_gene_audit.csv` — B104G analysis result or compact derived dataset: B104G_dhimmel_selected_gene_audit.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_duplicate_GO_candidate_full_universe_test.csv` — B104G analysis result or compact derived dataset: B104G_duplicate_GO_candidate_full_universe_test.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_former_equivalence_classes_resolved_by_node_order.csv` — B104G analysis result or compact derived dataset: B104G_former_equivalence_classes_resolved_by_node_order.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_full_121_label_validation.csv` — Machine-readable validation checks for the associated analysis or delivery.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_full_24_tissue_node_edge_verification.csv` — B104G analysis result or compact derived dataset: B104G_full_24_tissue_node_edge_verification.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_full_50_feature_validation.csv` — Machine-readable validation checks for the associated analysis or delivery.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_full_mapping_leakage_recalculation.csv` — B104G analysis result or compact derived dataset: B104G_full_mapping_leakage_recalculation.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_graph_selection_open_question_summary.json` — B104G analysis result or compact derived dataset: B104G_graph_selection_open_question_summary.json.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_mapping_components_touching_full_graph_gene_universe.csv` — Bipartite identifier-mapping components and ambiguity-resolution evidence.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_node_order_model_controls.csv` — B104G analysis result or compact derived dataset: B104G_node_order_model_controls.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/analysis/B104G_special_and_unmapped_gene_rows.csv` — B104G analysis result or compact derived dataset: B104G_special_and_unmapped_gene_rows.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/logs/analyze_B104G.stderr` — B104G assistant-generated investigation artifact: analyze_B104G.stderr.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/logs/analyze_B104G.stdout` — B104G assistant-generated investigation artifact: analyze_B104G.stdout.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/GO_annotations-9606-direct-allev.tsv.gz` — B104G compact normalized retained input used for later reproducibility: GO_annotations-9606-direct-allev.tsv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/GO_annotations-9606-direct-expev.tsv.gz` — B104G compact normalized retained input used for later reproducibility: GO_annotations-9606-direct-expev.tsv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/GO_annotations-9606-inferred-allev.tsv.gz` — B104G compact normalized retained input used for later reproducibility: GO_annotations-9606-inferred-allev.tsv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/GO_annotations-9606-inferred-expev.tsv.gz` — B104G compact normalized retained input used for later reproducibility: GO_annotations-9606-inferred-expev.tsv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz` — Row-preserving normalized derivative of GOA human GAF release 159.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz` — Row-preserving normalized derivative of GOA human GPI release 159.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/B104G_prior_compact_input_checksums.sha256` — B104G assistant-generated investigation artifact: B104G_prior_compact_input_checksums.sha256.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz` — B104G compact normalized retained input used for later reproducibility: B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz` — B104G compact normalized retained input used for later reproducibility: B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/core_verification_summary.json` — B104G compact normalized retained input used for later reproducibility: core_verification_summary.json.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/graphsage_row_to_entrez_topology_features.csv` — B104G compact normalized retained input used for later reproducibility: graphsage_row_to_entrez_topology_features.csv.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/retained_inputs/prior_compact_inputs/wl_residual_after_features.json` — B104G compact normalized retained input used for later reproducibility: wl_residual_after_features.json.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/scripts/analyze_B104G_full_row_mapping_and_dhimmel.py` — B104G analysis, download, validation, or packaging script: analyze_B104G_full_row_mapping_and_dhimmel.py.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/scripts/independent_validate_B104G_full_features.py` — B104G analysis, download, validation, or packaging script: independent_validate_B104G_full_features.py.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/scripts/independent_validate_B104G_full_labels.py` — B104G analysis, download, validation, or packaging script: independent_validate_B104G_full_labels.py.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/scripts/independent_validate_B104G_row_order.py` — B104G analysis, download, validation, or packaging script: independent_validate_B104G_row_order.py.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/B104G/scripts/prototype_test_node_order.py` — B104G analysis, download, validation, or packaging script: prototype_test_node_order.py.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/manifests/actual_input_file_manifest_through_B104G_20260829T151452Z.md` — Manifest of actual materialized inputs used by the analysis at that point in time.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/manifests/provenance_events_through_B104G_PREBUNDLE_20260829T151452Z.csv` — Append-only provenance-event history, including receipt, validation, and deletion states.
- `bundle_only/B104G_full_row_and_label_resolution_bundle_20260829T151452Z/manifests/source_ledger_through_B104G_PREBUNDLE_20260829T151452Z.csv` — Append-only source and artifact provenance ledger.
- `results/B104G_FINAL_DELIVERY_VALIDATION_20260829T151452Z.json` — Machine-readable validation checks for the associated analysis or delivery.
- `results/actual_input_file_manifest_through_B104G_20260829T151452Z.csv` — Manifest of actual materialized inputs used by the analysis at that point in time.
- `results/provenance_events_through_B104G_FINAL_20260829T151452Z.csv` — Append-only provenance-event history, including receipt, validation, and deletion states.
- `results/source_ledger_through_B104G_FINAL_20260829T151452Z.csv` — Append-only source and artifact provenance ledger.

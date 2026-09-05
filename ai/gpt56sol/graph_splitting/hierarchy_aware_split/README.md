# GraphSAGE PPI hierarchy-aware split investigation

This bundle contains the complete analysis supporting a hierarchy-aware train-validation-test split study for the GraphSAGE PPI / OhmNet tissue networks.

## Start here

- `EXECUTIVE_SUMMARY.md`: concise conclusions and recommended design.
- `REPORT.md`: methods, exact results, interpretation, caveats, and quality control.
- `RESULTS_SUMMARY.json`: machine-readable headline findings.
- `SOURCE_MANIFEST.csv`: input filenames, sizes, hashes, and roles.
- `ANALYSIS_MANIFEST.json`: exact-enumeration and Monte Carlo metadata.
- `results/CLAIMS.csv`: claim/evidence/confidence register.
- `results/discrepancy_register_additions.csv`: report-ready discrepancy entries.

## Recommended split resources

- `results/primary_leaf_split_tissues.tsv`: maximum-separation leaf-only stress-test split.
- `results/broad_coverage_leaf_split_tissues.tsv`: broad-training-branch sensitivity split.
- `results/method_definitions.csv`: precise construction rules.
- `results/primary_metric_sensitivity.csv`: agreement across hierarchy metrics.

## Comparison results

- `results/headline_comparison.csv`: main summary table.
- `results/paired_conditional_on_heldout_summary.csv`: hierarchy-aware minus matched-random paired effects.
- `results/leaf107_size_matched_conditional_summary.csv`: tight node/edge-size matching.
- `results/null_distribution_summary.csv`: one-million-split null summaries.
- `results/candidate_vs_null_percentiles.csv`: candidate locations in null distributions.
- `results/actual_split_metrics_recomputed.csv`: deposited validation/test metrics.

## Exact enumeration tables

Files matching:

    results/<universe>__<method>.tsv.gz

contain one row per feasible heldout quartet, including the deterministically selected training set and post hoc overlap metrics. The primary table is:

    results/leaf107__branch_distinct_node_stratified_minimax_wup.tsv.gz

It contains all 1,062 feasible primary leaf-only heldout quartets.

## Plots

- `plots/leaf107_mean_wup_ecdf.png`
- `plots/leaf107_row_overlap_ecdf.png`
- `plots/leaf107_lookup_f1_ecdf.png`
- `plots/leaf107_paired_wup_vs_f1_effect.png`
- `plots/leaf107_global_method_tradeoff.png`

## Reproducibility and verification

Core source code is under `scripts/`. Important checks are under `tests/`:

- `primary_leaf_exact_enumeration_independent_check_summary.csv`
- `independent_exact_optima_check.csv`
- `null_deterministic_rerun_sha256.txt`
- `null_first1000_independent_metric_check.csv`
- `leaf107_conditional_sample_coverage.csv`

The original 447 MB Monte Carlo binary streams are intentionally omitted from the release ZIP. They are deterministically regenerable from the source code and seeds in `ANALYSIS_MANIFEST.json`; compact summaries and the first 1,000 records of each stream are retained.

Run:

    python verify_bundle.py

from the extracted bundle directory to verify every released file listed in `SHA256SUMS.txt`.

## Interpretation boundary

Split optimization used hierarchy, graph size, thresholds, and leaf/internal status only. Gene identities and labels were used solely after split construction to measure overlap and the lookup diagnostic. The analysis does not claim to recover the authors' historical random split and does not estimate GraphSAGE neural-model performance on the proposed splits.

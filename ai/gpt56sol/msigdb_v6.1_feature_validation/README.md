# MSigDB v6.1 GraphSAGE PPI feature validation

This compact evidence package records the one-time test of MSigDB v6.1 C1/C3 as the canonical input for reconstructing the 50 GraphSAGE PPI features.

The result is **PASS**: all 50 selected memberships agree with the retained v6.0 screen, and all 2,847,200 released feature cells are reproduced exactly. See `REPORT.md` for interpretation and `validation_summary.json` for machine-readable details.

The original MSigDB and GraphSAGE inputs are not duplicated in this package. Their exact filenames, sizes, and SHA-256 hashes are recorded in `source_manifest.csv`.

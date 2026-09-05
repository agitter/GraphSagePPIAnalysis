# Attribution and provenance

## Original research and data

- The GraphSAGE manuscript and released PPI benchmark are attributed to the GraphSAGE authors. The exact supplied manuscript and data archives are identified by SHA-256 in `SOURCE_MANIFEST.csv`.
- The tissue networks, hierarchy, labels, and README are attributed to the OhmNet data authors. The exact supplied archives are identified by SHA-256 in `SOURCE_MANIFEST.csv`.

## User contribution

The user supplied the source archives, defined the forensic research questions, emphasized attribution and non-circular reproducibility, and required the manuscript/data discrepancies and the brain-midbrain containment result to be retained for the eventual report.

## Prior project work reused here

The prior GraphSAGE PPI investigation established the canonical node-to-Entrez mapping, exact label reconstruction, split membership, hierarchy-layer classification, and the original lookup diagnostic. This analysis used the fixed prior evidence bundle listed in `SOURCE_MANIFEST.csv`; it did not silently reconstruct or replace those results.

## Current analysis contribution

OpenAI GPT-5.6 Pro, under the user's direction, designed and executed the hierarchy-aware split enumeration, matched null sampling, paired and size-matched comparisons, metric sensitivity analysis, quality-control checks, plots, and this evidence bundle on 2026-09-04.

## Evidence classes

- **Source-derived:** statements read from the supplied GraphSAGE manuscript or OhmNet README.
- **Directly computed:** graph counts, hierarchy relationships, split metrics, exact enumerations, and Monte Carlo summaries calculated from the supplied bytes.
- **Prior-project evidence:** canonical node identities and reconstructed labels imported from the checksummed prior evidence bundle.
- **Methodological judgment:** the recommendation to use leaf-only branch blocking, size stratification, minimax hierarchy optimization, and paired nulls.
- **Unresolved history:** the authors' literal historical split mechanism and the reason the deposited split does not satisfy the current-release threshold interpretation.

## Non-circularity rule

Gene identities, labels, overlap measurements, and lookup F1 were not used to choose hierarchy-aware splits. They were evaluated only after a split was fixed from hierarchy, graph-size thresholds, leaf/internal status, and deterministic tie-breaks.

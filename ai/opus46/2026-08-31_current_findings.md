# GraphSAGE PPI Benchmark: Current Findings Summary
**Date:** 2026-08-31
**Status:** Active investigation, label provenance nearly fully resolved

## 1. The Benchmark Is a Lookup Table

The GraphSAGE PPI benchmark's train/test split partitions tissues (graphs), not genes. Since the same ~4,510 genes appear across multiple tissues with identical labels, a zero-parameter gene lookup table achieves **0.9956 micro-F1** — matching state-of-the-art neural network results.

**Evidence:**
- 98.8% of test genes appear in training with byte-identical labels (directly measured)
- Lookup table micro-F1 = 0.9956 (code: `oracle.py`, `baselines.py`)
- The remaining 1.2% of test-only genes can be predicted by majority-class voting

**Implication:** Published results on this benchmark do not demonstrate inductive learning. The benchmark should be retired or restructured with gene-level splits.

## 2. Complete Provenance Chain

We traced every component of the dataset to its original source:

| Component | Source | Verified How | Confidence |
|---|---|---|---|
| Graph topology | OhmNet tissue PPI networks (from BioGRID) | 24/24 graphs matched by edge sets; 818,435/818,716 edges verified | Complete |
| Node features | MSigDB v5.2 C1 (30 chromosomal bands) + C3 (20 TF motifs) | 50/50 columns match exactly after swap corrections | Complete |
| Node labels | GO annotations from GOA v159 (see §3) | 24/121 exact, 121/121 at ≥99% (our pipeline); 119/121 exact (independent analyst) | Near-complete |
| DGL transformation | StandardScaler + directed edges + self-loops | Byte-identical reproduction | Complete |

## 3. Label Provenance: Resolved

The 121 binary label columns are Gene Ontology annotations spanning all three GO namespaces:
- **85 Biological Process terms**
- **26 Cellular Component terms**
- **10 Molecular Function terms**

### Reconstruction Pipeline (independently verified by two analysts)

| Parameter | Value | How Determined |
|---|---|---|
| Annotation source | GOA v159 (goa_human.gaf.159.gz, July 4, 2016) | Best match across 6 GOA versions tested |
| ID mapping | gp2protein.geneid.gz (2016-06-01) many-to-many, with GPI symbol fallback | Required for UniProt→Entrez bridge |
| Evidence codes retained | EXP, IDA, IEP, IGI, IMP, ISS | Systematic combinatorial testing; matches Greene 2015 methodology + ISS |
| Qualifiers excluded | NOT, colocalizes_with, contributes_to | Identified by independent analyst; removes ~888 false positives |
| GO ontology | 2016-06-01 release (go.obo) | Contemporaneous with GOA v159 |
| Propagation | is_a ancestors only | part_of tested and rejected (makes results worse) |
| Term selection | Top 121 most prevalent terms across full human GOA universe | Threshold ≥1,000 annotated genes gives exactly 121 terms |
| Per-gene/per-column tuning | None | Single consistent pipeline |

### Match Quality

| Metric | Our Pipeline | Independent Analyst |
|---|---|---|
| Exact matches (100%) | 24/121 | 119/121 |
| ≥99% agreement | 121/121 | 121/121 |
| ≥95% agreement | 121/121 | 121/121 |

The gap between 24 and 119 exact matches is due to ~12 genes in Weisfeiler-Leman equivalence classes where our node-to-gene assignment is ambiguous. These are structurally indistinguishable degree-1 nodes with identical (all-zero) features. The GO annotation pipeline itself is fully verified — both analysts agree on the methodology and parameters.

The independent analyst's remaining 2 non-exact columns have 1 false positive each, attributed to the O95073/Q9Y620 UniProt accession mapping ambiguity affecting gene 25788 (RAD54B).

## 4. What the Papers Got Wrong

| Claim in Paper | Reality |
|---|---|
| Features include "immunological signatures" (C7) | Features are C1 + C3 only; no C7 |
| Labels from "gene ontology sets collected from MSigDB" | Labels are from GOA/GO directly; 0/18,026 MSigDB sets match |
| Labels are Biological Process only (implied) | Labels span BP (85), CC (26), and MF (10) |

## 5. Key Biological Insights

### Annotation Instability
GO annotations are not stable over time. Between March and September 2016:
- 7,443 TAS (Traceable Author Statement) annotations were removed from GOA
- 616 IDA (experimental) annotations were removed
- Evidence codes were reclassified between releases
- The same gene-GO association can have different evidence codes in GOA vs gene2go

This instability means the exact label reproduction depends on the precise annotation snapshot date, which we have narrowed to July 2016 (GOA v159) but cannot pin to a single day.

### UniProt Demerges
Some genes (calmodulins CALM1/2/3, histone H4 variants) share a single UniProt accession that was later demerged into separate entries (e.g., P62158 → P0DP23/P0DP24/P0DP25 in May 2017). The many-to-many ID mapping handles this correctly for the 2016 timeframe.

## 6. Open Questions

1. **Exact gene2go snapshot**: A mid-2016 gene2go.gz would likely achieve 121/121 exact. No public archive contains this file. The closest available snapshots are May 2, 2016 (from dhimmel/gene-ontology git history) and December 23, 2016 (Wayback Machine).

2. **WL equivalence class resolution**: ~12 genes cannot be disambiguated by graph structure or features alone. Resolving them requires either the GO labels as tiebreaker (valid given verified pipeline) or a different graph matching approach.

3. **Graph selection criteria**: Why these specific 24 of 144 OhmNet tissues? Why do 3 tissues (adipose, heart, lung) include small connected components while the other 21 use the largest connected component only?

4. **Reproducibility script**: No preprocessing code exists in the OhmNet or GraphSAGE repositories. The complete reconstruction pipeline should be published as a reproducible script.

## 7. Files and Reproducibility

All source data, generated artifacts, and analysis code are inventoried in `artifact_inventory.tsv` with SHA-256 hashes. Key artifacts:

- **`node2gene_corrected.json`**: Node-to-Entrez mapping (4,278 genes, 19 swap corrections)
- **`matches_full.json`**: Graph-to-tissue mapping (24/24)
- **`2016-05-02-gene2go_all_organisms.gz`**: Recovered May 2016 gene2go from dhimmel's git history
- **`artifact_inventory.tsv`**: Complete file inventory with provenance metadata

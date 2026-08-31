# Updated investigator handoff — full GraphSAGE row and GO-label resolution

## Correction to the earlier result

The first 121/121 statement covered 56,411 rows with independently recovered GeneIDs. The remaining 533 rows were validated only as label-vector multisets within 183 topology/feature-equivalence classes.

A new, label-independent row-order reconstruction now assigns all 56,944 rows.

## Complete row-order reconstruction

For each of the 24 matched OhmNet tissue edgelists:

1. scan edges in original line order;
2. keep endpoint tokens as strings;
3. insert them into the historical 64-bit unrandomized CPython 2 dictionary used by a NetworkX 1.x graph;
4. iterate occupied dictionary slots from low to high;
5. use this as the tissue-local node order.

This agrees with all 56,411 prior independent row identities and disagrees with none. It resolves all 533 former ambiguous rows and all 183 equivalence classes without using labels.

Simple alternatives fail: first-occurrence order matches 24 rows, sorted GeneID 30, reverse sorted 17, and a Python 2 integer-key dictionary 23; the Python 2 string-key dictionary matches 56,411/56,411.

Claim strength: **strongly supported complete reconstruction**, not source-code proof. The original preprocessing code has not been found.

## Complete external validation

Under the full mapping:

- 24/24 tissue node sets exact;
- 24/24 tissue edge sets exact;
- 818,716/818,716 GraphSAGE links verified;
- 56,944 × 50 feature matrix exact, 0/2,847,200 differences;
- 56,944 × 121 label matrix exact, 0/6,890,224 differences;
- independent implementations reproduce the row map, feature matrix, and label matrix.

## Fixed GO-label transformation

Use GOA human release 159, May-2016 `gp2protein.geneid`, GPI159 accession/primary-symbol metadata, evidence `EXP, IDA, IEP, IGI, IMP, ISS`, exclude `NOT`, retain ordinary relations `involved_in`, `part_of`, and `enables`, exclude `colocalizes_with` and `contributes_to`, canonicalize alternate GO IDs, and propagate only through `is_a`.

Preserve many-to-many identifier mappings and inspect complete components before restricting to graph genes. Do not transfer O95073/FSBP annotations to GeneID 25788/RAD54B. No mapping, evidence, or relation policy is tuned per gene or per column.

ATP6AP2/GeneID 10159 remains absent from the GPI projection, but its complete GraphSAGE label vector is zero, so this does not hide a mismatch.

## Remaining GO ambiguity

All 121 columns are exact, but only 115 identities are unique from membership. The following pairs remain identical on all 4,301 graph genes:

- GO:0043228 versus GO:0043232;
- GO:0006464 versus GO:0036211;
- GO:0043230 versus GO:1903561.

A Python 2 dictionary-order analysis strongly supports a provisional orientation for columns 24/71, 39/63, and 48/70, but the best plausible order has LCS 94/121 rather than a perfect sequence. Treat those six identities as strongly supported, not proven.

## Dhimmel Entrez-native controls

The four genuine commit-pinned human tables do not directly reproduce any complete GraphSAGE label column. The closest, inferred experimental evidence, reaches 74 columns at ≥95% and six at ≥99%, but zero exact columns. Its documented policy includes IPI and propagates through both `is_a` and `part_of`, unlike the exact GraphSAGE transformation.

## Other unresolved construction details

- exact original feature-selection source code and cutoff operator;
- MSigDB release, because v5.0–v6.0 produce the same 50 feature vectors;
- whether GO term selection was top 121 or ≥1,000, and whether genes or proteins were counted;
- exact label-column order and three duplicate-pair orientations;
- why the 24 tissues and 20/2/2 split were selected;
- whether the original GO membership source was GOA v159 itself or an equivalent Entrez-native/preprocessed product.

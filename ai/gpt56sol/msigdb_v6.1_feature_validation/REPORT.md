# MSigDB v6.1 acceptance test for GraphSAGE PPI feature reconstruction

## Conclusion

**MSigDB v6.1 C1 and C3 are suitable as the canonical source files for the central reproduction package.** Applying the recovered source-order selection rule to the user-supplied v6.1 Entrez GMT files selects the same ordered 50 membership vectors as v6.0 and reconstructs every one of the 2,847,200 cells in the released GraphSAGE PPI feature matrix exactly.

The selected features remain 30 C1 positional sets followed by 20 C3 regulatory-target sets. The corrected v6.1 names should therefore be used as the canonical display names, while earlier names should be retained as provenance aliases rather than discarded.

## Inputs and provenance

| Input | Size (bytes) | SHA-256 |
|---|---:|---|
| `c1.all.v6.1.entrez.gmt` | 208,819 | `c87eddfc173c34420cd082b0b04a8edc65d6cf79ec9d86388426a6bdd3a2de02` |
| `c3.all.v6.1.entrez.gmt` | 1,171,503 | `8f9d8218a0eaf9476ba009bb1a97f6161bc64a73152ce8de07e739eb7e91d30e` |
| `graphsage_ppi.zip` | 27,029,260 | `53aeb76e54fd41b645e7edb48b62929240b89839495396b048086fd212503fbd` |
| Complete 56,944-row GraphSAGE-to-Entrez map | 651,211 | `a21a6532e287bfe82102cdc9d6642a0c9aff77f9c7449997cbafb597d1145737` |
| Prior v6.0 selected-feature details | 53,484 | `634cd1de253aa087a1b69c8ff71f95461e7fef8c6fadab5939cea6e3e3dfc7fb` |

The prior v6.0 details were extracted from the repository snapshot `GraphSagePPIAnalysis-ai-31590e4.zip`, SHA-256 `3339dcd09332c148e5df92bd4ff7caa8e3a6990459c4d5c5dec75ea8bd9ff88b`, corresponding to commit `31590e4d2ab5da3b74a7ac74b585df6c67bca478`.

## Test performed

The primary validator independently parsed the two v6.1 GMT files in deposited row order. It then:

1. counted unique Entrez GeneIDs in each complete source gene set;
2. processed C1 before C3;
3. retained qualifying sets in source order;
4. stopped after a global total of 50 sets;
5. projected each selected membership set onto the complete reconstructed sequence of 56,944 GraphSAGE rows;
6. compared the reconstructed matrix with `ppi/ppi-feats.npy` from the independently supplied GraphSAGE archive;
7. compared every selected membership hash with the previously retained v6.0 selected-feature record;
8. repeated the row-level check using an independently implemented 50-bit signature validator.

No released feature values or feature names were used to choose individual source sets. The only selection policy was the previously recovered global rule.

## Source-file and selection results

| Property | C1 v6.1 | C3 v6.1 |
|---|---:|---:|
| GMT rows | 326 | 836 |
| Sets with at least 200 unique Entrez members | 30 | 487 |
| Sets with more than 200 unique Entrez members | 30 | 484 |
| Sets selected before the global cap | 30 | 20 |

The selected sequence under `>= 200` is identical to the selected sequence under `> 200`. Three C3 sets later in the source file have exactly 200 members—`CEBP_C`, `HFH4_01`, and `TGTGTGA_MIR377`—but all occur after the 50-column cap has already been reached. The released matrix therefore still cannot distinguish which of the two threshold operators was historically used.

## Exact comparison with v6.0 and GraphSAGE

| Check | Result |
|---|---:|
| Selected C1/C3 composition | 30 / 20 |
| v6.1 selected membership vectors matching v6.0 by column | 50 of 50 |
| Ordered membership-sequence SHA-256 | `41c5d821e1b706ec4c8dceb47ab25c5dbad689998483e34e9e41d08094448101` |
| GraphSAGE feature matrix shape | 56,944 x 50 |
| Feature cells compared | 2,847,200 |
| Mismatching cells | **0** |
| Exact feature columns | 50 of 50 |
| Distinct mapped Entrez genes checked | 4,301 |
| Conflicts among repeated occurrences of the same gene | 0 |
| Reconstructed and observed raw matrix SHA-256 | `dc2420e5fe912f9c14292a6a14ffc0f5ef0ef781d830a0aa3e05085ba53d365a` |

The previously observed all-zero feature remains column 10 in zero-based indexing (feature 11 in one-based indexing), corresponding to `chryq11`.

## Name correction relevant to the selected 50 features

Exactly one of the selected 50 source names differs between v6.0 and v6.1:

| GraphSAGE column | Collection | v6.0 name | v6.1 name | Membership status |
|---:|---|---|---|---|
| 35 zero-based / 36 one-based | C3 | `AAAYWAACM_HFH4_01` | `AACTTT_UNKNOWN` | Identical 1,890-member set |

Its membership SHA-256 is:

```text
31d4723cefaf528436136ef4ad8cdafb940e267a6fa4b6c8d6107396634fb404
```

All other selected v6.0 and v6.1 display names agree, and all 50 memberships agree regardless of name.

## Independent quality-control check

A second implementation represented the expected 50 features for each Entrez GeneID as one integer bit mask. It then converted every observed GraphSAGE row into the same representation and compared all 56,944 row signatures. This check found:

- 4,301 distinct genes checked;
- 56,944 row signatures checked;
- zero row-signature mismatches;
- exact agreement of all selected v6.1 and retained v6.0 membership hashes;
- exact agreement between the `>= 200` and `> 200` selected sequences.

This validator did not consume the matrix generated by the primary implementation.

## Recommendation for the reproduction package

Use the complete official v6.1 Entrez GMT files as the upstream inputs:

```text
c1.all.v6.1.entrez.gmt
c3.all.v6.1.entrez.gmt
```

The production feature specification should record, for each selected column:

- zero- and one-based GraphSAGE column numbers;
- collection;
- v6.1 source row number;
- canonical v6.1 name;
- prior/historical alias where different;
- unique source-member count;
- membership SHA-256;
- expected projected positive count.

The documentation should say that v6.1 is a **corrected canonical reproduction source**, not that it has been identified as the original historical GraphSAGE source. The feature matrix cannot distinguish among the already tested versions 5.0, 5.1, 5.2, 6.0, and now 6.1 because all yield the same ordered memberships under the recovered rule.

The main workflow should use only v6.1. The multi-version comparison remains a one-time forensic acceptance analysis outside the central reproduction workflow.

## Files in this evidence package

- `validation_summary.json`: primary machine-readable result.
- `independent_validation.json`: independent bit-mask validation.
- `selected_features_v6.1.csv`: all selected source rows, names, hashes, counts, and cell comparisons.
- `v6.0_to_v6.1_selected_name_changes.csv`: selected names changed by v6.1.
- `prior_v6.0_selected_features_used.json`: compact retained v6.0 comparator and its provenance.
- `source_manifest.csv`: exact input sizes and hashes.
- `scripts/`: the two validation implementations.
- `SHA256SUMS`: hashes of every retained evidence file.

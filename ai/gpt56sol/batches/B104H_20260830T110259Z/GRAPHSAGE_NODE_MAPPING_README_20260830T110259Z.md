# GraphSAGE PPI node-to-Entrez mapping

## Files

- `graphsage_ppi_node_to_entrez_20260830T110259Z.tsv.gz` is the compact mapping intended for ordinary use.
- `graphsage_ppi_node_to_entrez_evidence_20260830T110259Z.tsv.gz` adds row-level topology, feature, label, source-file, and legacy-Python ordering evidence.
- `B104H_MAPPING_VALIDATION_20260830T110259Z.json` records full-table counts, source hashes, matrix hashes, and validation totals.

## Compact mapping columns

| Column | Meaning |
|---|---|
| `graphsage_node_id` | Node ID in `ppi-G.json`; also the key in `ppi-id_map.json`. |
| `feature_label_row_index` | Row in `ppi-feats.npy` and the corresponding 121-vector in `ppi-class_map.json`. In this dataset it equals `graphsage_node_id`. |
| `graph_index_1based` | One of the 24 disjoint tissue graphs, in deposited order. |
| `tissue` | Recovered OhmNet tissue name. |
| `split` | `train`, `validation`, or `test`. |
| `local_node_index_0based` | Position within that tissue graph. |
| `entrez_gene_id` | Recovered biological node identifier from the original OhmNet edgelist. |

## Evidence tiers

The tiers describe how the row identity was established before the final global cross-checks.

- **A — topology-only independent (55,878 rows):** graph-structural refinement uniquely identified the row.
- **B — topology plus feature independent (533 rows):** topology left a small equivalence class and independently reconstructed MSigDB features resolved it.
- **C — legacy-Python source-order inference (533 rows):** topology and features did not distinguish the class. The row was assigned by replaying original OhmNet edgelist insertion with a 64-bit, unrandomized CPython 2 string-key dictionary. This mechanism agrees with every one of the 56,411 A/B assignments and then yields exact topology, feature, and GO-label reconstruction on the remaining rows.

Tier C is a strongly supported historical-mechanism inference rather than direct documentary proof from the missing preprocessing script.

## Global validation

The complete mapping has:

- 56,944 GraphSAGE rows;
- 4,301 distinct Entrez Gene IDs;
- 24 tissue graphs;
- 818,716 undirected deposited edges;
- exact per-row OhmNet neighbor sets for all 56,944 rows;
- exact MSigDB-derived feature vectors for all 56,944 rows and all 2,847,200 feature cells;
- exact fixed-GOA-derived label vectors for all 56,944 rows and all 6,890,224 label cells.

The feature and label reconstructions were performed globally. They were not tuned per row or per gene.

## Appropriate claim

A careful statement is:

> A deterministic row-to-Entrez mapping has been reconstructed for all 56,944 GraphSAGE PPI nodes. It is independently confirmed on 56,411 rows and exactly reproduces every deposited edge, feature value, and label value. The remaining 533 row identities follow from an inferred legacy CPython 2/NetworkX node-order mechanism that is perfectly consistent with all independently resolved rows and all downstream data, but the original preprocessing script has not been recovered.

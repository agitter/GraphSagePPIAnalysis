# PPI Benchmark: Complete Gene Identity Recovery and Leakage Measurement

## All 24 graphs identified

Using the colleague's insight — small connected components present in OhmNet were removed
during GraphSAGE preprocessing — all 24 graphs now match an OhmNet tissue, with degree
sequences verified on the largest connected component:

| Graph | Split | Tissue | Nodes |
|---|---|---|---|
| g1 | train | adipose_tissue | 1,767 |
| g2 | train | adrenal_cortex | 1,377 |
| g3 | train | adrenal_gland | 2,263 |
| g4 | train | amygdala | 2,339 |
| g5 | train | aorta | 1,578 |
| g6 | train | astrocyte | 1,021 |
| g7 | train | artery | 1,823 |
| g8 | train | basal_ganglion | 2,488 |
| g9 | train | basophil | 591 |
| g10 | train | blood | 3,312 |
| g11 | train | blood_plasma | 2,401 |
| g12 | train | blood_platelet | 1,878 |
| g13 | train | bone | 1,819 |
| g14 | train | brain | 3,480 |
| g15 | train | colon | 2,794 |
| g16 | train | eye | 2,326 |
| g17 | train | forebrain | 2,650 |
| g18 | train | large_intestine | 2,815 |
| g19 | train | liver | 3,163 |
| g20 | train | gastrointestinal_tract | 3,021 |
| g21 | val | heart | 3,230 |
| g22 | val | kidney | 3,284 |
| g23 | **test** | **lung** | 3,224 |
| g24 | **test** | **midbrain** | 2,300 |

## 99.12% of nodes have verified Entrez gene IDs

Three methods, applied in order:

| Method | Nodes resolved |
|---|---|
| WL colour refinement (structural fingerprints) | 55,874 |
| Feature-vector disambiguation (same WL colour, different features) | 261 |
| Process of elimination (last remaining candidate) | 24 |
| **Total mapped** | **56,159 / 56,658 LCC nodes = 99.12%** |
| Unmapped (WL-symmetric, same features) | 499 |

**816,527 mapped edges verified against OhmNet at 100.0000%.**

## The leakage

The gene universe across all OhmNet tissues is 4,510 proteins. Training maps to 4,136
distinct genes. Both test graphs draw from this same pool:

| Test graph | Mapped genes | In training | With identical labels | Not in training |
|---|---|---|---|---|
| lung (g23) | 3,167 | 3,132 (98.89%) | 3,132 (98.89%) | 35 (1.11%) |
| midbrain (g24) | 2,276 | 2,276 (100.00%) | 2,276 (100.00%) | 0 (0.00%) |

Every test gene that appears in training carries a byte-identical 121-label vector.

## Gene lookup baseline on the full test set

| | |
|---|---|
| Test nodes | 5,524 |
| Covered by gene lookup | 5,408 (97.90%) |
| **micro-F1 on full test set** | **0.9932** |
| micro-F1 on covered nodes | 1.0000 |
| Constant baseline (full test) | 0.3935 |
| Published Cluster-GCN | 0.9936 |

A zero-parameter lookup table — find the same gene in training, copy its labels —
scores **0.9932** on the full test set. This is within 0.0004 of Cluster-GCN's 0.9936.

The 116 uncovered test nodes (2.1%) fall into two groups: 35 lung genes not seen in
training, and 81 nodes unmapped due to WL symmetry. The fallback for those is the constant
predictor. With perfect gene identification, the score would be higher.

## What this means

The PPI benchmark does not measure generalization. A model that memorizes which gene each
node represents — using neighbourhood structure as a proxy for the stripped Entrez IDs — and
recalls that gene's label vector from training achieves state-of-the-art performance. No
learning about protein function, graph structure, or tissue biology is required.

The split partitions graphs (tissues), not genes. Since the same ~4,500 genes appear across
all tissues with identical features and identical labels, this is equivalent to training and
testing on the same data with shuffled indices.

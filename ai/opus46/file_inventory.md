# File Inventory for GraphSAGE PPI Benchmark Investigation

## Goal
Identify exactly how the GraphSAGE PPI dataset (features, labels, graph structure)
was generated, to understand the label leakage that makes the benchmark trivial.

## PRIMARY DATA (the thing we're investigating)

### GraphSAGE PPI dataset
- `/mnt/user-data/uploads/graphsage_ppi.zip` — canonical dataset from snap.stanford.edu/graphsage/ppi.zip
  Contains: ppi-G.json, ppi-feats.npy, ppi-class_map.json, ppi-id_map.json, ppi-walks.txt
  56,944 nodes, 818,716 edges, 50 binary features, 121 binary labels

### DGL PPI dataset
- `/mnt/user-data/uploads/dgl_ppi.zip` — pre-processed version from data.dgl.ai/dataset/ppi.zip
  Contains: {train,valid,test}_{graph.json,feats.npy,labels.npy,graph_id.npy}
  Same data but: z-scored features, directed edges, self-loops on all nodes, split into files
  **STATUS: FULLY REPRODUCED** — transformation from GraphSAGE to DGL confirmed byte-identical

## GRAPH TOPOLOGY SOURCES

### OhmNet tissue networks
- `/home/claude/ohm/bio-tissue-networks/` — 144 tissue-specific PPI edgelists from SNAP
  Binary edges with Entrez Gene IDs. Source: BioGRID physical interactions.
  **STATUS: All 24 GraphSAGE graphs matched to specific OhmNet tissues**

### OhmNet tissue labels
- `/home/claude/ohm/bio-tissue-labels/` — 503 tissue-specific GO annotation files
  **STATUS: NOT the source of GraphSAGE labels (0/121 match)**

### OhmNet code
- `/mnt/user-data/uploads/ohmnet-master.zip` — embedding algorithm only, no data preprocessing
  **STATUS: No label/feature generation code found**

### GIANT/HumanBase gold standards
- `/mnt/user-data/uploads/HumanBase-blood.dat` — C1-C4 classified gene pairs (4,730 genes)
- `/mnt/user-data/uploads/HuamnBase-kidney.dat` — same format, same 4,730-gene universe
  **STATUS: Confirmed OhmNet 4,510 genes ⊂ GIANT 4,730 genes. C1 co-annotation
  structure is consistent with labels but not independently informative.**

### GIANT network sample
- `/mnt/user-data/uploads/blood_sample_tsv.gz` — sampled top/bottom edges from GIANT blood
  25,689 genes, 34.6M weighted functional edges — completely different from OhmNet's
  3,326-gene, 54K-edge binary PPI network.
  **STATUS: GIANT functional networks are NOT the source of OhmNet graph topology**

## FEATURE PROVENANCE SOURCES

### MSigDB gene set collections
- `/mnt/user-data/uploads/msigdb_v5_2_files_to_download_locally.zip` — MSigDB v5.2 (Sep 2016)
- `/mnt/user-data/uploads/msigdb_v5_2_chip_files_to_download_locally.zip` — chip mappings
- `/mnt/user-data/uploads/msigdb_v5_1_files_to_download_locally.zip` — MSigDB v5.1
- `/mnt/user-data/uploads/msigdb_v6_0_files_to_download_locally.zip` — MSigDB v6.0
  **STATUS: 50/50 features CONFIRMED from MSigDB v5.2**
  - Columns 0-29: C1 (positional gene sets, 30 chromosomal bands)
  - Columns 30-49: C3 (motif gene sets, 20 TF binding motifs)
  - NO C7 features despite paper claiming "immunological signatures"
  - Column 10: all-zero (no gene in universe is a member)

## LABEL PROVENANCE SOURCES

### GOA (Gene Ontology Annotation) files from EBI
- `/mnt/user-data/uploads/gene_association_goa_human_155.gz` — March 2016
- `/mnt/user-data/uploads/gene_association_goa_human_156.gz` — April 2016
- `/mnt/user-data/uploads/gene_association_goa_human_157.gz` — May 2016
- `/mnt/user-data/uploads/goa_human_gaf_158.gz` — June 2016
- `/mnt/user-data/uploads/goa_human_gaf_159.gz` — July 2016 (BEST MATCH)
- `/mnt/user-data/uploads/goa_human_gaf_160.gz` — September 2016
  **STATUS: Best result 77/121 columns at ≥99% using v159**

### GOA support files
- `/mnt/user-data/uploads/2016-06-01-gp2protein_human.gz` — UniProt-to-UniProt (not useful)
- `/mnt/user-data/uploads/2016-06-01-gp2protein_geneid.gz` — Entrez-to-UniProt mapping
  **STATUS: Used to bridge GOA UniProt accessions to Entrez IDs**

### GO ontology structure
- `/mnt/user-data/uploads/2016-06-01-go.obo` — GO hierarchy for is_a propagation
- `/mnt/user-data/uploads/2016-06-01-go.owl` — OWL format (not used)
- `/mnt/user-data/uploads/2016-06-01-gene_ontology.obo` — alternative format
  **STATUS: Used for true-path-rule propagation (is_a only; part_of makes it worse)**

### NCBI gene2go
- `/mnt/user-data/uploads/2016-12-23-gene2go_human_tsv.gz` — Dec 2016 Wayback Machine snapshot
- `/mnt/user-data/uploads/2026-08-14-gene2go_human_tsv.gz` — current version
  **STATUS: Dec 2016 covers 4,504/4,510 OhmNet genes (vs GOA's 4,481) but
  is 4 months too late. Current version is 10 years too late. No mid-2016 snapshot exists.**

### Bioconductor org.Hs.eg.db (frozen gene2go snapshots)
- `/mnt/user-data/uploads/bioconductor-annotation-org_Hs_eg_db_3_0_0_tar.gz` — Oct 2014
- `/mnt/user-data/uploads/bioconductor-annotation-org_Hs_eg_db_3_1_2_tar.gz` — Apr 2015
- `/mnt/user-data/uploads/bioconductor-annotation-org_Hs_eg_db_3_3_0_tar.gz` — May 2016
- `/mnt/user-data/uploads/bioconductor-annotation-org_Hs_eg_db_3_4_0_tar.gz` — Oct 2016
  **STATUS: All worse than direct GOA matching. Earlier dates are worse, not better.**

### Greene 2015 supplementary files
- `/mnt/user-data/uploads/Greene2015_Table6.xlsx` — 973 expert-curated GO BP terms
- `/mnt/user-data/uploads/Greene2015_Table9.xlsx` — GO-to-tissue mappings
- `/mnt/user-data/uploads/Greene2015_sup.pdf` — supplementary figures
  **STATUS: The 973 Greene terms filtered to ≥15 annotated genes gives 117-118,
  not 121. The 121 labels are NOT restricted to Greene's curated set.**

### Papers
- `/mnt/user-data/uploads/OhmNet.pdf`
- `/mnt/user-data/uploads/Greene2015.pdf`

## GENERATED ARTIFACTS

### Gene identity recovery
- `/mnt/user-data/outputs/node2gene_corrected.json` — node-to-Entrez mapping for all 24 graphs
  4,278 genes mapped. 19 swaps corrected via MSigDB feature matching.
  **STATUS: 50/50 features match MSigDB v5.2 exactly after corrections**

### Reports
- `/mnt/user-data/outputs/ppi_lookup_table_report.md` — leakage measurement
- `/mnt/user-data/outputs/ppi_gene_identity_report.md` — gene recovery methodology
- `/mnt/user-data/outputs/ppi_final_report.md` — earlier comprehensive report

### Scripts
- `/mnt/user-data/outputs/sample_giant_network.sh` — local script to sample GIANT networks

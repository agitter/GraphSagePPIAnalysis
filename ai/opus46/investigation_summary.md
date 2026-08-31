# GraphSAGE PPI Benchmark: Data Provenance Investigation
## Snapshot Summary

### Motivation

The GraphSAGE PPI benchmark achieves ~99.5% micro-F1 with a zero-parameter
gene lookup table because 98.8% of test genes appear in training with
identical labels. This investigation traces how the dataset was constructed
to understand why this leakage exists.

The dataset has three components: graph topology (edges), node features
(50 binary columns), and node labels (121 binary columns). We trace
each independently.

---

## Input Files

### The dataset under investigation

| File | Source | Contents |
|---|---|---|
| `graphsage_ppi.zip` | `snap.stanford.edu/graphsage/ppi.zip` | 56,944 nodes, 818,716 edges, 50 features, 121 labels, 24 graphs |
| `dgl_ppi.zip` | `data.dgl.ai/dataset/ppi.zip` | Same data, transformed (z-scored features, directed edges, self-loops) |

### OhmNet (intermediate source)

| File | Source | Contents |
|---|---|---|
| `bio-tissue-networks/*.edgelist` | `snap.stanford.edu/ohmnet/` | 144 tissue-specific PPI networks, Entrez Gene IDs, from BioGRID |
| `bio-tissue-labels/*.lab` | `snap.stanford.edu/ohmnet/` | 503 tissue-specific GO annotation files |
| `ohmnet-master.zip` | `github.com/snap-stanford/ohmnet` | Embedding algorithm only; no data preprocessing code |
| `OhmNet.pdf` | Zitnik & Leskovec, Bioinformatics 2017 | Describes data sources; cites Greene 2015 for annotations |

### GIANT/HumanBase (upstream source for OhmNet)

| File | Source | Contents |
|---|---|---|
| `HumanBase-blood.dat` | `hb.flatironinstitute.org` | Gold standard: 4,730 gene pairs classified C1-C4 |
| `HumanBase-kidney.dat` | `hb.flatironinstitute.org` | Same format, same 4,730-gene universe |
| `blood_sample_tsv.gz` | Sampled from `HumanBase-blood_top.gz` | GIANT functional network: 25,689 genes, 34.6M weighted edges |
| `Greene2015.pdf` | Greene et al., Nature Genetics 2015 | Describes GIANT construction; defines gold standard methodology |
| `Greene2015_Table6.xlsx` | Supplementary Table 6 | 973 expert-curated GO BP terms used for gold standard |
| `Greene2015_Table9.xlsx` | Supplementary Table 9 | GO-to-tissue mappings |
| HumanBase download script | `hb.flatironinstitute.org/download` | Python script showing API and URL structure |
| HumanBase data sources page | `hb.flatironinstitute.org` | Lists annotation versions; GIANT uses UniProt-GOA for GO BP |

### Feature provenance (MSigDB)

| File | Source | Contents |
|---|---|---|
| `msigdb_v5.2_files_to_download_locally.zip` | MSigDB (Broad Institute) | Gene set collections; C1 (positional) and C3 (motif) matched features |
| `msigdb_v5.1_files_to_download_locally.zip` | MSigDB | Earlier version for comparison |
| `msigdb_v6.0_files_to_download_locally.zip` | MSigDB | Later version for comparison |

### Label provenance candidates

#### GOA (Gene Ontology Annotation, from EBI)

| File | Source | Date | Purpose |
|---|---|---|---|
| `gene_association.goa_human.155.gz` | EBI GOA archive | Mar 2016 | Monthly GOA snapshot |
| `gene_association.goa_human.156.gz` | EBI GOA archive | Apr 2016 | Monthly GOA snapshot |
| `gene_association.goa_human.157.gz` | EBI GOA archive | May 2016 | Monthly GOA snapshot |
| `goa_human.gaf.158.gz` | EBI GOA archive | Jun 2016 | Monthly GOA snapshot (filename format changed) |
| `goa_human.gaf.159.gz` | EBI GOA archive | Jul 2016 | Monthly GOA snapshot — best match |
| `goa_human.gaf.160.gz` | EBI GOA archive | Sep 2016 | Monthly GOA snapshot |

#### gene2go (NCBI Gene)

| File | Source | Date | Purpose |
|---|---|---|---|
| `gene2go_may2016.gz` | dhimmel/gene-ontology git history (commit 962a5e1) | May 2, 2016 | Raw NCBI gene2go, all organisms |
| `2016-12-23-gene2go_human_tsv.gz` | Wayback Machine snapshot | Dec 23, 2016 | Human-filtered gene2go |
| `2026-08-14-gene2go_human_tsv.gz` | Current NCBI FTP | Aug 2026 | Too recent; 10 years of annotation changes |

#### GO ontology structure

| File | Source | Date | Purpose |
|---|---|---|---|
| `2016-06-01-go.obo` | GO release archive | Jun 1, 2016 | is_a relationships for true-path-rule propagation |
| `2016-06-01-go.owl` | GO release archive | Jun 1, 2016 | OWL format (not used) |
| `2016-06-01-gene_ontology.obo` | GO release archive | Jun 1, 2016 | Alternative format |

#### ID mapping files

| File | Source | Date | Purpose |
|---|---|---|---|
| `2016-06-01-gp2protein_geneid.gz` | GO release archive | Jun 1, 2016 | Entrez GeneID ↔ UniProt accession mapping |
| `2016-06-01-gp2protein_human.gz` | GO release archive | Jun 1, 2016 | UniProt ↔ UniProt (not useful) |

#### Bioconductor annotation packages

| File | Source | Date | Purpose |
|---|---|---|---|
| `org.Hs.eg.db_3.0.0.tar.gz` | Bioconductor 3.0 | Oct 2014 | Bundled gene2go snapshot |
| `org.Hs.eg.db_3.1.2.tar.gz` | Bioconductor 3.1 | Apr 2015 | Bundled gene2go snapshot |
| `org.Hs.eg.db_3.3.0.tar.gz` | Bioconductor 3.3 | May 2016 | Bundled gene2go snapshot |
| `org.Hs.eg.db_3.4.0.tar.gz` | Bioconductor 3.4 | Oct 2016 | Bundled gene2go snapshot |

---

## Generated Artifacts

| File | Contents |
|---|---|
| `node2gene_corrected.json` | Node-to-Entrez mapping for all 24 graphs (4,278 genes mapped, 19 swaps corrected) |
| `label_matching_detail.txt` | Per-column best GO term match with FP/FN gene lists |
| `file_inventory.md` | Earlier version of this inventory |
| `sample_giant_network.sh` | Script to sample GIANT top-edge networks locally |
| `find_gene2go_2016.py` | Script to search GitHub for archived gene2go copies |

---

## RESOLVED QUESTIONS

### 1. Gene identity recovery

**Result:** 4,278 of ~4,510 genes identified by Entrez ID.

**Method:** Weisfeiler-Leman colour refinement matches GraphSAGE nodes to
OhmNet edgelist genes by graph structure. 98.6% of nodes get unique colours.
The remaining ~600 nodes in WL equivalence classes are disambiguated using
MSigDB v5.2 feature membership as ground truth.

**Verification:** 818,435 of 818,716 edges (99.97%) reconstructed from OhmNet
edgelists. 50/50 feature columns match MSigDB v5.2 exactly after 19 swap
corrections across 6 gene pairs.

**Confidence:** High. The edge reconstruction rate and feature match are
independent verification. The 232 unmapped genes are in non-LCC components
or truly indistinguishable WL equivalence classes (e.g., histone H4 paralogs
encoding identical proteins).

### 2. Graph topology provenance

**Result:** All 24 GraphSAGE graphs are OhmNet tissue-specific PPI networks
derived from BioGRID physical interactions.

**Verification:** 24/24 graphs matched to named OhmNet tissues by edge set
comparison. Edge verification: 818,435 matched edges at 100.0000%.

**What GIANT is NOT:** GIANT tissue networks are dense weighted functional
networks (25,689 genes, 34.6M edges for blood). OhmNet graphs are sparse
binary PPI networks (3,326 genes, 54K edges for blood). These are completely
different data. OhmNet used BioGRID PPI edges, not GIANT functional edges.

**What GIANT IS:** The 4,730-gene universe in GIANT's gold standard is a
strict superset of OhmNet's 4,510 genes. The GO annotation methodology
from Greene 2015 likely influenced how labels were generated, even though
the graph topology came from a different source.

**Confidence:** High for topology source. The GIANT relationship to labels
is established directionally (OhmNet ⊂ GIANT) but the exact mechanism is
unresolved.

### 3. Feature provenance

**Result:** 50 features = MSigDB v5.2 C1 (30 chromosomal bands, columns 0-29)
+ C3 (20 TF binding motifs, columns 30-49).

**Verification:** 50/50 columns match MSigDB v5.2 gene sets exactly after
swap corrections. Column 10 is entirely zero (no gene in the universe is
a member of that chromosomal band).

**What the paper says vs reality:** The GraphSAGE paper claims features
include "immunological signatures" (C7). This is incorrect — no C7 features
are present. The features are C1 + C3 only.

**Confidence:** Complete. 50/50 exact match leaves no ambiguity.

### 4. DGL transformation

**Result:** The DGL ppi.zip was created from GraphSAGE ppi.zip by:
1. Reorder nodes by graph_id (stable sort)
2. StandardScaler fit on train (float64, ddof=0), transform all, cast to float32
3. Convert to directed (double edges)
4. Add self-loop to every node
5. Split into per-split files

**Verification:** Byte-identical labels and graph_ids. Features within float32
rounding (max diff 5.37e-07). Created by Hao Zhang (sufeidechabei) for
DGL v0.2, PR #395 (Feb 17, 2019).

**Confidence:** Complete.

### 5. Leakage mechanism

**Result:** The 24 tissue graphs share ~4,510 genes. Train/test splits
partition graphs (tissues), not genes. 98.8% of test genes appear in
training with byte-identical labels. A zero-parameter gene lookup table
scores 0.9956 micro-F1.

**Confidence:** Complete. This is a direct measurement, not an inference.

---

## PARTIALLY RESOLVED QUESTIONS

### 6. Label provenance

**What we know:**
- The 121 labels are GO Biological Process annotations. This is supported by:
  - 78/121 columns match specific GO BP terms at ≥99% agreement (GOA v159)
  - 85/121 at ≥95% across all tested sources
  - The matching GO terms are biologically coherent (e.g., column 0 = GO:0050789
    "regulation of biological process")
- The labels are NOT from MSigDB (0/18,026 gene sets match across v5.1, v5.2, v6.0,
  all collections). The GraphSAGE paper's claim of MSigDB is incorrect.
- The labels are NOT the OhmNet tissue-specific labels (0/121 match)
- Propagation uses is_a relationships only (part_of makes results worse)
- The label-generating annotations were in Entrez Gene ID space (covers genes
  that UniProt-based sources miss, like histone H4 paralogs and calmodulins)

**What we cannot resolve:**

The exact gene2go snapshot and evidence code filter. We tested:
- 6 GOA versions (v155-v160, March-September 2016)
- 2 gene2go snapshots (May 2016, December 2016)
- 4 Bioconductor annotation packages (3.0-3.4, October 2014-October 2016)
- 16 evidence code filter combinations
- Intersection and union of May/December gene2go snapshots

No combination produces an exact match. The best results:
- GOA v159 (July 2016) with many-to-many UniProt mapping: 78/121 at ≥99%
- gene2go May 2016 directly: 37/121 at ≥99%
- gene2go December 2016 directly: 68/121 at ≥99%

**Why we're stuck:**

1. **No single evidence code filter works.** 1,908 of 4,060 labeled genes require
   at least one "excluded" evidence code (TAS, IBA, or NAS) to explain at least
   one of their labels. Including those codes improves coverage but adds false
   positives for other genes. This could mean:
   - The labels used a simple filter (e.g., exclude only IEA) from a gene2go
     snapshot where the evidence code assignments were different from our
     available snapshots
   - The labels used a different annotation source entirely
   - There is some other processing step we haven't identified

2. **Annotation instability.** Between March and September 2016, GOA underwent
   significant changes: 7,443 TAS annotations removed, 616 IDA annotations
   removed, 2,940 IEA annotations removed. gene2go shows similar patterns
   (549 IDA removed, 533 TAS removed between May and December 2016). Evidence
   codes are not stable identifiers — the same gene-GO association can change
   evidence codes between releases.

3. **No archived gene2go from the right date.** The labels were likely generated
   from a gene2go snapshot between May and September 2016. The May 2016 snapshot
   (from dhimmel's git history) is too early; the December 2016 snapshot (from
   Wayback Machine) is too late. No public archive contains a July-August 2016
   gene2go. A GitHub-wide search found no other copies from this period.

4. **ID mapping gaps when using GOA.** GOA uses UniProt accessions; the labels use
   Entrez Gene IDs. Our mapping (via gp2protein.geneid.gz + MSigDB symbols)
   covers 4,277 of 4,278 genes, but some mappings are incomplete due to UniProt
   demerges (e.g., calmodulin P62158 → P0DP23/P0DP24/P0DP25 in May 2017) and
   missing gp2protein entries (e.g., APOA4, LPA, PPP1R15B).

**What would resolve this:**
- The exact gene2go.gz file used to generate the labels (a mid-2016 snapshot
  that no longer exists in public archives)
- The preprocessing script from Marinka Zitnik or the OhmNet pipeline
- A gene2go snapshot from July-August 2016 found in a private archive

### 7. Graph selection and split assignment

**What we know:**
- 24 of OhmNet's 144 tissues were selected
- 20 assigned to training, 2 to validation, 2 to test
- The OhmNet paper says: ≥15K edges for training eligibility, "4 large networks"
  for val/test (≥35K edges)

**What we haven't verified:**
- Whether the selection was random or deliberate
- Whether the exact split assignment is reproducible from the stated criteria
- Why 3 tissues (adipose, heart, lung) include small connected components
  while the other 21 use LCC only

---

## OPEN QUESTIONS

### 8. The 121 GO term selection criterion

The OhmNet paper says "functions with at least 15 annotated proteins."
Using GOA v159 with our best evidence filter and the full OhmNet 4,510-gene
universe, a threshold of ≥500 gives exactly 121 GO BP terms. With the Greene
973-term restriction and experimental-only evidence, threshold ≥15 gives
117-118 terms (close but not 121). We have not determined the exact selection
criterion.

### 9. Alternative annotation sources

We have not tested whether the labels could come from:
- A pre-built annotation file distributed with GIANT/HumanBase
- The Sleipnir library's internal GO processing pipeline
- A different GO annotation database (e.g., Reactome, KEGG)
- A processed output from the Greene 2015 pipeline that combines
  multiple annotation sources
- An annotation set derived from the GIANT gold standard C1 edges

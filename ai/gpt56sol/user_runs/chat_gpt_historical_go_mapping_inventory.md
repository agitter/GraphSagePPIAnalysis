# Historical human GO annotation and Entrez–UniProt mapping inventory

**Research window:** primarily 2015–2016.  **Prepared:** 2026-08-23.

This inventory separates exact archived files from reconstruction sources and unresolved leads. It is intended as a staging document for later byte-level assembly and comparison.

## Bottom line

1. **Most likely GO annotation files:** EBI GOA `gene_association.goa_human.*.gz` (GAF) or, for a GPAD/GPI workflow, the matched `gp_association.goa_human.*.gz` plus `gp_information.goa_human.*.gz`. The names changed at release 158 in June 2016 to `goa_human.gaf/.gpa/.gpi`.
2. **Most likely Entrez↔UniProt mapping source:** UniProt `idmapping_selected.tab.2015_03.gz`, NCBI `gene2accession.gz`, or a bundled Bioconductor `org.Hs.eg.db` snapshot. The Bioconductor packages are especially useful because they preserve both `org.Hs.egUNIPROT` and `org.Hs.egGO` in one dated artifact.
3. **`human.xrefs.*.gz` is a legacy clue, not a 2015/2016 series:** the archived sequence ends at `human.xrefs.99.gz` in June 2011.
4. **GO `annotations/gp2protein` is secondary for a human workflow:** the historical files were chiefly model-organism ID→UniProt bridges (FlyBase, MGI, RGD, ZFIN, EcoCyc, GeneDB, etc.); no human-wide `gp2protein.human` file was found.

## Confidence-ranked sources

| Rank | Source | What it preserves | Status |
|---:|---|---|---|
| 1 | EBI GOA old/HUMAN | Exact monthly human GAF, GPAD, and GPI-family files | Confirmed and directly downloadable |
| 2 | Bioconductor `org.Hs.eg.db` 3.1–3.4 | Entrez↔UniProt and Entrez↔GO mappings in the same package | Confirmed exact package archives |
| 3 | UniProt `idmapping_selected.tab.2015_03.gz` | UniProt accession, GeneID, GO IDs, taxon, and other IDs | Confirmed exact March 2015 snapshot |
| 4 | NCBI Gene DATA | `gene2go`, `gene2accession`, `gene2refseq`, RefSeq↔UniProt collaboration map | Current rolling files confirmed; exact 2015/16 snapshots not versioned there |
| 5 | Ensembl releases 79–87 | Exact MySQL dumps containing external DB xrefs and ontology xrefs | Confirmed exact release snapshots; requires joins |
| 6 | HGNC complete-set archive | Entrez and UniProt columns in periodic complete-set snapshots | Archive confirmed; exact 2015/16 object names still unresolved |
| 7 | GO release `gp2protein` | Model-organism identifiers to UniProt accessions | Confirmed archive class, but weak human match |

## 1. EBI GOA human archive: exact 2015–2016 candidates

### Naming transition

- Releases **140–157** use `gene_association.goa_human`, `gp_association.goa_human`, and `gp_information.goa_human`.
- Release **158** (June 2016) begins `goa_human.gaf`, `goa_human.gpa`, and `goa_human.gpi`.
- The `goa_ref_human` variants are the reference-proteome subset.
- `complex`, `isoform`, and `rna` are auxiliary subsets introduced under the new naming scheme.

### `gene_association.goa_human` — 18 files

| Date | File | Size |
|---|---|---:|
| 2015-01-05 | [gene_association.goa_human.140.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.140.gz) | 5.8M |
| 2015-02-02 | [gene_association.goa_human.141.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.141.gz) | 5.8M |
| 2015-03-02 | [gene_association.goa_human.142.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.142.gz) | 6.1M |
| 2015-03-30 | [gene_association.goa_human.143.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.143.gz) | 5.7M |
| 2015-04-27 | [gene_association.goa_human.144.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.144.gz) | 6.1M |
| 2015-05-26 | [gene_association.goa_human.145.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.145.gz) | 6.1M |
| 2015-06-22 | [gene_association.goa_human.146.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.146.gz) | 6.2M |
| 2015-07-20 | [gene_association.goa_human.147.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.147.gz) | 6.2M |
| 2015-09-14 | [gene_association.goa_human.148.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.148.gz) | 6.3M |
| 2015-10-12 | [gene_association.goa_human.149.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.149.gz) | 6.4M |
| 2015-11-09 | [gene_association.goa_human.150.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.150.gz) | 6.4M |
| 2015-12-07 | [gene_association.goa_human.151.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.151.gz) | 6.1M |
| 2016-01-04 | [gene_association.goa_human.152.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.152.gz) | 6.2M |
| 2016-01-20 | [gene_association.goa_human.153.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.153.gz) | 6.2M |
| 2016-02-15 | [gene_association.goa_human.154.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.154.gz) | 6.2M |
| 2016-03-14 | [gene_association.goa_human.155.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.155.gz) | 6.2M |
| 2016-04-11 | [gene_association.goa_human.156.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.156.gz) | 6.2M |
| 2016-05-09 | [gene_association.goa_human.157.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.157.gz) | 6.1M |

### `gene_association.goa_ref_human` — 18 files

| Date | File | Size |
|---|---|---:|
| 2015-01-05 | [gene_association.goa_ref_human.140.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.140.gz) | 4.0M |
| 2015-02-02 | [gene_association.goa_ref_human.141.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.141.gz) | 4.0M |
| 2015-03-02 | [gene_association.goa_ref_human.142.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.142.gz) | 4.2M |
| 2015-03-30 | [gene_association.goa_ref_human.143.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.143.gz) | 4.3M |
| 2015-04-27 | [gene_association.goa_ref_human.144.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.144.gz) | 4.3M |
| 2015-05-26 | [gene_association.goa_ref_human.145.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.145.gz) | 4.4M |
| 2015-06-22 | [gene_association.goa_ref_human.146.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.146.gz) | 4.4M |
| 2015-07-20 | [gene_association.goa_ref_human.147.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.147.gz) | 4.4M |
| 2015-09-14 | [gene_association.goa_ref_human.148.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.148.gz) | 4.5M |
| 2015-10-12 | [gene_association.goa_ref_human.149.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.149.gz) | 4.6M |
| 2015-11-09 | [gene_association.goa_ref_human.150.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.150.gz) | 4.6M |
| 2015-12-07 | [gene_association.goa_ref_human.151.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.151.gz) | 4.6M |
| 2016-01-04 | [gene_association.goa_ref_human.152.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.152.gz) | 4.7M |
| 2016-01-20 | [gene_association.goa_ref_human.153.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.153.gz) | 4.7M |
| 2016-02-15 | [gene_association.goa_ref_human.154.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.154.gz) | 4.7M |
| 2016-03-14 | [gene_association.goa_ref_human.155.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.155.gz) | 4.7M |
| 2016-04-11 | [gene_association.goa_ref_human.156.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.156.gz) | 4.7M |
| 2016-05-09 | [gene_association.goa_ref_human.157.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_ref_human.157.gz) | 4.6M |

### `gp_association.goa_human` — 18 files

| Date | File | Size |
|---|---|---:|
| 2015-01-05 | [gp_association.goa_human.140.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.140.gz) | 3.6M |
| 2015-02-02 | [gp_association.goa_human.141.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.141.gz) | 3.7M |
| 2015-03-02 | [gp_association.goa_human.142.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.142.gz) | 3.9M |
| 2015-03-30 | [gp_association.goa_human.143.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.143.gz) | 3.7M |
| 2015-04-27 | [gp_association.goa_human.144.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.144.gz) | 3.9M |
| 2015-05-26 | [gp_association.goa_human.145.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.145.gz) | 3.9M |
| 2015-06-22 | [gp_association.goa_human.146.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.146.gz) | 4.0M |
| 2015-07-20 | [gp_association.goa_human.147.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.147.gz) | 4.0M |
| 2015-09-14 | [gp_association.goa_human.148.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.148.gz) | 4.1M |
| 2015-10-12 | [gp_association.goa_human.149.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.149.gz) | 4.2M |
| 2015-11-09 | [gp_association.goa_human.150.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.150.gz) | 4.1M |
| 2015-12-07 | [gp_association.goa_human.151.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.151.gz) | 4.0M |
| 2016-01-04 | [gp_association.goa_human.152.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.152.gz) | 4.0M |
| 2016-01-20 | [gp_association.goa_human.153.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.153.gz) | 4.0M |
| 2016-02-15 | [gp_association.goa_human.154.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.154.gz) | 4.1M |
| 2016-03-14 | [gp_association.goa_human.155.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.155.gz) | 4.1M |
| 2016-04-11 | [gp_association.goa_human.156.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.156.gz) | 4.0M |
| 2016-05-09 | [gp_association.goa_human.157.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.157.gz) | 4.0M |

### `gp_association.goa_ref_human` — 18 files

| Date | File | Size |
|---|---|---:|
| 2015-01-05 | [gp_association.goa_ref_human.140.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.140.gz) | 2.8M |
| 2015-02-02 | [gp_association.goa_ref_human.141.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.141.gz) | 2.8M |
| 2015-03-02 | [gp_association.goa_ref_human.142.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.142.gz) | 3.0M |
| 2015-03-30 | [gp_association.goa_ref_human.143.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.143.gz) | 3.1M |
| 2015-04-27 | [gp_association.goa_ref_human.144.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.144.gz) | 3.1M |
| 2015-05-26 | [gp_association.goa_ref_human.145.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.145.gz) | 3.1M |
| 2015-06-22 | [gp_association.goa_ref_human.146.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.146.gz) | 3.1M |
| 2015-07-20 | [gp_association.goa_ref_human.147.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.147.gz) | 3.2M |
| 2015-09-14 | [gp_association.goa_ref_human.148.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.148.gz) | 3.2M |
| 2015-10-12 | [gp_association.goa_ref_human.149.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.149.gz) | 3.3M |
| 2015-11-09 | [gp_association.goa_ref_human.150.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.150.gz) | 3.3M |
| 2015-12-07 | [gp_association.goa_ref_human.151.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.151.gz) | 3.3M |
| 2016-01-04 | [gp_association.goa_ref_human.152.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.152.gz) | 3.4M |
| 2016-01-20 | [gp_association.goa_ref_human.153.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.153.gz) | 3.4M |
| 2016-02-15 | [gp_association.goa_ref_human.154.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.154.gz) | 3.4M |
| 2016-03-14 | [gp_association.goa_ref_human.155.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.155.gz) | 3.4M |
| 2016-04-11 | [gp_association.goa_ref_human.156.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.156.gz) | 3.3M |
| 2016-05-09 | [gp_association.goa_ref_human.157.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_ref_human.157.gz) | 3.3M |

### `gp_information.goa_human` — 18 files

| Date | File | Size |
|---|---|---:|
| 2015-01-05 | [gp_information.goa_human.140.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.140.gz) | 1.5M |
| 2015-02-02 | [gp_information.goa_human.141.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.141.gz) | 1.5M |
| 2015-03-02 | [gp_information.goa_human.142.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.142.gz) | 1.8M |
| 2015-03-30 | [gp_information.goa_human.143.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.143.gz) | 2.2M |
| 2015-04-27 | [gp_information.goa_human.144.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.144.gz) | 2.2M |
| 2015-05-26 | [gp_information.goa_human.145.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.145.gz) | 2.2M |
| 2015-06-22 | [gp_information.goa_human.146.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.146.gz) | 2.2M |
| 2015-07-20 | [gp_information.goa_human.147.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.147.gz) | 2.2M |
| 2015-09-14 | [gp_information.goa_human.148.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.148.gz) | 2.2M |
| 2015-10-12 | [gp_information.goa_human.149.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.149.gz) | 2.2M |
| 2015-11-09 | [gp_information.goa_human.150.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.150.gz) | 2.2M |
| 2015-12-07 | [gp_information.goa_human.151.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.151.gz) | 2.2M |
| 2016-01-04 | [gp_information.goa_human.152.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.152.gz) | 2.2M |
| 2016-01-20 | [gp_information.goa_human.153.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.153.gz) | 2.2M |
| 2016-02-15 | [gp_information.goa_human.154.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.154.gz) | 2.2M |
| 2016-03-14 | [gp_information.goa_human.155.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.155.gz) | 2.2M |
| 2016-04-11 | [gp_information.goa_human.156.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.156.gz) | 3.2M |
| 2016-05-09 | [gp_information.goa_human.157.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.157.gz) | 3.2M |

### `gp_information.goa_ref_human` — 18 files

| Date | File | Size |
|---|---|---:|
| 2015-01-05 | [gp_information.goa_ref_human.140.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.140.gz) | 571K |
| 2015-02-02 | [gp_information.goa_ref_human.141.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.141.gz) | 571K |
| 2015-03-02 | [gp_information.goa_ref_human.142.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.142.gz) | 573K |
| 2015-03-30 | [gp_information.goa_ref_human.143.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.143.gz) | 580K |
| 2015-04-27 | [gp_information.goa_ref_human.144.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.144.gz) | 580K |
| 2015-05-26 | [gp_information.goa_ref_human.145.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.145.gz) | 581K |
| 2015-06-22 | [gp_information.goa_ref_human.146.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.146.gz) | 581K |
| 2015-07-20 | [gp_information.goa_ref_human.147.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.147.gz) | 581K |
| 2015-09-14 | [gp_information.goa_ref_human.148.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.148.gz) | 582K |
| 2015-10-12 | [gp_information.goa_ref_human.149.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.149.gz) | 582K |
| 2015-11-09 | [gp_information.goa_ref_human.150.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.150.gz) | 582K |
| 2015-12-07 | [gp_information.goa_ref_human.151.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.151.gz) | 582K |
| 2016-01-04 | [gp_information.goa_ref_human.152.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.152.gz) | 582K |
| 2016-01-20 | [gp_information.goa_ref_human.153.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.153.gz) | 582K |
| 2016-02-15 | [gp_information.goa_ref_human.154.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.154.gz) | 583K |
| 2016-03-14 | [gp_information.goa_ref_human.155.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.155.gz) | 583K |
| 2016-04-11 | [gp_information.goa_ref_human.156.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.156.gz) | 583K |
| 2016-05-09 | [gp_information.goa_ref_human.157.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_ref_human.157.gz) | 583K |

### `goa_human.gaf` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human.gaf.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.158.gz) | 4.6M |
| 2016-07-04 | [goa_human.gaf.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.159.gz) | 4.7M |
| 2016-09-14 | [goa_human.gaf.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.160.gz) | 4.4M |
| 2016-10-03 | [goa_human.gaf.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.161.gz) | 4.4M |
| 2016-10-31 | [goa_human.gaf.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.162.gz) | 4.7M |
| 2016-11-28 | [goa_human.gaf.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.163.gz) | 4.8M |

### `goa_human.gpa` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human.gpa.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.158.gz) | 3.5M |
| 2016-07-04 | [goa_human.gpa.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.159.gz) | 3.5M |
| 2016-09-14 | [goa_human.gpa.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.160.gz) | 3.2M |
| 2016-10-03 | [goa_human.gpa.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.161.gz) | 3.2M |
| 2016-10-31 | [goa_human.gpa.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.162.gz) | 3.5M |
| 2016-11-28 | [goa_human.gpa.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.163.gz) | 3.6M |

### `goa_human.gpi` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-06 | [goa_human.gpi.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.158.gz) | 588K |
| 2016-07-04 | [goa_human.gpi.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.159.gz) | 589K |
| 2016-09-14 | [goa_human.gpi.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.160.gz) | 677K |
| 2016-10-03 | [goa_human.gpi.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.161.gz) | 677K |
| 2016-10-31 | [goa_human.gpi.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.162.gz) | 677K |
| 2016-11-28 | [goa_human.gpi.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.163.gz) | 678K |

## 2. EBI GOA auxiliary subsets in 2016

### `goa_human_complex.gaf` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human_complex.gaf.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gaf.158.gz) | 20K |
| 2016-07-04 | [goa_human_complex.gaf.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gaf.159.gz) | 20K |
| 2016-09-14 | [goa_human_complex.gaf.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gaf.160.gz) | 21K |
| 2016-10-03 | [goa_human_complex.gaf.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gaf.161.gz) | 22K |
| 2016-10-31 | [goa_human_complex.gaf.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gaf.162.gz) | 23K |
| 2016-11-28 | [goa_human_complex.gaf.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gaf.163.gz) | 23K |

### `goa_human_complex.gpa` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human_complex.gpa.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpa.158.gz) | 25K |
| 2016-07-04 | [goa_human_complex.gpa.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpa.159.gz) | 26K |
| 2016-09-14 | [goa_human_complex.gpa.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpa.160.gz) | 28K |
| 2016-10-03 | [goa_human_complex.gpa.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpa.161.gz) | 29K |
| 2016-10-31 | [goa_human_complex.gpa.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpa.162.gz) | 30K |
| 2016-11-28 | [goa_human_complex.gpa.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpa.163.gz) | 30K |

### `goa_human_complex.gpi` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-06 | [goa_human_complex.gpi.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpi.158.gz) | 16K |
| 2016-07-04 | [goa_human_complex.gpi.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpi.159.gz) | 17K |
| 2016-09-14 | [goa_human_complex.gpi.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpi.160.gz) | 18K |
| 2016-10-03 | [goa_human_complex.gpi.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpi.161.gz) | 19K |
| 2016-10-31 | [goa_human_complex.gpi.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpi.162.gz) | 20K |
| 2016-11-28 | [goa_human_complex.gpi.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_complex.gpi.163.gz) | 20K |

### `goa_human_isoform.gaf` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human_isoform.gaf.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gaf.158.gz) | 1.1M |
| 2016-07-04 | [goa_human_isoform.gaf.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gaf.159.gz) | 1.1M |
| 2016-09-14 | [goa_human_isoform.gaf.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gaf.160.gz) | 1.1M |
| 2016-10-03 | [goa_human_isoform.gaf.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gaf.161.gz) | 1.1M |
| 2016-10-31 | [goa_human_isoform.gaf.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gaf.162.gz) | 1.2M |
| 2016-11-28 | [goa_human_isoform.gaf.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gaf.163.gz) | 1.2M |

### `goa_human_isoform.gpa` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human_isoform.gpa.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpa.158.gz) | 609K |
| 2016-07-04 | [goa_human_isoform.gpa.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpa.159.gz) | 616K |
| 2016-09-14 | [goa_human_isoform.gpa.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpa.160.gz) | 602K |
| 2016-10-03 | [goa_human_isoform.gpa.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpa.161.gz) | 612K |
| 2016-10-31 | [goa_human_isoform.gpa.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpa.162.gz) | 668K |
| 2016-11-28 | [goa_human_isoform.gpa.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpa.163.gz) | 676K |

### `goa_human_isoform.gpi` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-06 | [goa_human_isoform.gpi.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpi.158.gz) | 1.2M |
| 2016-07-04 | [goa_human_isoform.gpi.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpi.159.gz) | 1.2M |
| 2016-09-14 | [goa_human_isoform.gpi.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpi.160.gz) | 1.4M |
| 2016-10-03 | [goa_human_isoform.gpi.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpi.161.gz) | 1.4M |
| 2016-10-31 | [goa_human_isoform.gpi.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpi.162.gz) | 1.8M |
| 2016-11-28 | [goa_human_isoform.gpi.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_isoform.gpi.163.gz) | 1.8M |

### `goa_human_rna.gaf` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human_rna.gaf.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gaf.158.gz) | 18K |
| 2016-07-04 | [goa_human_rna.gaf.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gaf.159.gz) | 24K |
| 2016-09-14 | [goa_human_rna.gaf.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gaf.160.gz) | 27K |
| 2016-10-03 | [goa_human_rna.gaf.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gaf.161.gz) | 30K |
| 2016-10-31 | [goa_human_rna.gaf.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gaf.162.gz) | 30K |
| 2016-11-28 | [goa_human_rna.gaf.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gaf.163.gz) | 31K |

### `goa_human_rna.gpa` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-07 | [goa_human_rna.gpa.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpa.158.gz) | 34K |
| 2016-07-04 | [goa_human_rna.gpa.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpa.159.gz) | 45K |
| 2016-09-14 | [goa_human_rna.gpa.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpa.160.gz) | 51K |
| 2016-10-03 | [goa_human_rna.gpa.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpa.161.gz) | 55K |
| 2016-10-31 | [goa_human_rna.gpa.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpa.162.gz) | 55K |
| 2016-11-28 | [goa_human_rna.gpa.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpa.163.gz) | 57K |

### `goa_human_rna.gpi` — 6 files

| Date | File | Size |
|---|---|---:|
| 2016-06-06 | [goa_human_rna.gpi.158.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpi.158.gz) | 1.6M |
| 2016-07-04 | [goa_human_rna.gpi.159.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpi.159.gz) | 1.6M |
| 2016-09-14 | [goa_human_rna.gpi.160.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpi.160.gz) | 1.6M |
| 2016-10-03 | [goa_human_rna.gpi.161.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpi.161.gz) | 1.6M |
| 2016-10-31 | [goa_human_rna.gpi.162.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpi.162.gz) | 1.6M |
| 2016-11-28 | [goa_human_rna.gpi.163.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human_rna.gpi.163.gz) | 1.6M |

## 3. Legacy `human.xrefs` series

The EBI README describes this as a cross-reference table connecting UniProt/Swiss-Prot/TrEMBL and Ensembl entries with IPI, RefSeq, HGNC, and legacy NCBI LocusLink identifiers. The archive contains **61 surviving versions**, from release 24 to release 99, but the last is dated 2011-06-10.

| Date | File | Size |
|---|---|---:|
| 2004-10-22 | [human.xrefs.24.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.24.gz) | 1.0M |
| 2004-12-02 | [human.xrefs.25.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.25.gz) | 1.2M |
| 2005-01-04 | [human.xrefs.27.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.27.gz) | 1.3M |
| 2005-02-04 | [human.xrefs.28.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.28.gz) | 1.5M |
| 2005-02-04 | [human.xrefs.29.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.29.gz) | 1.5M |
| 2005-04-04 | [human.xrefs.30.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.30.gz) | 1.5M |
| 2005-06-03 | [human.xrefs.31.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.31.gz) | 1.9M |
| 2005-07-10 | [human.xrefs.32.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.32.gz) | 1.9M |
| 2005-07-31 | [human.xrefs.33.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.33.gz) | 1.9M |
| 2005-09-07 | [human.xrefs.34.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.34.gz) | 2.3M |
| 2005-10-07 | [human.xrefs.35.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.35.gz) | 2.3M |
| 2005-11-09 | [human.xrefs.36.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.36.gz) | 2.3M |
| 2005-12-05 | [human.xrefs.37.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.37.gz) | 2.3M |
| 2006-01-21 | [human.xrefs.38.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.38.gz) | 2.3M |
| 2006-02-17 | [human.xrefs.39.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.39.gz) | 2.4M |
| 2006-05-24 | [human.xrefs.40.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.40.gz) | 2.3M |
| 2006-06-07 | [human.xrefs.41.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.41.gz) | 2.2M |
| 2006-07-07 | [human.xrefs.42.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.42.gz) | 2.3M |
| 2006-08-22 | [human.xrefs.43.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.43.gz) | 2.3M |
| 2006-08-31 | [human.xrefs.44.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.44.gz) | 2.3M |
| 2006-11-01 | [human.xrefs.45.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.45.gz) | 2.3M |
| 2006-11-29 | [human.xrefs.46.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.46.gz) | 2.3M |
| 2007-01-12 | [human.xrefs.47.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.47.gz) | 2.3M |
| 2007-02-10 | [human.xrefs.48.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.48.gz) | 2.3M |
| 2007-03-27 | [human.xrefs.49.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.49.gz) | 2.3M |
| 2007-04-19 | [human.xrefs.50.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.50.gz) | 2.3M |
| 2007-05-18 | [human.xrefs.51.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.51.gz) | 2.3M |
| 2007-08-01 | [human.xrefs.52.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.52.gz) | 2.3M |
| 2007-09-12 | [human.xrefs.54.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.54.gz) | 2.3M |
| 2007-09-30 | [human.xrefs.55.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.55.gz) | 2.3M |
| 2007-10-18 | [human.xrefs.56.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.56.gz) | 2.3M |
| 2007-11-12 | [human.xrefs.57.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.57.gz) | 2.4M |
| 2007-12-03 | [human.xrefs.58.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.58.gz) | 2.4M |
| 2008-01-15 | [human.xrefs.59.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.59.gz) | 2.4M |
| 2008-02-06 | [human.xrefs.60.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.60.gz) | 2.4M |
| 2008-02-26 | [human.xrefs.61.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.61.gz) | 2.4M |
| 2008-03-30 | [human.xrefs.62.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.62.gz) | 2.5M |
| 2008-04-29 | [human.xrefs.63.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.63.gz) | 2.5M |
| 2008-05-16 | [human.xrefs.64.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.64.gz) | 2.4M |
| 2008-06-03 | [human.xrefs.65.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.65.gz) | 2.5M |
| 2008-06-27 | [human.xrefs.66.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.66.gz) | 2.4M |
| 2008-08-28 | [human.xrefs.67.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.67.gz) | 2.4M |
| 2008-09-23 | [human.xrefs.68.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.68.gz) | 2.5M |
| 2008-11-24 | [human.xrefs.69.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.69.gz) | 2.5M |
| 2009-02-11 | [human.xrefs.70.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.70.gz) | 2.5M |
| 2010-04-20 | [human.xrefs.84.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.84.gz) | 3.0M |
| 2010-04-20 | [human.xrefs.85.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.85.gz) | 3.0M |
| 2010-06-09 | [human.xrefs.86.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.86.gz) | 3.1M |
| 2010-07-06 | [human.xrefs.87.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.87.gz) | 3.1M |
| 2010-07-21 | [human.xrefs.88.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.88.gz) | 3.1M |
| 2010-07-21 | [human.xrefs.89.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.89.gz) | 3.1M |
| 2010-09-28 | [human.xrefs.90.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.90.gz) | 3.1M |
| 2010-11-02 | [human.xrefs.91.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.91.gz) | 3.1M |
| 2010-11-26 | [human.xrefs.92.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.92.gz) | 3.0M |
| 2010-11-26 | [human.xrefs.93.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.93.gz) | 3.0M |
| 2011-01-11 | [human.xrefs.94.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.94.gz) | 3.0M |
| 2011-02-06 | [human.xrefs.95.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.95.gz) | 3.0M |
| 2011-04-05 | [human.xrefs.96.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.96.gz) | 3.2M |
| 2011-04-05 | [human.xrefs.97.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.97.gz) | 3.2M |
| 2011-05-31 | [human.xrefs.98.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.98.gz) | 3.2M |
| 2011-06-10 | [human.xrefs.99.gz](https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/human.xrefs.99.gz) | 3.2M |

## 4. GO monthly release archive and `gp2protein`

The GO archive provides monthly release roots. For each date below, the relevant browser path is `annotations/gp2protein/`. The archive documentation says this class maps gene-product identifiers—usually model-organism database IDs—to UniProt accessions.

| Monthly release | Release root | gp2protein directory |
|---|---|---|
| 2015-01-01 | [https://release.geneontology.org/2015-01-01/](https://release.geneontology.org/2015-01-01/) | [https://release.geneontology.org/2015-01-01/annotations/gp2protein/](https://release.geneontology.org/2015-01-01/annotations/gp2protein/) |
| 2015-02-01 | [https://release.geneontology.org/2015-02-01/](https://release.geneontology.org/2015-02-01/) | [https://release.geneontology.org/2015-02-01/annotations/gp2protein/](https://release.geneontology.org/2015-02-01/annotations/gp2protein/) |
| 2015-03-01 | [https://release.geneontology.org/2015-03-01/](https://release.geneontology.org/2015-03-01/) | [https://release.geneontology.org/2015-03-01/annotations/gp2protein/](https://release.geneontology.org/2015-03-01/annotations/gp2protein/) |
| 2015-04-01 | [https://release.geneontology.org/2015-04-01/](https://release.geneontology.org/2015-04-01/) | [https://release.geneontology.org/2015-04-01/annotations/gp2protein/](https://release.geneontology.org/2015-04-01/annotations/gp2protein/) |
| 2015-05-01 | [https://release.geneontology.org/2015-05-01/](https://release.geneontology.org/2015-05-01/) | [https://release.geneontology.org/2015-05-01/annotations/gp2protein/](https://release.geneontology.org/2015-05-01/annotations/gp2protein/) |
| 2015-06-01 | [https://release.geneontology.org/2015-06-01/](https://release.geneontology.org/2015-06-01/) | [https://release.geneontology.org/2015-06-01/annotations/gp2protein/](https://release.geneontology.org/2015-06-01/annotations/gp2protein/) |
| 2015-07-01 | [https://release.geneontology.org/2015-07-01/](https://release.geneontology.org/2015-07-01/) | [https://release.geneontology.org/2015-07-01/annotations/gp2protein/](https://release.geneontology.org/2015-07-01/annotations/gp2protein/) |
| 2015-08-01 | [https://release.geneontology.org/2015-08-01/](https://release.geneontology.org/2015-08-01/) | [https://release.geneontology.org/2015-08-01/annotations/gp2protein/](https://release.geneontology.org/2015-08-01/annotations/gp2protein/) |
| 2015-09-01 | [https://release.geneontology.org/2015-09-01/](https://release.geneontology.org/2015-09-01/) | [https://release.geneontology.org/2015-09-01/annotations/gp2protein/](https://release.geneontology.org/2015-09-01/annotations/gp2protein/) |
| 2015-10-01 | [https://release.geneontology.org/2015-10-01/](https://release.geneontology.org/2015-10-01/) | [https://release.geneontology.org/2015-10-01/annotations/gp2protein/](https://release.geneontology.org/2015-10-01/annotations/gp2protein/) |
| 2015-11-01 | [https://release.geneontology.org/2015-11-01/](https://release.geneontology.org/2015-11-01/) | [https://release.geneontology.org/2015-11-01/annotations/gp2protein/](https://release.geneontology.org/2015-11-01/annotations/gp2protein/) |
| 2015-12-01 | [https://release.geneontology.org/2015-12-01/](https://release.geneontology.org/2015-12-01/) | [https://release.geneontology.org/2015-12-01/annotations/gp2protein/](https://release.geneontology.org/2015-12-01/annotations/gp2protein/) |
| 2016-01-01 | [https://release.geneontology.org/2016-01-01/](https://release.geneontology.org/2016-01-01/) | [https://release.geneontology.org/2016-01-01/annotations/gp2protein/](https://release.geneontology.org/2016-01-01/annotations/gp2protein/) |
| 2016-02-01 | [https://release.geneontology.org/2016-02-01/](https://release.geneontology.org/2016-02-01/) | [https://release.geneontology.org/2016-02-01/annotations/gp2protein/](https://release.geneontology.org/2016-02-01/annotations/gp2protein/) |
| 2016-03-01 | [https://release.geneontology.org/2016-03-01/](https://release.geneontology.org/2016-03-01/) | [https://release.geneontology.org/2016-03-01/annotations/gp2protein/](https://release.geneontology.org/2016-03-01/annotations/gp2protein/) |
| 2016-04-01 | [https://release.geneontology.org/2016-04-01/](https://release.geneontology.org/2016-04-01/) | [https://release.geneontology.org/2016-04-01/annotations/gp2protein/](https://release.geneontology.org/2016-04-01/annotations/gp2protein/) |
| 2016-05-01 | [https://release.geneontology.org/2016-05-01/](https://release.geneontology.org/2016-05-01/) | [https://release.geneontology.org/2016-05-01/annotations/gp2protein/](https://release.geneontology.org/2016-05-01/annotations/gp2protein/) |
| 2016-06-01 | [https://release.geneontology.org/2016-06-01/](https://release.geneontology.org/2016-06-01/) | [https://release.geneontology.org/2016-06-01/annotations/gp2protein/](https://release.geneontology.org/2016-06-01/annotations/gp2protein/) |
| 2016-07-01 | [https://release.geneontology.org/2016-07-01/](https://release.geneontology.org/2016-07-01/) | [https://release.geneontology.org/2016-07-01/annotations/gp2protein/](https://release.geneontology.org/2016-07-01/annotations/gp2protein/) |
| 2016-08-01 | [https://release.geneontology.org/2016-08-01/](https://release.geneontology.org/2016-08-01/) | [https://release.geneontology.org/2016-08-01/annotations/gp2protein/](https://release.geneontology.org/2016-08-01/annotations/gp2protein/) |
| 2016-09-01 | [https://release.geneontology.org/2016-09-01/](https://release.geneontology.org/2016-09-01/) | [https://release.geneontology.org/2016-09-01/annotations/gp2protein/](https://release.geneontology.org/2016-09-01/annotations/gp2protein/) |
| 2016-10-01 | [https://release.geneontology.org/2016-10-01/](https://release.geneontology.org/2016-10-01/) | [https://release.geneontology.org/2016-10-01/annotations/gp2protein/](https://release.geneontology.org/2016-10-01/annotations/gp2protein/) |
| 2016-11-01 | [https://release.geneontology.org/2016-11-01/](https://release.geneontology.org/2016-11-01/) | [https://release.geneontology.org/2016-11-01/annotations/gp2protein/](https://release.geneontology.org/2016-11-01/annotations/gp2protein/) |
| 2016-12-01 | [https://release.geneontology.org/2016-12-01/](https://release.geneontology.org/2016-12-01/) | [https://release.geneontology.org/2016-12-01/annotations/gp2protein/](https://release.geneontology.org/2016-12-01/annotations/gp2protein/) |

### Names recovered from the official GO SVN history

The archive-reconstruction SVN log contains 2015/16 updates for these families:

- `gp2protein.fb` — FlyBase releases.
- `gp2protein.mgi.gz` — Mouse Genome Informatics.
- `gp2protein.rgd` / `gp2protein.rgd.gz` — Rat Genome Database, submitted weekly.
- ZFIN gp2protein submissions.
- `gp2protein.ecocyc` — EcoCyc.
- DictyBase and GeneDB submissions for several pathogens.

No human-wide `gp2protein.human` was found in the archive history. For a human workflow, treat this branch as corroborating context unless the original pipeline used a model-organism ID namespace.

## 5. UniProt historical mapping sources

### Exact surviving 2015 mapping snapshots

| File | Purpose | Notes |
|---|---|---|
| [`idmapping_selected.tab.2015_03.gz`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.2015_03.gz) | Wide UniProt ID mapping table | GeneID is column 3; GO is column 7; NCBI taxon is column 13. Filter taxon 9606. |
| [`idmapping.dat.2015_03.gz`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.2015_03.gz) | Long-form mapping table | One source ID / target ID / target namespace per row. Filter for GeneID and human accessions. |

### Monthly UniProt release tarballs

The official previous-release area exposes monthly release directories. The knowledgebase tarballs can be used to rebuild a human GeneID↔UniProt crosswalk from UniProt flat-file `DR   GeneID;` records when a separate historical idmapping file is absent.

| Release | Directory |
|---|---|
| 2015-01 | [release-2015_01](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_01/) |
| 2015-02 | [release-2015_02](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_02/) |
| 2015-03 | [release-2015_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_03/) |
| 2015-04 | [release-2015_04](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_04/) |
| 2015-05 | [release-2015_05](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_05/) |
| 2015-06 | [release-2015_06](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_06/) |
| 2015-07 | [release-2015_07](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_07/) |
| 2015-08 | [release-2015_08](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_08/) |
| 2015-09 | [release-2015_09](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_09/) |
| 2015-10 | [release-2015_10](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_10/) |
| 2015-11 | [release-2015_11](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_11/) |
| 2015-12 | [release-2015_12](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2015_12/) |
| 2016-01 | [release-2016_01](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_01/) |
| 2016-02 | [release-2016_02](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_02/) |
| 2016-03 | [release-2016_03](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_03/) |
| 2016-04 | [release-2016_04](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_04/) |
| 2016-05 | [release-2016_05](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_05/) |
| 2016-06 | [release-2016_06](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_06/) |
| 2016-07 | [release-2016_07](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_07/) |
| 2016-08 | [release-2016_08](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_08/) |
| 2016-09 | [release-2016_09](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_09/) |
| 2016-10 | [release-2016_10](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_10/) |
| 2016-11 | [release-2016_11](https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2016_11/) |

## 6. NCBI Gene DATA candidates

These are strong candidates for the original file names, but the directory is rolling: most files are recalculated daily and NCBI does not retain ordinary historical snapshots in the same directory.

| File | Role | Historical issue |
|---|---|---|
| [`gene2go.gz`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz) | Entrez Gene ID ↔ GO term/evidence/category | Most likely Entrez-keyed annotation file; exact 2015/16 bytes not found in the official rolling area. |
| [`gene2accession.gz`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz) | GeneID ↔ genomic/RNA/protein accessions, including Swiss-Prot | Candidate direct bridge to UniProt/Swiss-Prot; rolling. |
| [`gene2refseq.gz`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz) | GeneID ↔ RefSeq accessions | Join to the collaboration map below. |
| [`gene_refseq_uniprotkb_collab.gz`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz) | RefSeq protein ↔ UniProtKB protein | Join through `gene2refseq` to obtain GeneID↔UniProt. |
| [`gene_info.gz`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz) | Gene-centric metadata and dbXrefs | Useful as a secondary bridge/check, but rolling. |
| [`ARCHIVE/`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/ARCHIVE/) | Retired NCBI Gene products | Not a chronological archive of daily `gene2go`/`gene2accession` snapshots. |

## 7. Bioconductor `org.Hs.eg.db`: exact all-in-one snapshots

Each package contains `org.Hs.egUNIPROT` (Entrez↔UniProt) and `org.Hs.egGO` (Entrez↔GO, with ontology and evidence). These are the best practical substitutes when the original standalone files are unknown.

| Bioconductor | Package | Entrez source date | GO source date | Download |
|---|---|---|---|---|
| 3.1 | 3.1.2 | 2015-03-17 | 2015-03-14 | [`org.Hs.eg.db_3.1.2.tar.gz`](https://bioconductor.statistik.tu-dortmund.de/packages/3.1/data/annotation/src/contrib/org.Hs.eg.db_3.1.2.tar.gz) |
| 3.2 | 3.2.3 | 2015-09-27 | 2015-09-19 | [`org.Hs.eg.db_3.2.3.tar.gz`](https://bioconductor.statistik.tu-dortmund.de/packages/3.2/data/annotation/src/contrib/org.Hs.eg.db_3.2.3.tar.gz) |
| 3.3 | 3.3.0 | 2016-03-14 | 2016-03-05 | [`org.Hs.eg.db_3.3.0.tar.gz`](https://bioconductor.statistik.tu-dortmund.de/packages/3.3/data/annotation/src/contrib/org.Hs.eg.db_3.3.0.tar.gz) |
| 3.4 | 3.4.0 | 2016-09-26 | 2016-09-21 | [`org.Hs.eg.db_3.4.0.tar.gz`](https://bioconductor.statistik.tu-dortmund.de/packages/3.4/data/annotation/src/contrib/org.Hs.eg.db_3.4.0.tar.gz) |

## 8. Ensembl exact-release reconstruction

Ensembl releases 79–87 span the target period. Their human core MySQL dumps preserve identifiers and ontology cross-references. Reconstruct with `external_db`, `xref`, `object_xref`, `ontology_xref`, and the gene/transcript/translation tables.

| Release | Approximate release point | Human core dump root |
|---:|---|---|
| 79 | 2015-03 | [https://ftp.ensembl.org/pub/release-79/mysql/homo_sapiens_core_79_38/](https://ftp.ensembl.org/pub/release-79/mysql/homo_sapiens_core_79_38/) |
| 80 | 2015-05 | [https://ftp.ensembl.org/pub/release-80/mysql/homo_sapiens_core_80_38/](https://ftp.ensembl.org/pub/release-80/mysql/homo_sapiens_core_80_38/) |
| 81 | 2015-07 | [https://ftp.ensembl.org/pub/release-81/mysql/homo_sapiens_core_81_38/](https://ftp.ensembl.org/pub/release-81/mysql/homo_sapiens_core_81_38/) |
| 82 | 2015-09 | [https://ftp.ensembl.org/pub/release-82/mysql/homo_sapiens_core_82_38/](https://ftp.ensembl.org/pub/release-82/mysql/homo_sapiens_core_82_38/) |
| 83 | 2015-12 | [https://ftp.ensembl.org/pub/release-83/mysql/homo_sapiens_core_83_38/](https://ftp.ensembl.org/pub/release-83/mysql/homo_sapiens_core_83_38/) |
| 84 | 2016-03 | [https://ftp.ensembl.org/pub/release-84/mysql/homo_sapiens_core_84_38/](https://ftp.ensembl.org/pub/release-84/mysql/homo_sapiens_core_84_38/) |
| 85 | 2016-07 | [https://ftp.ensembl.org/pub/release-85/mysql/homo_sapiens_core_85_38/](https://ftp.ensembl.org/pub/release-85/mysql/homo_sapiens_core_85_38/) |
| 86 | 2016-10 | [https://ftp.ensembl.org/pub/release-86/mysql/homo_sapiens_core_86_38/](https://ftp.ensembl.org/pub/release-86/mysql/homo_sapiens_core_86_38/) |
| 87 | 2016-12 | [https://ftp.ensembl.org/pub/release-87/mysql/homo_sapiens_core_87_38/](https://ftp.ensembl.org/pub/release-87/mysql/homo_sapiens_core_87_38/) |

Expected table files include:

- `external_db.txt.gz`
- `xref.txt.gz`
- `object_xref.txt.gz`
- `ontology_xref.txt.gz`
- `gene.txt.gz`, `transcript.txt.gz`, and `translation.txt.gz`

## 9. HGNC quarterly complete-set archive

HGNC documents an archive of monthly and quarterly complete-set files. Quarterly snapshots are retained and the complete set includes both `entrez_id` and `uniprot_ids` fields. This is a strong mapping source, but the exact object names for the 2015/16 quarterlies were not recoverable from the JavaScript-only archive browser during this search.

- Archive help: https://www.genenames.org/download/archive/
- Current complete-set format: https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt

## 10. Most likely original file pair by workflow style

### A. Workflow keyed by UniProt accession

- GO annotations: `gene_association.goa_human.<release>.gz` or `goa_human.gaf.<release>.gz`.
- Entrez mapping: UniProt `idmapping_selected.tab.2015_03.gz`, NCBI `gene2accession.gz`, or HGNC/Bioconductor mapping.

### B. Workflow keyed by Entrez Gene ID

- GO annotations: NCBI `gene2go.gz` or Bioconductor `org.Hs.egGO`.
- UniProt mapping: NCBI `gene2accession.gz`, UniProt idmapping, or Bioconductor `org.Hs.egUNIPROT`.

### C. GOA GPAD/GPI workflow

- Annotation rows: `gp_association.goa_human.<release>.gz` / `goa_human.gpa.<release>.gz`.
- Gene-product metadata: matching `gp_information.goa_human.<release>.gz` / `goa_human.gpi.<release>.gz`.
- Do not assume GPI alone is a complete Entrez↔UniProt crosswalk; inspect its xref column before using it as such.

## 11. Recommended assembly order

1. Download Bioconductor 3.1.2 and 3.4.0 first to bracket the target period with exact Entrez↔UniProt and Entrez↔GO tables.
2. Download matching GOA full GAF releases around the source dates: 142/143 (March 2015), 148 (September 2015), 155 (March 2016), and 160 (September 2016).
3. Add UniProt `idmapping_selected.tab.2015_03.gz` and compare the human GeneID↔UniProt pairs against Bioconductor 3.1.2.
4. If the workflow expects GPAD/GPI rather than GAF, obtain the same-numbered `gp_association` and `gp_information` files.
5. Use Ensembl or HGNC snapshots only to fill discrepancies or to reproduce a pipeline that explicitly cited those databases.

## 12. Gaps that remain

- No official byte-for-byte 2015/16 `gene2go.gz` or `gene2accession.gz` snapshot was located in the NCBI Gene DATA archive.
- No official historical human-specific `HUMAN_9606_idmapping*.gz` snapshot was located; the global March 2015 UniProt snapshots are available instead.
- Exact 2015/16 HGNC quarterly object names remain unresolved.
- Anonymous indexed searches did not surface convincing manuscript supplements, GitHub repositories, Zenodo, Figshare, Dryad, or OSF deposits containing the exact old NCBI or human-specific UniProt files. This does not prove they do not exist.

## Primary source roots

- EBI GOA old/HUMAN: https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/
- GO release archive: https://release.geneontology.org/
- UniProt idmapping: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
- UniProt previous releases: https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/
- NCBI Gene DATA: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
- Ensembl releases: https://ftp.ensembl.org/pub/
- HGNC archive help: https://www.genenames.org/download/archive/
- Bioconductor mirror: https://bioconductor.statistik.tu-dortmund.de/packages/

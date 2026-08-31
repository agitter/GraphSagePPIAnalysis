# B000 inventory receipt

- Inventory: `local_upload_inventory_20260827T145903Z.csv`
- SHA-256: `dae73b7f54d28a089917ba8cf5b7a7f62b2ea0496cba097f3c1e1bfd94e0e9df`
- Rows: 19
- Declared bytes: 337,849,631 (322.20 MiB)
- Status: parsed successfully; no B101+ raw candidate bytes have been uploaded.
- Deletion status: the uploaded B000 CSV may be removed after this receipt is saved; retain your local original inventory and all raw candidate files.

## Recommended upload sequence

| Batch | Files | Total MiB | Purpose |
|---|---|---:|---|
| B101 | `goa_human.gaf.159.gz`<br>`goa_human.gpa.159.gz`<br>`goa_human.gpi.159.gz` | 8.73 | July 2016 GOA/GPA/GPI core test; highest information gain. |
| B102 | `goa_human.gaf.158.gz`<br>`goa_human.gpa.158.gz`<br>`goa_human.gpi.158.gz` | 8.67 | June 2016 adjacent-release comparison. |
| B103 | `goa_human.gaf.160.gz`<br>`goa_human.gpa.160.gz`<br>`goa_human.gpi.160.gz` | 8.24 | September 2016 adjacent-release comparison. |
| B104 | `gene_association.goa_human.157.gz`<br>`gp_association.goa_human.157.gz` | 10.13 | May 2016 legacy-format boundary; incomplete until gp_information.157 is located. |
| B105 | `2016-06-01-go.obo` | 33.15 | June 2016 ontology propagation rules. |
| B106 | `2016-06-01-gp2protein.geneid.gz` | 36.75 | GO release Entrez-to-UniProt mapping; separate upload because of size. |
| B201 | `2016-12-23-gene2go_human.tsv.gz` | 3.13 | Compact human-filtered Dec 2016 gene2go comparator. |
| B202 | `2016-12-23-gene2go.gz` | 18.82 | Raw Dec 2016 gene2go; only if B201 cannot be validated adequately. |
| B107 | `gene_association.goa_human.156.gz` | 6.17 | April 2016 legacy GAF only; defer unless needed. |
| B108 | `gene_association.goa_human.155.gz` | 6.24 | March 2016 legacy GAF only; defer unless needed. |
| B401 | `HumanBase-blood.dat` | 8.01 | GIANT/HumanBase gold-standard candidate. |
| B402 | `HumanBase-blood_top.gz` | 174.15 | Large GIANT network; defer unless targeted analysis requires it. |

## Important inventory gaps

The following previously requested files were not reported by this B000 inventory. This only describes this inventory run; it does not prove the files are absent elsewhere on your machine.

- `gp_association.goa_human.155.gz`
- `gp_information.goa_human.155.gz`
- `gp_association.goa_human.156.gz`
- `gp_information.goa_human.156.gz`
- `gp_information.goa_human.157.gz`
- `idmapping_selected.tab.2015_03.gz`
- `idmapping.dat.2015_03.gz`
- `org.Hs.eg.db_3.1.2.tar.gz`
- `org.Hs.eg.db_3.2.3.tar.gz`
- `org.Hs.eg.db_3.3.0.tar.gz`
- `org.Hs.eg.db_3.4.0.tar.gz`
- `HumanBase-kidney.dat`
- `blood_sample_tsv.gz`

`gp_information.goa_human.157.gz` is the most immediately useful missing companion because B104 otherwise lacks the legacy gene-product information file. The Bioconductor and 2015 UniProt files remain useful later but do not block B101.

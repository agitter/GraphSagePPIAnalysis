# Source acquisition and checksum verification

Generated: 2026-08-24T16:10:11.351585+00:00

This file accompanies `actual_input_file_manifest.csv` and `source_ledger.csv`. The former contains only files actually supplied and used or inspected. The latter also records web references, historical candidates, and generated outputs.

## Commands

Verify the supplied files in place:

```bash
python scripts/download_or_verify_sources.py \
  --manifest results/actual_input_file_manifest.csv \
  --dest inputs \
  --verify-only \
  --log results/input_verification_log.csv
```

Download missing public direct-file inputs recorded in the manifest and verify any recorded checksums:

```bash
python scripts/download_or_verify_sources.py \
  --manifest results/actual_input_file_manifest.csv \
  --dest inputs \
  --download-missing \
  --log results/download_log.csv
```

Download selected historical candidates whose URLs are recorded in `source_ledger.csv`:

```bash
python scripts/download_or_verify_sources.py   --manifest results/source_ledger.csv   --dest historical_inputs   --download-missing   --include-historical   --artifact goa_human.gaf.159.gz   --artifact goa_human.gpa.159.gz   --artifact goa_human.gpi.159.gz   --log results/historical_download_log.csv
```

Historical downloads should be selected with one or more `--artifact` options. The ledger includes multi-gigabyte UniProt-wide mappings, so an unfiltered historical download can consume substantial storage and time.

MSigDB requires an authenticated account. The manifest therefore records its official download page and the SHA-256 of the supplied archives rather than fabricating a reusable direct URL.

## Actual supplied inputs

| Artifact | Status | SHA-256 | Direct/canonical URL | Source page | How obtained |
| --- | --- | --- | --- | --- | --- |
| `graphsage_ppi.zip` | present | `53aeb76e54fd41b645e7edb48b62929240b89839495396b048086fd212503fbd` | https://snap.stanford.edu/graphsage/ppi.zip | https://snap.stanford.edu/graphsage/ | supplied_by_user; not downloaded by corrected run |
| `dgl_ppi.zip` | present | `1f5b2b09ac0f897fa6aa1338c64ab75a5473674cbba89380120bede8cddb2a6a` | https://data.dgl.ai/dataset/ppi.zip | https://github.com/dmlc/dgl/blob/master/python/dgl/data/ppi.py | supplied_by_user; not downloaded by corrected run |
| `bio-tissue-networks.tar.gz` | present | `2c79e17f4a7c8680a7cbf8b20cef4acf356a7523c9a75fce586494153c0603d1` | https://snap.stanford.edu/ohmnet/bio-tissue-networks.tar.gz | https://snap.stanford.edu/ohmnet/ | supplied_by_user; not downloaded by corrected run |
| `bio-tissue-labels.tar.gz` | present | `6abf272940d2407849bd779e5f85c0377a2fb07c2351d1ebc82e3d06a46bc11d` | https://snap.stanford.edu/ohmnet/bio-tissue-labels.tar.gz | https://snap.stanford.edu/ohmnet/ | supplied_by_user; not downloaded by corrected run |
| `bio-tissue-hierarchy.tar.gz` | present | `c4568a68bb83319bff854eecf73a93f698fe2c41ed6e95639af974dd024ffef7` | https://snap.stanford.edu/ohmnet/bio-tissue-hierarchy.tar.gz | https://snap.stanford.edu/ohmnet/ | supplied_by_user; not downloaded by corrected run |
| `bio-tissue-readme.txt` | present | `7f3372f8ae3a90852951c73b18980386ceb4ad2f5d32d81366adf22fd75e2b20` | https://snap.stanford.edu/ohmnet/bio-tissue-readme.txt | https://snap.stanford.edu/ohmnet/ | supplied_by_user; not downloaded by corrected run |
| `msigdb_v5.1_files_to_download_locally.zip` | present | `5a8b3f10ea92f8e71eaaa0705ab9d3a5229d838864eec9699544b75884bc9e29` | https://www.gsea-msigdb.org/gsea/downloads.jsp | https://www.gsea-msigdb.org/gsea/downloads.jsp | supplied_by_user; not downloaded by corrected run |
| `msigdb_v5.2_files_to_download_locally.zip` | present | `a618c1c60b11570036034e6357e73e80ee43065ec7a57c1dbd238f205405fbdb` | https://www.gsea-msigdb.org/gsea/downloads.jsp | https://www.gsea-msigdb.org/gsea/downloads.jsp | supplied_by_user; not downloaded by corrected run |
| `msigdb_v5.2_chip_files_to_download_locally.zip` | present | `252befc853b5e01cfe99439ec50a7ab8de747cd0abe6224a93077f2d9b0b20fc` | https://www.gsea-msigdb.org/gsea/downloads.jsp | https://www.gsea-msigdb.org/gsea/downloads.jsp | supplied_by_user; not downloaded by corrected run |
| `msigdb_v6.0_files_to_download_locally.zip` | present | `39fa82c4cedc9183c532afb1c1431683536b5945e50fb8be4f5bcce3ac136edf` | https://www.gsea-msigdb.org/gsea/downloads.jsp | https://www.gsea-msigdb.org/gsea/downloads.jsp | supplied_by_user; not downloaded by corrected run |
| `Greene2015.pdf` | present | `15c734d37bf63dc586d9bfb95673612209a2f2d298a0a1dc84fa63a1d7a17ce2` | https://doi.org/10.1038/ng.3259 | https://www.nature.com/articles/ng.3259 | supplied_by_user; not downloaded by corrected run |
| `Greene2015_sup.pdf` | present | `89e84c545590a3d34890f24cf6543a336b59cace8926a83701f503afcd979ed9` | https://www.nature.com/articles/ng.3259#Sec23 | https://www.nature.com/articles/ng.3259 | supplied_by_user; not downloaded by corrected run |
| `Greene2015_Table6.xlsx` | present | `691b9d895ac6d0f6ed7abedb96d9b206965fe221e3ccdae940b4daa5db50533e` | https://www.nature.com/articles/ng.3259#Sec23 | https://www.nature.com/articles/ng.3259 | supplied_by_user; not downloaded by corrected run |
| `Greene2015_Table9.xlsx` | present | `18ae68f28d9b84f4b1cb7f7c7c1cc8eb76716414de2089a329f070a5aeca6cd5` | https://www.nature.com/articles/ng.3259#Sec23 | https://www.nature.com/articles/ng.3259 | supplied_by_user; not downloaded by corrected run |
| `OhmNet.pdf` | present | `e60daf8341d0e322ce58e7c6ad194f7e4b573df7c8aba1716ad78c98992b02fe` | https://doi.org/10.1093/bioinformatics/btx252 | https://academic.oup.com/bioinformatics/article/33/14/i190/3953967 | supplied_by_user; not downloaded by corrected run |
| `investigation_summary_2026_08_23.md` | present | `fe23f5d35c1c3a21bba13c2241e0c8c783c31a0a31026c8e5f1c0d7a8c320d16` | No external URL supplied |  | supplied_by_user; not downloaded by corrected run |
| `Pasted markdown(1).md` | present | `f856e6f29f02a17c83650724f2b8434383dc7952c03e5200d4bb9c47c1c8782e` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | supplied_by_user; not downloaded by corrected run |
| `historical_go_mapping_inventory.md` | present | `ad10c6ed2919409a3dda0f2e89f4a5f7709cb739f9e3129297b1fa9f497f4b0c` | No external URL supplied |  | supplied_by_user; not downloaded by corrected run |

## Historical GO and identifier candidates not materialized in this runtime

Their absence is explicit: none of the numerical GO-source conclusions in this corrected run uses these files.

| Artifact | Candidate URL | Source/index | Status/role |
| --- | --- | --- | --- |
| `gene_association.goa_human.155.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.155.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-03-14 |
| `gp_association.goa_human.155.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.155.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-03-14 |
| `gp_information.goa_human.155.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.155.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-03-14 |
| `gene_association.goa_human.156.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.156.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-04-11 |
| `gp_association.goa_human.156.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.156.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-04-11 |
| `gp_information.goa_human.156.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.156.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-04-11 |
| `gene_association.goa_human.157.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gene_association.goa_human.157.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-05-09 |
| `gp_association.goa_human.157.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_association.goa_human.157.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-05-09 |
| `gp_information.goa_human.157.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/gp_information.goa_human.157.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-05-09 |
| `goa_human.gaf.158.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.158.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-06-07 |
| `goa_human.gpa.158.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.158.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-06-07 |
| `goa_human.gpi.158.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.158.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-06-07 |
| `goa_human.gaf.159.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.159.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-07-04 |
| `goa_human.gpa.159.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.159.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-07-04 |
| `goa_human.gpi.159.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.159.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-07-04 |
| `goa_human.gaf.160.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.160.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-09-14 |
| `goa_human.gpa.160.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.160.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-09-14 |
| `goa_human.gpi.160.gz` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.160.gz | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | 2016-09-14 |
| `2016-06-01-go.obo` | https://release.geneontology.org/2016-06-01/ontology/go.obo | https://release.geneontology.org/2016-06-01/ |  |
| `2016-06-01-go-basic.obo` | https://release.geneontology.org/2016-06-01/ontology/go-basic.obo | https://release.geneontology.org/2016-06-01/ |  |
| `2016-06-01-gp2protein.geneid.gz` | https://release.geneontology.org/2016-06-01/annotations/gp2protein/gp2protein.geneid.gz | https://release.geneontology.org/2016-06-01/ |  |
| `org.Hs.eg.db_3.1.2.tar.gz` | https://bioconductor.statistik.tu-dortmund.de/packages/3.1/data/annotation/src/contrib/org.Hs.eg.db_3.1.2.tar.gz | https://bioconductor.org/packages/3.1/data/annotation/html/org.Hs.eg.db.html | Entrez 2015-03-17; GO 2015-03-14 |
| `org.Hs.eg.db_3.2.3.tar.gz` | https://bioconductor.statistik.tu-dortmund.de/packages/3.2/data/annotation/src/contrib/org.Hs.eg.db_3.2.3.tar.gz | https://bioconductor.org/packages/3.2/data/annotation/html/org.Hs.eg.db.html | Entrez 2015-09-27; GO 2015-09-19 |
| `org.Hs.eg.db_3.3.0.tar.gz` | https://bioconductor.statistik.tu-dortmund.de/packages/3.3/data/annotation/src/contrib/org.Hs.eg.db_3.3.0.tar.gz | https://bioconductor.org/packages/3.3/data/annotation/html/org.Hs.eg.db.html | Entrez 2016-03-14; GO 2016-03-05 |
| `org.Hs.eg.db_3.4.0.tar.gz` | https://bioconductor.statistik.tu-dortmund.de/packages/3.4/data/annotation/src/contrib/org.Hs.eg.db_3.4.0.tar.gz | https://bioconductor.org/packages/3.4/data/annotation/html/org.Hs.eg.db.html | Entrez 2016-09-26; GO 2016-09-21 |
| `idmapping_selected.tab.2015_03.gz` | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.2015_03.gz | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ |  |
| `idmapping.dat.2015_03.gz` | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.2015_03.gz | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ |  |
| `NCBI current gene2go.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz | https://ftp.ncbi.nlm.nih.gov/gene/DATA/ |  |
| `NCBI current gene2accession.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz | https://ftp.ncbi.nlm.nih.gov/gene/DATA/ |  |
| `NCBI current gene2refseq.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz | https://ftp.ncbi.nlm.nih.gov/gene/DATA/ |  |
| `NCBI current gene_refseq_uniprotkb_collab.gz` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz | https://ftp.ncbi.nlm.nih.gov/gene/DATA/ |  |

# Actual runtime input manifest through B101

This manifest contains **21 actual materialized inputs**. B101 adds the three GOA release-159 files. User-local files declared in B000A but not materialized in the runtime are tracked separately and are not mislabeled as actual inputs.

| Artifact | Status | Bytes | SHA-256 | Source |
|---|---|---:|---|---|
| `graphsage_ppi.zip` | present | 27,029,260 | `53aeb76e54fd41b645e7edb48b62929240b89839495396b048086fd212503fbd` | https://snap.stanford.edu/graphsage/ppi.zip |
| `dgl_ppi.zip` | present | 7,088,822 | `1f5b2b09ac0f897fa6aa1338c64ab75a5473674cbba89380120bede8cddb2a6a` | https://data.dgl.ai/dataset/ppi.zip |
| `bio-tissue-networks.tar.gz` | present | 15,568,869 | `2c79e17f4a7c8680a7cbf8b20cef4acf356a7523c9a75fce586494153c0603d1` | https://snap.stanford.edu/ohmnet/bio-tissue-networks.tar.gz |
| `bio-tissue-labels.tar.gz` | present | 1,793,328 | `6abf272940d2407849bd779e5f85c0377a2fb07c2351d1ebc82e3d06a46bc11d` | https://snap.stanford.edu/ohmnet/bio-tissue-labels.tar.gz |
| `bio-tissue-hierarchy.tar.gz` | present | 421,785 | `c4568a68bb83319bff854eecf73a93f698fe2c41ed6e95639af974dd024ffef7` | https://snap.stanford.edu/ohmnet/bio-tissue-hierarchy.tar.gz |
| `bio-tissue-readme.txt` | present | 1,465 | `7f3372f8ae3a90852951c73b18980386ceb4ad2f5d32d81366adf22fd75e2b20` | https://snap.stanford.edu/ohmnet/bio-tissue-readme.txt |
| `msigdb_v5.1_files_to_download_locally.zip` | present | 69,702,250 | `5a8b3f10ea92f8e71eaaa0705ab9d3a5229d838864eec9699544b75884bc9e29` | https://www.gsea-msigdb.org/gsea/downloads.jsp |
| `msigdb_v5.2_files_to_download_locally.zip` | present | 92,768,160 | `a618c1c60b11570036034e6357e73e80ee43065ec7a57c1dbd238f205405fbdb` | https://www.gsea-msigdb.org/gsea/downloads.jsp |
| `msigdb_v5.2_chip_files_to_download_locally.zip` | present | 67,891,640 | `252befc853b5e01cfe99439ec50a7ab8de747cd0abe6224a93077f2d9b0b20fc` | https://www.gsea-msigdb.org/gsea/downloads.jsp |
| `msigdb_v6.0_files_to_download_locally.zip` | present | 72,144,318 | `39fa82c4cedc9183c532afb1c1431683536b5945e50fb8be4f5bcce3ac136edf` | https://www.gsea-msigdb.org/gsea/downloads.jsp |
| `Greene2015.pdf` | present | 1,729,378 | `15c734d37bf63dc586d9bfb95673612209a2f2d298a0a1dc84fa63a1d7a17ce2` | https://doi.org/10.1038/ng.3259 |
| `Greene2015_sup.pdf` | present | 6,701,167 | `89e84c545590a3d34890f24cf6543a336b59cace8926a83701f503afcd979ed9` | https://www.nature.com/articles/ng.3259#Sec23 |
| `Greene2015_Table6.xlsx` | present | 37,992 | `691b9d895ac6d0f6ed7abedb96d9b206965fe221e3ccdae940b4daa5db50533e` | https://www.nature.com/articles/ng.3259#Sec23 |
| `Greene2015_Table9.xlsx` | present | 140,272 | `18ae68f28d9b84f4b1cb7f7c7c1cc8eb76716414de2089a329f070a5aeca6cd5` | https://www.nature.com/articles/ng.3259#Sec23 |
| `OhmNet.pdf` | present | 1,458,562 | `e60daf8341d0e322ce58e7c6ad194f7e4b573df7c8aba1716ad78c98992b02fe` | https://doi.org/10.1093/bioinformatics/btx252 |
| `investigation_summary_2026_08_23.md` | present | 13,819 | `fe23f5d35c1c3a21bba13c2241e0c8c783c31a0a31026c8e5f1c0d7a8c320d16` |  |
| `Pasted markdown(1).md` | present | 237,584 | `f856e6f29f02a17c83650724f2b8434383dc7952c03e5200d4bb9c47c1c8782e` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ |
| `historical_go_mapping_inventory.md` | present | 55,969 | `ad10c6ed2919409a3dda0f2e89f4a5f7709cb739f9e3129297b1fa9f497f4b0c` |  |
| `goa_human.gaf.159.gz` | present_at_B101_analysis_time; conversation_copy_deletion_clearance_issued | 4,890,788 | `b098509f9ce70fbe7e93e7fabb3c8a767b8201fe25662ff6a57959d6cee43838` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.159.gz |
| `goa_human.gpa.159.gz` | present_at_B101_analysis_time; conversation_copy_deletion_clearance_issued | 3,665,510 | `dbe9073e726e804f51b0fe80ecd2897798f22a558e75aeb1c8bfc229357b574a` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.159.gz |
| `goa_human.gpi.159.gz` | present_at_B101_analysis_time; conversation_copy_deletion_clearance_issued | 602,689 | `f783a8c71c464b3b98673eb4dc0ae836ad61bd6c72415ae62828d20d0e735612` | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.159.gz |

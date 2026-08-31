# Requested additional historical inputs

The corrected local analysis is complete for the supplied files. The remaining exact GO-label reconstruction requires historical bytes that were not materialized in this runtime. Candidate URLs are listed in `source_ledger.csv` and `SOURCE_ACQUISITION.md`.

Highest priority:

1. `goa_human.gaf.159.gz`
2. `goa_human.gpa.159.gz`
3. `goa_human.gpi.159.gz`
4. The analogous release 158 and 160 GAF/GPA/GPI files
5. Release 157 `gene_association`, `gp_association`, and `gp_information` files
6. `2016-06-01-go.obo` or `go-basic.obo`
7. `2016-06-01-gp2protein.geneid.gz`
8. `org.Hs.eg.db_3.3.0.tar.gz` and `org.Hs.eg.db_3.4.0.tar.gz`
9. Any human `gene2go` snapshot dated June–September 2016

No claim in the corrected report depends on these files until their bytes are present and hashed.

# Repackage the small UniProt audit outputs

The uploaded `uniprot_2016_mapping.zip` contains only the audit ledger. The ledger says all three large source archives were verified and deleted successfully, but it does not include the extracted record text needed to inspect `GN` and `DR   GeneID;` lines.

No large UniProt file needs to be downloaded again.

Run the packager against the directory that contains the `2016_04`, `2016_05`, and `2016_06` result subdirectories:

```bash
python package_uniprot_audit_outputs.py \
  --root /path/to/uniprot_audit_results \
  --ledger /path/to/uniprot_audit_results/uniprot_2016_mapping_audit_ledger.csv \
  --output uniprot_2016_mapping_complete.zip
```

The output filename receives a UTC datestamp automatically. The script:

1. normalizes Windows-style paths recorded in the ledger;
2. verifies every extracted `.dat` and `.tsv` SHA-256;
3. verifies each provenance JSON against the hash in the ledger notes;
4. excludes all large Swiss-Prot archives;
5. writes a package manifest and performs a ZIP integrity test.

Upload only the resulting small ZIP.

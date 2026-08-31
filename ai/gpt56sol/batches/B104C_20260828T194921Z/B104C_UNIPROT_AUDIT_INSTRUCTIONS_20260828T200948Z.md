# Sequential UniProt 2016 Mapping Audit

## Purpose

Audit the historical Swiss-Prot records for `O95073` and `Q9Y620` immediately before and after the May 9, 2016 `gp2protein.geneid` generation date without retaining all three approximately 1.5 GB archives at once.

The script processes releases `2016_04`, `2016_05`, and `2016_06` sequentially. For each release it:

1. downloads UniProt's official `RELEASE.metalink`;
2. downloads or resumes the reviewed-only Swiss-Prot archive;
3. verifies the official release identifier, exact archive size, and MD5;
4. computes a local SHA-256;
5. streams the entire `uniprot_sprot.dat` member and extracts complete records containing `O95073` or `Q9Y620`;
6. writes the complete records, a compact TSV summary, JSON provenance, and an append-only CSV ledger;
7. hashes all retained outputs; and
8. deletes the large archive only after every check succeeds.

On any failure, the archive or partial download is retained for diagnosis or resumption.

## Required free space

Allow approximately 1.7 GB of free space. Only one full archive is present at a time unless `--keep-archives` is used.

## Test first

```bash
python download_extract_uniprot_2016_mapping_audit.py --self-test
python download_extract_uniprot_2016_mapping_audit.py --dry-run
```

## Run all three releases sequentially

```bash
python download_extract_uniprot_2016_mapping_audit.py \
  --work-dir uniprot_audit_work \
  --output-dir uniprot_audit_results
```

After success, upload the entire small `uniprot_audit_results` directory as a ZIP. Do not upload the downloaded UniProt archives.

## Process one release only

```bash
python download_extract_uniprot_2016_mapping_audit.py \
  --releases 2016_05 \
  --work-dir uniprot_audit_work \
  --output-dir uniprot_audit_results
```

## Preserve the large archives intentionally

```bash
python download_extract_uniprot_2016_mapping_audit.py \
  --keep-archives \
  --work-dir uniprot_audit_work \
  --output-dir uniprot_audit_results
```

This is not recommended under current storage constraints.

## Official releases audited

| Release | Date | Archive bytes | Official MD5 |
|---|---:|---:|---|
| 2016_04 | 2016-04-13 | 1,516,525,310 | `e607b83de1ac87e6f63b13715c049a3f` |
| 2016_05 | 2016-05-11 | 1,504,161,063 | `fe9525832026b03ab34f0971b43c0c81` |
| 2016_06 | 2016-06-08 | 1,504,963,399 | `e3a5ac5a166efc95e9ad06465d5bd2c4` |

The script retrieves these values again from each official `RELEASE.metalink` and refuses to proceed if they disagree with the pinned audit values.

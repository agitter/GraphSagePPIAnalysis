# SAFE TO DELETE — BATCH B104E

Issued: `2026-08-29T12:27:34Z`

The following conversation attachments may be deleted after the final B104E bundle and validation files are available:

```text
uniprot_2016_mapping_audit_ledger.zip
Bytes: 30449
SHA-256: 4ac8cb1a900215ded9dc35a4fc44a4abaaff3e5e774cf62769592f1ab153a7b0

goa_date_screen_results.zip
Bytes: 54246
SHA-256: 502c8ffdb7b809c1665e82d31db52b72e7e855e7e9fcc06a8a2f46a64bc30de9
```

Retained for the UniProt package:

- three extracted complete DAT record files;
- three compact TSV summaries;
- three per-release provenance JSON files;
- three official RELEASE metalinks;
- the append-only audit ledger;
- runtime validation and scientific interpretation outputs.

Retained for the GOA package:

- one detailed JSON for every release 158 through 169;
- the release summary;
- run metadata;
- append-only source-integrity, analysis, and cleanup events;
- runtime validation and scientific summary outputs.

The user-side large UniProt and GOA parent downloads were already deleted only after successful validation by their respective scripts. Keep the local compact audit-result directories or the frozen B104E bundle.

After deleting the two conversation attachments, report:

```text
Deleted B104E
```

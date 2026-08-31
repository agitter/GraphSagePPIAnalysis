# Corrected GraphSAGE PPI provenance analysis bundle

This bundle supersedes the earlier incomplete deliverables.

Start with:

- `results/MASTER_REPRODUCTION_REPORT.md`
- `results/RUN_STATUS.md`
- `results/EXECUTION_DIAGNOSTICS.md`
- `results/actual_input_file_manifest.csv`
- `results/source_ledger.csv`
- `results/SOURCE_ACQUISITION.md`

The raw user-supplied datasets are not duplicated in this bundle. Their exact local paths, canonical source URLs, sizes, and SHA-256 checksums are recorded in the manifests. The scripts expect the supplied raw files under `/mnt/data` by default; paths can be changed through their command-line options.

The historical GOA, gene2go, ontology, mapping, and Bioconductor files listed in the source ledger were not materialized in the corrected runtime. Their absence is explicit in the report and no numerical conclusion is attributed to them.

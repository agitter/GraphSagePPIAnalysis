# B104G reproduction instructions

## External raw inputs

Place these existing user-held files at the indicated paths or substitute explicit paths in the command:

```text
/mnt/data/graphsage_ppi.zip
/mnt/data/bio-tissue-networks.tar.gz
```

Their expected SHA-256 values are:

```text
graphsage_ppi.zip
53aeb76e54fd41b645e7edb48b62929240b89839495396b048086fd212503fbd

bio-tissue-networks.tar.gz
2c79e17f4a7c8680a7cbf8b20cef4acf356a7523c9a75fce586494153c0603d1
```

## Compact retained inputs

The B104G bundle includes:

```text
retained_inputs/prior_compact_inputs/
retained_inputs/GO_annotations-9606-*.tsv.gz
```

Verify their hashes with:

```bash
sha256sum -c retained_inputs/prior_compact_inputs/B104G_prior_compact_input_checksums.sha256
```

The four dhimmel source tables are retained as deterministic gzip files. Decompress them to a temporary directory:

```bash
mkdir -p work/dhimmel
for f in retained_inputs/GO_annotations-9606-*.tsv.gz; do
  gzip -cd "$f" > "work/dhimmel/$(basename "${f%.gz}")"
done
```

## Primary analysis

From the B104G batch directory:

```bash
python scripts/analyze_B104G_full_row_mapping_and_dhimmel.py \
  --graphsage-zip /mnt/data/graphsage_ppi.zip \
  --ohmnet-tar /mnt/data/bio-tissue-networks.tar.gz \
  --core-summary retained_inputs/prior_compact_inputs/core_verification_summary.json \
  --known-row-map retained_inputs/prior_compact_inputs/graphsage_row_to_entrez_topology_features.csv \
  --residual-classes retained_inputs/prior_compact_inputs/wl_residual_after_features.json \
  --gaf retained_inputs/prior_compact_inputs/B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz \
  --gpi retained_inputs/prior_compact_inputs/B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz \
  --gp2protein retained_inputs/prior_compact_inputs/B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz \
  --go-terms retained_inputs/prior_compact_inputs/B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz \
  --go-edges retained_inputs/prior_compact_inputs/B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz \
  --symbol-map retained_inputs/prior_compact_inputs/B104F_MSigDB52_symbol_to_Entrez_relevant.tsv.gz \
  --column-map retained_inputs/prior_compact_inputs/B104C_inferred_unique_121_GO_column_order_20260828T194921Z.csv \
  --msigdb-normalized retained_inputs/prior_compact_inputs/B104C_msigdb_v5.0_normalized_entrez_gene_sets.tsv.gz \
  --feature-rule retained_inputs/prior_compact_inputs/B104E_exact_MSigDB52_feature_generation_rule_20260829T121535Z.csv \
  --dhimmel \
    work/dhimmel/GO_annotations-9606-direct-allev.tsv \
    work/dhimmel/GO_annotations-9606-direct-expev.tsv \
    work/dhimmel/GO_annotations-9606-inferred-allev.tsv \
    work/dhimmel/GO_annotations-9606-inferred-expev.tsv \
  --batch-dir work/B104G_rerun \
  --summary-json work/B104G_rerun/B104G_ANALYSIS_SUMMARY.json
```

The source script is deterministic for the same inputs.

## Independent validations

The independent scripts use the same raw GraphSAGE and OhmNet archives plus the primary output map and compact retained inputs. Their source code and previously accepted output JSONs are in the bundle. The accepted validations report:

```text
row map:      56,944 / 56,944 assignments identical
features:      0 / 2,847,200 mismatched cells
labels:        0 / 6,890,224 mismatched cells
```

## Expected matrix hashes

```text
full 56,944 × 50 feature matrix
274ecfee66596b7e9dfb19b71a7fb39a2611a3c140bc199c25062c6ca75bfca1

full 56,944 × 121 label matrix
677cc50459190ba22afc0762356e573ca56e422c80fd36204669a81694afa78d
```

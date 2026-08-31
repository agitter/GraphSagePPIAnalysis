# Sequential GOA date-range screen

## Purpose

Screen GOA releases 158 through 169 without retaining more than one GAF/GPI pair at a time. The range spans the last pre-OhmNet release control, the exact-membership release 159, all releases available before the `ppi-class_map.json` archive timestamp of 10 May 2017, and release 169 as a post-timestamp negative control.

The screen holds these factors fixed:

- May-2016 historical GeneID–UniProt mapping;
- ambiguity-preserving GraphSAGE projection, with `O95073 -> 25788` excluded;
- evidence codes `EXP, IDA, IEP, IGI, IMP, ISS`;
- ordinary aspect relations `involved_in`, `part_of`, and `enables`;
- `is_a` propagation only;
- archived 1 June 2016 ontology;
- 64-bit unrandomized CPython-2 dictionary-order model used in B104C/B104D.

Because the ontology is held fixed, this is a source-date screen, not a replacement for B105.

## Files

Place these two downloaded files in the same working directory:

```text
screen_goa_release_date_range.py
B104D_GOA_date_screen_reference_pack_20260829T010311Z.zip
```

## Run the likely source range and controls

```bash
python screen_goa_release_date_range.py \
  --reference-pack B104D_GOA_date_screen_reference_pack_20260829T010311Z.zip \
  --releases 158-169 \
  --cache-dir goa_date_screen_cache \
  --output-dir goa_date_screen_results
```

The script downloads one release pair, verifies gzip integrity, records SHA-256, analyzes it, writes compact results, and deletes files downloaded by the script after successful analysis. Files supplied through `--local-release-dir` are never deleted.

To reuse the release-158 and release-159 files already on your machine:

```bash
python screen_goa_release_date_range.py \
  --reference-pack B104D_GOA_date_screen_reference_pack_20260829T010311Z.zip \
  --releases 158-169 \
  --local-release-dir /path/to/your/existing/files \
  --cache-dir goa_date_screen_cache \
  --output-dir goa_date_screen_results
```

## Retained outputs

The result directory contains:

- one JSON report per release;
- one compact cross-release CSV;
- an append-only event log;
- run metadata recording all fixed assumptions and the reference-pack hash.

Zip only the result directory for upload:

```bash
python - <<'PY'
from pathlib import Path
import shutil
root = Path('goa_date_screen_results')
shutil.make_archive('goa_date_screen_results', 'zip', root)
print('Created goa_date_screen_results.zip')
PY
```

Do not upload the downloaded GAF/GPI files.

## Interpretation limits

A later release that matches exactly is a viable source snapshot under the fixed May-2016 mapping and June-2016 ontology. A later release that fails is ruled out under those fixed assumptions. New post-May-2016 UniProt accessions are handled only through unique primary-symbol fallbacks when possible, and each fallback count is reported.

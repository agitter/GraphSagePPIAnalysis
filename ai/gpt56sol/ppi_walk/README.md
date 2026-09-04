# GraphSAGE PPI walk reproduction investigation

**Investigation date:** 2026-09-04  
**Status:** exact byte-level replay achieved for both deposited walk files; no historical seed identified.

This bundle records a forensic investigation of the stochastic
`*-walks.txt` files distributed with GraphSAGE. It is deliberately separate
from the planned non-circular reproduction of the deterministic supervised
GraphSAGE/DGL data.

## tl;dr

We used the published walk files to infer one valid sequence of node-to-node choices that could have produced them, after reconstructing the program’s exact node and neighbor ordering. Replaying the original GraphSAGE walk algorithm with those inferred choices instead of random selections reproduces both files byte-for-byte, without knowing the authors’ original random seed or generator state. For a new graph, we could not reproduce the authors’ exact sequence of random walks without their original seed or full random-generator state.

## Main result

A compact, target-derived neighbor-index tape regenerates both deposited walk
files **byte for byte** when replayed through the archived five-step walk
control flow and the reconstructed NetworkX 1.11 / CPython 2 ordered graph.

| Case | Target pairs | Target bytes | Target SHA-256 | Compressed tape |
|---|---:|---:|---|---:|
| toy PPI | 1,895,817 | 18,523,412 | `4edfb102dad4a1100b992c16412e14a58e051f1d8a7f01c7373336d0d6e864b0` | 2,044,119 bytes |
| full PPI | 8,730,249 | 100,459,319 | `d0d5b5b3f2521d3727e88be3147029bac27a5409ee60c5106b9f12a58e75dbcc` | 9,138,914 bytes |

The tape is derived from the deposited output and is therefore a **forensic
replay resource**, not evidence that the original authors used the same tape,
seed, or decomposition.

## Verify from the original ZIP files

Python 3.10 or newer is sufficient; the verifier uses only the standard
library. It checks the three user-supplied archive hashes, extracts the toy and
full graph/target pairs into a temporary directory, regenerates both walk
files, and performs streaming byte comparisons plus SHA-256 checks.

```bash
python scripts/verify_bundle.py \
  --source-zip /path/to/GraphSAGE-a0fdef9.zip \
  --ppi-zip /path/to/graphsage_ppi.zip \
  --dgl-zip /path/to/dgl_ppi.zip \
  --output-json verification.json
```

Expected final field:

```json
"all_checks_pass": true
```

## Replay one case

```bash
python scripts/replay_choice_tape.py \
  --graph /path/to/ppi-G.json \
  --tape resources/full_choice_tape.u16le.gz \
  --output reproduced-ppi-walks.txt \
  --json replay.json
```

The output intentionally has no final newline, matching the archived writer.

## Why there was no `random.seed` call

The archived command-line utility does not need an explicit seed in order to
run. Importing Python's `random` module creates a process-global generator and
initializes it automatically. The GraphSAGE training scripts seed NumPy and
TensorFlow, but those calls do not seed Python's separate `random` generator.
Consequently, running the utility as supplied yields a valid but generally
new walk file on each fresh process. The fixed deposited walk file then acts as
an input artifact for later training.

The archive does not document whether omission of a Python-random seed was a
conscious policy or an oversight. This bundle does not infer intent.

## What was searched

The investigation used an early-abort oracle: each candidate generated only
as many choices as needed to reach its first mismatch. It searched:

- 51 common and repository/date-motivated seeds under Python-2 and Python-3
  choice behavior;
- every integer seed from `0` through `2^24 - 1` under both choice families;
- every integer seed from `0` through `2^32 - 1` under Python-2 float-based
  choice, for the toy and full graphs and five plausible neighbor orders;
- 51 common seeds after every prior-draw offset from 0 through 10,000,000;
- seeds 0, 1, 42, and 123 after every prior-draw offset through 100,000,000;
- historical clock-fallback-style values at 1/256-second resolution around the
  ZIP timestamps, under five graph-order variants per dataset.

No candidate reproduced either target. In the complete 32-bit screen, 136,241
seeds survived the first four-choice filter; the best two matched only eight
emitted full-PPI pairs and failed on the ninth.

See `REPORT.md`, `walk_reproduction_summary.json`, and
`evidence/seed_search/` for the precise search boundaries and retained
negative evidence.

## Bundle structure

```text
README.md
REPORT.md
ATTRIBUTION_AND_PROVENANCE.md
source_and_target_manifest.csv
archive_member_manifest.csv
walk_reproduction_summary.json
seed_search_register.csv
resources/
  toy_choice_tape.u16le.gz
  full_choice_tape.u16le.gz
  *_tape_metadata.json
scripts/
  py2_walk_machine.py
  replay_choice_tape.py
  verify_bundle.py
  derive_replay_tape*.py
  search_*.cpp / search_*.py
  nondeterminism_demo.py
evidence/
  decomposition/
  machines/
  seed_search/
tests/
  clean_verification.json
SHA256SUMS
```

## Interpretation boundary

The exact replay answers the practical question, “Can the deposited outputs be
regenerated by a documented deterministic algorithm?” with **yes**. It does
not answer, “What seed or entropy bytes did the original process use?” The
latter remains open and is unnecessary for the stated practical criterion.

# GraphSAGE PPI data reconstruction

This directory contains a small, auditable workflow for reconstructing the
published GraphSAGE protein-protein interaction (PPI) data from historical
upstream sources and deriving the corresponding DGL data representation.

The workflow is intentionally narrower than the surrounding forensic research.
It reconstructs the deterministic data needed by **supervised GraphSAGE** and by
the downstream DGL PPI dataset. It does not rerun the historical source screens,
tissue-split experiments, leakage analyses, literature review, or the stochastic
`ppi-walks.txt` generation.

## Current implementation status

The first implementation milestone is complete in this patch:

- acquisition and checksum verification for upstream and reference files;
- reconstruction of all 24 tissue graph blocks;
- reconstruction of the complete 56,944-row node-to-Entrez mapping;
- reconstruction of all 818,716 logical GraphSAGE edge records;
- reconstruction of all 50 MSigDB input features from public v6.1 C1/C3 files;
- target-independent invariant checks;
- independent validation against the released GraphSAGE topology, split flags,
  and feature matrix.

GO-label reconstruction, final GraphSAGE serialization, and the DGL conversion
remain later milestones. Commands and output names are already organized so
those stages can be added without changing the central workflow boundary.

## The non-circularity rule

Reconstruction and validation are separate commands:

```text
upstream sources -> reconstruction -> rebuilt data
                                      |
released targets --------------------+-> validation only
```

`pixi run reproduce` does not read `graphsage_ppi.zip` or `dgl_ppi.zip`.
It must finish when those reference archives are absent. Only
`pixi run validate` is permitted to read them.

This distinction matters because a workflow that reads the desired output while
claiming to regenerate it can accidentally encode the answer in the
reconstruction.

## Directory layout

```text
reproduction/
├── README.md                 # operational entry point
├── METHODS.md                # scientific transformations and evidence levels
├── THIRD_PARTY_NOTICES.md    # data licenses and attribution
├── Snakefile                 # one readable dependency graph
├── pixi.toml                 # environment and task definitions
├── pixi.lock                 # generated and committed after Pixi resolution
├── pyproject.toml            # installable Python package
├── spec/                     # accepted human-readable reconstruction specification
├── src/graphsage_ppi_repro/  # scientific implementation
├── tests/                    # focused synthetic and integration tests
├── build/                    # generated intermediates; ignored by Git
└── results/                  # generated manifests and validation; ignored by Git
```

The package deliberately avoids separate `config/`, `resources/`, `scripts/`,
and `docs/` trees. Small accepted specifications live together in `spec/`;
importable and tested executable code lives under `src/`; operational and
scientific explanations are consolidated here and in `METHODS.md`.

## Prerequisites

Install [Pixi](https://pixi.sh/). The package targets Linux x86-64 and native
Windows with Git Bash. The workflow uses Python for path and archive handling
and does not depend on the shell's current directory.

From `ai/reproduction/`, run:

```bash
pixi install
pixi run test
pixi run reproduce
pixi run validate
```

The canonical environment lock is `pixi.lock`. The initial implementation
patch includes `pixi.toml`; the lock must be generated with Pixi and committed
before the environment is considered frozen.

## Inputs and cache

External inputs are stored in the project-level ignored directory:

```text
ai/data/
```

`spec/sources.tsv` records, for every source:

- whether it is an upstream reconstruction input or a validation-only target;
- the expected filename, byte size, and SHA-256;
- a primary URL and an optional future mirror;
- a short license and provenance note.

A valid cached file is reused. A missing public file is downloaded to a
temporary `.part` path, checked, and moved into place only after its size,
SHA-256, and basic archive or text structure pass. Expected hashes are never
updated automatically from downloaded content.

The current topology-and-feature milestone uses:

- the OhmNet tissue-network archive;
- MSigDB v6.1 C1 Entrez gene sets;
- MSigDB v6.1 C3 Entrez gene sets.

Validation additionally uses the released GraphSAGE PPI ZIP. Later stages will
activate the historical GOA, GO, GeneID-UniProt, and DGL records already listed
in `spec/sources.tsv`.

## Commands

### Reconstruct without targets

```bash
pixi run reproduce
```

The current workflow performs these stages:

```text
verify upstream files
        |
reconstruct selected OhmNet graph blocks and biological row order
        |
select and project the 50 MSigDB feature sets
        |
check frozen target-independent invariants
        |
write a provenance and checksum manifest
```

Principal generated files are:

```text
build/topology/node_to_entrez.tsv.gz
build/topology/edges.tsv.gz
build/topology/summary.json
build/features/ppi-feats.npy
build/features/selected_features.tsv
build/features/summary.json
build/reconstruction_checks.json
results/reconstruction_manifest.json
```

### Validate against the released target

```bash
pixi run validate
```

For the current milestone this checks:

- consecutive node IDs and the identity ID map;
- exact training, validation, and test flags for every row;
- exact graph-wise undirected edge multisets;
- exact feature shape, dtype, and all 2,847,200 feature values;
- byte equality of the individual `ppi-feats.npy` file.

Validation writes:

```text
results/validation.json
results/validation.md
```

ZIP timestamps and compression bytes are not treated as scientific outputs.
The stochastic `ppi-walks.txt` member is explicitly excluded.

### Run fast tests

```bash
pixi run test
```

The tests use small committed synthetic inputs. They cover source verification,
CPython 2.7 dictionary behavior, graph reconstruction, MSigDB selection,
historical NumPy serialization, a target-free reconstruction, and a subsequent
independent validation step.

### Remove generated files

```bash
pixi run clean
```

Only `build/` and `results/` are removed. Downloaded data under `ai/data/` are
preserved.

## Specification files

- `sources.tsv` - external file identities and acquisition roles.
- `specification.yaml` - the accepted global transformation policy and compact
  expected invariants.
- `selected_graphs.tsv` - the 24 released tissue identities, graph order, and
  deposited 20/2/2 split. This is a recovered data-level specification; the
  historical selection algorithm remains unknown.
- `feature_columns.tsv` - the expected 30 C1 and 20 C3 results. The code derives
  these columns from GMT source order and checks them against the table.
- `label_columns.tsv` - the recovered 121 GO-column order for the later label
  milestone.
- `identifier_decisions.tsv` - explicit evidence-backed resolutions for
  historical many-to-many identifier components.

A specification table is not a substitute for reconstruction. For example,
`feature_columns.tsv` is used to verify which columns were derived; the matrix
is built from the GMT memberships, not copied from the table.

## Applying patches from elsewhere in the repository

Patches are generated relative to the Git repository root and therefore contain
paths beginning with `ai/`. From either `ai/` or `ai/reproduction/`:

```bash
PATCH=/path/to/GraphSagePPIAnalysis-reproduction-stage1.patch
ROOT=$(git rev-parse --show-toplevel)
git -C "$ROOT" apply --check "$PATCH"
git -C "$ROOT" apply "$PATCH"
git -C "$ROOT" status --short
```

## Evidence language

The project distinguishes:

- **data-level exact**: directly regenerated and compared exactly;
- **strongly inferred**: an implementation mechanism explains all observations
  but was not found in original preprocessing source;
- **documented**: stated by an official source;
- **open**: evidence does not distinguish among plausible histories.

In particular, the full row mapping is data-level exact under the implemented
mechanism. The historical explanation involving legacy NetworkX construction
and 64-bit unrandomized CPython 2.7 dictionary iteration remains strongly
inferred rather than source-code proven.

See [METHODS.md](METHODS.md) for the complete technical explanation.

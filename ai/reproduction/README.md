# GraphSAGE PPI data reconstruction

This directory contains a small, auditable workflow that reconstructs the
published GraphSAGE protein-protein interaction (PPI) data from historical
upstream sources and will derive the corresponding DGL representation.

The workflow is intentionally narrower than the surrounding forensic research.
It rebuilds deterministic data needed by **supervised GraphSAGE** and DGL. It
does not rerun the historical source screens, tissue-split experiments, leakage
analyses, literature review, or stochastic `ppi-walks.txt` generation.

## Current implementation status

The workflow currently reconstructs and validates:

- all 24 selected OhmNet tissue graph blocks;
- the complete 56,944-row GraphSAGE node order and Entrez identities;
- all 818,716 logical GraphSAGE edge records;
- all 50 MSigDB input features from public v6.1 C1/C3 files;
- all 121 GO labels from dated GOA, GeneID-UniProt, and GO sources;
- target-independent counts and content hashes;
- independent equality checks against released GraphSAGE topology, split flags,
  features, and labels.

Final GraphSAGE packaging and DGL conversion remain later stages.

## The non-circularity rule

Reconstruction and validation are separate commands:

```text
upstream sources -> reconstruction -> rebuilt data
                                      |
released targets --------------------+-> validation only
```

`pixi run reproduce` does not read `graphsage_ppi.zip` or `dgl_ppi.zip` and must
finish when those reference archives are absent. Only `pixi run validate` is
permitted to read the GraphSAGE reference. The DGL reference will be activated
when DGL reconstruction is implemented.

## Directory layout

```text
reproduction/
├── README.md                 # operational entry point
├── METHODS.md                # scientific transformations and evidence levels
├── THIRD_PARTY_NOTICES.md    # data licenses and attribution
├── Snakefile                 # one readable dependency graph
├── pixi.toml                 # environment and task definitions
├── pixi.lock                 # committed environment lock
├── pyproject.toml            # installable Python package
├── spec/                     # accepted human-readable specifications
├── src/graphsage_ppi_repro/  # scientific implementation
├── tests/                    # focused synthetic and integration tests
├── data/                     # immutable cached inputs; ignored by Git
├── build/                    # generated intermediates; ignored by Git
└── results/                  # generated manifests and validation; ignored by Git
```

The package deliberately avoids separate `config/`, `resources/`, `scripts/`,
and `docs/` trees. Accepted specifications live together in `spec/`; importable
and tested code lives under `src/`; explanations are consolidated here and in
`METHODS.md`.

## Prerequisites and commands

Install [Pixi](https://pixi.sh/), then run from `ai/reproduction/`:

```bash
pixi install
pixi run test
pixi run reproduce
pixi run validate
```

The package is tested on Linux and native Windows with Git Bash. Pixi tasks are
the supported interface. Direct Snakemake invocation is useful for debugging
but assumes all required cached files already exist.

### Reconstruct without released targets

```bash
pixi run reproduce
```

Pixi first acquires or verifies seven checksum-locked upstream files. Snakemake
then runs:

```text
verify immutable upstream sources
        |
reconstruct selected OhmNet graphs and biological row order
        |
        +------------------------+
        |                        |
reconstruct 50 features   reconstruct 121 GO labels
        |                        |
        +------------+-----------+
                     |
check target-independent invariants
                     |
write reconstruction manifest
```

Principal outputs are:

```text
build/topology/node_to_entrez.tsv.gz
build/topology/edges.tsv.gz
build/topology/summary.json
build/features/ppi-feats.npy
build/features/selected_features.tsv
build/features/summary.json
build/labels/ppi-labels.npy
build/labels/ppi-class_map.json
build/labels/selected_labels.tsv
build/labels/summary.json
build/reconstruction_checks.json
results/reconstruction_manifest.json
```

### Validate against released GraphSAGE data

```bash
pixi run validate
```

This independently verifies:

- consecutive node IDs, identity ID-map semantics, and split flags;
- exact graph-wise undirected edge multisets;
- exact feature shape, dtype, and all 2,847,200 feature cells;
- byte equality of the individual `ppi-feats.npy` member;
- exact class-map semantics and all 6,890,224 GO-label cells.

Validation writes:

```text
results/validation.json
results/validation.md
results/label_validation.json
results/label_validation.md
```

ZIP timestamps and compression bytes are not scientific outputs. The optional
stochastic `ppi-walks.txt` member is explicitly excluded.

### Run fast tests

```bash
pixi run test
```

The tests use small committed fixtures. They cover source-cache safety, legacy
CPython 2.7 ordering, graph reconstruction, MSigDB selection, historical NumPy
serialization, many-to-many identifier resolution, GAF filtering, alternate GO
IDs, `is_a` propagation, label projection, and target-independent versus
reference-based validation.

### Remove generated files

```bash
pixi run clean
```

Only `build/` and `results/` are removed. Files under `data/` are preserved.

## Inputs and cache

External inputs are stored in the package-local ignored directory:

```text
reproduction/data/
```

The workflow does not read or write the sibling `ai/data/` scratch directory.

`spec/sources.tsv` records every source's role, filename, byte size, SHA-256,
primary URL, optional mirror, license note, and description. Existing valid
files are verified in place without changing bytes or timestamps. Invalid
existing files cause a hard failure. Missing public files are downloaded to a
`.part` path and installed atomically only after all checks pass.

Pixi performs acquisition before Snakemake starts. Snakemake treats cached files
only as inputs; all rule outputs are confined to `build/` and `results/`.

Current upstream reconstruction inputs are:

- the OhmNet tissue-network archive;
- MSigDB v6.1 C1 and C3 Entrez GMT files;
- GOA human release-159 GAF and GPI files;
- the 2016-06-01 GeneID-UniProt mapping;
- the 2016-06-01 GO ontology.

Validation additionally uses the released GraphSAGE PPI ZIP.

## Specification files

- `sources.tsv` records external file identities and acquisition roles.
- `specification.yaml` records the accepted global transformation policy and
  compact expected invariants.
- `selected_graphs.tsv` records the 24 released tissue identities, graph order,
  and deposited 20/2/2 split. The historical selection algorithm remains open.
- `feature_columns.tsv` records the expected 30 C1 and 20 C3 results. The code
  derives them from the GMT files before checking this table.
- `label_columns.tsv` records the recovered 121-column GO order and explicitly
  marks three membership-indistinguishable term pairs as provisional.
- `identifier_decisions.tsv` records the small set of evidence-backed amendments
  needed beyond direct historical GeneID-UniProt edges.

These tables are reviewable reconstruction specifications, not hidden copies of
output matrices.

## Applying patches

Patches are generated relative to the Git repository root and contain paths
beginning with `ai/`. When a patch file is placed in the repository root:

```bash
ROOT=$(git rev-parse --show-toplevel)
PATCH="$ROOT/<patch-filename>.patch"

git -C "$ROOT" apply --check "$PATCH" &&
git -C "$ROOT" apply "$PATCH"
```

## Evidence language

The project distinguishes:

- **data-level exact**: regenerated and compared exactly;
- **strongly inferred**: an implementation mechanism explains all observations
  but is not documented in original preprocessing source;
- **documented**: stated by an official source;
- **open**: the evidence does not distinguish plausible histories.

The complete row mapping and reconstructed matrices are data-level exact under
the implemented workflow. Historical preprocessing mechanisms and the ordering
of three duplicate label-vector pairs remain explicitly qualified.

See [METHODS.md](METHODS.md) for the technical explanation.

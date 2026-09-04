# GraphSAGE PPI Core Reproduction Repository Outline

**Planning snapshot:** 2026-09-01 11:46:42 America/Chicago (2026-09-01T16:46:42Z)

**Status:** This document supersedes the broader planning outline dated 2026-09-01T16:29:03Z. It is intentionally limited to the clean, deterministic reproduction of the GraphSAGE PPI dataset and its downstream DGL representation.

## Purpose

The `reproduction/` subdirectory will contain a small, readable, and testable workflow that starts from independently acquired external source files and rebuilds the deterministic data products used by supervised GraphSAGE and DGL.

The workflow is intended for:

- the project investigators;
- an independent agent or human reviewer;
- external users who want to rerun the reconstruction;
- continuous integration after the local workflow is stable.

The package is not intended to preserve every exploratory analysis performed during the investigation. Historical searches, alternative hypotheses, literature analysis, and leakage experiments will remain elsewhere in the repository.

## Reproduction scope

The core workflow will:

1. Download public external sources and validate manually prepared restricted sources.
2. Verify file sizes, checksums, archive structure, and internal release metadata.
3. Reconstruct the selected OhmNet tissue graphs and the GraphSAGE node order.
4. Reconstruct the 50 GraphSAGE node-feature columns from a minimal MSigDB input.
5. Reconstruct the 121 GraphSAGE GO-label columns from historical GOA and mapping sources.
6. Assemble the deterministic GraphSAGE files required by supervised training.
7. Transform the reconstructed GraphSAGE data into the DGL representation.
8. Compare reconstructed intermediate and final artifacts with independently downloaded GraphSAGE and DGL reference files.
9. Write a run manifest, checksums, and concise machine- and human-readable validation summaries.

The core workflow will not rerun:

- historical GOA release screens;
- alternative `gene2go`, Bioconductor, or dhimmel annotation analyses;
- broad MSigDB version comparisons;
- label-column-order model searches;
- literature discrepancy analysis;
- dataset-impact searches;
- leakage experiments;
- PPI source-attribution audits beyond the inputs needed for reconstruction;
- stochastic unsupervised random-walk generation.

The deterministic reproduction target is therefore:

> Reconstruct every deterministic artifact needed by supervised GraphSAGE and DGL. Auxiliary stochastic files that are not used by that path are outside scope.

## Main technology choices

- **Environment:** Pixi, with `pixi.lock` as the only canonical environment lock.
- **Workflow:** One readable Snakemake `Snakefile` initially.
- **Scientific implementation:** A small installable Python package under `src/`.
- **Testing:** pytest for unit tests and full regression checks.
- **External-source inventory:** A flat `config/sources.tsv`.
- **Scientific reconstruction policy:** A compact `config/reconstruction.yaml`.
- **Generated intermediates:** `build/`, ignored by Git.
- **Generated final outputs:** `results/`, ignored by Git.
- **Canonical execution platform:** Linux x86-64; WSL2 is the recommended Windows route.

There will be no parallel `environment.yml`, Snakemake profiles, multipurpose audit framework, or large user-facing CLI at the start.

## Relationship to the repository root

```text
GraphSagePPIAnalysis/
|-- README.md
|-- data/                 # External source and reference files; ignored
|-- papers/               # Local manuscript copies; ignored
|-- reproduction/         # Clean deterministic reproduction described here
|-- reports/              # Scientific findings and narrative reports
|-- claims/               # Future machine-readable claims and evidence table
|-- gpt56sol/             # Investigation history and GPT-5.6 artifacts
`-- opus46/               # Independent-agent investigation artifacts
```

The `reproduction/` package should explain how to obtain inputs, rebuild outputs, and verify equality. Broader scientific interpretation belongs under repository-level `reports/` and `claims/`.

## Proposed `reproduction/` layout

```text
reproduction/
|-- README.md
|-- pixi.toml
|-- pixi.lock
|-- pyproject.toml
|-- Snakefile
|
|-- config/
|   |-- sources.tsv
|   `-- reconstruction.yaml
|
|-- resources/
|   |-- recovered/
|   |   |-- README.md
|   |   |-- selected_graphs.tsv
|   |   |-- label_columns.tsv
|   |   `-- identifier_decisions.tsv
|   |
|   `-- expected/
|       |-- expected.yaml
|       `-- feature_columns.tsv
|
|-- scripts/
|   `-- prepare_msigdb_input.py
|
|-- src/
|   `-- graphsage_ppi_repro/
|       |-- __init__.py
|       |-- sources.py
|       |-- manifest.py
|       |-- legacy_dict.py
|       |-- topology.py
|       |-- features.py
|       |-- labels.py
|       |-- graphsage.py
|       |-- dgl.py
|       `-- validate.py
|
|-- tests/
|   |-- data/
|   |-- conftest.py
|   |-- test_sources.py
|   |-- test_legacy_dict.py
|   |-- test_topology.py
|   |-- test_features.py
|   |-- test_labels.py
|   |-- test_graphsage.py
|   `-- test_dgl.py
|
|-- docs/
|   |-- architecture.md
|   |-- data-sources.md
|   |-- running-the-workflow.md
|   |-- recovered-decisions.md
|   `-- limitations.md
|
|-- build/                # Generated intermediates; ignored
`-- results/              # Generated outputs and validation; ignored
```

## Top-level files

### `README.md`

The user entry point. It should state:

- exactly what the workflow reconstructs;
- installation prerequisites;
- how to obtain or prepare each input;
- how to run the workflow and tests;
- where generated outputs appear;
- what equality guarantees are checked;
- which historical details remain unresolved but do not prevent reproduction.

### `pixi.toml`

Defines the software environment and a very small task interface:

```text
pixi run reproduce
pixi run test
pixi run clean
```

`reproduce` runs the complete Snakemake workflow. `test` runs pytest. `clean` removes generated `build/` and `results/` files but never removes downloaded source data.

### `pixi.lock`

The sole canonical environment lock. It records exact package builds for supported platforms.

### `pyproject.toml`

Defines the installable Python package and development-tool settings. Pixi remains responsible for the canonical environment.

### `Snakefile`

Describes the dependency graph and calls the Python modules. It should remain declarative and contain little or no scientific logic.

Likely rules are:

```text
acquire_sources
prepare_msigdb
reconstruct_node_mapping
reconstruct_topology
reconstruct_features
reconstruct_labels
assemble_graphsage
validate_graphsage
assemble_dgl
validate_dgl
write_summary
```

The workflow should remain in one Snakefile until that file becomes genuinely difficult to review.

## `config/`

Contains the two small, human-maintained files that define source acquisition and the canonical reconstruction policy.

### `config/sources.tsv`

A flat inventory of only the external files needed by the core workflow.

Proposed columns:

```text
source_id
role
filename
url
archive_url
sha256
size_bytes
access
description
```

The `role` field distinguishes upstream sources from released reference targets. The `access` field can use a small vocabulary such as:

```text
public
prepared_manual
```

Core entries are expected to include:

- the OhmNet tissue-network archive;
- the OhmNet README when needed to document file semantics;
- the released GraphSAGE PPI archive as a reference target;
- the released DGL PPI archive as a reference target;
- GOA human release-159 GAF and GPI files;
- the historical GeneID-to-UniProt mapping files;
- the historical GO ontology;
- the prepared minimal MSigDB input.

The workflow must fail if downloaded bytes do not match the expected size and SHA-256. It must never update expected hashes automatically.

### `config/reconstruction.yaml`

Records the scientific transformation policy used by the canonical workflow. It should remain short and readable.

Expected sections are:

```yaml
node_order:
  method: cpython27_string_dict
  word_size_bits: 64
  hash_randomization: false

features:
  collections: [C1, C3]
  minimum_source_members: 200
  maximum_columns: 50

labels:
  evidence_codes: [EXP, IDA, IEP, IGI, IMP, ISS]
  exclude_not: true
  biological_process_relation: involved_in
  cellular_component_relation: part_of
  molecular_function_relation: enables
  excluded_relations: [colocalizes_with, contributes_to]
  propagation_relations: [is_a]

resources:
  selected_graphs: resources/recovered/selected_graphs.tsv
  label_columns: resources/recovered/label_columns.tsv
  identifier_decisions: resources/recovered/identifier_decisions.tsv
```

Machine-specific paths, CPU counts, download URLs, and expected output hashes do not belong in this file.

## `resources/`

Contains small committed tables and expected values. These are part of the reproducible specification rather than generated outputs.

### `resources/recovered/`

Stores evidence-backed decisions recovered during the forensic investigation.

#### `selected_graphs.tsv`

Defines:

- the 24 OhmNet tissue graphs used by GraphSAGE;
- their graph order;
- their train, validation, or test assignment;
- the source edgelist associated with each graph.

The graph identities and assignments are known from the released data. The original rule used to select those 24 graphs from the larger OhmNet archive remains unknown, so the recovered table is explicit.

#### `label_columns.tsv`

Defines the 121 GraphSAGE label columns. It should include:

- column index;
- primary GO ID;
- GO namespace;
- alternative GO ID when another term has the same membership vector;
- an identity status such as `uniquely_determined` or `strongly_supported_provisional`.

The table makes the six duplicate-vector column identities visible rather than hiding them in code.

#### `identifier_decisions.tsv`

Documents identifier-mapping components that cannot safely be represented as a one-to-one dictionary. It records:

- all retained mappings in the relevant component;
- the selected annotation projection;
- the evidence for that decision;
- whether the decision is exact, strongly supported, or unresolved.

No special mapping decision should exist only as an unexplained Python conditional.

### `resources/expected/`

Stores compact values used to verify that the workflow still reproduces the known dataset.

#### `expected.yaml`

Contains stable invariants such as:

- graph, node, and link counts;
- feature and label matrix shapes;
- split counts;
- distinct Entrez GeneID count;
- canonical matrix or table hashes;
- exact comparison requirements for GraphSAGE and DGL outputs.

External source hashes remain in `sources.tsv`; they should not be duplicated here.

#### `feature_columns.tsv`

Lists the expected 50 feature columns, their MSigDB collection, source name, and output position. The reconstruction code should derive this list from the configured rule; this table is an expected result used for validation, not an input that bypasses feature selection.

## `scripts/`

Contains a small preparation utility for the one source that cannot be downloaded unattended from a public URL.

### `scripts/prepare_msigdb_input.py`

A licensed user supplies a privately downloaded MSigDB archive. The script writes the minimal, deterministic source needed by the reproduction:

- original-order C1 Entrez GMT rows;
- original-order C3 Entrez GMT rows;
- source archive version and SHA-256;
- output SHA-256 and row counts;
- a small provenance record.

The Snakemake workflow validates this prepared input but does not require MSigDB login credentials. Public redistribution of the prepared input should occur only after its licensing has been reviewed and permission confirmed where necessary.

## `src/graphsage_ppi_repro/`

Contains the scientific implementation. Modules are separated by responsibility so that acquisition, reconstruction, and target comparison remain distinct.

### `sources.py`

Reads `sources.tsv`, downloads public files, validates manually prepared files, checks file signatures, verifies sizes and hashes, and writes source-acquisition records.

### `manifest.py`

Writes a run manifest containing:

- Git revision and dirty-state flag;
- Pixi lock hash;
- software versions;
- command line;
- source hashes;
- generated artifact hashes;
- start and end times;
- validation outcomes.

### `legacy_dict.py`

Implements the required 64-bit, unrandomized CPython 2.7 string-dictionary behavior in modern Python. This module should be heavily commented and tested independently.

### `topology.py`

Reads the selected OhmNet edgelists, applies the recovered node ordering and graph assignments, and constructs the canonical 24-graph topology and node-to-Entrez table.

### `features.py`

Parses the prepared original-order MSigDB source, applies the configured source-level membership threshold and 50-column cap, projects memberships onto the ordered GraphSAGE rows, and writes the feature matrix and derived feature-column table.

### `labels.py`

Performs the complete label reconstruction:

1. Parse GOA GAF and GPI files, historical mapping files, and the GO ontology.
2. Preserve many-to-many identifier mappings.
3. Apply the documented component decisions.
4. Retain the configured evidence codes.
5. Exclude `NOT` annotations.
6. Retain only the configured ordinary relations.
7. Exclude qualified relations from ordinary binary membership.
8. Canonicalize alternate GO IDs.
9. Propagate through transitive `is_a` ancestors only.
10. Apply the recovered label-column specification.
11. Project labels onto all ordered GraphSAGE rows.

### `graphsage.py`

Writes the deterministic GraphSAGE files needed by supervised training and downstream conversion, including graph structure, ID map, class map, and feature array.

### `dgl.py`

Transforms the reconstructed GraphSAGE representation into the DGL graph, feature, label, graph-ID, split, directed-edge, and self-loop representation.

### `validate.py`

Compares reconstructed outputs with independently downloaded GraphSAGE and DGL targets. It is the only scientific module that should need access to those targets.

Comparisons should distinguish:

```text
byte_exact
array_exact
structurally_exact
numeric_tolerance
```

The reconstruction modules must not read reference targets.

## `tests/`

Contains fast tests for difficult transformations plus full regression checks.

### `tests/data/`

Contains tiny synthetic inputs for:

- GAF evidence, negation, and relation handling;
- OBO `is_a`, `part_of`, and alternate-ID behavior;
- many-to-many identifier mapping;
- CPython 2 dictionary collisions and resizing;
- source-order MSigDB filtering;
- repeated GeneIDs across graph splits;
- small GraphSAGE and DGL examples.

### `tests/conftest.py`

Defines reusable pytest fixtures, temporary directories, and test-data paths.

### Test modules

- `test_sources.py`: manifest parsing, downloading, hash verification, and invalid-content rejection.
- `test_legacy_dict.py`: Python 2 hashing, probing, resizing, and table iteration.
- `test_topology.py`: node order, graph construction, edge normalization, and split assignment.
- `test_features.py`: source-order selection, threshold behavior, column limit, and exact expected columns.
- `test_labels.py`: evidence filtering, relation handling, `is_a` closure, alternate IDs, and ambiguous mapping components.
- `test_graphsage.py`: deterministic GraphSAGE artifact structure and contents.
- `test_dgl.py`: scaling, graph grouping, directed edges, self-loops, dtypes, and split outputs.

The complete Snakemake run supplies the large end-to-end regression test against the released reference datasets.

## `docs/`

Contains only technical documentation for the reproduction package.

### `architecture.md`

Explains module boundaries and the separation between source acquisition, reconstruction, and target validation.

### `data-sources.md`

Documents each core source, its version, access method, licensing status, role in the reconstruction, and any required manual preparation.

### `running-the-workflow.md`

Explains Pixi installation, source preparation, workflow commands, expected disk use, generated outputs, and common failures.

### `recovered-decisions.md`

Explains why graph selection, label-column identities, and identifier-mapping decisions are committed resource tables instead of being rediscovered from the target data during every run.

### `limitations.md`

Records only uncertainties relevant to interpreting the reproduction, such as:

- the historical MSigDB version is not identifiable from the feature matrix;
- `>= 200` and `> 200` are observationally equivalent for the relevant source sets;
- the original rule that selected the 24 tissue graphs is unknown;
- six GO column identities remain strongly supported rather than directly documented;
- GOA v159 may be the literal source or may reproduce an equivalent historical intermediate;
- exact DGL container-byte equality may be unnecessary if logical content is exact.

Broader scientific interpretation remains outside this package.

## Generated directories

### `build/`

Ignored directory for reproducible intermediate files, including:

- extracted archives;
- normalized source tables;
- canonical edge lists;
- intermediate identifier maps;
- propagated GO memberships;
- reconstructed but not yet packaged artifacts.

Nothing in `build/` should be required as a committed input.

### `results/`

Ignored directory for final outputs from one run, including:

- complete node-to-Entrez mapping;
- reconstructed feature and label metadata;
- reconstructed deterministic GraphSAGE files;
- reconstructed DGL files;
- source and output checksums;
- run manifest;
- machine-readable validation JSON;
- concise human-readable validation summary.

Frozen results from formal releases can be attached to a GitHub Release or deposited in an external archive rather than committed repeatedly to Git history.

## Core dependency graph

```text
External sources
|
|-- OhmNet tissue edgelists
|      `-- selected graph specification + CPython 2 string-dictionary order
|              `-- node mapping and GraphSAGE topology
|
|-- Prepared MSigDB C1/C3 input
|      `-- source-order membership filtering + 50-column cap
|              `-- GraphSAGE features
|
|-- GOA v159 + GPI + GeneID/UniProt mapping + GO ontology
|      `-- identifier mapping + annotation filtering + is_a propagation
|              `-- GraphSAGE labels
|
`-- Recovered graph and label-column resources
       `-- deterministic GraphSAGE artifacts
               |-- validation against released GraphSAGE target
               `-- DGL transformation
                       `-- validation against released DGL target
```

## Non-circularity requirement

The reconstruction must complete when the downloaded GraphSAGE and DGL reference targets are unavailable. Only validation may read those targets.

This is the main safeguard against incorporating desired output values into the reconstruction while claiming to regenerate them.

## Initial implementation order

1. Create the Pixi project, Python package skeleton, source table, and minimal Snakefile.
2. Implement source acquisition, checksum verification, and run manifests.
3. Implement and test CPython 2.7 string-dictionary ordering.
4. Reconstruct node identities and topology.
5. Implement the private-to-minimal MSigDB preparation script and reconstruct features.
6. Reconstruct GO labels.
7. Assemble deterministic GraphSAGE artifacts and compare them with the released target.
8. Reconstruct DGL artifacts from the rebuilt GraphSAGE data and compare them with the released target.
9. Produce a concise validation summary and output checksum inventory.
10. Add GitHub Actions only after the complete local workflow succeeds from a clean checkout.

## Success criterion

The reproduction is successful when a clean environment can acquire or validate all required inputs, rebuild the deterministic supervised GraphSAGE PPI data products and DGL representation without reading the target outputs, and then pass the configured exact or explicitly tolerated comparisons against the independently acquired released targets.

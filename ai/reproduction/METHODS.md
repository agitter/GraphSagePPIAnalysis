# Reconstruction methods

## 1. Purpose and boundaries

The released GraphSAGE PPI dataset contains 24 tissue-specific graph components,
56,944 tissue-instance rows, 818,716 graph-link records, 50 binary input
features, and 121 binary GO labels. A row represents one occurrence of a gene in
one tissue graph; the same biological gene can occur in several graphs.

This package reconstructs the deterministic dataset used by supervised
GraphSAGE and will derive the DGL PPI representation. Exploratory source
screens, alternative hypotheses, split analyses, leakage experiments,
literature discrepancies, and stochastic walks remain elsewhere in the
repository.

The principal design constraint is non-circularity. Reconstruction may use only
upstream biological sources and committed evidence-backed specifications. Only
`graphsage_ppi_repro.validate` may open released GraphSAGE or DGL targets.

## 2. Evidence-backed specifications

Some facts are recoverable from the deposited dataset but not uniquely from
public prose. They are committed as small tables under `spec/` so that they are
visible and reviewable instead of hidden in Python constants.

### 2.1 Selected graph blocks

`selected_graphs.tsv` records the exact OhmNet tissue file corresponding to each
GraphSAGE component, its component order, and its train, validation, or test
role. These identities and assignments are data-level exact. The historical
rule that selected the 24 layers remains open, so the workflow reproduces the
released selection instead of claiming that manuscript thresholds rediscover it.

### 2.2 Feature and label columns

`feature_columns.tsv` records the expected feature order, source row, membership
hash, source membership count, and projected row count. The implementation still
parses the complete GMT sources and derives the columns before comparison.

`label_columns.tsv` records the 121 GO terms in released column order. Membership
uniquely identifies 115 columns. Six columns form three indistinguishable pairs;
their within-pair orientation is marked strongly supported and provisional.

### 2.3 Identifier decisions

`identifier_decisions.tsv` records 15 explicit include or exclude amendments to
direct historical UniProt-GeneID edges. These rows represent 13 protein
accessions in graph-relevant ambiguous or missing-edge cases. Each row includes
its rationale, evidence status, and repository evidence reference. No such
amendment exists only as a Python conditional.

## 3. Source acquisition and provenance

`spec/sources.tsv` is a flat checksum-locked inventory. An `upstream` source may
influence reconstruction; a `reference` source may be read only by validation.

Acquisition follows this sequence:

1. Verify and reuse an existing valid cache file without changing it.
2. Reject an existing invalid file without replacing it.
3. Download a missing public source to an adjacent temporary path.
4. Reject HTML error or login pages masquerading as data.
5. Verify exact size, SHA-256, and a lightweight format check.
6. Atomically install the verified file at its final path.

Pixi performs acquisition before Snakemake starts. Snakemake lists cached source
files only as inputs and confines generated outputs to `build/` and `results/`.
A future mirror may be added to `mirror_url`, but it must serve bytes with the
same expected SHA-256 rather than a different source version.

The reconstruction manifest records source verification, specification and
artifact hashes, software and platform details, and Git state.

## 4. Topology and biological row order

### 4.1 Immediate graph source

The selected components are reconstructed from OhmNet tissue-specific edgelists,
whose nodes are Entrez GeneIDs. For each selected member the workflow verifies
the decompressed SHA-256, retains source line order and self-loops, and rejects
duplicate undirected records.

### 4.2 Recovered legacy ordering mechanism

The anonymous GraphSAGE row order is reconstructed by one global mechanism:

1. Read each selected edgelist in source line order.
2. Retain Entrez identifiers as ASCII strings.
3. Insert both endpoints as edges are encountered.
4. Model a 64-bit, unrandomized CPython 2.7 dictionary.
5. Iterate occupied hash-table slots in table order.
6. Concatenate the 24 local sequences in frozen component order.

`legacy_order.py` implements only the required string hash, probing, resize,
reinsertion, lookup, and table-iteration behavior. A separate positive control
reproduces the irregular serialized key order of the released class map.

The complete 56,944-row mapping is data-level exact under this mechanism. The
claim that the original authors literally followed this preprocessing path is
strongly inferred because the original preprocessing source was not located.

### 4.3 Outputs

The topology stage writes a row-wise node-to-Entrez table, a source-provenance
edge table, and a summary containing graph boundaries, counts, archive-member
hashes, and stable content hashes. Compressed TSV outputs use deterministic gzip
metadata.

## 5. Feature reconstruction

### 5.1 Canonical source

The package uses public MSigDB v6.1 Entrez GMT files for C1 and C3. A separate
acceptance analysis established that v6.1 yields the same ordered 50 membership
vectors as tested versions 5.0 through 6.0 and exactly reproduces the released
matrix. Version 6.1 corrects the selected C3 name `AAAYWAACM_HFH4_01` to
`AACTTT_UNKNOWN` without changing membership.

This makes v6.1 a corrected canonical reproduction source; it does not identify
the historical version used by GraphSAGE.

### 5.2 Selection and projection

The accepted rule reads C1 and then C3 in deposited row order, retains sets with
at least 200 distinct source Entrez IDs, stops at a global maximum of 50 columns,
and only then projects memberships onto the reconstructed rows. It yields 30 C1
and 20 C3 columns.

The code uses `>= 200`; `> 200` is observationally equivalent for the selected
prefix. Selection before projection explains why `chryq11` is retained despite
becoming an all-zero GraphSAGE column.

### 5.3 Matrix representation

The 56,944 by 50 feature matrix is float64. It is written with the NumPy v1.0
header alignment found in the released member, permitting both array equality
and byte equality for `ppi-feats.npy`.

## 6. Identifier mapping and GO labels

### 6.1 Dated inputs

Label reconstruction uses:

- GOA human release-159 GAF;
- GOA human release-159 GPI;
- `2016-06-01-gp2protein.geneid.gz`;
- the 2016-06-01 GO OBO ontology.

The GPI file defines the relevant UniProt accessions. The workflow reads all
direct GeneID edges for those accessions from the historical mapping. It does
not collapse a many-to-many mapping by taking a first or arbitrary match.

### 6.2 Evidence-backed mapping amendments

Most accessions use their direct historical GeneID edges unchanged. The 15
rows in `identifier_decisions.tsv` add missing uniquely supported assignments or
remove cross-assignments in previously audited ambiguous components. This table
includes the O95073 case: O95073-to-25788 was a real historical cross-reference,
but FSBP annotations are not projected onto the GraphSAGE RAD54B node represented
by GeneID 25788.

The central package consumes the accepted decision table; it does not rerun the
broader forensic component and symbol searches that produced it, and it never
fits decisions against the released label matrix during reconstruction.

### 6.3 Annotation filtering and ontology propagation

One global policy is applied to every GAF row:

- retain evidence codes `EXP`, `IDA`, `IEP`, `IGI`, `IMP`, and `ISS`;
- exclude annotations qualified by `NOT`;
- use ordinary aspect relations `involved_in`, `part_of`, and `enables`;
- do not treat `colocalizes_with` or `contributes_to` as ordinary membership;
- replace alternate GO IDs with their primary IDs;
- propagate to the annotated term and transitive `is_a` ancestors only;
- do not propagate ontology `part_of` edges.

Term prevalence is counted once per mapped historical human GeneID. The top 121
terms must equal the set in `label_columns.tsv`; the table then supplies released
column order. The resulting namespaces are 85 Biological Process, 26 Cellular
Component, and 10 Molecular Function columns.

There are 121 columns but 118 distinct membership vectors. The three duplicate
pairs remain explicit in the specification rather than being silently assigned
names by output matching.

### 6.4 Label outputs

The label stage writes:

- a 56,944 by 121 uint8 NumPy matrix;
- a numeric-node-ordered logical GraphSAGE class map;
- a selected-term table with names, prevalence ranks, and counts;
- a summary with source hashes, filter counts, mapping decisions, and output
  hashes.

GeneID 10159 is absent from the accepted GPI-based mapping and receives an all-
zero vector, matching the released data. Repeated occurrences of every mapped
gene are checked for identical vectors.

## 7. GraphSAGE artifact assembly

The transparent intermediate tables and matrices are assembled into the four
files read by the supervised GraphSAGE loader:

```text
results/graphsage/ppi/ppi-G.json
results/graphsage/ppi/ppi-id_map.json
results/graphsage/ppi/ppi-class_map.json
results/graphsage/ppi/ppi-feats.npy
```

The graph contains nodes in ascending GraphSAGE ID order with the reconstructed
validation and test flags. Each undirected edge is written once as `(min, max)`,
and those pairs are sorted lexicographically. This is a simple canonical order
that is independent of dictionary iteration. It reproduces the released graph
structure exactly but does not claim to reproduce the historical ordering of
the JSON `links` list.

The identity map is written in ascending numeric node order. The class map uses
the independently recovered 64-bit CPython 2.7 string-dictionary key order;
this reproduces the released class-map bytes exactly. The feature NPY is copied
without reserialization because the upstream feature stage already reproduces
that file byte-for-byte. The final directory contains no `ppi-walks.txt`.

A checksum file at `results/graphsage/SHA256SUMS` records paths relative to its
own directory and can be checked with ordinary `sha256sum -c`.

## 8. Target-independent checks

`check-milestone` compares generated summaries with compact expectations in
`specification.yaml`. It verifies topology counts and hashes, feature dimensions
and hashes, label dimensions and positive-cell count, namespace counts, the 118
distinct label vectors, the one unmapped graph GeneID, exact term-set recovery,
the complete label-matrix data hash, and the canonical GraphSAGE artifact
hashes.

These checks use committed expectations but never open a released reference
archive. They catch regressions while preserving the reconstruction boundary.

## 9. Independent GraphSAGE validation

The GraphSAGE ZIP is acquired as a `reference` source only for
`pixi run validate`. Validation independently parses the four reconstructed
files and the corresponding released members and checks:

- exact graph metadata, nodes, split flags, and undirected edge multiset;
- a unique, sorted canonical edge list in the reconstructed graph;
- exact ID-map semantics and bytes;
- exact class-map values, all 6,890,224 label cells, and bytes;
- exact feature shape, dtype, all 2,847,200 cells, and NPY bytes;
- agreement of row counts across all four files;
- absence of the optional unsupervised walk file.

The released `ppi-G.json` stores the same undirected edges in a historical
NetworkX/Python-dependent order. Link-array order does not change the logical
node-link graph, so graph validation is structural. The validation report still
records graph byte equality explicitly; it is expected to be false under the
canonical ordering. The ID map, class map, and feature file are expected to be
byte-identical.

## 10. Remaining implementation stage

DGL outputs will next be derived solely from reconstructed GraphSAGE data. The
transformation will include graph grouping, training-only feature
standardization, directed edge expansion, one self-loop per node, labels, and
split packaging. DGL validation will require logical and data equivalence rather
than identical framework-container bytes.

## 11. Quality control, attribution, and uncertainty

- Scientific rules are centralized in `specification.yaml` and explained here.
- Recovered identities and exceptional decisions are visible in TSV files.
- Source hashes are never learned or updated automatically.
- Reconstruction and target validation are separate modules and workflow paths.
- Difficult historical behavior has focused synthetic tests and positive
  controls.
- No command rewrites expected results to accept current output.

The reconstruction builds on the human investigators, earlier Claude/Opus
analysis, and subsequent GPT-5.6 Sol analysis. The investigators framed the
problem, supplied and audited sources, required cautious biological mapping, and
corrected overclaims. Earlier Opus work established major parts of tissue
matching, structural node alignment, feature-family discovery, and the initial
DGL transformation. GPT-5.6 Sol independently extended and validated the work,
including the complete row-order mechanism, feature selection rule,
qualifier-aware GO policy, and full-row reconstruction.

The executable package records accepted results and evidence levels. Detailed
intellectual history and literature discrepancies remain in repository-level
reports rather than being duplicated here.

# Reconstruction methods

## 1. Purpose and boundaries

The released GraphSAGE PPI dataset contains 24 tissue-specific graph components,
56,944 tissue-instance rows, 818,716 graph-link records, 50 binary input
features, and 121 binary GO labels. The same biological gene can occur in more
than one tissue graph; a row is therefore a tissue-specific occurrence rather
than a globally unique protein.

This package rebuilds deterministic data products used by supervised GraphSAGE
and the DGL PPI derivative. It does not preserve the exploratory search process
that led to the accepted rules. Those investigations, alternative hypotheses,
and literature discrepancies remain elsewhere in the repository.

The principal design constraint is non-circularity: reconstruction code may use
only upstream biological sources and committed, evidence-backed specifications.
Only `graphsage_ppi_repro.validate` may open the released GraphSAGE and DGL
archives.

## 2. Evidence-backed specifications

Some facts are observable from the deposited data but cannot be rediscovered
from public prose alone. They are recorded in small tables under `spec/` so they
remain visible and reviewable rather than hidden in Python constants.

### 2.1 Selected graph blocks

`selected_graphs.tsv` records the exact OhmNet tissue file corresponding to each
of the 24 GraphSAGE graph components, its deposited component order, and its
training, validation, or test role.

The identities and assignments are data-level exact. The historical rule that
selected these 24 layers from the larger OhmNet release is open. The central
workflow therefore reproduces the released selection instead of pretending
that the manuscript's reported edge thresholds uniquely regenerate it.

### 2.2 Expected feature columns

`feature_columns.tsv` records the derived feature order, source row, membership
hash, source member count, and projected row count. It is an expected-result
table. The implementation still parses the complete C1 and C3 GMT files and
applies the global rule before comparing its derived columns with the table.

### 2.3 Label columns and identifier decisions

`label_columns.tsv` and `identifier_decisions.tsv` support the later GO-label
milestone. They expose recovered column identities, duplicate-vector ambiguity,
and historically complex identifier components that would otherwise become
unexplained special cases in code.

## 3. Source acquisition and provenance

`spec/sources.tsv` is a flat, checksum-locked inventory. Each record is either:

- `upstream`, meaning it may influence reconstruction; or
- `reference`, meaning validation alone may read it.

Acquisition uses the following sequence:

1. Reuse the local file when size, SHA-256, and structure are valid, without
   modifying the file or its timestamp.
2. Reject an existing invalid file without replacing it.
3. Download a missing public source to an adjacent temporary path.
4. Reject HTML error or login pages masquerading as data.
5. Verify exact size and SHA-256.
6. Perform a lightweight ZIP, TAR, gzip, GMT, or OBO check as appropriate.
7. Install the verified temporary file atomically at the final cache path.

Acquisition and reconstruction have separate ownership of files. A Pixi
preflight task ensures that the required immutable files exist under `data/`.
Snakemake lists those files only as inputs. It may generate source-verification
reports under `build/`, but no rule declares a cached source or reference archive
as an output. Consequently, cleaning or rerunning generated workflow products
cannot cause Snakemake to remove or update source files.

A future independent mirror can be added to the existing `mirror_url` column.
The mirror must serve the same checksum-identified bytes; it is not permitted to
be an alternate source version.

Each reconstruction run writes a manifest containing source verification,
specification hashes, generated artifact hashes, software and platform details,
and Git state.

## 4. Topology and biological row order

### 4.1 Immediate graph source

The selected graph components are reconstructed from OhmNet tissue-specific
edgelists. OhmNet represents nodes as Entrez GeneIDs, and the same ID in two
layers denotes the same gene.

For every selected archive member the implementation verifies the SHA-256 of
the decompressed edgelist before parsing it. It retains self-loops and source
line order and rejects duplicated undirected records.

### 4.2 Recovered legacy ordering mechanism

The anonymous GraphSAGE row order is reproduced by this global mechanism:

1. Read each selected OhmNet edgelist in source line order.
2. Retain Entrez identifiers as ASCII strings.
3. Insert both endpoints while edges are encountered.
4. Model an insertion-only, 64-bit, unrandomized CPython 2.7 dictionary.
5. Iterate occupied hash-table slots in table order.
6. Concatenate the 24 local row sequences in the frozen component order.

`legacy_order.py` implements only the required historical behavior: the Python
2 string hash, perturb probing, resizing, reinsertion during resize, lookup, and
occupied-slot iteration. It is not a general Python 2 emulator.

The resulting 56,944-row mapping agrees with all independently anchored rows and
exactly supports the released graph, feature, and label data. A separate
positive control inserts string keys `"0"` through `"56943"` and reproduces the
irregular serialized key order of the released class-map dictionary. These
facts make the mapping data-level exact under the mechanism, while the claim
that the original authors literally used this preprocessing path remains
strongly inferred because the original preprocessing source was not located.

### 4.3 Topology outputs

The milestone writes:

- a row-wise node-to-Entrez table;
- an edge table retaining source line and tissue provenance;
- graph boundaries, counts, source-member hashes, and stable content hashes.

Compressed TSV outputs use a fixed gzip timestamp and no embedded filename, so
identical logical rows produce identical gzip bytes.

## 5. Feature reconstruction

### 5.1 Canonical public source

The package uses MSigDB v6.1 Entrez GMT files for collections C1 and C3. A
separate acceptance analysis established that v6.1 produces the same ordered 50
membership vectors as tested releases 5.0 through 6.0 and reproduces every
released feature cell. Version 6.1 corrects the selected C3 name
`AAAYWAACM_HFH4_01` to `AACTTT_UNKNOWN` without changing that set's membership.

This makes v6.1 a corrected canonical reproduction source. It does not identify
which MSigDB version the GraphSAGE authors historically used.

### 5.2 Selection rule

The accepted global rule is:

1. Read C1 in deposited GMT order.
2. Keep sets with at least 200 distinct source Entrez IDs.
3. Continue with C3 in deposited GMT order.
4. Append qualifying sets until the global list reaches 50 columns.
5. Project the selected source memberships onto the reconstructed GraphSAGE row
   sequence.

This produces 30 C1 and 20 C3 columns. The code uses `>= 200`. The historical
choice between `>= 200` and `> 200` is not identifiable from the output because
both select the same first 50 columns in the tested releases.

Selection occurs before projection. The eleventh column in zero-based position
10, `chryq11`, has more than 200 source members but no members in the GraphSAGE
gene universe; it is consequently retained as an all-zero output column. This
is strong evidence for source-level filtering before projection and against a
later empty-column removal step.

### 5.3 Matrix representation

The reconstructed feature matrix is float64 with shape 56,944 by 50. It is
written using the NumPy v1.0 header form and 16-byte alignment found in the
released file. Modern NumPy normally uses a different header-padding convention;
the data array is unchanged, but reproducing the earlier header permits exact
byte comparison of `ppi-feats.npy` itself.

The package records independent hashes of the float64 data bytes, a uint8
binary representation, and the complete NPY file.

## 6. Target-independent checks

`check-milestone` compares generated summaries with compact expectations in
`specification.yaml`. It checks graph, row, gene, edge, and split counts; stable
topology content hashes; feature dimensions and collection counts; the all-zero
column position; and feature-data hashes.

These checks do not read the released GraphSAGE archive. They are useful for
catching regressions during reconstruction while preserving the non-circularity
boundary.

## 7. Independent target validation

The GraphSAGE reference ZIP is acquired as a `reference` source during the
Pixi validation preflight and appears only as an input to the validation branch
of the Snakemake graph. Snakemake writes a verification report under `build/`
but does not own the cached archive. For the current milestone, validation
independently parses the reconstructed and released representations and checks:

- node IDs and the identity ID map;
- per-row train/validation/test flags;
- graph-wise undirected edge multiplicities;
- feature shape, dtype, and exact values;
- byte equality of the deposited NPY member.

The topology comparison is structural because JSON object formatting and edge
record direction are not scientifically meaningful when the graph-wise
undirected multisets are exact. The feature NPY is checked both as data and as
bytes because its individual serialization is stable and reproducible.

The final GraphSAGE and DGL stages will extend this contract. For DGL, the
accepted success criterion is logical and data equivalence rather than identical
archive or framework-container bytes.

## 8. Later milestones

### 8.1 GO labels

The accepted label reconstruction will use GOA human release 159 GAF and GPI,
the historical 2016-06-01 GeneID-UniProt mapping, and the June 2016 GO ontology.
The global policy retains evidence codes EXP, IDA, IEP, IGI, IMP, and ISS;
excludes `NOT`; handles ordinary aspect relations separately from
`colocalizes_with` and `contributes_to`; canonicalizes alternate GO IDs; and
propagates only through transitive `is_a` ancestry.

Identifier mapping will preserve many-to-many components and apply the few
explicit decisions recorded in `identifier_decisions.tsv` before annotation
projection. The final matrix must reproduce all 56,944 by 121 binary cells.

### 8.2 GraphSAGE assembly

The package will write the deterministic graph JSON, identity map, class map,
and feature array needed by supervised GraphSAGE. The optional stochastic walk
file will remain excluded.

### 8.3 DGL derivation

DGL outputs will be derived from reconstructed GraphSAGE logical data, never
from the downloaded DGL target. The transformation will include graph grouping,
training-only feature standardization, directed edge expansion, one self-loop
per node, labels, and split packaging. Validation will compare graph membership,
row order, integer arrays, labels, directed edges, loops, dtypes, and feature
values under a fixed numerical tolerance.

## 9. Quality-control policy

- Scientific rules are centralized in `specification.yaml` and explained here.
- Recovered identities and exceptional decisions are visible in TSV files.
- Source hashes are never learned or updated automatically.
- Reconstruction and target validation are separate functions and workflow
  branches.
- Difficult historical behavior has focused unit tests and positive controls.
- Synthetic tests run quickly and do not require public downloads.
- Full-data validation reports the exact failed check and mismatching graph set.
- No “accept current output” command rewrites expected results.

## 10. Attribution and uncertainty

The reconstruction builds on contributions from the human investigators,
earlier Claude/Opus analysis, and subsequent GPT-5.6 Sol analysis. The human
investigators framed the provenance problem, supplied and audited sources,
required careful biological identifier handling, and repeatedly corrected
claims that exceeded the evidence. The earlier Opus work established major
parts of tissue matching, structural node alignment, feature-family discovery,
and the initial DGL transformation. GPT-5.6 Sol independently extended and
validated the work, including the complete legacy row-order mechanism, feature
selection rule, qualifier-aware GO policy, and full-row reconstruction.

The executable package records accepted results and evidence levels; the full
intellectual and chronological history remains in the surrounding project
reports rather than being duplicated here.

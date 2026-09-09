# Release and clean-clone checklist

This checklist turns the scientific reconstruction into a release process that
an independent reviewer can execute from a fresh repository checkout. It does
not add new reconstruction rules. It checks that the existing rules remain
readable, target-independent, deterministic, and fully validated.

## 1. What continuous integration checks

The repository workflow is `.github/workflows/reproduction.yml`.

For every push or pull request that changes `ai/reproduction/`, the fast job:

1. checks out a clean copy of the repository;
2. installs the environment strictly from the committed `pixi.lock`;
3. runs Ruff linting and formatting checks;
4. runs the focused pytest suite.

A full clean-clone job also runs:

- for pushes to the `ai` branch;
- when started manually with `workflow_dispatch`;
- on the first day of each month.

Pull requests receive the fast gate only. The full job downloads historical
source data and therefore runs only for trusted branch, manual, or scheduled
events.

The full job starts without a `data/` cache. It first reconstructs GraphSAGE and
DGL while the released GraphSAGE and DGL ZIP files are absent. It then repeats
the reconstruction and compares hashes of all 16 scientific output files. Only
after those checks does it download the two released targets and run independent
validation. Small manifests, summaries, validation reports, and checksum files
are retained as a GitHub Actions artifact for 30 days. Raw source data and the
large reconstructed datasets are not uploaded by the workflow.

The scheduled run is also a source-availability monitor. A failure to download a
checksum-locked historical input should be investigated as an upstream
availability problem, not worked around by changing its expected hash.

## 2. Manual release gate

Run the release gate from a new clone, not from a long-lived working directory.
Replace `<commit-or-tag>` with the revision being released.

```bash
git clone https://github.com/agitter/GraphSagePPIAnalysis.git
cd GraphSagePPIAnalysis
git checkout <commit-or-tag>
cd ai/reproduction

pixi install --frozen
pixi run test
```

Before reconstruction, confirm that neither validation target is present:

```bash
test ! -e data/graphsage_ppi.zip
test ! -e data/dgl_ppi.zip
```

Run the target-independent reconstruction:

```bash
pixi run reproduce

test ! -e data/graphsage_ppi.zip
test ! -e data/dgl_ppi.zip
```

This is the principal non-circularity check. Reconstruction is allowed to use
only the seven upstream sources and the committed specifications. The released
GraphSAGE and DGL archives may be opened only by validation.

Record the scientific artifact hashes, rebuild, and compare:

```bash
find results/graphsage/ppi results/dgl/ppi -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > /tmp/reproduction-run-1.sha256

pixi run clean
pixi run reproduce

find results/graphsage/ppi results/dgl/ppi -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > /tmp/reproduction-run-2.sha256

diff -u \
  /tmp/reproduction-run-1.sha256 \
  /tmp/reproduction-run-2.sha256
```

No `diff` output is expected. This comparison deliberately excludes the run
manifest and validation reports because they contain run-specific provenance.

Run independent validation and verify the generated checksum inventories:

```bash
pixi run validate

cat results/graphsage_validation.md
cat results/dgl_validation.md

(
  cd results/graphsage
  sha256sum -c SHA256SUMS
)

(
  cd results/dgl
  sha256sum -c SHA256SUMS
)
```

Both validation reports must say `Overall result: PASS`, and every checksum must
report `OK`.

Finally, confirm that the workflow has not changed a tracked file:

```bash
ROOT=$(git rev-parse --show-toplevel)
git -C "$ROOT" status --short
```

An empty result is expected. Cached inputs, generated intermediates, and results
are ignored; the workflow must not rewrite source code or committed
specifications.

## 3. Evidence to retain for a release

Retain these small files with the tagged release or its archival record:

```text
results/reconstruction_manifest.json
results/graphsage_validation.json
results/graphsage_validation.md
results/dgl_validation.json
results/dgl_validation.md
results/graphsage/SHA256SUMS
results/dgl/SHA256SUMS
/tmp/reproduction-run-1.sha256
/tmp/reproduction-run-2.sha256
```

Also record:

- the Git commit and tag;
- the SHA-256 of `pixi.lock`;
- the platform on which the full gate ran;
- the URL of the successful GitHub Actions run;
- any expected cross-platform floating-point differences already covered by the
  declared DGL tolerance.

Do not commit `data/`, `build/`, or `results/` merely to preserve a run. Attach
small evidence files to a release or deposit them in the project archive.

## 4. Expected non-byte-identical outputs

A release is valid even though these files are not byte-identical to their
historical counterparts:

- GraphSAGE `ppi-G.json`: the reconstructed graph has the exact nodes and
  undirected edge multiset but uses a clear canonical link order;
- DGL graph JSON files: directed links are canonically sorted;
- DGL feature NPY files: float64 values must satisfy the fixed absolute
  tolerance in `spec/specification.yaml`.

The GraphSAGE ID map, class map, and feature NPY file are expected to be
byte-identical. DGL graph-ID and label NPY files are also expected to be
byte-identical. `ppi-walks.txt` remains outside the package scope.

## 5. Independent input mirror policy

A durable mirror is important because several inputs are historical snapshots.
The mirror must preserve provenance rather than create a new, ambiguous source.

For each source independently:

1. Complete a redistribution and attribution review.
2. Mirror the exact original bytes without recompression, renaming, or format
   conversion.
3. Store the original filename, byte size, SHA-256, upstream URL, retrieval date,
   license or terms, and required citation beside the object.
4. Use an immutable, versioned URL from a durable research repository.
5. Add that URL to the existing `mirror_url` field in `spec/sources.tsv` without
   changing the expected size or SHA-256.
6. Test primary and mirror acquisition paths separately.
7. Treat a hash mismatch as an error; never update the accepted hash merely
   because a server returns different bytes.

The MSigDB v6.1 C1 and C3 inputs are identified in the project notices as
CC BY 4.0 and are the clearest initial mirror candidates. The OhmNet, historical
GOA, Gene Ontology, GraphSAGE, and DGL files remain `review_pending` in
`sources.tsv`; they should not be mirrored until their applicable terms and
notices have been recorded.

A source bundle may be offered as a convenience artifact, but it must not replace
the per-file identities in `sources.tsv`. The canonical identity of every input
remains its original filename, size, and SHA-256.

## 6. Publication gates not yet automated

Before a public versioned release, complete these human-review items:

- confirm the software license for the new reproduction code;
- finalize author and contributor attribution;
- add repository-level citation metadata;
- review third-party notices and mirror permissions;
- ensure the scientific report uses the same exact/strongly inferred/documented/
  open evidence language as the executable package;
- record the known manuscript/data discrepancies outside the central workflow;
- archive a successful clean-clone evidence set.

These items require investigator decisions and should not be inferred or silently
filled by the workflow.

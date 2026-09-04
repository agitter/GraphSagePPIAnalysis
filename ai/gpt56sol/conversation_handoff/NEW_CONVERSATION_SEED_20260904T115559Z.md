# Seed for a new ChatGPT conversation — GraphSAGE PPI reproduction project

I am continuing a long forensic/reproducibility investigation of the GraphSAGE PPI benchmark and its DGL derivative. Please treat the attached full handoff and critical artifacts as authoritative context, and ask me to resupply raw external inputs only when needed.

## Current core result

We now have a complete deterministic reconstruction of the released GraphSAGE PPI data at the data level:

- 24 tissue graphs
- 56,944 tissue-instance node rows
- 4,301 distinct Entrez Gene IDs under the complete node-order reconstruction
- 818,716 GraphSAGE link records
- 50 binary input features
- 121 binary GO labels

The exact released data are reproduced under a single global workflow; remaining uncertainty concerns historical implementation details rather than unexplained matrix mismatches.

### Topology and node identities

The 24 GraphSAGE graph components are exact matches to named OhmNet tissue-specific PPI edgelists.

The crucial node-order breakthrough is a strongly supported legacy Python mechanism:

1. read original OhmNet edgelists in line order;
2. keep Entrez IDs as strings;
3. build the graph under legacy NetworkX semantics;
4. emulate 64-bit unrandomized CPython 2.7 string dictionaries;
5. iterate dictionary table order.

This agrees with all 56,411 identities that had already been independently anchored by graph topology/features and resolves the remaining 533 rows, giving a complete 56,944-row mapping. A separate CPython 2 simulation reproduces the full irregular key order of `ppi-class_map.json`, providing an independent positive control.

Treat the complete row mapping as exact under this mechanism; describe the historical mechanism itself as strongly inferred, not source-code proven.

### Features

The 50 observed feature columns are exactly:

- 30 MSigDB C1 positional sets
- 20 MSigDB C3 motif/regulatory sets
- no C7 columns in the released matrix

A strongly supported exact reconstruction hypothesis is:

1. process C1 in source order;
2. then C3 in source order;
3. keep source sets around the 200-member boundary;
4. append qualifying sets until a global cap of 50;
5. project to GraphSAGE rows.

`chryq11` is important evidence: it has 204 members in the source set but no GraphSAGE genes, so its released column is all zero. This strongly supports source-level filtering before graph projection and no later removal of empty columns.

MSigDB 5.0, 5.1, 5.2, and 6.0 all yield the same 50 observed vectors. The historical version cannot be inferred from the matrix. The clean reproduction will likely use v6.0 C1/C3 for licensing reasons, with documentation that earlier tested versions give identical results.

### Labels

Exact label reconstruction uses:

- GOA human release 159 GAF
- GOA human release 159 GPI
- historical May/June 2016 `gp2protein.geneid`
- June 2016 GO ontology

Global policy:

- evidence: EXP, IDA, IEP, IGI, IMP, ISS
- exclude `NOT`
- retain ordinary relations:
  - BP: `involved_in`
  - CC: `part_of`
  - MF: `enables`
- do not treat `colocalizes_with` or `contributes_to` as ordinary binary membership
- canonicalize GO alt IDs
- propagate direct term + transitive `is_a` ancestors only
- do not propagate ontology `part_of`

Candidate terms span all three namespaces:

- 85 BP
- 26 CC
- 10 MF

Term selection on the full historical human universe is exactly the top 121 terms by prevalence and observationally equivalent to an approximately >=1000 human-gene/protein threshold in release 159.

A controlled GOA release screen 158-169 found release 159 uniquely exact under the fixed transformation.

The early "121/121 exact" result covered only individually resolved rows; that was later corrected. Under the final row map, all 56,944 x 121 = 6,890,224 label cells are exact with one fixed mapping/filter/propagation policy and no per-gene tuning.

### Identifier mapping

Preserve historical many-to-many UniProt-GeneID mappings as bipartite components. Do not force one-to-one mappings or filter to GraphSAGE genes before component resolution.

Important O95073/Q9Y620 case:

- historical O95073 / FSBP carried GeneIDs 100861412 and 25788
- Q9Y620 / RAD54B carried GeneID 25788
- O95073->25788 was a real historical UniProt cross-reference
- nevertheless, projecting O95073/FSBP GO annotations onto the GraphSAGE RAD54B/25788 node creates 13 false positives
- reconstruction semantically separates O95073/FSBP->100861412 and Q9Y620/RAD54B->25788 for annotation projection

Do not describe the historical cross-reference as corrupt.

### Label-column ambiguity

There are 121 candidate terms but 118 distinct binary vectors. Three duplicated-vector pairs are membership-indistinguishable:

- cols 24/71: GO:0043228 / GO:0043232
- cols 39/63: GO:0006464 / GO:0036211
- cols 48/70: GO:0043230 / GO:1903561

A strong CPython 2 dictionary-order fingerprint supports the provisional orientation:

- 24 -> GO:0043228; 71 -> GO:0043232
- 39 -> GO:0006464; 63 -> GO:0036211
- 48 -> GO:1903561; 70 -> GO:0043230

Call these strongly supported/provisional, not proven.

### PPI provenance

GraphSAGE's immediate source is OhmNet, not simply BioGRID.

A BioSNAP audit established:

- union of 144 supplied OhmNet edgelists = 70,338 unique undirected pairs
- exactly equals BioSNAP combined OhmNet network
- all are a subset of the BioSNAP 342,353-edge global human interactome

The GraphSAGE dataset webpage linking only to BioGRID is incomplete/misleading. OhmNet's global interactome is composite; BioGRID is one contributor.

### DGL

DGL is a deterministic downstream transformation of GraphSAGE. Broadly: graph/component ordering, training-only StandardScaler feature normalization, directed edge expansion, self-loops, split packaging. The clean reproduction should derive DGL from reconstructed GraphSAGE, never from the downloaded DGL target.

### Leakage

The graph split is not gene-disjoint:

- 3,599 / 4,301 genes appear in multiple splits
- 3,067 appear in train, validation, and test
- 5,490 / 5,524 test rows have GeneIDs seen in training
- simple GeneID lookup gives about 0.9971784 test micro-F1 under the full map

Later we will build randomized-label controls at the gene level to demonstrate that the performance can survive destruction of biological meaning if the repeated identity structure is preserved. This is deferred until after the clean reproduction package.

## Confidence categories

Use four high-level categories:

1. byte/data-level exact
2. strongly inferred
3. documented
4. open

Do not collapse exact output reconstruction and inferred historical implementation into the same claim.

## Reproduction package — agreed scope

Exact scope statement:

> Reconstruct every deterministic artifact needed by supervised GraphSAGE and DGL, but do not attempt to regenerate the original stochastic `ppi-walks.txt` byte-for-byte.

The clean package should omit most exploratory audits and only reproduce the core path.

Technology choices:

- Snakemake
- Pixi only; no `environment.yml`
- committed `pixi.lock`
- pytest
- one Snakefile initially
- no large central CLI
- user-facing commands such as `pixi run reproduce`, `pixi run test`, `pixi run clean`
- readable modern Python reimplementation of legacy Python-2 dictionary behavior; do not install Python 2

Agreed structure:

```text
reproduction/
├── README.md
├── pixi.toml
├── pixi.lock
├── pyproject.toml
├── Snakefile
├── config/
│   ├── sources.tsv
│   └── reconstruction.yaml
├── resources/
│   ├── selected_graphs.tsv
│   ├── label_columns.tsv
│   ├── identifier_decisions.tsv
│   └── expected.yaml
├── src/
│   └── graphsage_ppi_repro/
│       ├── __init__.py
│       ├── sources.py
│       ├── manifest.py
│       ├── legacy_dict.py
│       ├── topology.py
│       ├── features.py
│       ├── labels.py
│       ├── graphsage.py
│       ├── dgl.py
│       └── validate.py
├── tests/
│   ├── data/
│   ├── test_legacy_dict.py
│   ├── test_identifier_mapping.py
│   ├── test_go_labels.py
│   ├── test_features.py
│   └── test_reproduction.py
├── docs/
│   ├── architecture.md
│   ├── data-sources.md
│   └── limitations.md
├── build/      # ignored
└── results/    # ignored
```

`config/sources.tsv` should be flat and maintainable, approximately:

```text
source_id role filename url archive_url sha256 size_bytes access license_note description
```

`reconstruction.yaml` should record only the scientific reconstruction policy, not execution settings.

Keep scientific narrative/reporting outside `reproduction/` at repo level; later create `claims/claims.csv`.

Core circularity rule:

> Reconstruction must not read GraphSAGE or DGL reference targets. Only validation may compare reconstructed outputs to those targets.

## MSigDB plan

The login wall complicates public CI. Because v6.0 yields the same observed features and changed licensing, plan to:

1. use v6.0 C1/C3 as canonical reproduction input;
2. verify redistribution terms for the exact needed content;
3. if permitted, archive only the minimal original-order C1/C3 input publicly;
4. provide a script that transforms a private original MSigDB package into that minimal input;
5. document that 5.0-6.0 all reproduce the same observed feature matrix.

Avoid private CI if a legal public archival copy of the minimal v6.0 source can be used.

## `ppi-walks.txt` side investigation

This file is optional and only for unsupervised GraphSAGE, so it is outside the deterministic reproduction scope.

The original supplied `graphsage.utils.py` shows:

- `WALK_LEN=5`
- `N_WALKS=50`
- training nodes only
- training-induced subgraph
- Python `random.choice(G.neighbors(curr_node))`
- no explicit `random.seed`

The user nevertheless wants a side investigation into whether the original walk output can be recovered by exhaustively/pruned search using the original GraphSAGE code, likely Python/NetworkX versions, target output, common/arbitrary seeds, and early-abort prefix matching. This side analysis had not started before the previous conversation ended.

If we do that search, do not let it delay the core reproduction implementation.

## Repository organization outside reproduction

Current repo root is approximately:

```text
README.md
data/      # ignored external data
gpt56sol/  # prior GPT-5.6 artifacts
opus46/    # Claude/Opus artifacts
papers/    # ignored local papers
```

The user wants future scientific reports and `claims.csv` outside the software package.

## Important corrected/superseded claims

Do not regress to these old conclusions:

- not merely 4,268/4,278 genes mapped: complete mapping now has 4,301 distinct genes over all 56,944 rows
- labels are not BP-only: 85 BP, 26 CC, 10 MF
- feature matrix does not uniquely identify MSigDB v5.2
- PPI source is not simply BioGRID
- O95073->25788 was not a corrupt historical mapping
- early 121/121 label exactness did not cover every row; current final reconstruction does

## Immediate next step

The main next step is to implement the clean `reproduction/` package from scratch with maximal readability and comments, beginning with the repository skeleton/source manifest and then the topology/node-order stage.

Optionally, first perform the side `ppi-walks.txt` seed/order search if I ask for it.

Please use the accompanying `FULL_PROJECT_HANDOFF_20260904T115559Z.md` and `critical_artifacts/` for details rather than asking me to re-explain prior conclusions. I can resupply any omitted external input files when needed.

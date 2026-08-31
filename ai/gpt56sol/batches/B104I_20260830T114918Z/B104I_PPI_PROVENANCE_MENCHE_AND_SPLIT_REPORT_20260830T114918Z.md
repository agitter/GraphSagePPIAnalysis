# B104I — PPI provenance, Menche comparison, split semantics, and public-issue tracker

Generated: 20260830T114918Z UTC

## 1. Meaning of `split` in the node mapping

`split` is a single categorical value for each row/node instance: `train`, `validation`, or `test`. It is inherited from the tissue graph containing that row. No graph contains more than one split value.

Row counts:

- train: 44,906
- validation: 6,514
- test: 5,524

The same Entrez GeneID can occur in several tissue graphs, so a gene can occur in multiple splits. Among 4,301 distinct Entrez GeneIDs:

- 620 occur only in train;
- 69 occur only in validation;
- 13 occur only in test;
- 345 occur in train and validation;
- 166 occur in train and test;
- 21 occur in validation and test;
- 3,067 occur in all three.

Thus 3,599 genes occur in multiple splits. The node table is one row per tissue-specific node instance; the companion gene-level table uses a semicolon-delimited `splits_seen` field.

## 2. BioSNAP audit result

The uploaded audit passed its stated checks. The union of all 144 supplied OhmNet tissue edgelists contains 70,338 unique undirected pairs and is byte-level equivalent as an edge set to the BioSNAP combined OhmNet network. All 70,338 pairs occur in the 342,353-edge BioSNAP global physical interactome.

This establishes the external source chain:

```text
BioSNAP 342,353-edge global interactome
    contains all 70,338 unique OhmNet tissue-edge pairs
        exactly equals union of the 144 supplied tissue edgelists
            contains the 24 tissue graphs used by GraphSAGE
```

It does not provide per-edge attribution to BioGRID, IntAct/MIntAct, HPRD, CORUM, Rolland, or other constituent databases.

## 3. What Menche et al. contributes

Menche et al. 2015 reports an interactome of 13,460 proteins and 141,296 interactions. Its construction combines seven classes: regulatory interactions, binary yeast-two-hybrid/literature interactions, mostly low-throughput literature interactions, metabolic enzyme-coupled interactions, protein complexes, kinase-substrate pairs, and signaling interactions.

This is not numerically identical to the 21,557-node/342,353-edge OhmNet/BioSNAP global network. Therefore, the OhmNet citation to Menche should not be interpreted as saying that the exact Menche 2015 edge file was used unchanged. The best-supported interpretation is that OhmNet used a later or enlarged composite interactome in the same methodological/source lineage.

The terminology also needs care: the global resource is broader than a narrow binary-binding-only PPI map because its documented lineage includes regulatory, metabolic, kinase-substrate, signaling, and complex-derived edges. The safest description is **experimentally supported human physical/interactome network**, followed by the exact source details.

## 4. Literature differences to preserve in the final report

The accompanying evidence-versus-literature register records each claim as agreement, refinement, conflict/incompleteness, or open. The most important differences are:

1. GraphSAGE states C1, C3, and C7 feature families; the deposited 50 columns are exactly 30 C1 plus 20 C3 and contain no C7 column. The cap-based explanation is strongly supported but not documentary proof.
2. GraphSAGE says GO labels were collected from MSigDB; direct memberships from the tested MSigDB releases do not reproduce a single complete label column, while the fixed GOA release-159 transformation reproduces all deposited cells. MSigDB may have supplied term identities rather than memberships.
3. The paper's held-out graphs are genuinely unseen tissue graphs, but most Entrez gene identities recur across graph splits. The task is graph-inductive but not gene-identity-disjoint.
4. OhmNet's global-network size matches the BioSNAP 21,557/342,353 resource, not the Menche 2015 13,460/141,296 interactome.
5. “BioGRID network” is an inadequate shorthand; BioGRID is one contributor to a composite network.

## 5. GitHub issue follow-up

The issue tracker records #78, #86, #188, #190, #16, and #32. The first public replies should likely target #188 (node/protein correspondence), #190 (50-dimensional features), and #78 (overall source and preprocessing). Replies should wait until the report, source downloader, checksum manifest, and minimal reproduction script are frozen.

## 6. Evidence categories

- **Byte-level exact:** direct equality or complete reconstruction from retained files.
- **Strongly inferred:** a single historically plausible mechanism explains all available data and survives independent controls.
- **Documented:** explicitly stated in a manuscript, official page, archive metadata, or source code.
- **Open:** multiple histories remain observationally equivalent or an intermediate artifact is missing.

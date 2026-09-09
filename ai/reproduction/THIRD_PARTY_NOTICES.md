# Third-party data notices

This repository contains software for acquiring and transforming third-party
data. Except for small evidence-backed specifications and test fixtures, the
upstream datasets and released reference archives are not committed here. They
are downloaded into the ignored package-local `reproduction/data/` cache and
remain governed by their own terms.

## Molecular Signatures Database (MSigDB) v6.1

The feature reconstruction uses:

- `c1.all.v6.1.entrez.gmt`
- `c3.all.v6.1.entrez.gmt`

Source: Broad Institute / Gene Set Enrichment Analysis MSigDB release 6.1

- https://data.broadinstitute.org/gsea-msigdb/msigdb/release/6.1/
- https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp

MSigDB versions 6.0 through 7.5.1 are identified by the publisher as licensed
under the Creative Commons Attribution 4.0 International license. The additional
collection-specific restrictions listed by MSigDB concern specified C2 and M2
materials; this workflow uses C1 and C3 only.

License: https://creativecommons.org/licenses/by/4.0/

The workflow preserves the original GMT files in the local cache and records
their exact SHA-256 values. It selects source rows by documented rules and
projects their Entrez memberships onto the reconstructed GraphSAGE rows. The
resulting matrix and selected-column table are adapted outputs. The project
report should cite MSigDB and the appropriate underlying gene-set sources in
addition to retaining this license notice.

## OhmNet tissue networks

The topology reconstruction uses the OhmNet multi-layer human tissue network:

- https://snap.stanford.edu/ohmnet/

Suggested citation supplied with the data:

> Marinka Zitnik and Jure Leskovec. Predicting multicellular function through
> multi-layer tissue networks. Bioinformatics, 2017.

The supplied README states that each layer is a tissue-specific human PPI
network, that nodes are Entrez GeneIDs, and that the same ID in two layers is
the same gene. The redistribution terms for the complete raw archive should be
confirmed before creating the planned independent input mirror. Until then, the
workflow downloads the archive from its official source and does not commit it.

## Gene Ontology and UniProt-GOA

The label reconstruction uses dated Gene Ontology and UniProt-GOA files whose
exact source URLs, release identifiers, sizes, and hashes are recorded in
`spec/sources.tsv`:

- GOA human release-159 GAF and GPI files;
- the 2016-06-01 GO ontology;
- the associated historical GeneID-to-UniProt mapping.

Relevant project sites:

- https://geneontology.org/
- https://www.ebi.ac.uk/GOA/

Publications and redistributed source bundles should retain applicable upstream
notices and cite the Gene Ontology Consortium and UniProt-GOA project. The
workflow keeps these inputs unmodified in the ignored cache. Their redistribution
terms must be reviewed before they are added to the planned independent mirror.

## GraphSAGE reference data

The released GraphSAGE PPI ZIP is a validation target only:

- https://snap.stanford.edu/graphsage/
- https://github.com/williamleif/GraphSAGE

It is not used to construct topology, node identities, features, or labels. The
optional stochastic `ppi-walks.txt` member is outside this package's scope.

Suggested citation:

> William L. Hamilton, Rex Ying, and Jure Leskovec. Inductive Representation
> Learning on Large Graphs. NeurIPS, 2017.

## DGL reference data

The DGL PPI archive is used only to validate the downstream logical and data
transformation:

- https://data.dgl.ai/dataset/ppi.zip
- https://www.dgl.ai/

DGL outputs are reconstructed from the rebuilt GraphSAGE data rather than copied
from this reference archive. The reproduction does not redistribute the target
ZIP and does not require the DGL Python package to write its JSON/NumPy format.

## Software license versus data licenses

The repository's software license does not supersede any third-party data
license. Users who redistribute source files, mirrors, or derived materials are
responsible for retaining required notices, providing attribution, identifying
modifications, and complying with collection-specific restrictions.

# GraphSAGE PPI data-provenance and reproducibility issues

This tracker is for future public replies after the reconstruction report and code package are frozen. Claims should be limited to released artifacts and should link to reproducible outputs.

## #78 — A question about the ppi dataset

- **Status:** Open
- **Opened:** 2019-05-03
- **Question:** Exact raw dataset, node features, and source links; later comments ask how C1/C3/C7 became 50 dimensions and request preprocessing code.
- **Published Response:** A collaborator linked the raw graph source and MSigDB and described C1/C3/C7 as features and GO as labels.
- **Our Relevant Outputs:** Exact OhmNet source mapping; recovered node-to-Entrez table; exact feature reconstruction; exact GOA label reconstruction; evidence-versus-literature register.
- **Future Response Priority:** Highest
- **Url:** https://github.com/williamleif/GraphSAGE/issues/78

## #86 — the description of dataset

- **Status:** Open
- **Opened:** 2019-06-09
- **Question:** Requests a usable description of the dataset.
- **Published Response:** No substantive response visible.
- **Our Relevant Outputs:** Complete source chain, schemas, biological identifiers, tissue names, feature names, provisional GO label names, transformations, and caveats.
- **Future Response Priority:** High
- **Url:** https://github.com/williamleif/GraphSAGE/issues/86

## #188 — PPI ID-protein name correspondence?

- **Status:** Open
- **Opened:** 2022-07-14
- **Question:** Requests a mapping between PPI IDs and protein names.
- **Published Response:** No answer visible.
- **Our Relevant Outputs:** 56,944-row node-to-Entrez mapping plus evidence-rich validation table.
- **Future Response Priority:** Highest
- **Url:** https://github.com/williamleif/GraphSAGE/issues/188

## #190 — Question about PPI: how to process C1, C3 and C7 into 50-dim features?

- **Status:** Open
- **Opened:** 2022-08-17
- **Question:** Asks how large C1/C3/C7 collections were converted to 50 dimensions.
- **Published Response:** No answer visible.
- **Our Relevant Outputs:** Exact 50-column values and order from 30 C1 plus 20 C3 sets; strongly supported source-size/global-cap hypothesis; chryq11 zero-column evidence; explicit remaining uncertainty.
- **Future Response Priority:** Highest
- **Url:** https://github.com/williamleif/GraphSAGE/issues/190

## #16 — Is f1_score evaluation wrong in ppi_eval.py?

- **Status:** Closed
- **Opened:** 2017-11-14
- **Question:** Questions per-class micro-F1 computation in ppi_eval.py.
- **Published Response:** Issue is closed; the visible page does not show a substantive resolution.
- **Our Relevant Outputs:** Future leakage and evaluation report should distinguish paper metric, repository evaluation code, modern global micro-F1, and gene-lookup controls.
- **Future Response Priority:** Medium (evaluation phase)
- **Url:** https://github.com/williamleif/GraphSAGE/issues/16

## #32 — Replicating the results in the paper

- **Status:** Open
- **Opened:** 2018-03-30
- **Question:** Reports supervised PPI F1 below the paper and asks for exact hyperparameters.
- **Published Response:** No substantive response visible.
- **Our Relevant Outputs:** Future reproducibility package and leakage experiments; separate hyperparameter replication from dataset leakage.
- **Future Response Priority:** Medium (ML phase)
- **Url:** https://github.com/williamleif/GraphSAGE/issues/32


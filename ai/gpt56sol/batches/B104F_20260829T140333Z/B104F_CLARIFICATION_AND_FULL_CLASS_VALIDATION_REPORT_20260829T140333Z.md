# B104F — Clarifications, full-class GO-label validation, and input audit

Generated: 2026-08-29T14:03:33Z

## Executive corrections

1. The MSigDB feature-generation rule is a **strong, parsimonious reconstruction hypothesis**, not proof of the original preprocessing code. It exactly reproduces the observed 50-column feature matrix, and the empty `chryq11` column is strong support for source-level size filtering before intersecting with graph genes, but alternative code paths could select the same 50 sets.
2. The earlier statement that the GO label matrix was “completely reconstructed” needed a scope qualifier. The original cell-by-cell proof covered 56,411 rows corresponding to 4,268 GeneIDs independently resolved from topology and input features. It excluded 533 rows in 183 topology/feature equivalence classes.
3. This batch extends the validation without forcing an arbitrary row-to-gene permutation. Under one fixed global GOA/mapping policy, every one of the 183 unresolved classes has an **exact multiset match** between observed GraphSAGE 121-vectors and predicted GOA 121-vectors for its candidate GeneIDs. Therefore all 56,944 rows are reproduced exactly **up to permutations among nodes that are indistinguishable by graph topology and the 50 input features**.
4. Of the 533 unresolved rows, 438 have a unique predicted GO vector within their equivalence class. Ninety-five remain ambiguous because multiple candidate genes have identical predicted 121-vectors.
5. The four newly uploaded `dhimmel-...tsv` files are not TSV data. They are GitHub HTML blob pages saved with `.tsv` filenames and are rejected as scientific inputs.

## MSigDB feature rule: evidence and limits

The following deterministic rule reproduces all 50 feature columns on all 56,411 independently resolved rows:

- process C1, then C3, then C7;
- preserve GMT row order;
- retain sets having at least 200 or at least 201 unique Entrez IDs (the data cannot distinguish these operators);
- append qualifying sets until a global cap of 50 is reached.

The result is 30 C1 sets followed by 20 C3 sets. The cap is reached before C7.

### Strong supporting evidence

- Cell-level equality is exact for every resolved row and every feature column.
- A second implementation reproduces the same equality.
- The selected set sequence is stable across supplied MSigDB versions 5.0, 5.1, 5.2, and 6.0.
- `chryq11` has 204 genes in the source MSigDB set and therefore passes a source-level threshold, but has zero genes in the resolved graph universe, producing an all-zero deposited column. This is unlikely under a workflow that filtered or inspected columns after intersecting with graph genes.
- The first qualifying C3 set after the 50-column cap is not included, giving a natural stopping boundary.

### What remains unproven

- Whether the original code literally iterated C1, C3, and C7, or only C1 and C3.
- Whether the intended test was `>= 200` or `> 200`.
- Whether the global cap was implemented programmatically or an explicit 50-set list happened to be equivalent.
- Which MSigDB version was used; the feature matrix does not distinguish versions 5.0–6.0.

Recommended wording: **“An exact and strongly supported reconstruction of the feature-selection rule,”** not “the original procedure is proven.”

## Exact GO-label policy

The fixed policy is:

- GOA human release 159.
- GAF object IDs are UniProt accessions.
- Historical Entrez–UniProt edges come from the May 9, 2016 `gp2protein.geneid` file.
- GPI159 defines the accession/reference-proteome universe and supplies primary symbols for component resolution and unique fallback.
- Evidence codes: `EXP`, `IDA`, `IEP`, `IGI`, `IMP`, `ISS`.
- Exclude `NOT`.
- Retain only ordinary aspect relations: `involved_in`, `part_of`, `enables`.
- Do not collapse `colocalizes_with` or `contributes_to` into ordinary membership.
- Canonicalize alternate GO IDs.
- Propagate through ontology `is_a` only, including the direct term.
- Do not propagate through ontology `part_of`.
- Preserve many-to-many identifier mappings unless a complete component has a unique primary-symbol bijection.
- Use a unique primary-symbol fallback only when an accession has no historical GeneID edge.
- Do not project O95073/FSBP annotations to the GraphSAGE node for GeneID 25788/RAD54B.

No mapping edge, evidence filter, relation filter, or GO rule varies by label column or by gene.

## Full-row validation

### Independently resolved portion

- 56,411 GraphSAGE rows.
- 4,268 unique Entrez GeneIDs.
- Zero mismatched rows under the fixed policy.

### Topology/feature-equivalence portion

- 183 equivalence classes.
- 533 rows.
- 4,301 total candidate GeneIDs across resolved and unresolved portions.
- 183/183 classes have identical observed and predicted label-vector multisets.
- 438/533 unresolved rows have a unique GeneID candidate based on the fixed GOA vector within the class.
- 95/533 remain ambiguous because two or more candidate GeneIDs share the same vector.

This does not use an arbitrary permutation to score success: each class is compared as a multiset. Assigning the 438 uniquely matching rows afterward uses the deposited labels and therefore is not independent gene-identity evidence.

## The 12 genes reported by the independent verifier

All 12 are absent from the 4,268 independently resolved-gene table. They occur in the 533-row topology/feature equivalence classes. This explains why a forced row-order mapping can make them appear as annotation errors.

| GeneID | Symbol | GPI159 accession | Fixed GOA positives | Result within every relevant class |
|---:|---|---|---:|---|
| 3248 | HPGD | P15428 | 37 | one exact row |
| 3988 | LIPA | P38571 | 10 | one exact row in each of 4 classes |
| 8564 | KMO | O15229 | 32 | one exact row in each of 4 classes |
| 27201 | GPR78 | Q96P69 | 6 | one exact row |
| 30061 | SLC40A1 | Q9NP59 | 17 | one exact row in each of 4 classes |
| 51166 | AADAT | Q8N5Z0 | 14 | one exact row in each of 4 classes |
| 51312 | SLC25A37 | Q9NYZ2 | 0 | one exact all-zero row in each of 2 classes |
| 55471 | NDUFAF7 | Q7L592 | 14 | one exact row in each of 16 classes |
| 55801 | IL26 | Q9NPH9 | 51 | one exact row |
| 56994 | CHPT1 | Q8WUD6 | 14 | one exact row in each of 4 classes |
| 79017 | GGCT | O75223 | 21 | one exact row in each of 16 classes |
| 121599 | SPIC | Q8N5J4 | 0 | exact all-zero vector, shared with GeneID 8609 in the brain class |

Eleven of the twelve have a unique exact row match in every class in which they occur. GeneID 121599 cannot be distinguished from GeneID 8609 in its class using the deposited label vector because both are all zero.

## GPI fields and their role

### Used

- `DB_Object_ID`: UniProt accession key joining GAF/GPI to `gp2protein.geneid`.
- `DB_Object_Symbol`: primary-symbol consistency within mapping components and unique fallback when there is no historical edge.
- `Taxon`: audit/validation; the file is human-specific.

### Audit only

- `DB_Object_Name`.
- `DB_Object_Synonyms`; broad synonym joins were not used.
- `Properties`.

### Not usable for GeneID mapping in this release

- `DB_Xrefs`: empty in every GPI159 row.
- `Parent_Object_ID`: empty in every GPI159 row.

`gp2protein.geneid` supplies the actual GeneID–UniProt edges. GPI does not replace that file; it supplies the GOA object universe and semantic metadata for those accessions.

## Why another reproduction may obtain only 24/121 exact

The result is highly sensitive to several choices. Examples measured on the independently resolved 4,268 genes include:

- fixed exact policy: 121/121 exact;
- retain `colocalizes_with` and `contributes_to`: 90/121 exact after the O95073 correction;
- use default relations but retain the O95073→25788 projection: 108/121 exact;
- make both of those mistakes: 89/121 exact;
- use the common “experimental evidence” set `EXP, IDA, IEP, IGI, IMP, IPI` rather than replacing IPI with ISS: 3/121 exact under `is_a` propagation;
- propagate through ontology `part_of`: strongly degrades the result.

A 24/121 outcome is therefore evidence that at least one policy differs; it is not evidence that GOA159 cannot reproduce the labels.

## Audit of the four uploaded dhimmel files

All four uploaded files begin with `<!DOCTYPE html>` and contain GitHub page HTML. Their sizes are approximately 147–151 KB, not the multi-megabyte annotation tables. They must not be parsed as TSV.

They are recorded as `rejected_wrong_media_type`. A commit-pinned downloader is supplied in this batch. It rejects HTML, requires the `go_id` TSV header, requires a plausible minimum size, writes atomically, and records URL, retrieval time, byte size, SHA-256, and line count.

The genuine dhimmel products remain useful controls, but their published processing policy differs from the exact GraphSAGE policy: the website describes “experimental” as including IPI but not ISS, and inferred annotations as propagating through both `is_a` and `part_of`. Therefore they should not be expected to reproduce the exact matrix without reprocessing.

## Claim-strength recommendations

### Proven by exact computation on supplied data

- The fixed release-159 transformation reproduces all 56,411 independently resolved rows exactly.
- Every unresolved topology/feature class has an exact observed-versus-predicted label-vector multiset.
- The 12 reported genes are members of unresolved topology/feature classes, not failures of the fixed GOA transformation.

### Strongly supported, not source-code proven

- The parsimonious MSigDB feature-selection rule.
- Python 2 dictionary iteration as the source of label-column order.
- The orientation of the three duplicate-vector GO-term pairs.

### Still open

- Exact original preprocessing code.
- Exact row-to-GeneID identity for 95 rows that remain indistinguishable even after considering the 121 labels.
- Whether GOA was read directly or an Entrez-native product carried the same associations.
- The exact source of the label-column ordering.

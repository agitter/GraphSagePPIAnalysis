# Response to independent verifier

## 1. Which source/version produced 121/121?

The exact source combination was not a `gene2go` table. It was:

- GOA human GAF release 159 (`goa_human.gaf.159.gz`, generated 2016-07-04);
- GPI release 159 for UniProt object IDs and primary symbols;
- the May 9, 2016 `gp2protein.geneid` crosswalk;
- the June 1, 2016 `go.obo` archive ontology, using only `is_a` edges.

The exact filtering and propagation policy is:

```text
Evidence: EXP, IDA, IEP, IGI, IMP, ISS
Exclude: NOT
GAF relations retained: involved_in, part_of, enables
GAF relations excluded: colocalizes_with, contributes_to
Ontology propagation: is_a only, including the direct term
Ontology part_of propagation: no
```

Alternate GO IDs are canonicalized.

The original 121/121 cell-level comparison was on 4,268 GeneIDs independently resolved from graph topology and the 50 input features. A subsequent class-level validation covers all 56,944 rows: the remaining 533 rows form 183 topology/feature equivalence classes, and all 183 classes have identical observed and predicted label-vector multisets under the same global policy.

If you get 24/121, first check the six-code evidence set, qualifier relations, ontology `part_of`, and identifier permutation inside unresolved classes.

## 2. What happened with GeneIDs 3248, 3988, 8564, 27201, 30061, 51166, 51312, 55471, 55801, 56994, 79017, and 121599?

They were not included in the initial 4,268 independently resolved-gene comparison. Every one occurs in a topology/feature equivalence class where multiple GraphSAGE rows and multiple OhmNet GeneIDs are indistinguishable by the graph and the 50 input features.

Their historical accession mappings are ordinary and fixed:

```text
3248   HPGD      P15428
3988   LIPA      P38571
8564   KMO       O15229
27201  GPR78     Q96P69
30061  SLC40A1   Q9NP59
51166  AADAT     Q8N5Z0
51312  SLC25A37  Q9NYZ2
55471  NDUFAF7   Q7L592
55801  IL26      Q9NPH9
56994  CHPT1     Q8WUD6
79017  GGCT      O75223
121599 SPIC      Q8N5J4
```

Under the fixed GOA transformation, eleven of the twelve have exactly one matching GraphSAGE label row in every equivalence class in which they occur. GeneID 121599 has an all-zero vector shared with GeneID 8609 in one brain class, so those two remain interchangeable.

Do not assign these genes to rows by sorted order. Compare the class as a multiset or use the supplied row-match table. A forced permutation is the likely source of your reported errors.

## 3. Was the mapping tuned per gene or per label?

No. One mapping and one annotation policy are used for every gene and every label column.

Mapping policy:

1. Preserve all May-2016 `gp2protein.geneid` edges.
2. Resolve a complete many-to-many component only when primary symbols establish a unique bijection.
3. Use a unique primary-symbol fallback only if an accession has no historical GeneID edge.
4. Do not use broad synonyms.
5. Do not project O95073/FSBP annotations to GeneID 25788/RAD54B; retain Q9Y620 as the RAD54B protein and O95073 as FSBP for annotation projection.

The O95073 rule is global and was decided from the complete mapping component and symbols, not by optimizing individual labels. Date-matched UniProt records confirm that O95073 carried both GeneID cross-references in 2016, so the exact original upstream mechanism may instead have been an Entrez-native annotation source that already separated FSBP from RAD54B.

For the 533 unresolved rows, no row-to-gene permutation is selected when scoring the full validation: observed and predicted vectors are compared as class-level multisets.

## 4. Which GPI fields were used, and how is GPI different from `gp2protein.geneid`?

Used from GPI159:

- `DB_Object_ID`: UniProt accession; joins GAF annotations to the mapping.
- `DB_Object_Symbol`: primary symbol; used only for component consistency and unique fallback.
- `Taxon`: validation.

Audit only:

- `DB_Object_Name`.
- `DB_Object_Synonyms`; not used for mapping.
- `Properties`.

Not used because empty throughout GPI159:

- `DB_Xrefs`.
- `Parent_Object_ID`.

`gp2protein.geneid` is the actual two-column GeneID↔UniProt edge source. GPI defines the release-159 GOA object/reference-proteome universe and supplies semantic metadata, but it does not contain the GeneID crosswalk.

## Reproduction diagnostics

On the 4,268 independently resolved genes:

```text
Exact policy                                                121 exact columns
Include colocalizes_with and contributes_to                  90 exact columns
Keep O95073 -> 25788 but otherwise use exact relation policy 108 exact columns
Make both mistakes                                           89 exact columns
Use EXP/IDA/IEP/IGI/IMP/IPI instead of .../ISS                3 exact columns
```

The complete scripts and per-class/per-gene outputs are in the B104F bundle.

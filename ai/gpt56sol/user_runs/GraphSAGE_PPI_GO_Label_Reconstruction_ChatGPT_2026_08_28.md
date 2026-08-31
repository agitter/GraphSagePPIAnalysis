# GraphSAGE PPI GO Label Reconstruction — Investigator Handoff

## Executive result

We can exactly reproduce all 121 GraphSAGE PPI label columns on the 4,268 independently resolved Entrez Gene IDs after collapsing repeated tissue instances.

The compared matrix contains:

- 4,268 resolved genes
- 121 label columns
- 516,428 binary cells
- 0 false positives
- 0 false negatives

A second, independently implemented set-based pipeline reproduces the same result.

This equality applies to every GraphSAGE node row whose biological identity was independently recovered. It does not by itself assign GeneIDs to the remaining unresolved anonymous nodes.

## Required historical inputs

The exact reconstruction uses:

- GraphSAGE `ppi.zip`
- OhmNet tissue-network edgelists
- GOA human release 159 GAF and GPI
- May 2016 `gp2protein.geneid`
- the June 2016 historical human UniProt accession set
- the June 1, 2016 GO ontology
- the independently recovered GraphSAGE-node-to-Entrez mapping

## Gene-identity handling

### Preserve many-to-many mappings

Treat the Entrez–UniProt mapping as a bipartite graph, not a dictionary.

Do not:

- arbitrarily select one accession per GeneID;
- arbitrarily select one GeneID per accession;
- discard every ambiguous component;
- restrict the mapping graph to GraphSAGE genes before resolving components.

Retain legitimate one-to-many and many-to-one mappings.

Examples include:

- P69905 → GeneIDs 3039 and 3040
- P62158 → GeneIDs 801, 805, and 808
- P62805 → multiple histone H4 GeneIDs
- P0DMV8/P0DMV9 ↔ HSPA1A/HSPA1B

Use a primary-symbol match only when an accession has no usable historical GeneID edge and the symbol match is unique. Resolve a square many-to-many component only when the symbols provide a unique bijection.

### O95073/Q9Y620 correction

The complete historical mapping component includes:

- Q9Y620 ↔ GeneID 25788, symbol RAD54B
- O95073 ↔ GeneID 100861412, symbol FSBP
- an additional historical O95073 ↔ 25788 edge

Resolve this component before restricting to GraphSAGE genes.

Retain:

```text
Q9Y620 → 25788
O95073 → 100861412
```

Do not project O95073/FSBP annotations to GeneID 25788/RAD54B.

Filtering to the GraphSAGE universe too early hides GeneID 100861412 and incorrectly assigns O95073 to 25788. That error creates 13 false-positive label cells.

The correction is based on the complete mapping component and concordant gene symbols, not on the observed labels.

## GOA annotation policy

For every GOA release-159 GAF row:

1. Exclude `NOT` annotations.
2. Retain only these evidence codes:

```text
EXP
IDA
IEP
IGI
IMP
ISS
```

3. Retain only the default relation for each ontology aspect:

```text
Biological Process:   involved_in
Cellular Component:   part_of
Molecular Function:   enables
```

4. Do not collapse these qualified relations into ordinary binary membership:

```text
colocalizes_with
contributes_to
```

5. Canonicalize alternate GO IDs.
6. Propagate each retained annotation to the term and all transitive `is_a` ancestors.
7. Do not propagate through ontology `part_of` edges.

## Why the relation policy matters

An earlier reconstruction treated every positive GAF qualifier as ordinary membership and produced 901 excess gene-label assignments.

Their decomposition was:

- 501 assignments dependent on `colocalizes_with`
- 387 assignments dependent on `contributes_to`
- 13 assignments caused by the O95073→25788 mapping error

After excluding the two qualified relations and correcting the mapping component:

```text
false positives = 0
false negatives = 0
exact columns   = 121/121
```

## The labels are not all Biological Process

The exact GO candidates contain:

- 85 Biological Process terms
- 26 Cellular Component terms
- 10 Molecular Function terms

Do not describe the dataset as containing 121 Biological Process labels.

## 121 columns versus 118 distinct vectors

The deposited matrix has 121 columns but only 118 distinct membership vectors.

Three duplicated vector pairs each admit two exact GO terms:

```text
columns 24 and 71:
    GO:0043228
    GO:0043232

columns 39 and 63:
    GO:0006464
    GO:0036211

columns 48 and 70:
    GO:0043230
    GO:1903561
```

The union contains 121 distinct exact GO-term candidates.

The matrix alone cannot identify which candidate was intended for each duplicated column. A source term list or preprocessing script is needed.

## Term-selection reconstruction

Do not calculate prevalence only inside the GraphSAGE gene universe.

Apply the exact transformation above to the full historical human GOA release-159 universe.

This gives:

```text
19,056 historical human GeneIDs
16,338 nonempty propagated GO terms
```

The 121 exact GraphSAGE GO candidates are exactly:

```text
the 121 most prevalent full-human propagated GO terms
```

They are also exactly:

```text
all propagated GO terms annotating at least 1,000
distinct historical human genes or proteins
```

The boundary is:

```text
rank 121: 1,007 genes
rank 122:   997 genes
```

The result is reproduced independently by a standard-library/set implementation.

Counting distinct Swiss-Prot/GPI accessions instead of GeneIDs gives the same 121-term set, so the original preprocessing may have described the threshold in terms of proteins.

A 24-configuration grid varied evidence filters, relation handling, and propagation. Only the exact six-code/default-relation/`is_a` policy recovered either the exact top-121 set or the exact ≥1,000 set.

Release 158 produces the same top-121 set but only 120 terms at the ≥1,000 threshold. Membership vectors strongly favor release 159.

## MSigDB result

Direct Entrez memberships from MSigDB versions 5.1, 5.2, and 6.0 do not reproduce any GraphSAGE label column exactly.

This remains true when searching:

- C5 alone
- every MSigDB collection

Closest mismatches:

```text
MSigDB 5.1: 307 genes
MSigDB 5.2:  10 genes
MSigDB 6.0:  10 genes
```

MSigDB may have supplied term identities or metadata, while memberships were regenerated from GOA. MSigDB 5.0 remains untested.

## Important cautions

Do not:

- force mappings to one-to-one;
- discard all ambiguous accessions;
- filter mappings to GraphSAGE genes before component resolution;
- treat `colocalizes_with` as ordinary Cellular Component membership;
- treat `contributes_to` as ordinary Molecular Function membership;
- propagate ontology `part_of`;
- call all 121 labels Biological Process terms;
- equate 118 distinct vectors with 118 unambiguously identified source terms;
- claim that direct MSigDB memberships reproduce the labels.

## Highest-priority remaining questions

1. Was the source rule explicitly “at least 1,000 genes,” “at least 1,000 proteins,” or “top 121 terms”?
2. What determined the label-column order?
3. Which GO term was intended for each duplicated-vector column?
4. Did MSigDB provide only a term list?
5. What caused the historical O95073→25788 edge?
6. Does the exact late-June ontology preserve the exact result?
7. Can the remaining anonymous graph genes be identified without consulting labels?

## Next tests

1. B105: late-June or July 1, 2016 `go.obo`.
2. MSigDB v5.0.
3. `org.Hs.eg.db` 3.3.0 and 3.4.0.
4. UniProt Swiss-Prot releases 2016_04, 2016_05, and 2016_06 for O95073 and Q9Y620.
5. Search preprocessing code and intermediate files for constants such as `1000` or `121`, and for the column-order source.
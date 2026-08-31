# GraphSAGE PPI GO labels: MSigDB tests, term-selection reconstruction, and remaining provenance questions

## Scope and terminology

The exact label equality established in B104A concerns the collapsed matrix for the 4,268 independently resolved Entrez Gene IDs. It contains 4,268 × 121 = 516,428 binary cells. Re-expanding by the recovered row-to-gene map reproduces every GraphSAGE node row whose biological identity was independently resolved. This does not assign Entrez IDs to the small set of still-unresolved graph nodes.

The 121 deposited columns contain 118 distinct binary membership vectors. Three vector pairs are duplicated. Each duplicated vector has two distinct GO terms whose memberships are identical within the resolved GraphSAGE gene universe:

| GraphSAGE columns | Exact GO candidates |
|---|---|
| 24 and 71 | GO:0043228 and GO:0043232 |
| 39 and 63 | GO:0006464 and GO:0036211 |
| 48 and 70 | GO:0043230 and GO:1903561 |

Consequently, the matrix alone identifies 121 distinct exact GO-term candidates, but it does not uniquely assign one of the two candidates to each member of the three duplicate column pairs.

## Exact membership reconstruction retained from B104A

The exact reconstruction uses:

1. GOA human release 159.
2. The May 2016 `gp2protein.geneid` Entrez–UniProt mapping.
3. A bipartite, ambiguity-preserving identifier model rather than a forced one-to-one map.
4. Resolution of mapping components before restricting to GraphSAGE genes.
5. Evidence codes `EXP`, `IDA`, `IEP`, `IGI`, `IMP`, and `ISS`.
6. Only the default aspect relations: Biological Process `involved_in`, Cellular Component `part_of`, and Molecular Function `enables`.
7. Exclusion of `NOT`, `colocalizes_with`, and `contributes_to` from ordinary binary membership.
8. Canonicalization of GO alternate IDs.
9. Propagation through `is_a` edges only.

Under this policy, all 121 deposited columns match exactly on all 4,268 resolved genes.

## Can direct MSigDB memberships reproduce the labels?

No direct match was found in any uploaded version.

| MSigDB version | Scope searched | Sets tested | Exact columns | Closest column/set mismatch | Columns at least 99% | Columns at least 95% |
|---|---:|---:|---:|---:|---:|---:|
| 5.1 | C5 | 1,454 | 0 | 307 genes | 0 | 0 |
| 5.1 | all collections | 13,311 | 0 | 307 genes | 0 | 0 |
| 5.2 | C5 | 6,166 | 0 | 10 genes | 1 | 2 |
| 5.2 | all collections | 18,890 | 0 | 10 genes | 1 | 2 |
| 6.0 | C5 | 5,917 | 0 | 10 genes | 1 | 2 |
| 6.0 | all collections | 18,643 | 0 | 10 genes | 1 | 2 |

The test used Entrez memberships directly from each MSigDB archive and compared them with the same 4,268-gene universe. Thus, the uploaded MSigDB 5.1, 5.2, and 6.0 archives cannot be the direct source of the deposited membership vectors.

This does not rule out a weaker role for MSigDB. It may have provided term names, a candidate list, or provenance metadata while memberships were regenerated from GOA. The GraphSAGE paper says the 121 GO sets were collected from MSigDB, and a later GraphSAGE GitHub response likewise identifies MSigDB as the feature/label information source. The exact membership evidence instead points to GOA-derived memberships.

MSigDB 5.0 remains untested because its archive is present only on the user's machine. It should be analyzed as a separate batch. No additional chip file is needed for this direct Entrez-membership test; the v5.2 chip archive is already available.

## New exact reconstruction of the term-selection set

The term selection can now be reconstructed exactly when prevalence is calculated in the full historical human annotation universe rather than only among GraphSAGE genes.

### Procedure

1. Retain all GPI159-linked historical GeneID–UniProt edges, preserving many-to-many mappings.
2. Remove the cross-symbol `O95073 → 25788` edge while retaining `Q9Y620 → 25788` and `O95073 → 100861412`.
3. Apply the exact B104A evidence and relation policy to GAF159.
4. Propagate every accepted annotation through `is_a` only.
5. Count distinct historical human GeneIDs for every propagated GO term across the full GPI159-linked human universe.

The full universe contains 19,056 historical human GeneIDs linked to GPI159 accessions. There are 16,338 propagated GO terms with nonzero membership under the accepted policy.

### Result

The 121 most prevalent propagated GO terms are exactly the 121 exact GO-term candidates identified from the GraphSAGE label matrix:

- true positives: 121
- extra terms: 0
- missing terms: 0

The same set is also obtained by the natural round cutoff:

```text
at least 1,000 distinct historical human GeneIDs
```

The boundary is sharp:

| Rank | GO ID | Name | Full-human GeneIDs |
|---:|---|---|---:|
| 120 | GO:0051254 | positive regulation of RNA metabolic process | 1,016 |
| 121 | GO:0031399 | regulation of protein modification process | 1,007 |
| 122 | GO:0043167 | ion binding | 997 |
| 123 | GO:0051049 | regulation of transport | 991 |

Any threshold from 998 through 1,007 produces the same set; 1,000 is the most plausible round threshold.

This result was independently reproduced by a second implementation using only Python standard-library CSV parsing and ordinary sets. It did not use pandas, NumPy, or the integer-bitset implementation used for the first analysis.

### Counting-unit sensitivity

The same exact 121-term set is obtained by either:

- distinct historical GeneIDs; or
- distinct GPI159/Swiss-Prot accessions.

For accessions, the rank-121 term has 1,006 proteins and rank 122 has 990. Therefore the data do not yet distinguish whether the original code counted genes or proteins.

### Release sensitivity

Using release 158 annotations with the same mapping, relation, evidence, and ontology logic:

- the top 121 terms are still the same exact set;
- a threshold of at least 1,000 selects only 120 terms because GO:0031399 has 995 genes.

Thus, “top 121 by full-human prevalence” is stable across releases 158 and 159, whereas the exact round threshold of 1,000 points more specifically to release 159. Label memberships themselves strongly favor release 159.

### Aspect composition

The selected set contains:

- 85 Biological Process terms;
- 26 Cellular Component terms;
- 10 Molecular Function terms.

Selecting the top 85, 26, and 10 terms within those three aspects also recovers the exact set, but this is less parsimonious than the single global prevalence rule because the aspect counts emerge automatically from the global ranking.

Direct-annotation prevalence does not work: the top 121 terms by direct membership overlap only 13 of the 121 candidates. Ontology propagation is therefore essential to both membership reconstruction and term selection.

## Policy-robustness test of the selection rule

To test whether the exact 121-term prevalence result was merely a generic consequence of broad GO terms, I repeated the full-human selection analysis across 24 processing combinations:

- six evidence policies;
- default-only versus all positive GAF relations;
- no propagation versus `is_a` propagation.

Only one combination reproduced either the exact top-121 set or the exact ≥1,000 set:

```text
Evidence: EXP, IDA, IEP, IGI, IMP, ISS
Relations: involved_in, part_of, enables
Propagation: is_a
```

Every alternative evidence, relation, or propagation policy introduced at least one extra or missing term. This is independent support for the same transformation that produced the exact gene-by-label membership matrix; the term-selection match is not simply a robust property of many plausible GO processing choices.

## Current best hypothesis for why 121 terms were selected

The strongest hypothesis is:

> Generate full-human GO memberships from GOA release 159 using the six evidence codes, default aspect relations, and `is_a` propagation; then retain every GO term annotating at least 1,000 human genes/proteins.

This yields exactly 121 terms. An equivalent implementation could simply take the 121 most prevalent terms.

The remaining distinction between these formulations can be tested by:

1. repeating the count on GOA release 160;
2. repeating it with the exact late-June ontology used by GAF159;
3. testing older GOA releases to locate when the boundary crossed 1,000;
4. searching preprocessing code or intermediate term lists for an explicit constant such as `1000` or `121`;
5. examining MSigDB 5.0 to see whether a source list independently contains all 121 term identities;
6. testing column order against source-file order, GO ID order, prevalence rank, and historical language/runtime iteration order.

The deposited column order is not explained by prevalence rank, GO ID order, ontology depth, or the tested MSigDB order fields. The low order correlations make a random or implementation-dependent ordering plausible, but that is not yet established.

## O95073 / GeneID 25788: files needed for date-matched confirmation

The exact label reconstruction requires not projecting O95073/FSBP annotations onto GeneID 25788/RAD54B. Current UniProt and NCBI records support:

- `Q9Y620 = RAD54B = GeneID 25788`
- `O95073 = FSBP = GeneID 100861412`

The historical `gp2protein.geneid` row connecting O95073 to 25788 should nevertheless be audited in date-matched sources rather than dismissed solely using current records.

### Highest-value official UniProt snapshots

Download the Swiss-Prot-only archives for:

- UniProt 2016_04, released 13 April 2016;
- UniProt 2016_05, released 11 May 2016;
- UniProt 2016_06, released 8 June 2016.

These bracket the May 9 `gp2protein` generation date and the July 2016 GOA release. The archives are approximately 1.5 GB each and contain reviewed entries, which is sufficient because O95073 and Q9Y620 are Swiss-Prot records. Inspect the `GN` and `DR   GeneID;` lines for both accessions.

### Smaller date-matched alternatives already available locally

Analyze these before transferring multi-gigabyte UniProt archives:

- `bioconductor-annotation-org.Hs.eg.db_3.3.0.tar.gz` — May 2016-era Entrez/UniProt maps;
- `bioconductor-annotation-org.Hs.eg.db_3.4.0.tar.gz` — later 2016 comparison.

Also search for NCBI snapshots around April–July 2016 containing:

- `gene_info.gz` or `Homo_sapiens.gene_info.gz`;
- `gene2refseq.gz`;
- `gene_refseq_uniprotkb_collab.gz`.

The latter two can be joined through RefSeq protein accessions to test the Entrez–UniProt relationship independently.

## B105 ontology inputs

First try the exact inferred OBO counterpart of the ontology date in the GAF159 header:

```bash
curl -fL \
  -o 2016-06-29-go.obo \
  http://purl.obolibrary.org/obo/go/releases/2016-06-29/go.obo
```

If that path does not resolve, download the official nearest monthly archive:

```bash
curl -fL \
  -o 2016-07-01-go.obo \
  https://release.geneontology.org/2016-07-01/ontology/go.obo
```

B105 is now a robustness test. The May-31 ontology already gives exact membership and term-selection matches, so a different late-June closure is no longer needed to explain the deposited matrix.

## Remaining open questions

1. Was the source rule explicitly “at least 1,000 proteins,” “at least 1,000 genes,” or “top 121 terms”?
2. What source determined the 121-column order?
3. Which GO ID was intended for each member of the three duplicate-vector column pairs?
4. Did MSigDB supply only a term list, or is the paper's MSigDB statement simply inaccurate?
5. What caused the historical O95073→25788 cross-reference, and in which source/release did it first appear or disappear?
6. Can the exact result be reproduced with the late-June ontology and GOA release 160?
7. Can the remaining unresolved graph nodes be biologically identified without consulting the deposited labels?

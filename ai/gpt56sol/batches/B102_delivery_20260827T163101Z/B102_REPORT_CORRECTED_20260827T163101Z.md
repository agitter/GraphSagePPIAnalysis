# Batch B102 — June 2016 gp2protein mapping verification

Corrected delivery generated: `2026-08-27T16:31:01.375910+00:00`  
Canonical scientific run: `20260827T162132Z` (`exit 0`, 63.795 seconds)

## Correction to the initial report

The three uploaded GO data files matched the size and SHA-256 values recorded in the 65-row full inventory. The inventory CSV itself cannot contain its own stable hash and therefore did **not** “match an inventory row.” Its received copy was independently hashed as `4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4`. No scientific result changes because the accepted parser already treated self-verification as not applicable.

## Input integrity and provenance

| file | bytes | sha256 | inventory comparison | compression test |
| --- | --- | --- | --- | --- |
| 2016-06-01-annotations-README | 1,999 | 32134f9555d6710a9bb488fe652fef04cd84facf516f669939487d137f8bcc88 | exact match | not applicable |
| 2016-06-01-gp2protein.geneid.gz | 38,540,366 | f3a2d329ada32f03e4c3ec60c55ef77cfe3626c3c875832e59612e5f316504e7 | exact match | gzip -t passed |
| 2016-06-01-gp2protein.human.gz | 378,274 | 5a62823541d718c212b61efe741b61f67b10c9a8fb71114f9cd3e33f9cc501dd | exact match | gzip -t passed |
| local_upload_inventory_full_20260827T160408Z.csv | 12,067 | 4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4 | not applicable; inventory cannot self-list | not applicable |

The official historical archive paths are recorded in the updated manifests. Remote bytes were not re-downloaded in the container, so the report does not claim remote-byte identity.

## Full local inventory

The unfiltered user inventory contains **65 files** totaling **2,680,734,828 bytes**. Every row has an exact SHA-256 declaration. The exact spelling `HuamnBase-kidney.dat` remains preserved. This supersedes the earlier pattern-filtered inventory for statements about what the user has locally.

## Concrete file semantics

- `gp2protein.human.gz`: **70,625 rows**, all UniProtKB self-maps; it defines a historical human-accession set but contains no Entrez IDs.
- `gp2protein.geneid.gz`: generated on `2016-05-09`, with **7,296,170 valid data rows** and **0 malformed rows**. It uses one GeneID–UniProt pair per line.
- Human filtering by the self-map yields **23,046 rows / 23,046 unique GeneID–UniProt pairs**.
- Of **70,625** historical human accessions, **22,677** have a GeneID link and **47,948** do not.

## Relationship to GOA release 159

- GPI objects: **21,002**; **21,001** occur in the historical human accession set. The sole exception is `P0DO97`.
- GPI objects with any full historical GeneID link: **18,869**.
- Annotated GAF objects with any full historical GeneID link: **18,094 / 19,194**.
- Nonempty GAF `Gene_Product_Form_ID` values: **0**.

## GraphSAGE gene coverage

Using all historical links from GPI objects covers **4,259 / 4,268** independently resolved GraphSAGE genes. The nine uncovered Entrez IDs are `176, 337, 3108, 4018, 7957, 10159, 29901, 55125, 84919`.

Five GPI accessions map to more than one resolved GraphSAGE GeneID: `P69905`, `P0DMV8`, `P0DMV9`, `P62158`, and `P62805`. The retained missing-gene table distinguishes absent GeneID rows, non-reference accession mappings, and reference-accession gaps.

## Direct GO-label reconstruction

The direct, non-propagated screen tested **6 mapping strategies × 9 evidence filters × 5 term scopes × 2 comparison scopes × 121 label columns**, yielding **65,340** per-label result rows.

The best full-universe unrestricted-BP configuration used `gp2protein_all_graphsage_links` with `all_except_IEA_ND_NAS`. Its median best agreement was **81.6776%**, its closest column still differed at **280 genes**, and it produced **0/121** columns at 95% agreement and **0/121** exact columns.

The closest result anywhere in the direct grid was label column **72** versus **`GO:0045944`**, with **273 mismatched genes** (93.6036% agreement; TP=387, FP=7, FN=266).

**Conclusion:** the historical May-2016 gp2protein mapping improves provenance and explains almost all resolved GraphSAGE genes, but it does not make direct GOA v159 annotations reproduce the 121 labels. The next discriminating transformation is ontology propagation, followed by release-specific term selection or Entrez-native annotation products if propagation remains insufficient.

## Validation and retained evidence

- A separate brute-force matcher reproduced all 121 optimized best-term choices and distances for the main historical configuration.
- The canonical derivatives were re-read successfully: 70,625 human self-map rows; 25,983 retained relevant GeneID rows; 21,002 GPI mapping rows.
- A later duplicate successful run produced a decompressed grid identical to the accepted grid: `True`; its decompressed GPI mapping table was also identical: `True`.
- Complete failed-attempt and duplicate-run diagnostics are included rather than referenced without being present.

## Next analysis

Batch B103 should contain `2016-06-01-go.obo`. It will test no propagation, `is_a` propagation, `is_a + part_of`, alternative relations, obsolete/alternate identifiers, and the 121-term selection criterion. The GAF v159 header identifies a later June 2016 ontology date, so the June 1 ontology is a near-date test rather than an assumed exact companion; any remaining near-match will trigger acquisition of the exact header-specified or next monthly ontology.

# Label Provenance Analysis — Current Status

## Gene coverage breakdown (verified, single consistent accounting)

Starting point: 4,278 mapped genes with known Entrez IDs.

| Category | Count | Description |
|---|---|---|
| Have GOA BP annotations (after evidence filter) | 3,706 | Can be matched against GraphSAGE labels |
| In GOA but only excluded evidence codes (IEA/IBA/NAS/TAS/ND) | 549 | Have BP annotations but only from excluded codes |
| Not in GOA at all | 23 | No UniProt accession in GOA v159 maps to these Entrez IDs |
| **Total** | **4,278** | |

The 23 genes not in GOA include 13 histone H4 variants (8294, 8359-8368, 121504, 554313),
2 calmodulins (801, 805), and 8 other paralogs/pseudogenes. These are genes where NCBI
assigns separate Entrez IDs but UniProt maps them to a canonical family member.

## Best label match: GOA v159, evidence filter {IEA,IBA,NAS,TAS,ND} excluded

| Match quality | Columns |
|---|---|
| ≥99% agreement | 77 of 121 |
| 95-99% agreement | 8 of 121 |
| <95% agreement | 36 of 121 |

## Concrete example: Column 0 = GO:0050789 (regulation of biological process)

Agreement: 99.23% (4,245 of 4,278 genes agree)

### 12 false positives (GOA says annotated, GraphSAGE says label=0):
43, 624, 1215, 3248, 3848, 9702, 10928, 27201, 54576, 54658, 55801, 79849

These genes gained this annotation between the data freeze and GOA v159.

### 21 false negatives (GraphSAGE says label=1, GOA says not annotated):
801, 805, 1139, 3039, 3303, 4831, 8294, 8359, 8360, 8361, 8362, 8363,
8364, 8365, 8366, 8367, 8368, 30061, 121504, 121599, 554313

- 19 are from the 23 genes NOT IN GOA AT ALL (paralogs/histones)
- 2 are in GOA but lack this specific GO term (30061, 121599)

## Why 36 columns remain below 95%

These columns have 200-1400 false negatives. The false negatives come from two sources:
1. The 23 genes not in GOA (contributes ~20 FN per column for columns where they have labels)
2. The 549 genes with only excluded evidence codes (these genes have NO BP annotations
   in our reconstruction, so if GraphSAGE gives them label=1, we produce a false negative)

The 549 genes are the main bottleneck, not the 23 paralogs. These are genes that have
GO BP annotations in gene2go (NCBI's native Entrez-based file) but whose annotations
come exclusively from evidence codes we're filtering out (IEA, IBA, NAS, TAS, ND).

This is a CONTRADICTION: if GraphSAGE labels were generated with the same evidence
filter, these 549 genes should have label=0. But they have non-zero labels. This means
either:
  (a) The evidence filter is wrong — the labels include some codes we're excluding
  (b) The labels come from gene2go where the same gene has DIFFERENT evidence codes
      than in GOA (NCBI and EBI annotate independently)
  (c) Our UniProt-to-Entrez mapping is assigning annotations to the wrong gene

Option (b) is most likely: gene2go and GOA are independent databases. A gene might have
experimental evidence in gene2go but only IEA in GOA, or vice versa.

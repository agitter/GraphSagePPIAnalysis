# Identifier mapping watchlist

This table separates direct 2016 mapping evidence, GPI membership, current UniProt state, and contextual symbol fallbacks. It deliberately does not collapse many-to-many relationships or infer replacement merely because accessions differ.

| GeneID | Symbol | Accession | Historical evidence | GPI158/159 | Current interpretation | Decision |
|---:|---|---|---|---|---|---|
| 7957 | EPM2A | `B3EWF7` | none | 1/1 | Distinct reviewed EPM2A isoform entry; coexists with O95278 | Retain as a historically contextual fallback. It yields the exact 17-label row under the accepted model; do not call it a replacement for O95278. |
| 7957 | EPM2A | `O95278` | gp2protein | 0/0 | Canonical laforin entry; coexists with B3EWF7 | Historical GeneID edge is direct but the accession is absent from GPI158/159, so it contributes no GOA rows in these reference-proteome files. |
| 7957 | EPM2A | `H0UI04` | gp2protein | 0/0 |  | Keep recorded; absent from historical human self-map and GPI158/159. |
| 29901 | SAC3D1 | `A6NKF1` | none | 1/1 | Reviewed SAC3D1 entry; current page lists F8WC89/A0A6I8PRW4/H9KVA8 as mapped potential isoforms | Retain as historically contextual fallback. No accepted-evidence GAF159 rows survive, and the observed label row is all zero. |
| 29901 | SAC3D1 | `F8WC89` | gp2protein | 0/0 | Potential isoform mapped to A6NKF1 on current UniProt page | Track without treating as a replacement relationship; absent from GPI158/159 and irrelevant to current zero label row. |
| 29901 | SAC3D1 | `A0A6I8PRW4` | none | 0/0 | Potential isoform mapped to A6NKF1 on current UniProt page | Track without treating as a replacement relationship; absent from GPI158/159 and irrelevant to current zero label row. |
| 29901 | SAC3D1 | `H9KVA8` | none | 0/0 | Potential isoform mapped to A6NKF1 on current UniProt page | Track without treating as a replacement relationship; absent from GPI158/159 and irrelevant to current zero label row. |
| 10159 | ATP6AP2 | `O75787` | gp2protein | 0/0 | Current primary accession; GeneID 10159; MANE NM_005765.3 / NP_005756.2 | Canonical anchor. Absent from GPI158/159. The observed GraphSAGE row is all zero, so absence causes no false negative. |
| 10159 | ATP6AP2 | `B7Z9I3` | none | 0/0 | Current secondary accession of O75787; it is not one of the eight accessions named in the 2005 replacement event | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q5QTQ7` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q6T7F5` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q8NBP3` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q8NG15` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q96FV6` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q96LB5` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q9H2P8` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `Q9UG89` | none | 0/0 | Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10 | Do not treat as an independent current reviewed mapping. It is absent from GPI158/159. |
| 10159 | ATP6AP2 | `A0A1C7CYW4` | none | 0/0 | Current alternative ATP6AP2 product | Absent from GPI158/159; not evidence for a 2016 reference-proteome mapping. |
| 10159 | ATP6AP2 | `Q8NBJ9` | none | 1/1 | GPI158/159 symbol SIDT2 includes synonym PSEC0072 | Never map ATP6AP2 via PSEC0072 alone; that synonym would collide with SIDT2/Q8NBJ9 in the historical GPI. |

## ATP6AP2 conclusion

O75787 is the defensible accession anchor. The many accessions shown by NCBI are not all independent primary Swiss-Prot records: UniProt records most of the Q-accessions as secondary or replaced accessions of O75787. None of the listed ATP6AP2 accessions appears as a GPI158 or GPI159 object. Because GeneID 10159 has an all-zero GraphSAGE label row, leaving it unmapped does not explain any residual false negative. A targeted historical reference-proteome provenance check remains useful, but it is not on the critical path for the 121-label reconstruction.

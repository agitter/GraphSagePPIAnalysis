#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, json
from pathlib import Path

STAMP='20260828T030759Z'
ROOT=Path('/mnt/data/ppi_repro_corrected/batches/B104_20260828T030759Z')
ANA=ROOT/'analysis'
GPI158=ROOT/f'derived/B104_goa_human_gpi158_normalized_{STAMP}.tsv.gz'
GPI159=ROOT/'retained_inputs/B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz'
GAF159=ROOT/'retained_inputs/B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
HIST=ROOT/'retained_inputs/B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
CUR=ROOT/f'derived/B104_repaired_B103_current_UniProt_idmapping_2026-08-27_{STAMP}.tsv.gz'
LABELS=ROOT/'retained_inputs/collapsed_gene_labels_topology_features.csv'
BEST={'EXP','IDA','IEP','IGI','IMP','ISS'}

# Load local facts

def load_gpi(p):
 d={}
 with gzip.open(p,'rt') as f:
  for r in csv.DictReader(f,delimiter='\t'): d[r['DB_Object_ID']]=r
 return d
g158=load_gpi(GPI158); g159=load_gpi(GPI159)
hist={}
with gzip.open(HIST,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  hist.setdefault(int(r['GeneID']),[]).append(r)
cur={}
with gzip.open(CUR,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'): cur.setdefault(int(r['From']),[]).append(r)
labels={}
with LABELS.open() as f:
 for r in csv.DictReader(f): labels[int(r['entrez_gene_id'])]=[int(r[f'label_{j}']) for j in range(121)]
gaf={}
with gzip.open(GAF159,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['Is_NOT']=='1' or 'NOT' in (r['Qualifier'] or '').split('|'): continue
  gaf.setdefault(r['DB_Object_ID'],[]).append(r)

rows=[]
def add(gene,sym,acc,role,current_status='',current_relation='',note='',source_url=''):
 anns=[r for r in gaf.get(acc,[]) if r['Evidence_Code'] in BEST]
 rows.append({
  'GeneID':gene,'gene_symbol':sym,'accession':acc,'role_in_audit':role,
  'historical_gp2protein_edge':int(any(r['UniProtKB_accession']==acc for r in hist.get(gene,[]))),
  'historical_human_self_map':next((r['in_human_self_map'] for r in hist.get(gene,[]) if r['UniProtKB_accession']==acc),''),
  'present_in_GPI158':int(acc in g158),'GPI158_symbol':g158.get(acc,{}).get('DB_Object_Symbol',''),
  'present_in_GPI159':int(acc in g159),'GPI159_symbol':g159.get(acc,{}).get('DB_Object_Symbol',''),
  'current_query_returned_for_GeneID_2026_08_27':int(any(r['Entry']==acc for r in cur.get(gene,[]))),
  'current_query_reviewed_field':next((r.get('Reviewed','') for r in cur.get(gene,[]) if r['Entry']==acc),''),
  'current_status_from_official_UniProt_web':current_status,
  'current_relation_or_history':current_relation,
  'best_evidence_GAF159_rows':len(anns),
  'best_evidence_GAF159_direct_GO_IDs':'|'.join(sorted({r['GO_ID'] for r in anns})),
  'GraphSAGE_positive_label_count':sum(labels.get(gene,[])),
  'mapping_decision_or_note':note,
  'official_source_url':source_url
 })

add(7957,'EPM2A','B3EWF7','GPI158/159 reference-proteome accession and unique primary-symbol fallback','reviewed','Distinct reviewed EPM2A isoform entry; coexists with O95278','Retain as a historically contextual fallback. It yields the exact 17-label row under the accepted model; do not call it a replacement for O95278.','https://www.uniprot.org/uniprotkb/B3EWF7')
add(7957,'EPM2A','O95278','May-2016 gp2protein and current canonical EPM2A accession','reviewed','Canonical laforin entry; coexists with B3EWF7','Historical GeneID edge is direct but the accession is absent from GPI158/159, so it contributes no GOA rows in these reference-proteome files.','https://www.uniprot.org/uniprotkb/O95278')
add(7957,'EPM2A','H0UI04','May-2016 non-human-self-map historical edge','','','Keep recorded; absent from historical human self-map and GPI158/159.','')
add(29901,'SAC3D1','A6NKF1','GPI158/159 reference-proteome accession and unique primary-symbol fallback','reviewed','Reviewed SAC3D1 entry; current page lists F8WC89/A0A6I8PRW4/H9KVA8 as mapped potential isoforms','Retain as historically contextual fallback. No accepted-evidence GAF159 rows survive, and the observed label row is all zero.','https://www.uniprot.org/uniprotkb/A6NKF1')
for a in ['F8WC89','A0A6I8PRW4','H9KVA8']:
 add(29901,'SAC3D1',a,'May-2016 or current alternative SAC3D1 accession',('unreviewed' if a!='A6NKF1' else ''),'Potential isoform mapped to A6NKF1 on current UniProt page','Track without treating as a replacement relationship; absent from GPI158/159 and irrelevant to current zero label row.',('https://www.uniprot.org/uniprotkb/'+a))
# ATP6AP2 current primary, secondary accessions, and unreviewed product reported by NCBI
add(10159,'ATP6AP2','O75787','May-2016 gp2protein and current primary reviewed accession','reviewed','Current primary accession; GeneID 10159; MANE NM_005765.3 / NP_005756.2','Canonical anchor. Absent from GPI158/159. The observed GraphSAGE row is all zero, so absence causes no false negative.','https://www.uniprot.org/uniprotkb/O75787')
add(10159,'ATP6AP2','B7Z9I3','Current O75787 secondary accession','','Current secondary accession of O75787; it is not one of the eight accessions named in the 2005 replacement event','Do not treat as an independent current reviewed mapping. It is absent from GPI158/159.','https://www.uniprot.org/uniprotkb/O75787')
for a in ['Q5QTQ7','Q6T7F5','Q8NBP3','Q8NG15','Q96FV6','Q96LB5','Q9H2P8','Q9UG89']:
 add(10159,'ATP6AP2',a,'O75787 secondary/replaced accession','','Explicitly replaced by O75787 in UniProt release 5.0 on 2005-05-10','Do not treat as an independent current reviewed mapping. It is absent from GPI158/159.', 'https://www.uniprot.org/uniprotkb/O75787/history')
add(10159,'ATP6AP2','A0A1C7CYW4','Current NCBI-listed UniProtKB/TrEMBL product','unreviewed','Current alternative ATP6AP2 product','Absent from GPI158/159; not evidence for a 2016 reference-proteome mapping.','https://www.uniprot.org/uniprotkb/A0A1C7CYW4')
add(10159,'ATP6AP2','Q8NBJ9','Unsafe historical synonym collision: SIDT2','','GPI158/159 symbol SIDT2 includes synonym PSEC0072','Never map ATP6AP2 via PSEC0072 alone; that synonym would collide with SIDT2/Q8NBJ9 in the historical GPI.','https://www.uniprot.org/uniprotkb/Q8NBJ9')

out=ANA/f'B104_identifier_mapping_watchlist_{STAMP}.csv'
with out.open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

md=ANA/f'B104_identifier_mapping_watchlist_{STAMP}.md'
with md.open('w',encoding='utf-8') as f:
 f.write('# Identifier mapping watchlist\n\n')
 f.write('This table separates direct 2016 mapping evidence, GPI membership, current UniProt state, and contextual symbol fallbacks. It deliberately does not collapse many-to-many relationships or infer replacement merely because accessions differ.\n\n')
 f.write('| GeneID | Symbol | Accession | Historical evidence | GPI158/159 | Current interpretation | Decision |\n|---:|---|---|---|---|---|---|\n')
 for r in rows:
  histtxt='gp2protein' if r['historical_gp2protein_edge'] else 'none'
  gpitxt=f"{r['present_in_GPI158']}/{r['present_in_GPI159']}"
  f.write(f"| {r['GeneID']} | {r['gene_symbol']} | `{r['accession']}` | {histtxt} | {gpitxt} | {r['current_relation_or_history']} | {r['mapping_decision_or_note']} |\n")
 f.write('\n## ATP6AP2 conclusion\n\n')
 f.write('O75787 is the defensible accession anchor. The many accessions shown by NCBI are not all independent primary Swiss-Prot records: UniProt records most of the Q-accessions as secondary or replaced accessions of O75787. None of the listed ATP6AP2 accessions appears as a GPI158 or GPI159 object. Because GeneID 10159 has an all-zero GraphSAGE label row, leaving it unmapped does not explain any residual false negative. A targeted historical reference-proteome provenance check remains useful, but it is not on the critical path for the 121-label reconstruction.\n')
print(out)
print(md)

#!/usr/bin/env python3
import argparse,csv,json,re
from pathlib import Path

def parse_record(rec):
 out={'id_line':'','primary_accession':'','all_accessions':[],'entry_version_date':'','entry_version':'','gene_name':'','geneids':[],'hgnc':[],'miscellaneous':[]}
 for line in rec.splitlines():
  if line.startswith('ID   '):out['id_line']=line[5:].strip()
  elif line.startswith('AC   '):
   vals=[x.strip() for x in line[5:].split(';') if x.strip()];out['all_accessions']+=vals
  elif line.startswith('DT   ') and 'entry version' in line:
   m=re.search(r'^(\d{2}-[A-Z]{3}-\d{4}), entry version (\d+)\.',line[5:].strip())
   if m:out['entry_version_date']=m.group(1);out['entry_version']=int(m.group(2))
  elif line.startswith('GN   Name='):
   out['gene_name']=line.split('Name=',1)[1].split(';',1)[0]
  elif line.startswith('DR   GeneID;'):
   out['geneids'].append(line.split(';')[1].strip())
  elif line.startswith('DR   HGNC;'):
   out['hgnc'].append(';'.join(x.strip() for x in line[5:].split(';')[:2]))
  elif line.startswith('CC   -!- MISCELLANEOUS:'):
   out['miscellaneous'].append(line.split(':',1)[1].strip())
 out['primary_accession']=out['all_accessions'][0] if out['all_accessions'] else ''
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input-root',type=Path,required=True);ap.add_argument('--output-csv',type=Path,required=True);ap.add_argument('--output-json',type=Path,required=True);a=ap.parse_args()
 rows=[]
 for rel in ('2016_04','2016_05','2016_06'):
  p=a.input_root/rel/f'O95073_Q9Y620_{rel}.dat';text=p.read_text()
  records=[r+'//\n' for r in text.split('//\n') if r.strip()]
  for rec in records:
   x=parse_record(rec);x['release']=rel;x['record_path']=str(p);x['geneids']='|'.join(x['geneids']);x['all_accessions']='|'.join(x['all_accessions']);x['hgnc']='|'.join(x['hgnc']);x['miscellaneous']='|'.join(x['miscellaneous']);rows.append(x)
 if len(rows)!=6:raise SystemExit(len(rows))
 a.output_csv.parent.mkdir(parents=True,exist_ok=True)
 with a.output_csv.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 byacc={}
 for r in rows:byacc.setdefault(r['primary_accession'],[]).append(r)
 out={'records':rows,'stable_geneid_crossreferences':{acc:len(set(r['geneids'] for r in rr))==1 for acc,rr in byacc.items()},'stable_gene_names':{acc:len(set(r['gene_name'] for r in rr))==1 for acc,rr in byacc.items()},'interpretation_constraints':[
  'The official reviewed UniProt records themselves linked O95073 to both GeneID 100861412 and GeneID 25788 in all three releases.',
  'Therefore the O95073-to-25788 edge is not a parser artifact and should not be described as disproven by date-matched UniProt.',
  'The exact GraphSAGE label reconstruction instead shows that the label-generating mapping behaved as if FSBP annotations were not assigned to RAD54B/GeneID 25788.',
  'This behavior could arise from symbol-aware component resolution or from an Entrez-native annotation source; the exact source remains unresolved.'
 ]}
 a.output_json.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()

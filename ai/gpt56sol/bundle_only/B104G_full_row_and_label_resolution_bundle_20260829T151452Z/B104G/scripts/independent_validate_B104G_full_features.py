#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np

def open_text(p):return gzip.open(p,'rt',newline='',encoding='utf-8') if str(p).endswith('.gz') else open(p,newline='',encoding='utf-8')
def main():
 p=argparse.ArgumentParser()
 for x in ['graphsage_zip','row_map','msigdb_normalized','feature_rule','output']:p.add_argument('--'+x.replace('_','-'),type=Path,required=True)
 a=p.parse_args()
 with open_text(a.row_map) as f: row={int(r['graphsage_row']):int(r['entrez_gene_id']) for r in csv.DictReader(f)}
 with zipfile.ZipFile(a.graphsage_zip) as z:obs=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')),allow_pickle=False).astype(np.uint8)
 sets={}
 with gzip.open(a.msigdb_normalized,'rt',newline='',encoding='utf-8') as f:
  for r in csv.DictReader(f,delimiter='\t'):sets[(r['category'],r['standard_name'])]={int(x) for x in r['member_Entrez_IDs'].split('|') if x}
 rule=[]
 with open(a.feature_rule,newline='',encoding='utf-8') as f:
  for r in csv.DictReader(f):rule.append((int(r['graphsage_feature_column']),r['collection'],r['set_name']))
 rule.sort(); gene_sets=[sets[(c,n)] for _,c,n in rule]
 exp=np.zeros_like(obs)
 for r,g in row.items():exp[r]=[g in s for s in gene_sets]
 d=np.argwhere(exp!=obs)
 result={'independent_implementation':True,'rows':len(row),'columns':50,'cells':int(exp.size),'mismatched_cells':int(len(d)),
  'exact_columns':int(sum(np.array_equal(exp[:,j],obs[:,j]) for j in range(50))),
  'observed_matrix_sha256':hashlib.sha256(obs.tobytes()).hexdigest(),'expected_matrix_sha256':hashlib.sha256(exp.tobytes()).hexdigest(),
  'first_mismatches':[{'row':int(r),'gene':row[int(r)],'column':int(c)} for r,c in d[:20]],'all_checks_pass':not len(d)}
 a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
 if len(d):raise SystemExit(1)
if __name__=='__main__':main()

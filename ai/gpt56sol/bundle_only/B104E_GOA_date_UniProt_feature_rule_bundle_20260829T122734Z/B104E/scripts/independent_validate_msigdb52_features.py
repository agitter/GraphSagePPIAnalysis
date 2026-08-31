#!/usr/bin/env python3
"""Independent loop-based validation of the exact MSigDB feature rule."""
import argparse,csv,io,json,zipfile
from pathlib import Path
import numpy as np

def main():
 ap=argparse.ArgumentParser();
 ap.add_argument('--graphsage-zip',type=Path,required=True);ap.add_argument('--msigdb-zip',type=Path,required=True);ap.add_argument('--row-map',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
 a=ap.parse_args()
 rows=[]
 with a.row_map.open(newline='') as f:
  for r in csv.DictReader(f): rows.append((int(r['graphsage_row']),int(r['entrez_gene_id'])))
 with zipfile.ZipFile(a.graphsage_zip) as z: x=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')),allow_pickle=False)
 members={
 'C1':'msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c1.all.v5.2.entrez.gmt',
 'C3':'msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c3.all.v5.2.entrez.gmt',
 'C7':'msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c7.all.v5.2.entrez.gmt'}
 selected=[]
 with zipfile.ZipFile(a.msigdb_zip) as z:
  for c in ('C1','C3','C7'):
   for line_number,line in enumerate(z.read(members[c]).decode().splitlines(),start=1):
    p=line.split('\t'); g=set()
    for value in p[2:]:
     try:g.add(int(value))
     except ValueError:pass
    if len(g)>=200 and len(selected)<50:selected.append((c,line_number,p[0],g))
   if len(selected)>=50:break
 if len(selected)!=50:raise SystemExit(len(selected))
 mismatch_examples=[]; mismatches=0; positives_observed=[0]*50; positives_expected=[0]*50
 for row,gene in rows:
  for j,(_,_,_,geneset) in enumerate(selected):
   observed=int(x[row,j]); expected=int(gene in geneset)
   positives_observed[j]+=observed; positives_expected[j]+=expected
   if observed!=expected:
    mismatches+=1
    if len(mismatch_examples)<20:mismatch_examples.append({'row':row,'gene':gene,'column':j,'observed':observed,'expected':expected})
 result={'implementation':'independent nested Python loops','resolved_rows':len(rows),'cells_checked':len(rows)*50,'mismatches':mismatches,'mismatch_examples':mismatch_examples,
         'selected_names':[r[2] for r in selected], 'selected_collections':[r[0] for r in selected],
         'observed_positive_counts':positives_observed,'expected_positive_counts':positives_expected,
         'selected_c1':sum(r[0]=='C1' for r in selected),'selected_c3':sum(r[0]=='C3' for r in selected),'selected_c7':sum(r[0]=='C7' for r in selected)}
 a.output.write_text(json.dumps(result,indent=2)+'\n')
 if mismatches or positives_observed!=positives_expected:raise SystemExit('validation failed')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()

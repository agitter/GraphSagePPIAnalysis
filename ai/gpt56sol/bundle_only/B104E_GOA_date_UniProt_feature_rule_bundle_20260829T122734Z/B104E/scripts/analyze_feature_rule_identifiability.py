#!/usr/bin/env python3
import argparse,csv,json,zipfile
from pathlib import Path

def parse(z,m,c):
 out=[]
 for i,l in enumerate(z.read(m).decode().splitlines()):
  p=l.split('\t');s=set()
  for x in p[2:]:
   try:s.add(int(x))
   except:pass
  out.append((c,p[0],len(s),i))
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--zip',type=Path,required=True);ap.add_argument('--out-csv',type=Path,required=True);ap.add_argument('--out-json',type=Path,required=True);a=ap.parse_args()
 mem={'C1':'msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c1.all.v5.2.entrez.gmt','C3':'msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c3.all.v5.2.entrez.gmt','C7':'msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c7.all.v5.2.entrez.gmt'}
 with zipfile.ZipFile(a.zip) as z:rows={c:parse(z,m,c) for c,m in mem.items()}
 def sel(t):
  out=[]
  for c in ('C1','C3','C7'):
   for r in rows[c]:
    if r[2]>=t and len(out)<50:out.append(r)
   if len(out)>=50:break
  return out
 baseline=sel(200);bn=[(x[0],x[1]) for x in baseline]
 tests=[]
 for t in range(1,1001):
  s=sel(t);tests.append({'threshold_ge':t,'selected_count':len(s),'same_ordered_50_as_threshold_200':int([(x[0],x[1]) for x in s]==bn),'selected_C1':sum(x[0]=='C1' for x in s),'selected_C3':sum(x[0]=='C3' for x in s),'selected_C7':sum(x[0]=='C7' for x in s)})
 with a.out_csv.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(tests[0]));w.writeheader();w.writerows(tests)
 valid=[r['threshold_ge'] for r in tests if r['same_ordered_50_as_threshold_200']]
 result={'thresholds_producing_same_ordered_selected_sets':valid,'minimum_selected_set_size':min(x[2] for x in baseline),'first_unselected_qualifying_after_cap':{'collection':sel(200)[-1][0],'last_selected_name':sel(200)[-1][1]},'interpretation':'The exact feature matrix does not distinguish >=200 from >=201 (equivalently >200), because both produce the same selected sequence. A cutoff described as 200 is natural but not source-code proven.'}
 a.out_json.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()

#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np

def sha(path):
 h=hashlib.sha256();
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def load_map(path):
 with path.open(newline='') as f:return {int(r['graphsage_row']):int(r['entrez_gene_id']) for r in csv.DictReader(f)}
def membership_hash(s):return hashlib.sha256(('\n'.join(map(str,sorted(s)))+'\n').encode()).hexdigest()
def select_from_normalized_v50(path):
 cats={'C1':[],'C3':[],'C7':[]}
 with gzip.open(path,'rt',newline='') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   c=r['category']
   if c not in cats:continue
   s={int(x) for x in r['member_Entrez_IDs'].split('|') if x}
   cats[c].append((int(r['order_all']),r['standard_name'],s))
 for c in cats:cats[c].sort()
 stream=[]
 for c in ('C1','C3','C7'):
  stream += [(c,n,s,o) for o,n,s in cats[c] if len(s)>=200]
 return stream[:50],{c:sum(len(s)>=200 for _,_,s in cats[c]) for c in cats}
def select_from_zip(path):
 with zipfile.ZipFile(path) as z:
  names=z.namelist();stream=[];qcounts={}
  for c in ('C1','C3','C7'):
   hits=[n for n in names if f'/{c.lower()}.all.' in n.lower() and n.lower().endswith('.entrez.gmt')]
   if len(hits)!=1:raise ValueError((c,hits))
   q=[]
   for i,l in enumerate(z.read(hits[0]).decode().splitlines()):
    p=l.split('\t');s=set()
    for x in p[2:]:
     try:s.add(int(x))
     except:pass
    if len(s)>=200:q.append((c,p[0],s,i))
   qcounts[c]=len(q);stream+=q
  return stream[:50],qcounts
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--graphsage-zip',type=Path,required=True);ap.add_argument('--row-map',type=Path,required=True);ap.add_argument('--v50-normalized',type=Path,required=True);ap.add_argument('--version-zip',action='append',nargs=2,required=True);ap.add_argument('--output-csv',type=Path,required=True);ap.add_argument('--output-json',type=Path,required=True);a=ap.parse_args()
 mp=load_map(a.row_map);rr=sorted(mp);genes=[mp[i] for i in rr]
 with zipfile.ZipFile(a.graphsage_zip) as z:x=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')),allow_pickle=False)[rr]
 inputs=[('5.0',a.v50_normalized,'retained normalized derivative')]+[(v,Path(p),'prior user-supplied raw archive') for v,p in a.version_zip]
 out=[];details={}
 for ver,path,source_type in inputs:
  selected,q=select_from_normalized_v50(path) if ver=='5.0' else select_from_zip(path)
  mism=[]
  for j,(_,_,s,_) in enumerate(selected):mism.append(sum(int(x[k,j])!=int(g in s) for k,g in enumerate(genes)))
  seqhash=hashlib.sha256(('\n'.join(membership_hash(s) for _,_,s,_ in selected)+'\n').encode()).hexdigest()
  out.append({'version':ver,'input_path':str(path),'input_type':source_type,'input_sha256':sha(path),'selected_C1':sum(c=='C1' for c,_,_,_ in selected),'selected_C3':sum(c=='C3' for c,_,_,_ in selected),'selected_C7':sum(c=='C7' for c,_,_,_ in selected),'qualifying_C1':q['C1'],'qualifying_C3':q['C3'],'qualifying_C7':q['C7'],'exact_columns':sum(m==0 for m in mism),'total_mismatches':sum(mism),'membership_sequence_sha256':seqhash,'first_name':selected[0][1],'last_name':selected[-1][1]})
  details[ver]=[{'column':j,'collection':c,'name':n,'unique_entrez_count':len(s),'membership_sha256':membership_hash(s),'mismatches':mism[j]} for j,(c,n,s,_) in enumerate(selected)]
 with a.output_csv.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 result={'rows':out,'details':details,'all_exact':all(r['exact_columns']==50 and r['total_mismatches']==0 for r in out),'all_selected_membership_sequences_identical':len({r['membership_sequence_sha256'] for r in out})==1,'provenance_note':'Version 5.0 was evaluated from the hash-verified complete normalized derivative retained before the user-confirmed raw deletion. It was not read from a residual raw mount.'}
 a.output_json.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'rows':out,'all_exact':result['all_exact'],'all_sequences_identical':result['all_selected_membership_sequences_identical']},indent=2))
 if not result['all_exact'] or not result['all_selected_membership_sequences_identical']:raise SystemExit(1)
if __name__=='__main__':main()

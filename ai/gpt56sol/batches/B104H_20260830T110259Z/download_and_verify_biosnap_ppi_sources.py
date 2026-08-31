#!/usr/bin/env python3
"""Download and verify the two small BioSNAP files most relevant to GraphSAGE PPI edge provenance.

The script keeps the source files by default. Pass --delete-raw-after-package to delete them
only after all integrity and comparison checks succeed and a compact results ZIP has been made.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,tarfile,urllib.request,zipfile
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path

URLS={
 'ohmnet_combined':'https://snap.stanford.edu/biodata/datasets/10013/files/PPT-Ohmnet_tissues-combined.edgelist.gz',
 'global_ppi':'https://snap.stanford.edu/biodata/datasets/10000/files/PP-Pathways_ppi.csv.gz',
}

def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def download(url,path):
 tmp=path.with_suffix(path.suffix+'.part')
 req=urllib.request.Request(url,headers={'User-Agent':'ppi-provenance-audit/1.0'})
 with urllib.request.urlopen(req,timeout=120) as r,open(tmp,'wb') as w:
  while True:
   b=r.read(1<<20)
   if not b:break
   w.write(b)
 tmp.replace(path)

def parse_pairs_gz(path):
 pairs=set(); malformed=[]
 with gzip.open(path,'rt',encoding='utf-8',errors='replace') as f:
  for n,line in enumerate(f,1):
   s=line.strip()
   if not s or s.startswith('#'):continue
   parts=s.replace(',','\t').split()
   ints=[]
   for x in parts:
    try:ints.append(int(x))
    except ValueError:continue
    if len(ints)==2:break
   if len(ints)<2:
    if n<=20:malformed.append({'line':n,'text':s[:200]})
    continue
   a,b=ints[:2];pairs.add((a,b) if a<=b else (b,a))
 return pairs,malformed

def local_ohmnet_union(tar_path):
 pairs=set(); triples=0; tissues=0
 with tarfile.open(tar_path,'r:gz') as tf:
  for m in tf.getmembers():
   if not m.isfile() or not m.name.endswith('.edgelist'):continue
   tissues+=1
   for line in tf.extractfile(m):
    if not line.strip():continue
    a,b=map(int,line.split()[:2]);pairs.add((a,b) if a<=b else (b,a));triples+=1
 return pairs,triples,tissues

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--ohmnet-tar',type=Path,required=True)
 ap.add_argument('--graphsage-zip',type=Path,required=True)
 ap.add_argument('--output-dir',type=Path,default=Path('biosnap_ppi_audit_results'))
 ap.add_argument('--delete-raw-after-package',action='store_true')
 a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 now=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
 paths={k:a.output_dir/Path(u).name for k,u in URLS.items()}
 for k,u in URLS.items():
  if not paths[k].exists():download(u,paths[k])
  with gzip.open(paths[k],'rb') as f:
   while f.read(1<<20):pass
 local_union,local_rows,local_tissues=local_ohmnet_union(a.ohmnet_tar)
 combined,combined_bad=parse_pairs_gz(paths['ohmnet_combined'])
 global_pairs,global_bad=parse_pairs_gz(paths['global_ppi'])
 with zipfile.ZipFile(a.graphsage_zip) as z:
  G=json.loads(z.read('ppi/ppi-G.json'))
 graph_pairs=set()
 # GraphSAGE row IDs are anonymous; this check is completed after joining the published row map.
 # Here record only raw edge count. Full pair comparison can be rerun with the mapping table.
 result={
  'generated_at_utc':now,'sources':{},
  'local_ohmnet_tissues':local_tissues,'local_ohmnet_edge_rows_across_tissues':local_rows,
  'local_ohmnet_unique_edge_pairs':len(local_union),
  'biosnap_combined_unique_edge_pairs':len(combined),
  'biosnap_combined_equals_local_ohmnet_union':combined==local_union,
  'biosnap_combined_missing_from_local_count':len(combined-local_union),
  'local_missing_from_biosnap_combined_count':len(local_union-combined),
  'biosnap_global_unique_edge_pairs':len(global_pairs),
  'local_ohmnet_union_is_subset_of_biosnap_global':local_union<=global_pairs,
  'local_ohmnet_edges_absent_from_biosnap_global_count':len(local_union-global_pairs),
  'graphsage_raw_link_count':len(G['links']),
  'malformed_preview':{'ohmnet_combined':combined_bad,'global_ppi':global_bad},
 }
 for k,p in paths.items():result['sources'][k]={'url':URLS[k],'local_path':str(p),'size_bytes':p.stat().st_size,'sha256':sha256(p)}
 outjson=a.output_dir/f'biosnap_ppi_provenance_audit_{now}.json';outjson.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 readme=a.output_dir/f'README_{now}.md';readme.write_text(
  '# BioSNAP PPI provenance audit\n\n'
  'This audit compares the union of the 144 supplied OhmNet tissue edgelists with the BioSNAP OhmNet combined network, then tests whether that union is a subset of the 21,557-node/342,353-edge BioSNAP global physical interactome.\n\n'
  f'- Results: `{outjson.name}`\n- Raw source hashes are embedded in the JSON.\n',encoding='utf-8')
 zipout=a.output_dir/f'biosnap_ppi_provenance_audit_{now}.zip'
 with zipfile.ZipFile(zipout,'w',compression=zipfile.ZIP_DEFLATED) as z:
  z.write(outjson,outjson.name);z.write(readme,readme.name)
 if not result['biosnap_combined_equals_local_ohmnet_union'] or not result['local_ohmnet_union_is_subset_of_biosnap_global']:
  raise SystemExit('Comparison did not pass; raw files retained for diagnosis')
 if a.delete_raw_after_package:
  for p in paths.values():p.unlink()
 print(json.dumps({'results_zip':str(zipout),'results_zip_sha256':sha256(zipout),'checks_pass':True},indent=2))
if __name__=='__main__':main()

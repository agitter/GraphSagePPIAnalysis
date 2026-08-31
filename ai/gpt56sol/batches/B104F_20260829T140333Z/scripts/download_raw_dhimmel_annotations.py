#!/usr/bin/env python3
"""Download four commit-pinned dhimmel/gene-ontology TSVs and reject HTML error pages.
Writes atomically and records URL, UTC retrieval time, size, SHA-256, and line count.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,os,tempfile,urllib.request
from pathlib import Path
COMMIT='962a5e12f8590400c2891cde93fd6a783b26e02e'
FILES=[
'GO_annotations-9606-direct-allev.tsv',
'GO_annotations-9606-direct-expev.tsv',
'GO_annotations-9606-inferred-allev.tsv',
'GO_annotations-9606-inferred-expev.tsv',
]
BASE=f'https://raw.githubusercontent.com/dhimmel/gene-ontology/{COMMIT}/annotations/taxid_9606'

def digest(path):
 h=hashlib.sha256();lines=0
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b);lines+=b.count(b'\n')
 return h.hexdigest(),lines

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--manifest',type=Path);a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
 manifest=a.manifest or a.output_dir/f'dhimmel_GO_annotations_{COMMIT[:7]}_manifest.csv';rows=[]
 for name in FILES:
  url=f'{BASE}/{name}';dest=a.output_dir/name
  fd,tmp=tempfile.mkstemp(prefix=name+'.',suffix='.part',dir=a.output_dir);os.close(fd);tmp=Path(tmp)
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'GraphSAGE-PPI-reproduction/1.0'})
   with urllib.request.urlopen(req,timeout=120) as src,tmp.open('wb') as out:
    while True:
     b=src.read(1<<20)
     if not b:break
     out.write(b)
   head=tmp.read_bytes()[:4096]
   if b'<!DOCTYPE html' in head or b'<html' in head.lower():raise RuntimeError(f'{url} returned HTML, not TSV')
   first=next((x for x in tmp.read_text(errors='replace').splitlines() if x.strip()),'')
   if not first.startswith('go_id\t'):raise RuntimeError(f'unexpected first nonempty line in {name}: {first[:120]!r}')
   if tmp.stat().st_size < 1_000_000:raise RuntimeError(f'{name} is unexpectedly small: {tmp.stat().st_size} bytes')
   os.replace(tmp,dest);sha,lines=digest(dest)
   rows.append({'filename':name,'commit':COMMIT,'source_url':url,'retrieved_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'size_bytes':dest.stat().st_size,'sha256':sha,'line_count':lines,'validation':'TSV header valid; HTML rejected; minimum size passed'})
   print(f'OK {name} {dest.stat().st_size} {sha}')
  except Exception:
   tmp.unlink(missing_ok=True);raise
 with manifest.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 print(f'Manifest: {manifest}')
if __name__=='__main__':main()

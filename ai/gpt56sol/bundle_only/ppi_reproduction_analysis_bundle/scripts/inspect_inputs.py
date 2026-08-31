#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json, os, sys, tarfile, zipfile
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd

ROOT=Path('/mnt/data')
OUT=ROOT/'ppi_repro'/'results'
OUT.mkdir(parents=True,exist_ok=True)
FILES=[
'investigation_summary_2026_08_23.md','Greene2015_Table9.xlsx','Greene2015_Table6.xlsx',
'Greene2015_sup.pdf','Greene2015.pdf','OhmNet.pdf',
'msigdb_v6.0_files_to_download_locally.zip','msigdb_v5.1_files_to_download_locally.zip',
'msigdb_v5.2_chip_files_to_download_locally.zip','msigdb_v5.2_files_to_download_locally.zip',
'graphsage_ppi.zip','dgl_ppi.zip','bio-tissue-readme.txt','bio-tissue-labels.tar.gz',
'bio-tissue-hierarchy.tar.gz','bio-tissue-networks.tar.gz','historical_go_mapping_inventory.md',
'Pasted markdown(1).md']

def sha256(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def zip_member_summary(z:zipfile.ZipFile,n:str)->dict[str,Any]:
 info=z.getinfo(n); d={'name':n,'size':info.file_size,'compressed':info.compress_size}
 low=n.lower()
 try:
  if low.endswith('.npy'):
   with z.open(n) as f:
    a=np.load(f,allow_pickle=False)
   d.update({'shape':list(a.shape),'dtype':str(a.dtype),'min':float(np.nanmin(a)) if a.size else None,'max':float(np.nanmax(a)) if a.size else None})
  elif low.endswith('.json') and info.file_size<100_000_000:
   with z.open(n) as f: obj=json.load(f)
   if isinstance(obj,dict):
    d['json_type']='dict'; d['keys']=list(obj)[:20]; d['length']=len(obj)
    if 'nodes' in obj and 'links' in obj: d.update({'nodes':len(obj['nodes']),'links':len(obj['links'])})
   elif isinstance(obj,list): d.update({'json_type':'list','length':len(obj)})
  elif low.endswith(('.txt','.tsv','.csv','.gmt','.chip','.edgelist','.lab')) and info.file_size<5_000_000:
   with z.open(n) as f: raw=f.read(4096).decode('utf-8','replace')
   d['head']=raw.splitlines()[:5]
 except Exception as e:
  d['inspection_error']=repr(e)
 return d

manifest=[]
for name in FILES:
 p=ROOT/name
 if not p.exists():
  manifest.append({'name':name,'exists':False}); continue
 rec={'name':name,'exists':True,'path':str(p),'size':p.stat().st_size,'sha256':sha256(p)}
 try:
  if zipfile.is_zipfile(p):
   with zipfile.ZipFile(p) as z:
    names=z.namelist(); rec['archive']='zip'; rec['member_count']=len(names)
    rec['members']=[zip_member_summary(z,n) for n in names]
  elif tarfile.is_tarfile(p):
   with tarfile.open(p,'r:*') as t:
    mem=[m for m in t.getmembers() if m.isfile()]
    rec['archive']='tar'; rec['member_count']=len(mem)
    rec['members']=[{'name':m.name,'size':m.size} for m in mem]
  elif p.suffix.lower()=='.xlsx':
   xl=pd.ExcelFile(p,engine='openpyxl'); rec['sheets']=xl.sheet_names; rec['sheet_summaries']={}
   for s in xl.sheet_names:
    df=pd.read_excel(p,sheet_name=s,engine='openpyxl')
    rec['sheet_summaries'][s]={'shape':list(df.shape),'columns':[str(x) for x in df.columns], 'head':df.head(8).fillna('').astype(str).to_dict('records')}
  elif p.suffix.lower() in ('.md','.txt'):
   rec['head']=p.read_text(errors='replace').splitlines()[:30]
 except Exception as e:
  rec['inspection_error']=repr(e)
 manifest.append(rec)

(OUT/'input_inventory.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
lines=['# Input inventory','', 'Generated from the exact mounted files. SHA-256 values are the provenance anchors.','']
for r in manifest:
 lines.append(f"## `{r['name']}`")
 if not r.get('exists'):
  lines.append('- **Missing**'); lines.append(''); continue
 lines += [f"- Path: `{r['path']}`",f"- Size: {r['size']:,} bytes",f"- SHA-256: `{r['sha256']}`"]
 if r.get('archive'):
  lines.append(f"- Archive: {r['archive']}, {r['member_count']:,} file members")
  lines.append('')
  lines.append('| Member | Size | Shape / type |')
  lines.append('|---|---:|---|')
  for m in r.get('members',[]):
   extra=''
   if 'shape' in m: extra=f"{m['shape']} {m.get('dtype','')}"
   elif 'nodes' in m: extra=f"JSON graph: {m['nodes']:,} nodes, {m['links']:,} links"
   elif 'json_type' in m: extra=f"JSON {m['json_type']}, n={m.get('length')}"
   lines.append(f"| `{m['name']}` | {m['size']:,} | {extra} |")
 elif 'sheets' in r:
  for s,d in r['sheet_summaries'].items(): lines.append(f"- Sheet `{s}`: {tuple(d['shape'])}; columns: {', '.join(d['columns'])}")
 lines.append('')
(OUT/'input_inventory.md').write_text('\n'.join(lines))
print(json.dumps({'status':'ok','files':len(manifest),'out':str(OUT/'input_inventory.md')}))

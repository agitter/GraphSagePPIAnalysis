#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures,datetime as dt,gzip,hashlib,json,os,re,shutil,ssl,time,urllib.error,urllib.parse,urllib.request
from pathlib import Path
import pandas as pd
ROOT=Path('/mnt/data');BASE=ROOT/'ppi_repro';DL=BASE/'downloads';OUT=BASE/'results';DL.mkdir(parents=True,exist_ok=True)
UA='Mozilla/5.0 reproducibility-research/1.0'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def fetch(item):
 url,dest,label=item;dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists() and dest.stat().st_size>100:return {'label':label,'url':url,'path':str(dest),'status':'existing','size':dest.stat().st_size,'sha256':sha(dest)}
 tmp=dest.with_suffix(dest.suffix+'.part');err=''
 for attempt in range(3):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':UA})
   with urllib.request.urlopen(req,timeout=180,context=ssl.create_default_context()) as r:
    clen=r.headers.get('Content-Length')
    if clen and int(clen)>800_000_000:raise RuntimeError(f'refusing file >800MB: {clen}')
    with tmp.open('wb') as f:shutil.copyfileobj(r,f,1<<20)
    final=r.geturl();ctype=r.headers.get('Content-Type');code=getattr(r,'status',200)
   if tmp.stat().st_size<100:raise RuntimeError('download too small')
   tmp.replace(dest)
   return {'label':label,'url':url,'final_url':final,'path':str(dest),'status':'downloaded','http_status':code,'content_type':ctype,'size':dest.stat().st_size,'sha256':sha(dest)}
  except Exception as e:
   err=repr(e)
   try:tmp.unlink()
   except:pass
   time.sleep(2**attempt)
 return {'label':label,'url':url,'path':str(dest),'status':'error','error':err}
items=[]
base='https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN'
for r in range(140,158):
 for stem in ['gene_association.goa_human','gp_information.goa_human','gene_association.goa_ref_human','gp_information.goa_ref_human']:
  fn=f'{stem}.{r}.gz';items.append((f'{base}/{fn}',DL/'ebi_goa'/fn,f'EBI {fn}'))
for r in range(158,164):
 for stem in ['goa_human.gaf','goa_human.gpi','goa_human.gpa']:
  fn=f'{stem}.{r}.gz';items.append((f'{base}/{fn}',DL/'ebi_goa'/fn,f'EBI {fn}'))
items.append((f'{base}/README',DL/'ebi_goa'/'README','EBI GOA README'))
for y in [2015,2016]:
 for m in range(1,13):
  date=f'{y:04d}-{m:02d}-01';rb=f'https://release.geneontology.org/{date}'
  items.append((f'{rb}/ontology/go-basic.obo',DL/'go_releases'/date/'go-basic.obo',f'GO basic {date}'))
  items.append((f'{rb}/annotations/gp2protein/gp2protein.geneid.gz',DL/'go_releases'/date/'gp2protein.geneid.gz',f'gp2protein GeneID {date}'))
# Bioconductor one destination per version: primary then mirror only if primary fails, handled after parallel batch.
bio=[('3.0','3.0.0'),('3.1','3.1.2'),('3.2','3.2.3'),('3.3','3.3.0'),('3.4','3.4.0')]
for branch,ver in bio:
 fn=f'org.Hs.eg.db_{ver}.tar.gz';items.append((f'https://bioconductor.org/packages/{branch}/data/annotation/src/contrib/{fn}',DL/'bioconductor'/fn,f'Bioconductor {ver} primary'))
# Human-specific UniProt files already discovered by broad downloader.
sm=OUT/'source_manifest.csv'
if sm.exists():
 try:
  d=pd.read_csv(sm).fillna('')
  for _,r in d.iterrows():
   u=str(r.get('url',''))
   if r.get('source_type')=='UNIPROT' and re.search(r'HUMAN|9606',u,re.I):
    rel=str(r.get('relative_path',''));dest=DL/rel if not rel.startswith('/') else Path(rel)
    items.append((u,dest,str(r.get('label','UniProt human mapping'))))
   if r.get('source_type')=='NCBI_WAYBACK':
    rel=str(r.get('relative_path',''));dest=DL/rel
    items.append((u,dest,str(r.get('label','NCBI Wayback'))))
 except:pass
# dedupe URL/destination
seen=set();items2=[]
for x in items:
 key=(x[0],str(x[1]))
 if key not in seen:seen.add(key);items2.append(x)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:rows=list(ex.map(fetch,items2))
# mirror fallback for Bioconductor failures
for branch,ver in bio:
 fn=f'org.Hs.eg.db_{ver}.tar.gz';dest=DL/'bioconductor'/fn
 if not dest.exists() or dest.stat().st_size<100:
  rows.append(fetch((f'https://bioconductor.statistik.tu-dortmund.de/packages/{branch}/data/annotation/src/contrib/{fn}',dest,f'Bioconductor {ver} mirror fallback')))
now=dt.datetime.now(dt.timezone.utc).isoformat()
for r in rows:r['recorded_at_utc']=now
pd.DataFrame(rows).to_csv(OUT/'priority_source_manifest.csv',index=False)
(OUT/'priority_source_manifest.json').write_text(json.dumps(rows,indent=2))
counts=pd.Series([r['status'] for r in rows]).value_counts().to_dict()
lines=['# Priority historical source acquisition','',f'Status counts: `{counts}`','', '| Label | Status | Bytes | SHA-256 | URL |','|---|---|---:|---|---|']
for r in rows:lines.append(f"| {r['label']} | {r['status']} | {r.get('size','')} | `{r.get('sha256','')}` | {r['url']} |")
(OUT/'priority_source_manifest.md').write_text('\n'.join(lines))
print(json.dumps(counts))

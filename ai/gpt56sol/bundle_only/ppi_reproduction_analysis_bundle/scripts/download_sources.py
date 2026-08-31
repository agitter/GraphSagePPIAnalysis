#!/usr/bin/env python3
from __future__ import annotations
import csv, datetime as dt, gzip, hashlib, json, os, re, shutil, ssl, subprocess, sys, time, traceback, urllib.error, urllib.parse, urllib.request
from pathlib import Path

ROOT=Path('/mnt/data'); BASE=ROOT/'ppi_repro'; DL=BASE/'downloads'; OUT=BASE/'results';DL.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
UA='Mozilla/5.0 reproducibility-research/1.0'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def fetch(url,dest,timeout=90,retries=2):
 dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists() and dest.stat().st_size>0:return {'status':'existing','size':dest.stat().st_size,'sha256':sha(dest)}
 tmp=dest.with_suffix(dest.suffix+'.part')
 err=None
 for k in range(retries+1):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':UA})
   with urllib.request.urlopen(req,timeout=timeout,context=ssl.create_default_context()) as r, tmp.open('wb') as f:
    shutil.copyfileobj(r,f,1<<20); code=getattr(r,'status',200); final=r.geturl();ctype=r.headers.get('Content-Type')
   tmp.replace(dest)
   return {'status':'downloaded','http_status':code,'final_url':final,'content_type':ctype,'size':dest.stat().st_size,'sha256':sha(dest)}
  except Exception as e:
   err=repr(e)
   try:tmp.unlink()
   except:pass
   time.sleep(2*(k+1))
 return {'status':'error','error':err}

def get_text(url,timeout=60):
 req=urllib.request.Request(url,headers={'User-Agent':UA})
 with urllib.request.urlopen(req,timeout=timeout,context=ssl.create_default_context()) as r:return r.read().decode('utf-8','replace')

entries=[]
def add(source_type,label,url,rel,notes=''):
 entries.append({'source_type':source_type,'label':label,'url':url,'relative_path':rel,'notes':notes})

# Uploaded inputs: exact local provenance.
uploads={
'graphsage_ppi.zip':'SNAP GraphSAGE PPI archive supplied by user','dgl_ppi.zip':'DGL PPI archive supplied by user',
'bio-tissue-networks.tar.gz':'OhmNet network archive supplied by user','bio-tissue-labels.tar.gz':'OhmNet label archive supplied by user',
'bio-tissue-hierarchy.tar.gz':'OhmNet hierarchy archive supplied by user','bio-tissue-readme.txt':'OhmNet README supplied by user',
'msigdb_v5.1_files_to_download_locally.zip':'MSigDB v5.1 supplied by user','msigdb_v5.2_files_to_download_locally.zip':'MSigDB v5.2 supplied by user',
'msigdb_v5.2_chip_files_to_download_locally.zip':'MSigDB v5.2 chip files supplied by user','msigdb_v6.0_files_to_download_locally.zip':'MSigDB v6.0 supplied by user',
'Greene2015_Table6.xlsx':'Greene et al. Supplementary Table 6 supplied by user','Greene2015_Table9.xlsx':'Greene et al. Supplementary Table 9 supplied by user',
'Greene2015.pdf':'Greene et al. paper supplied by user','Greene2015_sup.pdf':'Greene et al. supplement supplied by user','OhmNet.pdf':'OhmNet paper supplied by user',
'investigation_summary_2026_08_23.md':'Prior investigation summary supplied by user'}
manifest=[]
for fn,note in uploads.items():
 p=ROOT/fn
 manifest.append({'source_type':'uploaded','label':fn,'url':'','relative_path':str(p.relative_to(ROOT)) if p.exists() else fn,'notes':note,
                  'status':'present' if p.exists() else 'missing','size':p.stat().st_size if p.exists() else None,'sha256':sha(p) if p.exists() else None})

# EBI GOA historical releases in target period.
base='https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN'
for r in range(140,158):
 for stem,typ in [('gene_association.goa_human','GOA full human GAF'),('gene_association.goa_ref_human','GOA reference-proteome GAF'),
                  ('gp_information.goa_human','GOA full human GPI predecessor'),('gp_information.goa_ref_human','GOA reference GPI predecessor'),
                  ('gp_association.goa_human','GOA full human GPAD predecessor')]:
  fn=f'{stem}.{r}.gz';add('EBI_GOA',f'{typ} release {r}',f'{base}/{fn}',f'ebi_goa/{fn}')
for r in range(158,164):
 for stem,typ in [('goa_human.gaf','GOA full human GAF'),('goa_human.gpi','GOA full human GPI'),('goa_human.gpa','GOA full human GPAD')]:
  fn=f'{stem}.{r}.gz';add('EBI_GOA',f'{typ} release {r}',f'{base}/{fn}',f'ebi_goa/{fn}')
add('EBI_GOA','legacy human xrefs 99',f'{base}/human.xrefs.99.gz','ebi_goa/human.xrefs.99.gz','Last legacy human.xrefs series visible in archive')
add('EBI_GOA','archive README',f'{base}/README','ebi_goa/README')

# GO monthly release products. Use several possible annotation names; downloader records 404s explicitly.
months=[]
for y in (2015,2016):
 for m in range(1,13):months.append(f'{y:04d}-{m:02d}-01')
for date in months:
 rb=f'https://release.geneontology.org/{date}'
 add('GO_RELEASE',f'go-basic ontology {date}',f'{rb}/ontology/go-basic.obo',f'go_releases/{date}/go-basic.obo')
 add('GO_RELEASE',f'go ontology {date}',f'{rb}/ontology/go.obo',f'go_releases/{date}/go.obo')
 add('GO_RELEASE',f'GOA human annotation {date}',f'{rb}/annotations/goa_human.gaf.gz',f'go_releases/{date}/goa_human.gaf.gz')
 add('GO_RELEASE',f'gene association human {date}',f'{rb}/annotations/gene_association.goa_human.gz',f'go_releases/{date}/gene_association.goa_human.gz')
 add('GO_RELEASE',f'gp2protein GeneID {date}',f'{rb}/annotations/gp2protein/gp2protein.geneid.gz',f'go_releases/{date}/gp2protein.geneid.gz')

# Historical Bioconductor org.Hs.eg.db packages. Try primary and mirror aliases separately.
bio=[('3.0','3.0.0'),('3.1','3.1.2'),('3.2','3.2.3'),('3.3','3.3.0'),('3.4','3.4.0')]
for branch,ver in bio:
 fn=f'org.Hs.eg.db_{ver}.tar.gz'
 add('BIOCONDUCTOR',f'org.Hs.eg.db {ver}',f'https://bioconductor.org/packages/{branch}/data/annotation/src/contrib/{fn}',f'bioconductor/{fn}')
 add('BIOCONDUCTOR_MIRROR',f'org.Hs.eg.db {ver} mirror',f'https://bioconductor.statistik.tu-dortmund.de/packages/{branch}/data/annotation/src/contrib/{fn}',f'bioconductor_mirror/{fn}')

# UniProt release directory discovery; candidates are added after parsing listings.
uniprot_releases=['2015_03','2015_06','2015_09','2015_12','2016_03','2016_06','2016_09','2016_11']
for rel in uniprot_releases:
 root=f'https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-{rel}/knowledgebase/idmapping/'
 try:
  html=get_text(root)
  listing_path=DL/f'uniprot/release-{rel}-idmapping-index.html';listing_path.parent.mkdir(parents=True,exist_ok=True);listing_path.write_text(html)
  manifest.append({'source_type':'UNIPROT_INDEX','label':f'UniProt {rel} idmapping listing','url':root,'relative_path':str(listing_path.relative_to(BASE)),
                   'notes':'directory listing used for discovery','status':'downloaded','size':listing_path.stat().st_size,'sha256':sha(listing_path)})
  hrefs=re.findall(r'href=["\']([^"\']+)["\']',html,re.I)
  for href in hrefs:
   if href.startswith('?') or href.startswith('../'):continue
   if re.search(r'(HUMAN|9606|idmapping_selected|idmapping\.dat)',href,re.I) and not href.endswith('/'):
    url=urllib.parse.urljoin(root,href);add('UNIPROT',f'UniProt {rel} {href}',url,f'uniprot/release-{rel}/{Path(href).name}')
   if href.rstrip('/').endswith('by_organism'):
    sub=urllib.parse.urljoin(root,href)
    try:
     h2=get_text(sub);sp=DL/f'uniprot/release-{rel}-by_organism-index.html';sp.write_text(h2)
     manifest.append({'source_type':'UNIPROT_INDEX','label':f'UniProt {rel} by-organism listing','url':sub,'relative_path':str(sp.relative_to(BASE)),
                      'notes':'directory listing used for discovery','status':'downloaded','size':sp.stat().st_size,'sha256':sha(sp)})
     for h in re.findall(r'href=["\']([^"\']+)["\']',h2,re.I):
      if re.search(r'(HUMAN|9606)',h,re.I) and not h.endswith('/'):
       add('UNIPROT',f'UniProt {rel} human mapping {h}',urllib.parse.urljoin(sub,h),f'uniprot/release-{rel}/{Path(h).name}')
    except Exception as e:
     manifest.append({'source_type':'UNIPROT_INDEX','label':f'UniProt {rel} by-organism listing','url':sub,'relative_path':'','notes':'discovery failed','status':'error','error':repr(e)})
 except Exception as e:
  manifest.append({'source_type':'UNIPROT_INDEX','label':f'UniProt {rel} idmapping listing','url':root,'relative_path':'','notes':'discovery failed','status':'error','error':repr(e)})

# Wayback CDX discovery for historical NCBI gene2go and gene2accession.
for target,label in [('https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz','gene2go'),('ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz','gene2go_ftp'),
                     ('https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz','gene2accession')]:
 cdx='https://web.archive.org/cdx/search/cdx?'+urllib.parse.urlencode({'url':target,'from':'2015','to':'2017','output':'json','fl':'timestamp,original,statuscode,mimetype,digest,length','filter':'statuscode:200','collapse':'digest'})
 try:
  txt=get_text(cdx,120);cp=DL/f'wayback/{label}_cdx.json';cp.parent.mkdir(parents=True,exist_ok=True);cp.write_text(txt)
  manifest.append({'source_type':'WAYBACK_CDX','label':label+' captures','url':cdx,'relative_path':str(cp.relative_to(BASE)),'notes':'capture discovery',
                   'status':'downloaded','size':cp.stat().st_size,'sha256':sha(cp)})
  data=json.loads(txt)
  for row in data[1:]:
   ts,orig,status,mime,digest,length=row
   if ts[:4] not in ('2015','2016','2017'):continue
   url=f'https://web.archive.org/web/{ts}id_/{orig}'
   add('NCBI_WAYBACK',f'{label} capture {ts}',url,f'ncbi_wayback/{label}_{ts}.gz',f'CDX digest {digest}; length {length}')
 except Exception as e:
  manifest.append({'source_type':'WAYBACK_CDX','label':label+' captures','url':cdx,'relative_path':'','notes':'capture discovery failed','status':'error','error':repr(e)})

# Extract any explicit http(s) URLs from prior investigation/inventory as a tracked lead, but do not auto-download arbitrary pages.
for fn in ['investigation_summary_2026_08_23.md','historical_go_mapping_inventory.md']:
 p=ROOT/fn
 if p.exists():
  for u in sorted(set(re.findall(r'https?://[^\s)\]<>`"\']+',p.read_text(errors='replace')))):
   manifest.append({'source_type':'PRIOR_LEAD','label':fn,'url':u.rstrip('.,;'),'relative_path':'','notes':'URL extracted from prior notes; not automatically trusted or downloaded','status':'lead_only'})

# Download declared entries. Skip duplicate URLs and duplicate rel paths after first successful one.
seen=set(); path_success=set()
for e in entries:
 if e['url'] in seen:continue
 seen.add(e['url'])
 dest=DL/e['relative_path']
 if e['relative_path'] in path_success:
  continue
 r=fetch(e['url'],dest,timeout=120,retries=1)
 rec=e|r
 manifest.append(rec)
 if r['status'] in ('downloaded','existing'):path_success.add(e['relative_path'])

# Add local original path and timestamp.
now=dt.datetime.now(dt.timezone.utc).isoformat()
for r in manifest:r['recorded_at_utc']=now
(OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2))
fields=sorted(set().union(*(r.keys() for r in manifest)))
with (OUT/'source_manifest.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(manifest)
# Human-readable status report.
counts={}
for r in manifest:counts[r.get('status','unknown')]=counts.get(r.get('status','unknown'),0)+1
lines=['# Source acquisition manifest','',f'Generated: `{now}`','',f'Status counts: `{counts}`','',
       'Every successfully acquired binary has its SHA-256 in `source_manifest.csv` and `source_manifest.json`. Failed URLs are retained rather than silently discarded.','',
       '| Type | Label | Status | Bytes | SHA-256 | URL |','|---|---|---|---:|---|---|']
for r in manifest:
 lines.append(f"| {r.get('source_type','')} | {r.get('label','')} | {r.get('status','')} | {r.get('size','') or ''} | `{r.get('sha256','') or ''}` | {r.get('url','')} |")
(OUT/'source_manifest.md').write_text('\n'.join(lines))
print(json.dumps({'counts':counts,'records':len(manifest)}))

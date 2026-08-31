#!/usr/bin/env python3
from pathlib import Path
import hashlib,csv,json,subprocess,datetime,re
ROOT=Path('/mnt/data');BASE=ROOT/'ppi_repro';DL=BASE/'downloads';OUT=BASE/'results'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def infer_url(p):
 rel=str(p.relative_to(DL))
 if rel.startswith('ebi_goa/'):return 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/'+p.name
 m=re.match(r'go_releases/(20\d{2}-\d{2}-\d{2})/(.+)',rel)
 if m:
  date,name=m.groups()
  if name=='go-basic.obo':return f'https://release.geneontology.org/{date}/ontology/go-basic.obo'
  if name=='go.obo':return f'https://release.geneontology.org/{date}/ontology/go.obo'
  if name=='gp2protein.geneid.gz':return f'https://release.geneontology.org/{date}/annotations/gp2protein/gp2protein.geneid.gz'
  if name.endswith('.gaf.gz'):return f'https://release.geneontology.org/{date}/annotations/{name}'
 if rel.startswith('bioconductor/'):
  return 'historical Bioconductor package; exact attempted URLs recorded in download scripts/manifests'
 if rel.startswith('ncbi_wayback/'):return 'Wayback capture; exact URL recorded in source manifests when available'
 if rel.startswith('uniprot/'):return 'https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/ (discovered historical path)'
 return ''
rows=[];now=datetime.datetime.now(datetime.timezone.utc).isoformat()
# Supplied files
for p in sorted(ROOT.iterdir()):
 if p.is_file() and p.name not in {'ppi_reproduction_analysis_bundle.zip','test_output.png'}:
  rows.append({'kind':'uploaded_or_prior_generated','relative_path':str(p.relative_to(ROOT)),'absolute_path':str(p),'size':p.stat().st_size,'sha256':sha(p),'source_url':'','recorded_at_utc':now})
# Actual downloaded files only, no .part files or git internals.
if DL.exists():
 for p in sorted(DL.rglob('*')):
  if not p.is_file() or '.git' in p.parts or p.suffix=='.part':continue
  rows.append({'kind':'downloaded_or_cloned','relative_path':str(p.relative_to(BASE)),'absolute_path':str(p),'size':p.stat().st_size,'sha256':sha(p),'source_url':infer_url(p),'recorded_at_utc':now})
# Git repository commit provenance
repos=[]
for d in DL.rglob('.git') if DL.exists() else []:
 repo=d.parent
 try:
  commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip()
  remote=subprocess.check_output(['git','-C',str(repo),'remote','get-url','origin'],text=True).strip()
  repos.append({'path':str(repo),'remote':remote,'commit':commit})
 except:pass
with (OUT/'actual_input_file_manifest.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['kind']);w.writeheader();w.writerows(rows)
(OUT/'actual_input_file_manifest.json').write_text(json.dumps({'files':rows,'repositories':repos},indent=2))
lines=['# Actual input-file manifest','',f'Generated: `{now}`','',f'Files present and hashed: **{len(rows)}**.','',
       '| Kind | Relative path | Bytes | SHA-256 | Source URL / origin |','|---|---|---:|---|---|']
for r in rows:lines.append(f"| {r['kind']} | `{r['relative_path']}` | {r['size']} | `{r['sha256']}` | {r['source_url']} |")
lines += ['', '## Git repositories']
for r in repos:lines.append(f"- `{r['remote']}` at commit `{r['commit']}` in `{r['path']}`")
(OUT/'actual_input_file_manifest.md').write_text('\n'.join(lines))

#!/usr/bin/env python3
from pathlib import Path
import subprocess,json,re,os,datetime as dt,hashlib
BASE=Path('/mnt/data/ppi_repro');roots=[BASE/'downloads'/'entrez-gene',BASE/'downloads'/'repos'/'entrez-gene',BASE/'downloads'/'repos'/'ncbi-gene'];OUT=BASE/'downloads'/'gene2go_git_history';RES=BASE/'results';OUT.mkdir(parents=True,exist_ok=True)
def run(args,cwd=None):
 return subprocess.check_output(args,cwd=cwd,text=True,stderr=subprocess.STDOUT)
def sha(p):
 h=hashlib.sha256();h.update(p.read_bytes());return h.hexdigest()
records=[]
for repo in roots:
 if not (repo/'.git').exists():continue
 try:subprocess.run(['git','fetch','--all','--tags'],cwd=repo,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=300)
 except:pass
 try:objs=run(['git','rev-list','--objects','--all'],repo)
 except Exception as e:records.append({'repo':str(repo),'error':repr(e)});continue
 paths=sorted({line.split(' ',1)[1] for line in objs.splitlines() if ' ' in line and re.search(r'gene2go',line,re.I)})
 for path in paths:
  try:logs=run(['git','log','--all','--format=%H\t%cI','--',path],repo)
  except:continue
  seen_month=set()
  for line in logs.splitlines():
   if '\t' not in line:continue
   commit,date=line.split('\t',1);month=date[:7]
   if not ('2014-01'<=month<='2017-12') or month in seen_month:continue
   seen_month.add(month)
   try:data=subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=repo,stderr=subprocess.DEVNULL)
   except:continue
   ext='gz' if data[:2]==b'\x1f\x8b' else Path(path).suffix.lstrip('.') or 'txt'
   dest=OUT/f'{repo.name}_{month}_{commit[:10]}_{Path(path).name}'
   dest.write_bytes(data)
   records.append({'repo':str(repo),'path':path,'commit':commit,'commit_date':date,'export':str(dest),'size':len(data),'sha256':sha(dest)})
(RES/'gene2go_git_history.json').write_text(json.dumps(records,indent=2))
lines=['# Git-history gene2go recovery','']
for r in records:
 if 'export' in r:lines.append(f"- `{r['commit_date']}` `{r['commit']}` `{r['path']}` → `{r['export']}` ({r['size']:,} bytes, `{r['sha256']}`)")
 else:lines.append(f"- Error: `{r}`")
(RES/'gene2go_git_history.md').write_text('\n'.join(lines))

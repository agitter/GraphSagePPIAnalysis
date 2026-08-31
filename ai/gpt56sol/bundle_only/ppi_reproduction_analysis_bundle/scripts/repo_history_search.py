#!/usr/bin/env python3
from pathlib import Path
import subprocess,json,re,os
root=Path('/mnt/data/ppi_repro/downloads/repos');out=[]
def run(args,cwd=None,timeout=120):
 try:return subprocess.check_output(args,cwd=cwd,text=True,stderr=subprocess.STDOUT,timeout=timeout)
 except Exception as e:return 'ERROR '+repr(e)
for d in sorted(root.iterdir()) if root.exists() else []:
 if not (d/'.git').exists():continue
 # Fetch full history for the smaller directly relevant repos; fetch DGL PR 395 explicitly.
 if d.name in {'GraphSAGE','ohmnet'}:
  subprocess.run(['git','fetch','--unshallow'],cwd=d,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=180)
 if d.name=='dgl':
  subprocess.run(['git','fetch','origin','pull/395/head:pr-395'],cwd=d,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=180)
 rec={'repo':d.name,'head':run(['git','rev-parse','HEAD'],d).strip(),'searches':{}}
 for term in ['ppi.zip','ppi-G.json','class_map','MSigDB','gene2go','bio-tissue','GO:0050789','StandardScaler']:
  rec['searches'][term]=run(['git','log','--all','--oneline','-S'+term,'--'],d)[:20000]
 rec['ppi_files']=run(['bash','-lc',"git log --all --name-only --pretty=format: | sort -u | grep -Ei 'ppi|protein|gene.?ontology|msigdb|gene2go' | head -1000"],d)[:50000]
 # Grep working tree and all reachable commit snapshots for small text references.
 rec['grep_worktree']=run(['bash','-lc',"grep -RInE 'ppi\\.zip|ppi-G\\.json|MSigDB|gene2go|bio-tissue|GO:0050789|StandardScaler' . --exclude-dir=.git --exclude='*.min.js' | head -2000"],d)[:100000]
 if d.name=='dgl':
  rec['pr395_show']=run(['git','show','--stat','--oneline','pr-395'],d)[:50000]
  rec['pr395_files']=run(['git','show','--name-only','--pretty=format:','pr-395'],d)[:50000]
  rec['pr395_patch']=run(['git','show','--format=fuller','pr-395'],d)[:250000]
 out.append(rec)
Path('/mnt/data/ppi_repro/results/repository_history_search.json').write_text(json.dumps(out,indent=2))
lines=['# Repository and history search','']
for r in out:
 lines += [f"## {r['repo']}",f"Head: `{r['head']}`",'', '### Working-tree hits','```',r.get('grep_worktree',''),'```','', '### Historical searches']
 for k,v in r['searches'].items():lines += [f'#### `{k}`','```',v,'```']
 if r.get('pr395_patch'):lines += ['### DGL PR 395 patch','```diff',r['pr395_patch'],'```']
Path('/mnt/data/ppi_repro/results/repository_history_search.md').write_text('\n'.join(lines))

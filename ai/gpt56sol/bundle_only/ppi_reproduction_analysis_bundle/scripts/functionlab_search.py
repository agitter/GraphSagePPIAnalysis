#!/usr/bin/env python3
from pathlib import Path
import subprocess,json,re,os
root=Path('/mnt/data/ppi_repro/downloads/repos');res=[]
for name in ['sleipnir','humanbase','GIANT']:
 d=root/name
 if not d.exists():continue
 hits=[]
 for p in d.rglob('*'):
  if not p.is_file() or '.git' in p.parts or p.stat().st_size>20_000_000:continue
  rel=str(p.relative_to(d))
  if re.search(r'go|ontology|gold|annotation|fkt|function',rel,re.I):
   hits.append({'path':rel,'size':p.stat().st_size})
  try:t=p.read_text(errors='ignore')
  except:continue
  if re.search(r'EXP.*IDA.*IPI|gene2go|UniProt.?GOA|gold standard|ontology pruning|functional knowledge transfer|FKT',t,re.I|re.S):
   for m in list(re.finditer(r'EXP.*?IDA.*?IPI|gene2go|UniProt.?GOA|gold standard|ontology pruning|functional knowledge transfer|FKT',t,re.I|re.S))[:5]:
    hits.append({'path':rel,'context':t[max(0,m.start()-500):m.end()+1500]})
 res.append({'repo':name,'hits':hits[:3000]})
Path('/mnt/data/ppi_repro/results/functionlab_search.json').write_text(json.dumps(res,indent=2))
lines=['# FunctionLab repository search','']
for r in res:
 lines += [f"## {r['repo']}",'']
 for h in r['hits']:
  lines.append(f"### `{h['path']}`")
  if 'context' in h:lines += ['```',h['context'],'```']
  else:lines.append(f"Size: {h['size']}")
Path('/mnt/data/ppi_repro/results/functionlab_search.md').write_text('\n'.join(lines))

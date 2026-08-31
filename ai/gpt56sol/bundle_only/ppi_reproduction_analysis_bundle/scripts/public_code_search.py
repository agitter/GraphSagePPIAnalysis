#!/usr/bin/env python3
from pathlib import Path
import json,urllib.parse,urllib.request,time,re,html
OUT=Path('/mnt/data/ppi_repro/results');OUT.mkdir(exist_ok=True,parents=True)
queries=['ppi-class_map.json','ppi-G.json','snap.stanford.edu/graphsage/ppi.zip','bio-tissue-networks','bio-tissue-labels','GO:0050789','ppi-feats.npy','graph_id StandardScaler ppi','GraphSAGE PPI gene ontology 121','min_count 500 gene ontology ppi']
allres=[]
for q in queries:
 url='https://api.grep.app/v1/search?'+urllib.parse.urlencode({'q':q,'regexp':'false','case':'false'})
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0 reproducibility-research'})
  with urllib.request.urlopen(req,timeout=60) as r:data=json.load(r)
  allres.append({'query':q,'url':url,'response':data})
 except Exception as e:allres.append({'query':q,'url':url,'error':repr(e)})
 time.sleep(1)
(OUT/'public_code_search.json').write_text(json.dumps(allres,indent=2))
lines=['# Public code search','', 'Search provider: grep.app public code index.','']
for item in allres:
 lines += [f"## `{item['query']}`",'']
 if 'error' in item:lines += [f"Error: `{item['error']}`",''];continue
 hits=item['response'].get('hits',{}).get('hits',[])
 lines.append(f"Total indexed hits: **{item['response'].get('hits',{}).get('total',0)}**")
 for h in hits[:50]:
  repo=h.get('repo',{});repo=repo.get('raw',repo) if isinstance(repo,dict) else repo
  path=h.get('path',{});path=path.get('raw',path) if isinstance(path,dict) else path
  content=h.get('content',{}).get('snippet','') if isinstance(h.get('content'),dict) else ''
  content=re.sub('<[^>]+>','',content);content=html.unescape(content)
  lines += [f"- `{repo}/{path}`",'```',content[:2000],'```']
 lines.append('')
(OUT/'public_code_search.md').write_text('\n'.join(lines))

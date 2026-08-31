#!/usr/bin/env python3
from pathlib import Path
import urllib.request,urllib.parse,json,re,ssl,html,hashlib,shutil,time
BASE=Path('/mnt/data/ppi_repro');DL=BASE/'downloads'/'humanbase';OUT=BASE/'results';DL.mkdir(parents=True,exist_ok=True)
UA='Mozilla/5.0 reproducibility-research'
def get(url,binary=False,timeout=90):
 req=urllib.request.Request(url,headers={'User-Agent':UA})
 with urllib.request.urlopen(req,timeout=timeout,context=ssl.create_default_context()) as r:return r.read() if binary else r.read().decode('utf-8','replace'),r.geturl(),dict(r.headers)
def sha(p):
 h=hashlib.sha256();h.update(p.read_bytes());return h.hexdigest()
records=[];pages=['https://hb.flatironinstitute.org/download','https://hb.flatironinstitute.org/','https://humanbase.io/download','https://humanbase.io/']
urls=set()
for u in pages:
 try:
  txt,final,heads=get(u);p=DL/(re.sub(r'[^A-Za-z0-9]+','_',u).strip('_')+'.html');p.write_text(txt)
  records.append({'type':'page','url':u,'final_url':final,'status':'downloaded','path':str(p),'size':p.stat().st_size,'sha256':sha(p)})
  for x in re.findall(r'(?:href|src)=["\']([^"\']+)["\']',txt,re.I):urls.add(urllib.parse.urljoin(final,html.unescape(x)))
  for x in re.findall(r'https?://[^\s"\'<>]+',txt):urls.add(html.unescape(x.rstrip(');,]')))
 except Exception as e:records.append({'type':'page','url':u,'status':'error','error':repr(e)})
# Fetch JS and search endpoint strings.
for u in list(urls):
 if re.search(r'\.js(?:\?|$)',u):
  try:
   txt,final,heads=get(u);p=DL/('js_'+hashlib.sha1(u.encode()).hexdigest()+'.js');p.write_text(txt)
   records.append({'type':'javascript','url':u,'status':'downloaded','path':str(p),'size':p.stat().st_size,'sha256':sha(p)})
   for x in re.findall(r'https?://[^\s"\'<>]+',txt):urls.add(html.unescape(x.rstrip(');,]')))
   for x in re.findall(r'["\']([^"\']*(?:download|api|gold|network|tissue)[^"\']*)["\']',txt,re.I):
    if len(x)<500:urls.add(urllib.parse.urljoin(final,x))
  except Exception as e:records.append({'type':'javascript','url':u,'status':'error','error':repr(e)})
# Wayback captures of exact desired names and download page.
for target in ['https://hb.flatironinstitute.org/download','http://hb.flatironinstitute.org/download','https://hb.flatironinstitute.org/HumanBase-blood.dat','https://hb.flatironinstitute.org/HumanBase-kidney.dat','https://hb.flatironinstitute.org/blood_sample_tsv.gz']:
 cdx='https://web.archive.org/cdx/search/cdx?'+urllib.parse.urlencode({'url':target,'from':'2015','to':'2026','output':'json','fl':'timestamp,original,statuscode,mimetype,digest,length','filter':'statuscode:200','collapse':'digest'})
 try:
  txt,_,_=get(cdx);data=json.loads(txt);records.append({'type':'wayback_cdx','url':cdx,'status':'downloaded','captures':len(data)-1})
  for row in data[1:]:
   ts,orig,status,mime,digest,length=row;records.append({'type':'wayback_capture','timestamp':ts,'original':orig,'url':f'https://web.archive.org/web/{ts}id_/{orig}','digest':digest,'length':length,'status':'lead'})
 except Exception as e:records.append({'type':'wayback_cdx','url':cdx,'status':'error','error':repr(e)})
# grep.app exact filenames.
for q in ['HumanBase-blood.dat','HumanBase-kidney.dat','blood_sample_tsv.gz','HumanBase-blood_top.gz']:
 u='https://api.grep.app/v1/search?'+urllib.parse.urlencode({'q':q})
 try:
  txt,_,_=get(u);data=json.loads(txt);records.append({'type':'grep_app','query':q,'url':u,'status':'downloaded','result':data})
 except Exception as e:records.append({'type':'grep_app','query':q,'url':u,'status':'error','error':repr(e)})
(OUT/'humanbase_discovery.json').write_text(json.dumps(records,indent=2))
lines=['# HumanBase/GIANT source discovery','']
for r in records:
 if r['type'] in ('page','javascript'):lines.append(f"- {r['type']} `{r['url']}` — {r['status']} — {r.get('path','')} {r.get('error','')}")
 elif r['type']=='wayback_capture':lines.append(f"- Wayback `{r['timestamp']}` `{r['original']}` → {r['url']} ({r['length']} bytes)")
 elif r['type']=='grep_app':lines.append(f"- grep.app `{r['query']}`: {r.get('result',{}).get('hits',{}).get('total','error')} hits")
(OUT/'humanbase_discovery.md').write_text('\n'.join(lines))

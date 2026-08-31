#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,mimetypes,re
from pathlib import Path
NAMES=[
'dhimmel-gene-ontology-962a5e1-GO_annotations-9606-direct-allev.tsv',
'dhimmel-gene-ontology-962a5e1-GO_annotations-9606-direct-expev.tsv',
'dhimmel-gene-ontology-962a5e1-GO_annotations-9606-inferred-allev.tsv',
'dhimmel-gene-ontology-962a5e1-GO_annotations-9606-inferred-expev.tsv',
]
COMMIT='962a5e12f8590400c2891cde93fd6a783b26e02e'
def sha(p):
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();rows=[]
 for local_name in NAMES:
  p=a.input_dir/local_name;data=p.read_bytes();text=data[:5000].decode('utf-8','replace');base=local_name.split('962a5e1-',1)[-1]
  title='';m=re.search(r'<title>(.*?)</title>',p.read_text(errors='replace'),re.S)
  if m:title=' '.join(m.group(1).split())
  first_nonempty=next((x for x in p.read_text(errors='replace').splitlines() if x.strip()),'')
  is_html='<!DOCTYPE html' in text or '<html' in text.lower()
  tsv_header=first_nonempty.startswith('go_id\t')
  raw_url=f'https://raw.githubusercontent.com/dhimmel/gene-ontology/{COMMIT}/annotations/taxid_9606/{base}'
  rows.append({'artifact_name':local_name,'local_path':str(p),'size_bytes':p.stat().st_size,'sha256':sha(p),'detected_content':'GitHub HTML blob page' if is_html else 'unknown_or_data','first_nonempty_line':first_nonempty[:200],'html_title':title,'valid_annotation_tsv':int(tsv_header and not is_html),'scientific_input_status':'rejected_wrong_media_type' if is_html else 'requires_further_validation','commit_pinned_raw_url':raw_url,'notes':'Filename ends in .tsv, but bytes are an HTML document and must not be parsed as annotation data.' if is_html else ''})
 a.output.parent.mkdir(parents=True,exist_ok=True)
 with a.output.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
 print('\n'.join(f"{r['artifact_name']}: {r['scientific_input_status']} {r['sha256']}" for r in rows))
if __name__=='__main__':main()

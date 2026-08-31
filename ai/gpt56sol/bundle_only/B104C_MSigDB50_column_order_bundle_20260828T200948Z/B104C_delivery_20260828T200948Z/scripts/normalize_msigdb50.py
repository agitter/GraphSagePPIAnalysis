#!/usr/bin/env python3
from __future__ import annotations
import argparse
import gzip
import html
import re
import zipfile
from pathlib import Path

ATTR_RE = re.compile(r'([A-Z_]+)="([^"]*)"')
GO_RE = re.compile(r'GO:\d{7}')
FIELDS = ['order_all','order_C5','standard_name','systematic_name','category','subcategory','GO_ID','chip','external_details_url','member_Entrez_IDs']


def normalize(source: Path, output: Path) -> dict:
    rows=[]; c5=0; membership_sum=0; build_date=''
    with zipfile.ZipFile(source) as zf:
        xml=[n for n in zf.namelist() if n.lower().endswith('.xml')]
        if len(xml)!=1: raise RuntimeError(f'expected one XML member, got {xml}')
        with zf.open(xml[0]) as fh:
            for raw in fh:
                text=raw.decode('utf-8','replace')
                if '<MSIGDB ' in text:
                    attrs={k:html.unescape(v) for k,v in ATTR_RE.findall(text)}
                    build_date=attrs.get('BUILD_DATE','')
                if '<GENESET ' not in text: continue
                attrs={k:html.unescape(v) for k,v in ATTR_RE.findall(text)}
                category=attrs.get('CATEGORY_CODE','')
                order_c5=c5 if category=='C5' else ''
                if category=='C5': c5 += 1
                members=sorted({int(x) for x in attrs.get('MEMBERS_EZID','').split(',') if x.isdigit()})
                membership_sum += len(members)
                match=GO_RE.search(attrs.get('EXTERNAL_DETAILS_URL','')+' '+attrs.get('EXACT_SOURCE','')+' '+attrs.get('DESCRIPTION_BRIEF',''))
                rows.append([
                    len(rows),order_c5,attrs.get('STANDARD_NAME',''),attrs.get('SYSTEMATIC_NAME',''),category,
                    attrs.get('SUB_CATEGORY_CODE',''),match.group(0) if match else '',attrs.get('CHIP',''),
                    attrs.get('EXTERNAL_DETAILS_URL',''),'|'.join(map(str,members))])
    output.parent.mkdir(parents=True,exist_ok=True)
    # gzip.open does not expose mtime; use GzipFile directly for deterministic bytes.
    with output.open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
            gz.write(('\t'.join(FIELDS)+'\n').encode())
            for row in rows:
                def esc(x): return str(x).replace('\t',' ').replace('\r',' ').replace('\n',' ')
                gz.write(('\t'.join(esc(x) for x in row)+'\n').encode())
            gz.flush()
    return {'rows':len(rows),'C5_rows':c5,'unique_memberships_sum':membership_sum,'XML_build_date':build_date,'xml_member':xml[0]}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source',type=Path,default=Path('/mnt/data/msigdb_v5.0_files_to_download_locally.zip'))
    ap.add_argument('--output',type=Path,default=Path('/mnt/data/ppi_repro_corrected/batches/B104C_20260828T194921Z/retained_inputs/B104C_msigdb_v5.0_normalized_entrez_gene_sets_20260828T194921Z.tsv.gz'))
    args=ap.parse_args()
    print(normalize(args.source,args.output))
if __name__=='__main__': main()

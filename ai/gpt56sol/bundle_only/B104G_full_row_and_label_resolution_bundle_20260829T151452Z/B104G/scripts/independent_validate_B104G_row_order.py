#!/usr/bin/env python3
"""Independent validation of the full GraphSAGE row-to-GeneID ordering.

This implementation deliberately does not import or call the primary analysis
code. It uses a separately written CPython-2-style open-addressing table,
parses the raw tar/ZIP again, and checks all node sets, all edge sets, and
within-gene feature/label invariance.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, tarfile, zipfile
from pathlib import Path
from collections import defaultdict
import numpy as np


def h2(text: str) -> int:
    raw=text.encode('ascii')
    if not raw:return 0
    mask=(1<<64)-1
    value=(raw[0]<<7)&mask
    for ch in raw:
        value=((value*1000003)^ch)&mask
    value=(value^len(raw))&mask
    if value >= (1<<63):value-=1<<64
    return -2 if value==-1 else value


def legacy_keys(tokens):
    size=8; table=[None]*size; used=0
    def place(key, hv, target):
        mask=len(target)-1
        i=hv&mask; perturb=hv&((1<<64)-1)
        while target[i] is not None and target[i][0]!=key:
            i=(5*i+1+perturb)&mask; perturb>>=5
        if target[i] is None:
            target[i]=(key,hv); return 1
        return 0
    for key in tokens:
        hv=h2(key)
        added=place(key,hv,table)
        used+=added
        if added and used*3 >= len(table)*2:
            minimum=4*used if used<=50000 else 2*used
            newsize=8
            while newsize<=minimum:newsize*=2
            old=[x for x in table if x is not None]
            table=[None]*newsize
            for k,hv0 in old:place(k,hv0,table)
    return [x[0] for x in table if x is not None], len(table)


def canonical(u,v):return (u,v) if u<=v else (v,u)

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        while True:
            b=f.read(1<<20)
            if not b:break
            h.update(b)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--graphsage-zip',type=Path,required=True)
    ap.add_argument('--ohmnet-tar',type=Path,required=True)
    ap.add_argument('--core-summary',type=Path,required=True)
    ap.add_argument('--primary-map',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    core=json.load(open(a.core_summary)); tissues=core['partition']['tissues']; bounds=core['partition']['bounds']
    with zipfile.ZipFile(a.graphsage_zip) as z:
        G=json.loads(z.read('ppi/ppi-G.json'))
        F=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')),allow_pickle=False)
        ids=json.loads(z.read('ppi/ppi-id_map.json')); C=json.loads(z.read('ppi/ppi-class_map.json'))
    Y=np.zeros((len(ids),121),dtype=np.uint8)
    for node,row in ids.items():Y[int(row)]=C[node]
    expected={}
    with gzip_open_text(a.primary_map) as f:
        for r in csv.DictReader(f):expected[int(r['graphsage_row'])]=int(r['entrez_gene_id'])
    tar_data={}
    with tarfile.open(a.ohmnet_tar,'r:gz') as tf:
        for tissue in tissues:
            member=f'bio-tissue-networks/{tissue}.edgelist'
            fh=tf.extractfile(member); assert fh is not None
            toks=[]; seen=set(); edges=set(); nodes=set()
            for b in fh:
                p=b.decode().split();
                if len(p)<2:continue
                u0,v0=p[:2]
                for x in (u0,v0):
                    if x not in seen:seen.add(x);toks.append(x)
                u,v=int(u0),int(v0);nodes|={u,v};edges.add(canonical(u,v))
            tar_data[tissue]=(toks,nodes,edges)
    rowmap={}; table_sizes={}
    for i,tissue in enumerate(tissues):
        order,ts=legacy_keys(tar_data[tissue][0]); table_sizes[tissue]=ts
        assert len(order)==bounds[i+1]-bounds[i]
        for j,g in enumerate(order):rowmap[bounds[i]+j]=int(g)
    map_mismatch=[r for r in range(len(rowmap)) if rowmap[r]!=expected[r]]
    block=np.zeros(len(rowmap),dtype=np.int16)
    for i in range(24):block[bounds[i]:bounds[i+1]]=i
    gs=[set() for _ in range(24)]
    for e in G['links']:
        u,v=int(e['source']),int(e['target']);assert block[u]==block[v]
        gs[int(block[u])].add(canonical(rowmap[u],rowmap[v]))
    edge_bad=[];node_bad=[]
    for i,t in enumerate(tissues):
        _,nodes,edges=tar_data[t]
        mapped={rowmap[r] for r in range(bounds[i],bounds[i+1])}
        if mapped!=nodes:node_bad.append(t)
        if gs[i]!=edges:edge_bad.append(t)
    by_gene=defaultdict(list)
    for r,g in rowmap.items():by_gene[g].append(r)
    feature_conflicts=[];label_conflicts=[]
    for g,rows in by_gene.items():
        if any(not np.array_equal(F[rows[0]],F[r]) for r in rows[1:]):feature_conflicts.append(g)
        if any(not np.array_equal(Y[rows[0]],Y[r]) for r in rows[1:]):label_conflicts.append(g)
    result={
      'independent_implementation':True,
      'rows':len(rowmap),'unique_GeneIDs':len(by_gene),
      'primary_map_rows_compared':len(expected),'primary_map_mismatches':len(map_mismatch),
      'node_set_mismatch_tissues':node_bad,'edge_set_mismatch_tissues':edge_bad,
      'feature_vector_conflict_GeneIDs':feature_conflicts,'label_vector_conflict_GeneIDs':label_conflicts,
      'all_checks_pass':not map_mismatch and not node_bad and not edge_bad and not feature_conflicts and not label_conflicts,
      'table_sizes_by_tissue':table_sizes,
      'input_sha256':{'graphsage':sha(a.graphsage_zip),'ohmnet':sha(a.ohmnet_tar),'primary_map':sha(a.primary_map)},
    }
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))
    if not result['all_checks_pass']:raise SystemExit(1)

def gzip_open_text(path):
    import gzip
    if str(path).endswith('.gz'):return gzip.open(path,'rt',newline='',encoding='utf-8')
    return open(path,newline='',encoding='utf-8')

if __name__=='__main__':main()

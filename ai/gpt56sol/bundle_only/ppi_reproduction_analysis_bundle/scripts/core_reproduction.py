#!/usr/bin/env python3
from __future__ import annotations
import collections, csv, gzip, hashlib, json, math, os, re, shutil, sys, tarfile, traceback, zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import networkx as nx
from scipy.optimize import linear_sum_assignment
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import f1_score

ROOT=Path('/mnt/data')
BASE=ROOT/'ppi_repro'
EXT=BASE/'extracted'
OUT=BASE/'results'
for d in (EXT,OUT): d.mkdir(parents=True,exist_ok=True)
LOG=[]
def log(s):
    LOG.append(str(s)); print(s,flush=True)

def sha256_file(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def safe_extract_zip(src:Path,dst:Path):
    marker=dst/'.complete'
    if marker.exists(): return
    dst.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(src) as z:
        for info in z.infolist():
            q=(dst/info.filename).resolve()
            if not str(q).startswith(str(dst.resolve())): raise ValueError(f'unsafe zip member {info.filename}')
        z.extractall(dst)
    marker.write_text(sha256_file(src))

def safe_extract_tar(src:Path,dst:Path):
    marker=dst/'.complete'
    if marker.exists(): return
    dst.mkdir(parents=True,exist_ok=True)
    with tarfile.open(src,'r:*') as t:
        for m in t.getmembers():
            q=(dst/m.name).resolve()
            if not str(q).startswith(str(dst.resolve())): raise ValueError(f'unsafe tar member {m.name}')
        t.extractall(dst)
    marker.write_text(sha256_file(src))

def recursively_extract_zips(root:Path):
    changed=True
    while changed:
        changed=False
        for p in list(root.rglob('*.zip')):
            dst=p.with_suffix('')
            if not (dst/'.complete').exists():
                try: safe_extract_zip(p,dst); changed=True
                except zipfile.BadZipFile: pass

def find_one(root:Path,patterns:list[str])->Path:
    candidates=[]
    for pat in patterns: candidates.extend(root.rglob(pat))
    candidates=[p for p in candidates if p.is_file()]
    if not candidates: raise FileNotFoundError(f'none of {patterns} under {root}')
    candidates.sort(key=lambda p:(len(str(p)),str(p)))
    return candidates[0]

@dataclass
class GSData:
    G:nx.Graph
    feats:np.ndarray
    labels:np.ndarray
    graph_id:np.ndarray
    split:np.ndarray
    ext_ids:list[str]
    node_records:list[dict]
    paths:dict[str,str]

def load_graphsage(src:Path)->GSData:
    dst=EXT/'graphsage_ppi'; safe_extract_zip(src,dst)
    gpath=find_one(dst,['*G.json','*graph*.json'])
    idpath=find_one(dst,['*id_map.json','*idmap*.json'])
    cpath=find_one(dst,['*class_map.json','*classmap*.json'])
    fpath=find_one(dst,['*feats.npy','*feat*.npy'])
    obj=json.loads(gpath.read_text())
    id_map=json.loads(idpath.read_text())
    class_map=json.loads(cpath.read_text())
    feats=np.load(fpath,allow_pickle=False)
    n=feats.shape[0]
    row_to_ext=[None]*n
    for k,v in id_map.items(): row_to_ext[int(v)]=str(k)
    if any(x is None for x in row_to_ext): raise ValueError('id_map not dense')
    ext_to_row={x:i for i,x in enumerate(row_to_ext)}
    # node records indexed by row
    recs=[{} for _ in range(n)]
    for nd in obj.get('nodes',[]):
        eid=str(nd.get('id'))
        if eid not in ext_to_row and isinstance(nd.get('id'),float) and nd['id'].is_integer(): eid=str(int(nd['id']))
        if eid not in ext_to_row: continue
        recs[ext_to_row[eid]]=dict(nd)
    G=nx.Graph(); G.add_nodes_from(range(n))
    links=obj.get('links',obj.get('edges',[]))
    def torow(x):
        sx=str(x)
        if sx in ext_to_row:return ext_to_row[sx]
        if isinstance(x,float) and x.is_integer() and str(int(x)) in ext_to_row:return ext_to_row[str(int(x))]
        if isinstance(x,int) and 0<=x<n and row_to_ext[x]==str(x):return x
        raise KeyError(x)
    for e in links:
        if isinstance(e,dict): a,b=e['source'],e['target']
        else:a,b=e[:2]
        G.add_edge(torow(a),torow(b))
    labels=[]
    for eid in row_to_ext:
        v=class_map[eid] if eid in class_map else class_map.get(str(int(eid)) if eid.isdigit() else eid)
        if isinstance(v,list): labels.append(v)
        else: labels.append([v])
    labels=np.asarray(labels)
    # graph id
    keys=set().union(*(r.keys() for r in recs if r))
    gidkey=next((k for k in ['graph_id','graphid','graph','component'] if k in keys),None)
    if gidkey:
        graph_id=np.array([int(r[gidkey]) for r in recs],dtype=int)
    else:
        graph_id=np.full(n,-1,dtype=int)
        for i,cc in enumerate(nx.connected_components(G)):
            for u in cc:graph_id[u]=i
    split=np.array(['test' if bool(r.get('test',False)) else 'valid' if bool(r.get('val',False)) else 'train' for r in recs],dtype=object)
    # Validate graph-level split
    for gid in np.unique(graph_id):
        ss=set(split[graph_id==gid]);
        if len(ss)!=1: log(f'WARNING graph {gid} has split flags {ss}')
    return GSData(G,feats,labels,graph_id,split,row_to_ext,recs,{'graph':str(gpath),'id_map':str(idpath),'class_map':str(cpath),'feats':str(fpath)})

def load_edgelist_int(p:Path)->nx.Graph:
    G=nx.Graph()
    opener=gzip.open if p.suffix=='.gz' else open
    with opener(p,'rt',errors='replace') as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith(('#','!')):continue
            a=line.strip().split()
            if len(a)<2:continue
            try:u,v=int(float(a[0])),int(float(a[1]))
            except:continue
            if u!=v:G.add_edge(u,v)
            else:G.add_node(u)
    return G

def load_ohmnet_networks(src:Path)->dict[str,nx.Graph]:
    dst=EXT/'ohmnet_networks';safe_extract_tar(src,dst)
    out={}
    for p in sorted(dst.rglob('*.edgelist')):
        tissue=p.name[:-len('.edgelist')]
        out[tissue]=load_edgelist_int(p)
    if not out:
        for p in sorted(dst.rglob('*')):
            if p.is_file() and p.name!='.complete':
                g=load_edgelist_int(p)
                if g.number_of_edges():out[p.stem]=g
    return out

def graph_stats(G:nx.Graph)->dict[str,Any]:
    deg=sorted(d for _,d in G.degree())
    comps=sorted((len(c) for c in nx.connected_components(G)),reverse=True) if G.number_of_nodes() else []
    return {'n':G.number_of_nodes(),'m':G.number_of_edges(),'components':comps,'degree':deg,
            'deg_hash':hashlib.sha256(','.join(map(str,deg)).encode()).hexdigest()}

def degree_l1(a:list[int],b:list[int])->int:
    n=max(len(a),len(b)); aa=np.pad(np.asarray(a,dtype=np.int64),(0,n-len(a)));bb=np.pad(np.asarray(b,dtype=np.int64),(0,n-len(b)))
    return int(np.abs(aa-bb).sum())

def split_gs_graphs(gs:GSData)->dict[int,nx.Graph]:
    out={}
    for gid in sorted(map(int,np.unique(gs.graph_id))):
        nodes=np.where(gs.graph_id==gid)[0].tolist();out[gid]=gs.G.subgraph(nodes).copy()
    return out

def match_tissues(ggraphs:dict[int,nx.Graph],nets:dict[str,nx.Graph]):
    gids=sorted(ggraphs); tissues=sorted(nets)
    variants={}
    for t,G in nets.items():
        variants[(t,'full')]=G
        if G.number_of_nodes():
            cc=max(nx.connected_components(G),key=len);variants[(t,'lcc')]=G.subgraph(cc).copy()
    gstats={gid:graph_stats(ggraphs[gid]) for gid in gids}
    vstats={k:graph_stats(v) for k,v in variants.items()}
    cost=np.zeros((len(gids),len(tissues)),dtype=float); bestvar={}
    for i,gid in enumerate(gids):
        a=gstats[gid]
        for j,t in enumerate(tissues):
            opts=[]
            for typ in ('full','lcc'):
                b=vstats[(t,typ)]
                c=abs(a['n']-b['n'])*1e9+abs(a['m']-b['m'])*1e6+degree_l1(a['degree'],b['degree'])
                c+=abs(len(a['components'])-len(b['components']))*1e3
                opts.append((c,typ))
            c,typ=min(opts);cost[i,j]=c;bestvar[(gid,t)]=typ
    ri,ci=linear_sum_assignment(cost)
    matches=[]
    for i,j in zip(ri,ci):
        gid,t=gids[i],tissues[j];typ=bestvar[(gid,t)]; row=cost[i]
        order=np.argsort(row);second=float(row[order[1]]) if len(order)>1 else math.nan
        matches.append({'graph_id':gid,'tissue':t,'variant':typ,'cost':float(cost[i,j]),'second_best_cost':second,
                        'gs_stats':gstats[gid],'ohm_stats':vstats[(t,typ)]})
    return matches,variants

def wl_colors(G:nx.Graph,max_iter=50)->dict[Any,str]:
    colors={u:f'd{G.degree(u)}' for u in G}
    prev_classes=0
    for it in range(max_iter):
        sig={}
        for u in G:
            s=colors[u]+'|'+','.join(sorted(colors[v] for v in G.neighbors(u)))
            sig[u]=hashlib.sha256(s.encode()).hexdigest()
        classes=len(set(sig.values()))
        colors=sig
        if classes==prev_classes: break
        prev_classes=classes
    return colors

def initial_wl_mapping(A:nx.Graph,B:nx.Graph):
    ca,cb=wl_colors(A),wl_colors(B)
    da=collections.defaultdict(list);db=collections.defaultdict(list)
    for u,c in ca.items():da[c].append(u)
    for u,c in cb.items():db[c].append(u)
    mapping={};amb=[];bad=[]
    for c,aa in da.items():
        bb=db.get(c,[])
        if len(aa)!=len(bb):bad.append((c,len(aa),len(bb)));continue
        if len(aa)==1:mapping[aa[0]]=bb[0]
        elif aa:amb.append((c,aa,bb))
    return mapping,amb,bad,ca,cb

def extract_msigdb(src:Path,label:str):
    dst=EXT/label;safe_extract_zip(src,dst);recursively_extract_zips(dst);return dst

def parse_gmts(root:Path,version:str)->list[dict[str,Any]]:
    records=[]
    for p in sorted(root.rglob('*.gmt')):
        low=p.name.lower();
        with p.open('rt',errors='replace') as f:
            for line_no,line in enumerate(f,1):
                a=line.rstrip('\n\r').split('\t')
                if len(a)<3:continue
                name=a[0]; toks=[x.strip() for x in a[2:] if x.strip()]
                numeric=[int(x) for x in toks if re.fullmatch(r'\d+',x)]
                idtype='entrez' if ('entrez' in low or (toks and len(numeric)/len(toks)>.98)) else 'symbol'
                if idtype!='entrez':continue
                col='unknown'
                m=re.search(r'(^|[._-])(c\d(?:\.[a-z0-9_]+)?)([._-]|$)',low)
                if m:col=m.group(2)
                records.append({'version':version,'collection':col,'name':name,'description':a[1],
                                'genes':frozenset(numeric),'source':str(p),'line':line_no})
    return records

def feature_count_signatures(gs:GSData,gids:list[int],matches:list[dict],variants:dict,sets:list[dict]):
    col_sigs={j:tuple(int(gs.feats[gs.graph_id==gid,j].sum()) for gid in gids) for j in range(gs.feats.shape[1])}
    match_by_gid={x['graph_id']:x for x in matches}
    genes_by_gid={gid:set(variants[(match_by_gid[gid]['tissue'],match_by_gid[gid]['variant'])].nodes()) for gid in gids}
    sig_index=collections.defaultdict(list)
    for r in sets:
        sig=tuple(len(r['genes'] & genes_by_gid[gid]) for gid in gids)
        sig_index[sig].append(r)
    candidates={j:sig_index.get(sig,[]) for j,sig in col_sigs.items()}
    return candidates,col_sigs,genes_by_gid

def choose_feature_sets(candidates,gs,wl_maps,sets_all):
    # Evaluate count-signature candidates using uniquely WL-mapped nodes.
    chosen={}; evalrows=[]
    mapped_pairs=[]
    for gid,mp in wl_maps.items():
        for row,gene in mp.items():mapped_pairs.append((row,gene))
    for j,cands in candidates.items():
        rows=[]
        for r in cands:
            mism=sum(int(gs.feats[row,j]>0.5)!=(gene in r['genes']) for row,gene in mapped_pairs)
            rows.append((mism,r))
            evalrows.append({'column':j,'version':r['version'],'collection':r['collection'],'set_name':r['name'],
                             'source':r['source'],'mismatches_unique_wl':mism,'mapped_rows':len(mapped_pairs),
                             'set_size':len(r['genes'])})
        if rows:
            best=min(x[0] for x in rows);ties=[r for m,r in rows if m==best]
            # Prefer c1/c3 and lexicographically stable source/name; retain equivalent ties in report.
            ties.sort(key=lambda r:(0 if r['collection'].startswith(('c1','c3')) else 1,r['collection'],r['name'],r['source']))
            chosen[j]={'record':ties[0],'best_mismatch':best,'ties':ties}
    return chosen,evalrows

def resolve_ambiguous(gs:GSData,amb_by_gid,wl_maps,chosen):
    selected={j:x['record'] for j,x in chosen.items() if x['best_mismatch']==0}
    cols=sorted(selected)
    unresolved=[];assigned=[]
    for gid,classes in amb_by_gid.items():
        mp=wl_maps[gid]
        for color,rows,genes in classes:
            if len(rows)!=len(genes):
                unresolved.append({'graph_id':gid,'color':color,'rows':rows,'genes':genes,'reason':'class size mismatch'});continue
            if not cols:
                unresolved.append({'graph_id':gid,'color':color,'rows':rows,'genes':genes,'reason':'no selected feature sets'});continue
            C=np.zeros((len(rows),len(genes)),dtype=int)
            for i,row in enumerate(rows):
                obs=(gs.feats[row,cols]>0.5).astype(np.int8)
                for k,gene in enumerate(genes):
                    exp=np.array([gene in selected[j]['genes'] for j in cols],dtype=np.int8)
                    C[i,k]=int(np.abs(obs-exp).sum())
            ii,jj=linear_sum_assignment(C)
            # Assign only zero-cost; retain ambiguous identical rows/genes as unresolved but deterministic map is useful.
            for i,k in zip(ii,jj):
                if C[i,k]==0:
                    mp[rows[i]]=genes[k];assigned.append((gid,rows[i],genes[k]))
            badrows=[rows[i] for i in range(len(rows)) if rows[i] not in mp]
            if badrows:
                unresolved.append({'graph_id':gid,'color':color,'rows':badrows,'genes':[g for g in genes if g not in mp.values()],
                                   'reason':'nonzero or feature-indistinguishable assignment','min_cost':int(C[ii,jj].sum())})
    return unresolved,assigned

def verify_mapping(gs:GSData,ggraphs,matches,variants,wl_maps):
    mbg={x['graph_id']:x for x in matches}
    total=both=matched=0; rows=[]
    for gid,A in ggraphs.items():
        B=variants[(mbg[gid]['tissue'],mbg[gid]['variant'])];mp=wl_maps.get(gid,{})
        t=b=m=0
        for u,v in A.edges():
            t+=1
            if u in mp and v in mp:
                b+=1
                if B.has_edge(mp[u],mp[v]):m+=1
        total+=t;both+=b;matched+=m
        rows.append({'graph_id':gid,'tissue':mbg[gid]['tissue'],'variant':mbg[gid]['variant'],'edges':t,
                     'mapped_endpoint_edges':b,'matched_edges':m,'mapped_nodes':len(mp),'nodes':A.number_of_nodes()})
    return {'total_edges':total,'edges_with_both_endpoints_mapped':both,'matched_edges':matched,'per_graph':rows}

def load_dgl(src:Path):
    dst=EXT/'dgl_ppi';safe_extract_zip(src,dst)
    files=list(dst.rglob('*'))
    arrays={p.name:np.load(p,allow_pickle=False) for p in files if p.is_file() and p.suffix=='.npy'}
    jsons={p.name:json.loads(p.read_text()) for p in files if p.is_file() and p.suffix=='.json'}
    return dst,arrays,jsons

def pick_array(arrays,split,kind):
    aliases={'valid':['valid','val'],'train':['train'],'test':['test']}[split]
    kinds={'feat':['feat'],'label':['label'],'graph_id':['graph_id','graphid']}[kind]
    cand=[]
    for n,a in arrays.items():
        low=n.lower()
        if any(x in low for x in aliases) and any(x in low for x in kinds):cand.append((n,a))
    if not cand:return None,None
    cand.sort(key=lambda x:len(x[0]));return cand[0]

def json_edges(obj):
    if isinstance(obj,dict) and 'links' in obj:
        return [(int(e['source']),int(e['target'])) for e in obj['links']],len(obj.get('nodes',[])),bool(obj.get('directed',False))
    if isinstance(obj,dict) and 'edges' in obj:
        es=[]
        for e in obj['edges']:
            es.append((int(e.get('src',e.get('source'))),int(e.get('dst',e.get('target')))))
        return es,int(obj.get('num_nodes',obj.get('num_nodes_per_type',[0])[0] if isinstance(obj.get('num_nodes_per_type'),list) else 0)),True
    return [],0,False

def verify_dgl(gs:GSData,arrays,jsons):
    # Stable sort by graph_id within each split; original row order retained.
    scaler=StandardScaler(copy=True,with_mean=True,with_std=True).fit(gs.feats[gs.split=='train'].astype(np.float64))
    z=scaler.transform(gs.feats.astype(np.float64)).astype(np.float32)
    results={'scaler_mean':scaler.mean_.tolist(),'scaler_scale':scaler.scale_.tolist(),'splits':{},'array_files':list(arrays),'json_files':list(jsons)}
    for sp in ('train','valid','test'):
        mask=(gs.split==sp); rows=np.where(mask)[0]
        order=np.argsort(gs.graph_id[rows],kind='stable');rows=rows[order]
        fn,fa=pick_array(arrays,sp,'feat');ln,la=pick_array(arrays,sp,'label');gn,ga=pick_array(arrays,sp,'graph_id')
        rec={'n':len(rows),'expected_graph_ids':gs.graph_id[rows].tolist(),'files':{'feat':fn,'label':ln,'graph_id':gn}}
        if fa is not None:
            rec['feature_shape']=list(fa.shape);rec['feature_max_abs_diff']=float(np.max(np.abs(fa.astype(np.float64)-z[rows].astype(np.float64))))
            rec['feature_array_equal']=bool(np.array_equal(fa,z[rows]))
        if la is not None:
            rec['label_shape']=list(la.shape);rec['labels_equal']=bool(np.array_equal(la,gs.labels[rows]))
        if ga is not None:
            # allow DGL graph ids renumbered per split
            exp=gs.graph_id[rows]; uniq={g:i for i,g in enumerate(sorted(set(exp.tolist())))};ren=np.array([uniq[x] for x in exp])
            rec['graph_id_equal_global']=bool(np.array_equal(ga,exp));rec['graph_id_equal_renumbered']=bool(np.array_equal(ga,ren))
        # choose JSON matching split
        jc=[(n,o) for n,o in jsons.items() if any(x in n.lower() for x in ({'valid':['valid','val'],'train':['train'],'test':['test']}[sp]))]
        if jc:
            jn,jo=sorted(jc,key=lambda x:len(x[0]))[0];edges,nj,directed=json_edges(jo)
            local={r:i for i,r in enumerate(rows)}
            und=[(u,v) for u,v in gs.G.edges() if mask[u] and mask[v]]
            expected=collections.Counter()
            for u,v in und:
                a,b=local[u],local[v];expected[(a,b)]+=1;expected[(b,a)]+=1
            for i in range(len(rows)):expected[(i,i)]+=1
            got=collections.Counter(edges)
            rec['graph_json']=jn;rec['json_nodes']=nj;rec['json_directed']=directed;rec['json_edges']=len(edges)
            rec['expected_directed_selfloop_edges']=sum(expected.values());rec['edge_multiset_equal']=got==expected
            rec['missing_expected_edges']=sum((expected-got).values());rec['extra_edges']=sum((got-expected).values())
        results['splits'][sp]=rec
    return results

def leakage_analysis(gs:GSData,wl_maps):
    row_gene={r:g for mp in wl_maps.values() for r,g in mp.items()}
    train_rows=np.where(gs.split=='train')[0];test_rows=np.where(gs.split=='test')[0];valid_rows=np.where(gs.split=='valid')[0]
    train_map={}
    conflicts=0
    for r in train_rows:
        if r not in row_gene:continue
        g=row_gene[r];v=gs.labels[r]
        if g in train_map and not np.array_equal(train_map[g],v):conflicts+=1
        else:train_map[g]=v.copy()
    pred=np.zeros_like(gs.labels[test_rows]); known=ident=0
    for i,r in enumerate(test_rows):
        g=row_gene.get(r)
        if g in train_map:
            known+=1;pred[i]=train_map[g]
            if np.array_equal(pred[i],gs.labels[r]):ident+=1
    return {'mapped_unique_genes':len(set(row_gene.values())),'mapped_rows':len(row_gene),
            'train_nodes':len(train_rows),'valid_nodes':len(valid_rows),'test_nodes':len(test_rows),
            'test_genes_seen_in_train':known,'test_seen_fraction':known/len(test_rows) if len(test_rows) else None,
            'seen_test_labels_identical':ident,'seen_identical_fraction':ident/known if known else None,
            'lookup_micro_f1':float(f1_score(gs.labels[test_rows].ravel(),pred.ravel(),zero_division=0)),
            'train_gene_label_conflicts':conflicts}

def parse_ohmnet_labels(src:Path):
    dst=EXT/'ohmnet_labels';safe_extract_tar(src,dst)
    records=[]
    rgx=re.compile(r'^(.*?)_(GO:\d{7})\.lab$')
    for p in dst.rglob('*.lab'):
        m=rgx.match(p.name)
        if not m:continue
        tissue,term=m.groups();pos=set();n=0
        with p.open(errors='replace') as f:
            for line in f:
                a=line.strip().split()
                if not a:continue
                n+=1
                if len(a)>=2:
                    try:g=int(float(a[0]));y=float(a[-1]);
                    except:continue
                    if y>0.5:pos.add(g)
                else:
                    # no gene ID: retain only aggregate
                    try:
                        if float(a[0])>0.5:pos.add(n-1)
                    except:pass
        records.append({'tissue':tissue,'go_id':term,'positives':pos,'rows':n,'path':str(p)})
    return records

def compare_ohmnet_labels(gs:GSData,wl_maps,matches,label_records):
    by={(r['tissue'],r['go_id']):r for r in label_records}; mbg={x['graph_id']:x for x in matches}
    best=[]
    for j in range(gs.labels.shape[1]):
        bestj=None
        for (t,go),rec in by.items():
            gids=[gid for gid,x in mbg.items() if x['tissue']==t]
            if not gids:continue
            gid=gids[0];mp=wl_maps[gid]; rows=list(mp)
            if not rows:continue
            y=gs.labels[rows,j].astype(bool);p=np.array([mp[r] in rec['positives'] for r in rows])
            mism=int(np.sum(y!=p));agree=1-mism/len(rows)
            if bestj is None or agree>bestj['agreement']:
                bestj={'column':j,'tissue':t,'go_id':go,'agreement':agree,'mismatches':mism,'n':len(rows)}
        if bestj:best.append(bestj)
    return best

def workbook_go_ids(path:Path):
    out=[]
    xl=pd.ExcelFile(path,engine='openpyxl')
    for s in xl.sheet_names:
        df=pd.read_excel(path,sheet_name=s,header=None,engine='openpyxl')
        for i,row in df.iterrows():
            text=' | '.join('' if pd.isna(x) else str(x) for x in row)
            gos=re.findall(r'GO:\d{7}',text);btos=re.findall(r'BTO:\d{7}',text)
            if gos or btos:out.append({'file':path.name,'sheet':s,'excel_row':int(i)+1,'go_ids':';'.join(sorted(set(gos))),
                                      'bto_ids':';'.join(sorted(set(btos))),'text':text[:1000]})
    return out

def write_df(rows,path):
    pd.DataFrame(rows).to_csv(path,index=False)

def main():
    gs=load_graphsage(ROOT/'graphsage_ppi.zip');log(f'GraphSAGE n={len(gs.ext_ids)} m={gs.G.number_of_edges()} feats={gs.feats.shape} labels={gs.labels.shape} graph_ids={len(np.unique(gs.graph_id))}')
    ggraphs=split_gs_graphs(gs)
    nets=load_ohmnet_networks(ROOT/'bio-tissue-networks.tar.gz');log(f'OhmNet networks={len(nets)}')
    matches,variants=match_tissues(ggraphs,nets);write_df([{k:v for k,v in x.items() if not isinstance(v,dict)}|{'gs_n':x['gs_stats']['n'],'gs_m':x['gs_stats']['m'],'ohm_n':x['ohm_stats']['n'],'ohm_m':x['ohm_stats']['m']} for x in matches],OUT/'tissue_matches.csv')
    mbg={x['graph_id']:x for x in matches}
    wl_maps={};amb_by_gid={};wl_report=[]
    for gid,A in ggraphs.items():
        x=mbg[gid];B=variants[(x['tissue'],x['variant'])]
        mp,amb,bad,ca,cb=initial_wl_mapping(A,B);wl_maps[gid]=mp;amb_by_gid[gid]=amb
        wl_report.append({'graph_id':gid,'tissue':x['tissue'],'variant':x['variant'],'nodes':A.number_of_nodes(),
                          'unique_wl_mapped':len(mp),'ambiguous_classes':len(amb),'ambiguous_nodes':sum(len(a) for _,a,b in amb),
                          'bad_color_classes':len(bad)})
    write_df(wl_report,OUT/'wl_mapping_initial.csv')
    msroots={}
    for ver,fn in [('5.1','msigdb_v5.1_files_to_download_locally.zip'),('5.2','msigdb_v5.2_files_to_download_locally.zip'),('6.0','msigdb_v6.0_files_to_download_locally.zip')]:
        msroots[ver]=extract_msigdb(ROOT/fn,'msigdb_'+ver.replace('.','_'))
    sets_by_ver={v:parse_gmts(r,v) for v,r in msroots.items()};log('MSigDB sets '+str({v:len(s) for v,s in sets_by_ver.items()}))
    gids=sorted(ggraphs)
    cand52,col_sigs,genes_by_gid=feature_count_signatures(gs,gids,matches,variants,sets_by_ver['5.2'])
    chosen,evalrows=choose_feature_sets(cand52,gs,wl_maps,sets_by_ver['5.2']);write_df(evalrows,OUT/'feature_candidates_v5_2.csv')
    feature_summary=[]
    for j in range(gs.feats.shape[1]):
        ch=chosen.get(j); cands=cand52.get(j,[])
        feature_summary.append({'column':j,'ones_all_rows':int(gs.feats[:,j].sum()),'count_signature_candidates':len(cands),
                                'chosen_name':ch['record']['name'] if ch else None,'chosen_collection':ch['record']['collection'] if ch else None,
                                'chosen_source':ch['record']['source'] if ch else None,'mismatches_unique_wl':ch['best_mismatch'] if ch else None,
                                'best_ties':len(ch['ties']) if ch else None,'all_tied_names':'|'.join(r['name'] for r in ch['ties'][:50]) if ch else ''})
    write_df(feature_summary,OUT/'feature_column_mapping.csv')
    unresolved,assigned=resolve_ambiguous(gs,amb_by_gid,wl_maps,chosen);(OUT/'wl_unresolved.json').write_text(json.dumps(unresolved,indent=2))
    verify=verify_mapping(gs,ggraphs,matches,variants,wl_maps);(OUT/'topology_verification.json').write_text(json.dumps(verify,indent=2))
    write_df(verify['per_graph'],OUT/'topology_per_graph.csv')
    # Evaluate feature matches across all versions on collapsed recovered genes.
    # This is equivalent to row-wise exactness when features are gene-consistent, and
    # avoids an O(sets * columns * tissue-node-instances) loop.
    row_gene={r:g for mp in wl_maps.values() for r,g in mp.items()}; mapped_rows=sorted(row_gene)
    gene_rows=collections.defaultdict(list)
    for row,gene in row_gene.items(): gene_rows[gene].append(row)
    observed_pos={}; feature_gene_conflicts=[]
    gene_weight={g:len(rows) for g,rows in gene_rows.items()}
    for j in range(gs.feats.shape[1]):
        pos=set()
        for g,rows in gene_rows.items():
            vals=(gs.feats[rows,j]>0.5)
            if not np.all(vals==vals[0]): feature_gene_conflicts.append({'gene':g,'column':j,'rows':rows})
            if bool(vals[0]): pos.add(g)
        observed_pos[j]=pos
    (OUT/'feature_gene_conflicts.json').write_text(json.dumps(feature_gene_conflicts,indent=2))
    all_feature_matches=[]; mapped_gene_set=set(gene_rows)
    for ver,sets in sets_by_ver.items():
        projected=[]
        for r in sets:
            s=set(r['genes']) & mapped_gene_set
            projected.append((r,s))
        for j in range(gs.feats.shape[1]):
            obs=observed_pos[j];best=[];bestm=None
            for r,s in projected:
                diff=obs ^ s
                # Weight by number of tissue-node occurrences to retain the original row objective.
                m=sum(gene_weight[g] for g in diff)
                if bestm is None or m<bestm:bestm=m;best=[r]
                elif m==bestm:best.append(r)
            for r in best[:100]:
                all_feature_matches.append({'version':ver,'column':j,'mismatches':int(bestm or 0),'mapped_rows':len(mapped_rows),'mapped_genes':len(mapped_gene_set),'collection':r['collection'],'set_name':r['name'],'source':r['source'],'ties':len(best)})
    write_df(all_feature_matches,OUT/'feature_best_matches_all_versions.csv')
    # DGL
    _,arrays,jsons=load_dgl(ROOT/'dgl_ppi.zip');dgl=verify_dgl(gs,arrays,jsons);(OUT/'dgl_verification.json').write_text(json.dumps(dgl,indent=2))
    # Leakage
    leak=leakage_analysis(gs,wl_maps);(OUT/'leakage_verification.json').write_text(json.dumps(leak,indent=2))
    # OhmNet labels direct comparison
    labs=parse_ohmnet_labels(ROOT/'bio-tissue-labels.tar.gz');write_df([{k:v for k,v in r.items() if k!='positives'}|{'positive_count':len(r['positives'])} for r in labs],OUT/'ohmnet_label_inventory.csv')
    olbest=compare_ohmnet_labels(gs,wl_maps,matches,labs);write_df(olbest,OUT/'ohmnet_label_best_column_matches.csv')
    # Greene workbooks
    wbrows=workbook_go_ids(ROOT/'Greene2015_Table6.xlsx')+workbook_go_ids(ROOT/'Greene2015_Table9.xlsx');write_df(wbrows,OUT/'greene_workbook_go_bto_rows.csv')
    # Gene-label collapsed data for downstream historical GO testing.
    collapsed=[];bygene=collections.defaultdict(list)
    for r,g in row_gene.items():bygene[g].append(r)
    conflicts=[]
    for g,rows in sorted(bygene.items()):
        vals=gs.labels[rows];same=bool(np.all(vals==vals[0]))
        if not same:conflicts.append({'gene_id':g,'rows':rows})
        collapsed.append([g]+vals[0].astype(int).tolist()+[len(rows),same])
    cols=['entrez_gene_id']+[f'label_{i}' for i in range(gs.labels.shape[1])]+['occurrences','labels_identical_across_occurrences']
    pd.DataFrame(collapsed,columns=cols).to_csv(OUT/'collapsed_gene_labels.csv',index=False)
    (OUT/'gene_label_conflicts.json').write_text(json.dumps(conflicts,indent=2))
    pd.DataFrame([{'row_index':r,'graph_id':int(gs.graph_id[r]),'split':str(gs.split[r]),'entrez_gene_id':g} for r,g in sorted(row_gene.items())]).to_csv(OUT/'graphsage_row_to_entrez.csv',index=False)
    # Core summary JSON
    summary={'graphsage':{'nodes':len(gs.ext_ids),'edges':gs.G.number_of_edges(),'features':list(gs.feats.shape),'labels':list(gs.labels.shape),
                          'graphs':len(ggraphs),'splits':dict(collections.Counter(gs.split.tolist()))},
             'ohmnet_network_count':len(nets),'tissue_matches':matches,'wl_initial':wl_report,
             'mapping':{'mapped_rows':len(row_gene),'unique_entrez':len(set(row_gene.values())),'unresolved_classes':len(unresolved)},
             'topology':verify,'feature_columns':feature_summary,'dgl':dgl,'leakage':leak,
             'ohmnet_labels':{'files':len(labs),'best_exact_columns':sum(x['mismatches']==0 for x in olbest),
                              'best_99pct_columns':sum(x['agreement']>=.99 for x in olbest)},
             'greene':{'rows_with_go_or_bto':len(wbrows),'unique_go':len(set(g for r in wbrows for g in r['go_ids'].split(';') if g)),
                       'unique_bto':len(set(g for r in wbrows for g in r['bto_ids'].split(';') if g))}}
    (OUT/'core_summary.json').write_text(json.dumps(summary,indent=2))
    # Markdown report generated from actual numbers.
    s=summary; lines=['# Independent core reconstruction','',
    'This report was generated directly from the uploaded archives. No numerical claim from the prior investigation summary is used as an input.', '',
    '## Input-independent checks','',
    f"- GraphSAGE archive: **{s['graphsage']['nodes']:,} nodes**, **{s['graphsage']['edges']:,} undirected simple edges**, feature shape `{tuple(s['graphsage']['features'])}`, label shape `{tuple(s['graphsage']['labels'])}`, and **{s['graphsage']['graphs']} graph IDs**.",
    f"- Split node counts: `{s['graphsage']['splits']}`.",
    f"- OhmNet archive: **{s['ohmnet_network_count']} tissue network files**.", '',
    '## Graph-level provenance','',
    '| Graph ID | Split | OhmNet tissue | Variant | GS n | GS m | Ohm n | Ohm m | assignment cost |',
    '|---:|---|---|---|---:|---:|---:|---:|---:|']
    for x in matches:
        gid=x['graph_id'];sp=str(gs.split[np.where(gs.graph_id==gid)[0][0]])
        lines.append(f"| {gid} | {sp} | `{x['tissue']}` | {x['variant']} | {x['gs_stats']['n']} | {x['gs_stats']['m']} | {x['ohm_stats']['n']} | {x['ohm_stats']['m']} | {x['cost']:.0f} |")
    lines += ['', '## Node identity recovery', '',
              f"- Unique-WL stage mapped **{sum(x['unique_wl_mapped'] for x in wl_report):,}** tissue-node instances.",
              f"- After MSigDB-feature disambiguation, **{s['mapping']['mapped_rows']:,}/{s['graphsage']['nodes']:,}** tissue-node instances map to **{s['mapping']['unique_entrez']:,} unique Entrez Gene IDs**.",
              f"- Remaining unresolved class records: **{s['mapping']['unresolved_classes']}** (see `wl_unresolved.json`).",
              f"- Of {verify['total_edges']:,} GraphSAGE edges, {verify['edges_with_both_endpoints_mapped']:,} have both endpoints mapped and {verify['matched_edges']:,} are present in the matched OhmNet network.", '',
              '## Feature provenance', '',
              '| Col | Ones | Best v5.2 set | Collection | Mismatches on unique-WL nodes | Equivalent best ties |',
              '|---:|---:|---|---|---:|---:|']
    for x in feature_summary:
        lines.append(f"| {x['column']} | {x['ones_all_rows']} | `{x['chosen_name'] or ''}` | {x['chosen_collection'] or ''} | {x['mismatches_unique_wl'] if x['mismatches_unique_wl'] is not None else ''} | {x['best_ties'] if x['best_ties'] is not None else ''} |")
    lines += ['', '## DGL transformation', '']
    for sp,r in dgl['splits'].items():
        lines.append(f"- **{sp}**: n={r['n']}; labels equal={r.get('labels_equal')}; max feature difference={r.get('feature_max_abs_diff')}; graph IDs global/renumbered={r.get('graph_id_equal_global')}/{r.get('graph_id_equal_renumbered')}; directed+self-loop edge multiset equal={r.get('edge_multiset_equal')}; missing={r.get('missing_expected_edges')}, extra={r.get('extra_edges')}.")
    lines += ['', '## Leakage measurement', '',
              f"- Mapped unique genes: **{leak['mapped_unique_genes']:,}**.",
              f"- Test nodes whose Entrez gene occurs in training: **{leak['test_genes_seen_in_train']:,}/{leak['test_nodes']:,} ({leak['test_seen_fraction']:.6%})**.",
              f"- Of those seen test nodes, byte-identical label vectors: **{leak['seen_test_labels_identical']:,}/{leak['test_genes_seen_in_train']:,} ({leak['seen_identical_fraction']:.6%})**.",
              f"- Zero-parameter Entrez lookup micro-F1 (unseen genes predicted all-zero): **{leak['lookup_micro_f1']:.8f}**.", '',
              '## Direct OhmNet-label test', '',
              f"- Parsed **{len(labs)}** OhmNet `.lab` files. Best direct tissue-specific file comparison gives **{s['ohmnet_labels']['best_exact_columns']} exact** GraphSAGE columns and **{s['ohmnet_labels']['best_99pct_columns']} at ≥99%**. Per-column results are in `ohmnet_label_best_column_matches.csv`.", '',
              '## Outputs', '',
              '- `graphsage_row_to_entrez.csv`: independently recovered row→Entrez map.',
              '- `collapsed_gene_labels.csv`: one 121-bit label vector per recovered Entrez gene.',
              '- `feature_column_mapping.csv` and `feature_best_matches_all_versions.csv`: feature provenance.',
              '- `dgl_verification.json`: full transformation checks.',
              '- `tissue_matches.csv`, `wl_mapping_initial.csv`, `topology_per_graph.csv`: graph provenance.',
              '- `greene_workbook_go_bto_rows.csv`: all GO/BTO IDs extracted from Tables 6 and 9.', '']
    (OUT/'core_reconstruction_report.md').write_text('\n'.join(lines))
    log('DONE')

if __name__=='__main__':
    try: main()
    except Exception:
        err=traceback.format_exc();print(err,file=sys.stderr);(OUT/'core_reproduction_ERROR.txt').write_text(err);raise
    finally:
        (OUT/'core_reproduction.log').write_text('\n'.join(LOG))

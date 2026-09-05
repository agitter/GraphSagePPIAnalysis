#!/usr/bin/env python3
"""Exact enumeration of hierarchy-aware GraphSAGE/OhmNet split constructions.

Split construction uses only tissue hierarchy, graph-size thresholds, node/edge counts,
and leaf/internal status. Gene sets and labels are used only after each split has been
fixed, to evaluate overlap and the GeneID-lookup diagnostic.
"""
from __future__ import annotations
import argparse, csv, gzip, itertools, json, math, re, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

LCM = 12252240  # lcm(1..18), covering every WUP denominator in this hierarchy


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True, help='output bundle root')
    p.add_argument('--networks-dir', type=Path, required=True)
    p.add_argument('--hierarchy-edges', type=Path, required=True)
    p.add_argument('--layer-classification', type=Path, required=True)
    p.add_argument('--method', required=True, choices=[
        'conditional_mean_wup','conditional_minimax_wup',
        'branch_distinct_mean_wup','branch_distinct_coverall_mean_wup',
        'branch_distinct_node_stratified_mean_wup',
        'ancestor_blocked_mean_wup'
    ])
    p.add_argument('--universe', required=True, choices=['all144','leaf107'])
    p.add_argument('--output', type=Path, required=True)
    return p.parse_args()


def load_data(a):
    meta=pd.read_csv(a.root/'prepared/ohmnet_network_metadata.tsv',sep='\t')
    cls=pd.read_csv(a.layer_classification)
    meta=meta.merge(cls[['tissue','hierarchy_node','hierarchy_layer_type']],on='tissue',how='left',validate='one_to_one')
    names=meta.tissue.tolist(); idx={n:i for i,n in enumerate(names)}; n=len(names)
    # Gene universe and label weights.
    gdf=pd.read_csv(a.root/'prepared/all_ohmnet_genes.tsv',sep='\t')
    genes=gdf.entrez_gene_id.astype(int).tolist(); gidx={g:i for i,g in enumerate(genes)}
    weights=gdf.positive_label_count.astype(int).to_numpy()
    masks=[]; node_sets=[]
    for t in names:
        s=set()
        with open(a.networks_dir/f'{t}.edgelist') as f:
            for line in f:
                if line.strip():
                    x,y=map(int,line.split()[:2]);s.add(x);s.add(y)
        node_sets.append(s)
        m=0
        for g in s:m|=1<<gidx[g]
        masks.append(m)
    pos_total=np.array([sum(int(weights[gidx[g]]) for g in s) for s in node_sets],dtype=np.int64)
    # Hierarchy.
    parent={};nodes={'Root'}
    with open(a.hierarchy_edges) as f:
        for line in f:
            if line.strip():
                c,p=line.rstrip('\n').split('\t');parent[c]=p;nodes|={c,p}
    def chain(x):
        o=[x]
        while x!='Root':x=parent[x];o.append(x)
        return o
    anc={x:chain(x) for x in nodes};depth={x:len(anc[x])-1 for x in nodes}
    def lca(x,y):
        sy=set(anc[y]);return next(z for z in anc[x] if z in sy)
    def branch(x):return anc[x][-2]
    W=np.zeros((n,n),dtype=np.int32);D=np.zeros((n,n),dtype=np.int16);REL=np.zeros((n,n),dtype=np.uint8)
    for i,x in enumerate(meta.hierarchy_node):
        for j,y in enumerate(meta.hierarchy_node):
            z=lca(x,y);D[i,j]=depth[x]+depth[y]-2*depth[z]
            num=2*(depth[z]+1);den=(depth[x]+1)+(depth[y]+1)
            assert LCM%den==0
            W[i,j]=num*(LCM//den)
            REL[i,j]=int(x in anc[y] or y in anc[x])
    branches=np.array([branch(x) for x in meta.hierarchy_node],object)
    return dict(meta=meta,names=names,idx=idx,weights=weights,masks=masks,node_sets=node_sets,
                pos_total=pos_total,W=W,D=D,REL=REL,branches=branches)


def weight_sum(mask:int,weights:np.ndarray)->int:
    s=0
    while mask:
        bit=(mask & -mask).bit_length()-1;s+=int(weights[bit]);mask&=mask-1
    return s


def eval_overlap(T,H,d):
    tu=0
    for t in T:tu|=d['masks'][t]
    unseen_rows=0;fn=0;u=0;total_rows=0;total_pos=0
    for h in H:
        miss=d['masks'][h]&~tu
        unseen_rows+=miss.bit_count();fn+=weight_sum(miss,d['weights']);u|=d['masks'][h]
        total_rows+=len(d['node_sets'][h]);total_pos+=int(d['pos_total'][h])
    unseen_unique=(u&~tu).bit_count();unique=u.bit_count();tp=total_pos-fn
    return dict(total_rows=total_rows,unseen_rows=unseen_rows,row_overlap=(total_rows-unseen_rows)/total_rows,
                unique_genes=unique,unseen_unique=unseen_unique,unique_overlap=(unique-unseen_unique)/unique,
                total_positive_labels=total_pos,false_negative_labels=fn,positive_coverage=tp/total_pos,
                lookup_micro_f1=2*tp/(2*tp+fn))


def eval_hierarchy(T,H,d):
    W=d['W'][np.ix_(T,H)];D=d['D'][np.ix_(T,H)]
    return dict(mean_wup=float(W.mean()/LCM),max_wup=float(W.max()/LCM),
                mean_heldout_best_wup=float(W.max(axis=0).mean()/LCM),
                mean_distance=float(D.mean()),min_distance=int(D.min()),
                mean_heldout_nearest_distance=float(D.min(axis=0).mean()),
                ancestor_descendant_pairs=int(d['REL'][np.ix_(T,H)].sum()),
                train_branch_count=len(set(d['branches'][T])),heldout_branch_count=len(set(d['branches'][H])))


def stratification(E,values):
    E=np.asarray(E,dtype=int)
    order=E[np.argsort(values[E],kind='mergesort')]
    groups=np.array_split(order,5)
    lookup={int(x):b for b,g in enumerate(groups) for x in g}
    return groups,lookup


def select_training(H,E,d,method):
    H=np.asarray(H,dtype=int);hs=set(map(int,H));cand=np.array([x for x in E if int(x) not in hs],dtype=int)
    distinct='branch_distinct' in method
    if distinct and len(set(d['branches'][H]))<4:return None
    if method.startswith('branch_distinct'):
        cand=cand[~np.isin(d['branches'][cand],d['branches'][H])]
    if method.startswith('ancestor_blocked'):
        cand=cand[d['REL'][np.ix_(cand,H)].max(axis=1)==0]
    if len(cand)<20:return None
    W=d['W'][np.ix_(cand,H)];D=d['D'][np.ix_(cand,H)]
    total_w=W.sum(axis=1);max_w=W.max(axis=1);total_d=D.sum(axis=1);min_d=D.min(axis=1)
    def order(ids):
        ids=np.asarray(ids,dtype=int)
        return ids[np.lexsort((cand[ids],-total_d[ids],max_w[ids],total_w[ids]))]
    if method=='conditional_minimax_wup':
        z=int(np.partition(max_w,19)[19]);allowed=np.flatnonzero(max_w<=z);sel=order(allowed)[:20]
    elif method=='branch_distinct_coverall_mean_wup':
        mandatory=[]
        for br in sorted(set(d['branches'][cand])):
            ids=np.flatnonzero(d['branches'][cand]==br);mandatory.append(int(order(ids)[0]))
        if len(mandatory)>20:return None
        used=set(mandatory);rem=np.array([i for i in range(len(cand)) if i not in used],dtype=int)
        sel=np.array(mandatory+order(rem)[:20-len(mandatory)].tolist(),dtype=int)
    elif method=='branch_distinct_node_stratified_mean_wup':
        _,smap=stratification(E,d['meta'].node_count.to_numpy())
        out=[]
        for b in range(5):
            ids=np.array([i for i,x in enumerate(cand) if smap[int(x)]==b],dtype=int)
            if len(ids)<4:return None
            out.extend(order(ids)[:4].tolist())
        sel=np.array(out,dtype=int)
    else:
        sel=order(np.arange(len(cand)))[:20]
    T=cand[sel]
    obj=(int(total_w[sel].sum()),int(max_w[sel].max()),-int(total_d[sel].sum()))
    if method=='conditional_minimax_wup':obj=(int(max_w[sel].max()),int(total_w[sel].sum()),-int(total_d[sel].sum()))
    return T,obj


def partition_roles(T,H,d):
    H=tuple(map(int,H));rows=[]
    # Three unordered 2+2 partitions. Orientation is intentionally left symmetric.
    seen=set()
    for A in itertools.combinations(H,2):
        B=tuple(h for h in H if h not in A)
        key=tuple(sorted((tuple(sorted(A)),tuple(sorted(B)))))
        if key in seen:continue
        seen.add(key)
        ha=eval_hierarchy(np.asarray(T),np.asarray(A),d);hb=eval_hierarchy(np.asarray(T),np.asarray(B),d)
        cross=d['W'][np.ix_(A,B)]
        # Primary: minimize the worse train-to-role mean similarity. Secondary: balance the two roles.
        # Tertiary: minimize validation-test similarity (negative distance equivalent is not used here).
        score=(max(ha['mean_wup'],hb['mean_wup']),abs(ha['mean_wup']-hb['mean_wup']),float(cross.mean()/LCM),key)
        rows.append((score,A,B,ha,hb,eval_overlap(T,A,d),eval_overlap(T,B,d)))
    rows.sort(key=lambda x:x[0])
    best=rows[0]
    all_f1=[x[5]['lookup_micro_f1'] for x in rows]+[x[6]['lookup_micro_f1'] for x in rows]
    all_ro=[x[5]['row_overlap'] for x in rows]+[x[6]['row_overlap'] for x in rows]
    return dict(role_pair_A='|'.join(d['names'][i] for i in best[1]),role_pair_B='|'.join(d['names'][i] for i in best[2]),
                role_A_mean_wup=best[3]['mean_wup'],role_B_mean_wup=best[4]['mean_wup'],
                role_A_lookup_f1=best[5]['lookup_micro_f1'],role_B_lookup_f1=best[6]['lookup_micro_f1'],
                role_A_row_overlap=best[5]['row_overlap'],role_B_row_overlap=best[6]['row_overlap'],
                any_role_lookup_f1_min=min(all_f1),any_role_lookup_f1_max=max(all_f1),
                any_role_row_overlap_min=min(all_ro),any_role_row_overlap_max=max(all_ro))


def main():
    a=parse_args();d=load_data(a);m=d['meta']
    leaf=(m.hierarchy_layer_type=='leaf_layer').to_numpy()
    E=np.flatnonzero(m.meets_raw_15000.to_numpy()==1);L=np.flatnonzero(m.meets_raw_35000.to_numpy()==1)
    if a.universe=='leaf107':E=E[leaf[E]];L=L[leaf[L]]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    rows=[];best=None;best_count=0;start=time.time();feasible=0
    for hi,H in enumerate(itertools.combinations(map(int,L),4),1):
        if hi%10000==0: print(f'progress {hi}/{math.comb(len(L),4)} rows={len(rows)} elapsed={time.time()-start:.1f}s', flush=True)
        x=select_training(H,E,d,a.method)
        if x is None:continue
        feasible+=1;T,obj=x
        ho=eval_hierarchy(T,np.asarray(H),d);ov=eval_overlap(T,H,d)
        r=dict(universe=a.universe,method=a.method,heldout_indices='|'.join(map(str,H)),heldout_tissues='|'.join(d['names'][i] for i in H),
               training_indices='|'.join(map(str,map(int,T))),training_tissues='|'.join(d['names'][i] for i in T),
               objective_primary=obj[0],objective_secondary=obj[1],objective_tertiary=obj[2],
               training_node_total=int(m.iloc[T].node_count.sum()),training_edge_total=int(m.iloc[T].edge_count.sum()),
               heldout_node_total=int(m.iloc[list(H)].node_count.sum()),heldout_edge_total=int(m.iloc[list(H)].edge_count.sum()),
               training_internal_layers=int((m.iloc[T].hierarchy_layer_type=='internal_hierarchy_layer').sum()),
               heldout_internal_layers=int((m.iloc[list(H)].hierarchy_layer_type=='internal_hierarchy_layer').sum()),
               **ho,**ov)
        rows.append(r)
        if best is None or obj<best:best=obj;best_count=1
        elif obj==best:best_count+=1
    print(f'building dataframe rows={len(rows)} elapsed={time.time()-start:.1f}s',flush=True)
    df=pd.DataFrame(rows)
    print(f'writing {a.output} elapsed={time.time()-start:.1f}s',flush=True)
    df.to_csv(a.output,index=False,compression={'method':'gzip','compresslevel':9,'mtime':0})
    summary={'universe':a.universe,'method':a.method,'eligible_training_networks':len(E),'eligible_heldout_networks':len(L),
             'heldout_combinations_total':math.comb(len(L),4),'heldout_combinations_feasible':feasible,
             'best_objective':best,'best_objective_tie_count':best_count,'rows':len(df),'seconds':time.time()-start,
             'output':str(a.output),'best_rows':[]}
    best_df=df.loc[(df.objective_primary==best[0])&(df.objective_secondary==best[1])&(df.objective_tertiary==best[2])]
    for _,br in best_df.iterrows():
        H=tuple(map(int,str(br.heldout_indices).split('|')));T=np.array(list(map(int,str(br.training_indices).split('|'))),dtype=int)
        z=br.to_dict();z.update(partition_roles(T,H,d));summary['best_rows'].append(z)
    with open(a.output.with_suffix('').with_suffix('.summary.json'),'w') as f:json.dump(summary,f,indent=2,sort_keys=True)
    print(json.dumps({k:v for k,v in summary.items() if k!='best_rows'},indent=2))

if __name__=='__main__':main()

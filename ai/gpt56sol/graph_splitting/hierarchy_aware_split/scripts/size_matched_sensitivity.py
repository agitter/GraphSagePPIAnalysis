#!/usr/bin/env python3
import struct
from pathlib import Path
import numpy as np, pandas as pd
ROOT=Path('/mnt/data/GraphSAGE_PPI_hierarchy_aware_split_20260904');RES=ROOT/'results';WORK=ROOT/'work'
DT=np.dtype([('h0','<u2'),('h1','<u2'),('h2','<u2'),('h3','<u2'),('a0','<u2'),('a1','<u2'),('b0','<u2'),('b1','<u2'),('mean_wup','<f4'),('max_wup','<f4'),('mean_best_wup','<f4'),('mean_dist','<f4'),('min_dist','<f4'),('mean_nearest_dist','<f4'),('row_overlap','<f4'),('unique_overlap','<f4'),('pos_coverage','<f4'),('f1','<f4'),('roleA_mean_wup','<f4'),('roleB_mean_wup','<f4'),('role_worst_mean_wup','<f4'),('role_absdiff_wup','<f4'),('role_mean_f1','<f4'),('role_min_f1','<f4'),('role_max_f1','<f4'),('role_mean_row_overlap','<f4'),('role_min_row_overlap','<f4'),('role_max_row_overlap','<f4'),('train_nodes','<u4'),('train_edges','<u4'),('held_nodes','<u4'),('held_edges','<u4'),('train_branches','u1'),('held_branches','u1'),('train_internal','u1'),('held_internal','u1'),('related_pairs','u1')])
p=WORK/'leaf107__matched_stratified.bin'
with open(p,'rb') as f:f.read(8);N=struct.unpack('<Q',f.read(8))[0];f.read(20)
x=np.memmap(p,mode='r',dtype=DT,offset=32,shape=(N,))
k=((x['h0'].astype(np.int64)*144+x['h1'])*144+x['h2'])*144+x['h3']
order=np.argsort(k);ks=k[order];uk,start,cnt=np.unique(ks,return_index=True,return_counts=True);lookup={int(v):i for i,v in enumerate(uk)}
exact=pd.read_csv(RES/'leaf107__branch_distinct_node_stratified_minimax_wup.tsv.gz',sep='\t')
sn=float(np.std(x['train_nodes']));se=float(np.std(x['train_edges']))
rows=[]
for _,r in exact.iterrows():
 h=list(map(int,r.heldout_indices.split('|')));key=((h[0]*144+h[1])*144+h[2])*144+h[3];j=lookup[key];ix=order[start[j]:start[j]+cnt[j]]
 dist=((x['train_nodes'][ix].astype(float)-r.training_node_total)/sn)**2+((x['train_edges'][ix].astype(float)-r.training_edge_total)/se)**2
 ii=np.argsort(dist)
 for K in [25,50,100,200]:
  z=ix[ii[:min(K,len(ii))]]
  row={'heldout_tissues':r.heldout_tissues,'heldout_indices':r.heldout_indices,'K':min(K,len(z)),'available':len(ix),'max_standardized_squared_distance':float(np.max(dist[ii[:min(K,len(ii))]]))}
  for em,nm in [('mean_wup','mean_wup'),('row_overlap','row_overlap'),('unique_overlap','unique_overlap'),('lookup_micro_f1','f1'),('training_node_total','train_nodes'),('training_edge_total','train_edges')]:
   rv=float(np.mean(x[nm][z]));ov=float(r[em]);row[f'optimized_{em}']=ov;row[f'matched_random_mean_{em}']=rv;row[f'difference_{em}']=ov-rv
  rows.append(row)
df=pd.DataFrame(rows);df.to_csv(RES/'leaf107_size_matched_conditional_effects.csv',index=False)
s=[]
for K,g in df.groupby('K'):
 for m in ['mean_wup','row_overlap','unique_overlap','lookup_micro_f1','training_node_total','training_edge_total']:
  a=g[f'difference_{m}'].to_numpy();s.append({'K':K,'metric':m,'n_heldout_sets':len(a),'mean_difference':a.mean(),'median_difference':np.median(a),'q05':np.quantile(a,.05),'q95':np.quantile(a,.95),'fraction_optimized_lower':np.mean(a<0)})
pd.DataFrame(s).to_csv(RES/'leaf107_size_matched_conditional_summary.csv',index=False)
pd.DataFrame([{'min_samples_per_H':int(cnt.min()),'median_samples_per_H':float(np.median(cnt)),'max_samples_per_H':int(cnt.max()),'node_sd':sn,'edge_sd':se}]).to_csv(ROOT/'tests/leaf107_conditional_sample_coverage.csv',index=False)
print(pd.DataFrame(s).to_string(index=False))

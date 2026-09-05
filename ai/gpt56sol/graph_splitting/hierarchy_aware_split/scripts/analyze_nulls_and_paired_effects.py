#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, json, math, struct
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path('/mnt/data/GraphSAGE_PPI_hierarchy_aware_split_20260904'); RES=ROOT/'results'; WORK=ROOT/'work'
DT=np.dtype([
('h0','<u2'),('h1','<u2'),('h2','<u2'),('h3','<u2'),('a0','<u2'),('a1','<u2'),('b0','<u2'),('b1','<u2'),
('mean_wup','<f4'),('max_wup','<f4'),('mean_best_wup','<f4'),('mean_dist','<f4'),('min_dist','<f4'),('mean_nearest_dist','<f4'),
('row_overlap','<f4'),('unique_overlap','<f4'),('pos_coverage','<f4'),('f1','<f4'),
('roleA_mean_wup','<f4'),('roleB_mean_wup','<f4'),('role_worst_mean_wup','<f4'),('role_absdiff_wup','<f4'),
('role_mean_f1','<f4'),('role_min_f1','<f4'),('role_max_f1','<f4'),('role_mean_row_overlap','<f4'),('role_min_row_overlap','<f4'),('role_max_row_overlap','<f4'),
('train_nodes','<u4'),('train_edges','<u4'),('held_nodes','<u4'),('held_edges','<u4'),
('train_branches','u1'),('held_branches','u1'),('train_internal','u1'),('held_internal','u1'),('related_pairs','u1')])
assert DT.itemsize==117

def load_bin(p):
    with open(p,'rb') as f:
        magic=f.read(8); N=struct.unpack('<Q',f.read(8))[0];seed=struct.unpack('<Q',f.read(8))[0];rs=struct.unpack('<I',f.read(4))[0];nh=struct.unpack('<I',f.read(4))[0]
    assert magic==b'GSHNULL1' and rs==DT.itemsize
    x=np.memmap(p,mode='r',dtype=DT,offset=32,shape=(N,))
    return x,dict(N=N,seed=seed,record_size=rs,Hpool=nh)

qs=[0,0.0001,0.001,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,0.999,0.9999,1]
metrics=['mean_wup','max_wup','mean_best_wup','mean_dist','min_dist','mean_nearest_dist','row_overlap','unique_overlap','pos_coverage','f1','role_worst_mean_wup','role_absdiff_wup','role_mean_f1','role_min_f1','role_max_f1','role_mean_row_overlap','role_min_row_overlap','role_max_row_overlap','train_nodes','train_edges','held_nodes','held_edges','train_branches','held_branches','train_internal','held_internal','related_pairs']
summary=[];metadata={};conv=[];uniformity=[]
nulls={}
for u in ['all144','leaf107']:
  for mode in ['uniform','matched_stratified']:
    x,md=load_bin(WORK/f'{u}__{mode}.bin'); nulls[(u,mode)]=x; metadata[f'{u}__{mode}']=md
    for m in metrics:
      a=np.asarray(x[m],dtype=float)
      rec={'universe':u,'null':mode,'metric':m,'N':len(a),'mean':float(a.mean()),'sd':float(a.std(ddof=1)),'se_mean':float(a.std(ddof=1)/math.sqrt(len(a)))}
      for q,v in zip(qs,np.quantile(a,qs)):rec[f'q_{q:g}']=float(v)
      summary.append(rec)
    # block convergence for load-bearing means
    B=20;n=len(x)//B
    for b in range(B):
      sl=slice(b*n,(b+1)*n)
      for m in ['mean_wup','row_overlap','f1','role_min_f1','train_nodes','train_edges']:
        conv.append({'universe':u,'null':mode,'block':b,'start':b*n,'end':(b+1)*n,'metric':m,'mean':float(np.asarray(x[m][sl],float).mean())})
    # Heldout tissue marginal z scores; expected uniform because H is sampled uniformly over enumerated H sets.
    ids=np.concatenate([np.asarray(x[f'h{i}']) for i in range(4)])
    vals,counts=np.unique(ids,return_counts=True)
    expected=len(ids)/len(vals); z=(counts-expected)/math.sqrt(expected*(1-1/len(vals)))
    uniformity.append({'universe':u,'null':mode,'categories_observed':len(vals),'expected_count':expected,'min_count':int(counts.min()),'max_count':int(counts.max()),'max_abs_z':float(np.max(np.abs(z))),'chi_square':float(np.sum((counts-expected)**2/expected))})

pd.DataFrame(summary).to_csv(RES/'null_distribution_summary.csv',index=False)
pd.DataFrame(conv).to_csv(RES/'null_block_convergence.csv',index=False)
pd.DataFrame(uniformity).to_csv(ROOT/'tests/null_heldout_marginal_uniformity.csv',index=False)
(RES/'null_sampling_metadata.json').write_text(json.dumps(metadata,indent=2,sort_keys=True))

# Load reference evaluator and calculate actual split + exact optimum role partitions.
spec=importlib.util.spec_from_file_location('enum',ROOT/'scripts/enumerate_hierarchy_aware_splits.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class A:pass
a=A();a.root=ROOT;a.networks_dir=Path('/mnt/data/_ohmnet/bio-tissue-networks');a.hierarchy_edges=Path('/mnt/data/_ohmnet/bio-tissue-hierarchy/tissue.edges');a.layer_classification=ROOT/'prepared/ohmnet_144_layer_classification.csv';d=mod.load_data(a);meta=d['meta'];idx=d['idx']
actual_T=[idx[t] for t in meta.loc[meta.graphsage_split=='train','tissue']]
actual_V=[idx[t] for t in meta.loc[meta.graphsage_split=='validation','tissue']]
actual_S=[idx[t] for t in meta.loc[meta.graphsage_split=='test','tissue']]
actual_H=actual_V+actual_S
actual={}
for label,H in [('actual_combined_heldout',actual_H),('actual_validation',actual_V),('actual_test',actual_S)]:
    z={'label':label,'training_tissues':'|'.join(meta.loc[actual_T,'tissue']),'heldout_tissues':'|'.join(meta.loc[H,'tissue'])}
    z.update(mod.eval_hierarchy(np.array(actual_T,dtype=int),np.array(H,dtype=int),d));z.update(mod.eval_overlap(np.array(actual_T,dtype=int),H,d));actual[label]=z
pd.DataFrame(actual.values()).to_csv(RES/'actual_split_metrics_recomputed.csv',index=False)

complete_methods={
'all144':['conditional_mean_wup','conditional_minimax_wup','ancestor_blocked_mean_wup','branch_distinct_mean_wup','branch_distinct_minimax_wup','branch_distinct_coverall_mean_wup','branch_distinct_node_stratified_mean_wup','branch_distinct_node_stratified_minimax_wup'],
'leaf107':['conditional_mean_wup','conditional_minimax_wup','ancestor_blocked_mean_wup','branch_distinct_mean_wup','branch_distinct_minimax_wup','branch_distinct_coverall_mean_wup','branch_distinct_node_stratified_mean_wup','branch_distinct_node_stratified_minimax_wup','branch_distinct_coverall_node_stratified_mean_wup']}
opt=[]
for u,methods in complete_methods.items():
  for method in methods:
    p=RES/f'{u}__{method}.tsv.gz';df=pd.read_csv(p,sep='\t')
    keys=['objective_primary','objective_secondary','objective_tertiary'];best=min(map(tuple,df[keys].to_numpy()));b=df[(df[keys[0]]==best[0])&(df[keys[1]]==best[1])&(df[keys[2]]==best[2])]
    for _,r in b.iterrows():
      H=tuple(map(int,r.heldout_indices.split('|')));T=np.array(list(map(int,r.training_indices.split('|'))));z=r.to_dict();z.update(mod.partition_roles(T,H,d));z['best_tie_count']=len(b);opt.append(z)
opt_df=pd.DataFrame(opt);opt_df.to_csv(RES/'exact_global_optima_and_balanced_roles.csv',index=False)

# Empirical percentiles against relevant nulls.
def midrank(a,x):
    a=np.asarray(a);return float((np.count_nonzero(a<x)+0.5*np.count_nonzero(a==x))/len(a))
comparisons=[]
for _,r in opt_df.iterrows():
  u=r.universe
  modes=['uniform']
  if str(r.method).startswith('branch_distinct_node_stratified'):modes.append('matched_stratified')
  if str(r.method).startswith('branch_distinct_coverall_node_stratified'):modes.append('matched_stratified')
  for mode in modes:
    x=nulls[(u,mode)]
    metric_map={'mean_wup':'mean_wup','max_wup':'max_wup','row_overlap':'row_overlap','unique_overlap':'unique_overlap','lookup_micro_f1':'f1','role_A_lookup_f1':'role_min_f1','role_B_lookup_f1':'role_max_f1','training_node_total':'train_nodes','training_edge_total':'train_edges'}
    for src,dst in metric_map.items():
      comparisons.append({'candidate_universe':u,'candidate_method':r.method,'heldout_tissues':r.heldout_tissues,'null':mode,'candidate_metric':src,'null_metric':dst,'value':float(r[src]),'midrank_percentile':midrank(x[dst],float(r[src]))})
# actual only against all144 uniform; external because invalid under current thresholds
ar=actual['actual_combined_heldout']
for src,dst in [('mean_wup','mean_wup'),('max_wup','max_wup'),('row_overlap','row_overlap'),('unique_overlap','unique_overlap'),('lookup_micro_f1','f1')]:
    comparisons.append({'candidate_universe':'actual_external','candidate_method':'released_graphsage_split','heldout_tissues':ar['heldout_tissues'],'null':'all144_uniform','candidate_metric':src,'null_metric':dst,'value':float(ar[src]),'midrank_percentile':midrank(nulls[('all144','uniform')][dst],float(ar[src]))})
pd.DataFrame(comparisons).to_csv(RES/'candidate_vs_null_percentiles.csv',index=False)

# Paired conditional-on-H comparison for branch-distinct + node-stratified minimax.
paired=[];paired_summary=[]
for u in ['all144','leaf107']:
  exact=pd.read_csv(RES/f'{u}__branch_distinct_node_stratified_minimax_wup.tsv.gz',sep='\t')
  x=nulls[(u,'matched_stratified')]
  # integer key from sorted H IDs, base 144.
  k=((np.asarray(x['h0'],np.int64)*144+np.asarray(x['h1'],np.int64))*144+np.asarray(x['h2'],np.int64))*144+np.asarray(x['h3'],np.int64)
  keys=np.unique(k)
  # DataFrame group means via sorting/bincount
  order=np.argsort(k);ks=k[order];uk,start,cnt=np.unique(ks,return_index=True,return_counts=True)
  means={}
  for m in ['mean_wup','max_wup','row_overlap','unique_overlap','f1','train_nodes','train_edges']:
      vals=np.asarray(x[m],float)[order];means[m]=np.add.reduceat(vals,start)/cnt
  lookup={int(key):i for i,key in enumerate(uk)}
  for _,r in exact.iterrows():
      h=list(map(int,r.heldout_indices.split('|')));key=((h[0]*144+h[1])*144+h[2])*144+h[3];j=lookup[key]
      z={'universe':u,'heldout_tissues':r.heldout_tissues,'heldout_indices':r.heldout_indices,'random_samples_conditional':int(cnt[j])}
      for em,nm in [('mean_wup','mean_wup'),('max_wup','max_wup'),('row_overlap','row_overlap'),('unique_overlap','unique_overlap'),('lookup_micro_f1','f1'),('training_node_total','train_nodes'),('training_edge_total','train_edges')]:
          z[f'optimized_{em}']=float(r[em]);z[f'random_mean_{em}']=float(means[nm][j]);z[f'difference_{em}']=float(r[em]-means[nm][j])
      paired.append(z)
  pdf=pd.DataFrame([z for z in paired if z['universe']==u])
  for metric in ['mean_wup','max_wup','row_overlap','unique_overlap','lookup_micro_f1','training_node_total','training_edge_total']:
      a=pdf[f'difference_{metric}'].to_numpy()
      paired_summary.append({'universe':u,'metric':metric,'n_heldout_sets':len(a),'mean_difference':float(a.mean()),'median_difference':float(np.median(a)),'q05':float(np.quantile(a,.05)),'q95':float(np.quantile(a,.95)),'fraction_optimized_lower':float(np.mean(a<0)),'fraction_optimized_higher':float(np.mean(a>0))})
  rho,p=spearmanr(pdf['difference_mean_wup'],pdf['difference_lookup_micro_f1'])
  paired_summary.append({'universe':u,'metric':'spearman_delta_wup_vs_delta_f1','n_heldout_sets':len(pdf),'mean_difference':float(rho),'median_difference':float(p),'q05':np.nan,'q95':np.nan,'fraction_optimized_lower':np.nan,'fraction_optimized_higher':np.nan})
pd.DataFrame(paired).to_csv(RES/'paired_conditional_on_heldout_effects.csv',index=False)
pd.DataFrame(paired_summary).to_csv(RES/'paired_conditional_on_heldout_summary.csv',index=False)

print('NULL HEADLINES')
sdf=pd.DataFrame(summary)
for u in ['all144','leaf107']:
  for mode in ['uniform','matched_stratified']:
    print('\n',u,mode)
    for m in ['mean_wup','row_overlap','f1','role_min_f1','train_nodes','train_edges']:
      r=sdf[(sdf.universe==u)&(sdf['null']==mode)&(sdf.metric==m)].iloc[0]
      print(m,'mean',r['mean'],'q05',r['q_0.05'],'median',r['q_0.5'],'q95',r['q_0.95'])
print('\nPAIRED')
print(pd.DataFrame(paired_summary).to_string(index=False))
print('\nACTUAL',json.dumps(actual,indent=2))

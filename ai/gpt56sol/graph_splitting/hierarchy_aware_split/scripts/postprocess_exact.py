#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, itertools, json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path('/mnt/data/GraphSAGE_PPI_hierarchy_aware_split_20260904')
RES=ROOT/'results'
METHODS_BY_UNIVERSE={
'all144':['conditional_mean_wup','conditional_minimax_wup','ancestor_blocked_mean_wup','branch_distinct_mean_wup','branch_distinct_minimax_wup','branch_distinct_coverall_mean_wup','branch_distinct_node_stratified_mean_wup','branch_distinct_node_stratified_minimax_wup'],
'leaf107':['conditional_mean_wup','conditional_minimax_wup','ancestor_blocked_mean_wup','branch_distinct_mean_wup','branch_distinct_minimax_wup','branch_distinct_coverall_mean_wup','branch_distinct_node_stratified_mean_wup','branch_distinct_node_stratified_minimax_wup','branch_distinct_coverall_node_stratified_mean_wup']}
UNIVERSES=['all144','leaf107']

# Load code module for independent reference evaluation / role partitioning.
import importlib.util
spec=importlib.util.spec_from_file_location('enum',ROOT/'scripts/enumerate_hierarchy_aware_splits.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class A: pass
a=A();a.root=ROOT;a.networks_dir=Path('/mnt/data/_ohmnet/bio-tissue-networks');a.hierarchy_edges=Path('/mnt/data/_ohmnet/bio-tissue-hierarchy/tissue.edges');a.layer_classification=ROOT/'prepared/ohmnet_144_layer_classification.csv'
d=mod.load_data(a)

summary=[]; best_rows=[]; checks=[]
qs=[0,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99,1]
for u in UNIVERSES:
  for meth in METHODS_BY_UNIVERSE[u]:
    p=RES/f'{u}__{meth}.tsv.gz'
    df=pd.read_csv(p,sep='\t')
    # Best is lexicographic from exact objective.
    key=['objective_primary','objective_secondary','objective_tertiary']
    best_tuple=min(map(tuple,df[key].to_numpy()))
    bm=(df[key[0]]==best_tuple[0])&(df[key[1]]==best_tuple[1])&(df[key[2]]==best_tuple[2])
    bdf=df[bm].copy()
    for metric in ['mean_wup','max_wup','mean_heldout_best_wup','mean_distance','min_distance','mean_heldout_nearest_distance','ancestor_descendant_pairs','train_branch_count','heldout_branch_count','training_node_total','training_edge_total','heldout_node_total','heldout_edge_total','training_internal_layers','heldout_internal_layers','row_overlap','unique_overlap','positive_coverage','lookup_micro_f1']:
      arr=df[metric].to_numpy(float)
      rec={'universe':u,'method':meth,'metric':metric,'n':len(arr),'mean':float(np.mean(arr)),'sd':float(np.std(arr,ddof=1))}
      for q,v in zip(qs,np.quantile(arr,qs)):rec[f'q{q:g}']=float(v)
      summary.append(rec)
    for _,r in bdf.iterrows():
      H=tuple(map(int,str(r.heldout_indices).split('|')));T=np.array(list(map(int,str(r.training_indices).split('|'))),dtype=int)
      role=mod.partition_roles(T,H,d)
      z=r.to_dict();z.update(role);z['best_tie_count']=len(bdf);best_rows.append(z)
      # independent evaluator check
      eh=mod.eval_hierarchy(T,np.array(H),d);eo=mod.eval_overlap(T,H,d)
      fields={'mean_wup':eh['mean_wup'],'max_wup':eh['max_wup'],'mean_heldout_best_wup':eh['mean_heldout_best_wup'],'mean_distance':eh['mean_distance'],'min_distance':eh['min_distance'],'mean_heldout_nearest_distance':eh['mean_heldout_nearest_distance'],'ancestor_descendant_pairs':eh['ancestor_descendant_pairs'],'row_overlap':eo['row_overlap'],'unique_overlap':eo['unique_overlap'],'positive_coverage':eo['positive_coverage'],'lookup_micro_f1':eo['lookup_micro_f1']}
      maxerr=max(abs(float(r[k])-float(v)) for k,v in fields.items())
      checks.append({'universe':u,'method':meth,'heldout_tissues':r.heldout_tissues,'max_abs_numeric_difference':maxerr,'pass':maxerr<1e-10})

pd.DataFrame(summary).to_csv(RES/'exact_method_distribution_summary.csv',index=False)
pd.DataFrame(best_rows).to_csv(RES/'exact_global_optima_and_balanced_roles.csv',index=False)
pd.DataFrame(checks).to_csv(ROOT/'tests/independent_exact_optima_check.csv',index=False)
print(pd.DataFrame(best_rows)[['universe','method','heldout_tissues','training_tissues','mean_wup','max_wup','row_overlap','unique_overlap','lookup_micro_f1','role_pair_A','role_pair_B','role_A_mean_wup','role_B_mean_wup','role_A_lookup_f1','role_B_lookup_f1','training_node_total','training_edge_total','train_branch_count']].to_string(index=False))
print('\nchecks all pass',all(x['pass'] for x in checks),'n',len(checks))

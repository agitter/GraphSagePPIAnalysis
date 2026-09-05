#!/usr/bin/env python3
import importlib.util
from pathlib import Path
import numpy as np, pandas as pd
ROOT=Path('/mnt/data/GraphSAGE_PPI_hierarchy_aware_split_20260904')
spec=importlib.util.spec_from_file_location('enum',ROOT/'scripts/enumerate_hierarchy_aware_splits.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class A:pass
a=A();a.root=ROOT;a.networks_dir=Path('/mnt/data/_ohmnet/bio-tissue-networks');a.hierarchy_edges=Path('/mnt/data/_ohmnet/bio-tissue-hierarchy/tissue.edges');a.layer_classification=ROOT/'prepared/ohmnet_144_layer_classification.csv';d=m.load_data(a)
rows=[]
for u in ['all144','leaf107']:
 for mode in ['uniform','matched_stratified']:
  p=ROOT/'results'/f'{u}__{mode}__first1000.tsv';df=pd.read_csv(p,sep='\t')
  mx=0;bad=0
  for _,r in df.iterrows():
   H=np.array(list(map(int,r.heldout.split('|'))),dtype=int);T=np.array(list(map(int,r.training.split('|'))),dtype=int)
   eh=m.eval_hierarchy(T,H,d);eo=m.eval_overlap(T,H,d);pr=m.partition_roles(T,H,d)
   vals={'mean_wup':eh['mean_wup'],'max_wup':eh['max_wup'],'row_overlap':eo['row_overlap'],'unique_overlap':eo['unique_overlap'],'f1':eo['lookup_micro_f1'],'role_min_f1':min(pr['role_A_lookup_f1'],pr['role_B_lookup_f1'])}
   err=max(abs(float(r[k])-v) for k,v in vals.items());mx=max(mx,err);bad+=err>2e-6
  rows.append({'universe':u,'null':mode,'records_checked':len(df),'max_abs_difference':mx,'failures_tolerance_2e-6':bad,'pass':bad==0})
out=pd.DataFrame(rows);out.to_csv(ROOT/'tests/null_first1000_independent_metric_check.csv',index=False);print(out.to_string(index=False))

#!/usr/bin/env python3
import importlib.util,itertools
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path('/mnt/data/GraphSAGE_PPI_hierarchy_aware_split_20260904')
spec=importlib.util.spec_from_file_location('enum',ROOT/'scripts/enumerate_hierarchy_aware_splits.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class A:pass
a=A();a.root=ROOT;a.networks_dir=Path('/mnt/data/_ohmnet/bio-tissue-networks');a.hierarchy_edges=Path('/mnt/data/_ohmnet/bio-tissue-hierarchy/tissue.edges');a.layer_classification=ROOT/'prepared/ohmnet_144_layer_classification.csv';d=m.load_data(a);meta=d['meta']
leaf=(meta.hierarchy_layer_type=='leaf_layer').to_numpy();E=np.flatnonzero((meta.meets_raw_15000.to_numpy()==1)&leaf);L=np.flatnonzero((meta.meets_raw_35000.to_numpy()==1)&leaf)
# exact strata as production C++
o=E[np.lexsort((E,meta.node_count.to_numpy()[E]))]; groups=np.array_split(o,5);smap={int(x):b for b,g in enumerate(groups) for x in g}
prod=pd.read_csv(ROOT/'results/leaf107__branch_distinct_node_stratified_minimax_wup.tsv.gz',sep='\t');pmap={r.heldout_indices:r for _,r in prod.iterrows()}
fail=[];n=0
for H in itertools.combinations(map(int,L),4):
 if len(set(d['branches'][list(H)]))<4:continue
 hbr=set(d['branches'][list(H)]);cand=np.array([x for x in E if d['branches'][x] not in hbr],dtype=int)
 pools=[];zs=[]
 for b in range(5):
  c=cand[[smap[int(x)]==b for x in cand]]
  if len(c)<4:break
  maxw=d['W'][np.ix_(c,H)].max(axis=1);zs.append(int(np.sort(maxw)[3]));pools.append(c)
 else:
  z=max(zs);sel=[]
  for c in pools:
   W=d['W'][np.ix_(c,H)];D=d['D'][np.ix_(c,H)];sumw=W.sum(1);maxw=W.max(1);sumd=D.sum(1);ok=np.flatnonzero(maxw<=z)
   order=ok[np.lexsort((c[ok],-sumd[ok],maxw[ok],sumw[ok]))];sel.extend(map(int,c[order[:4]]))
  sel=sorted(sel);key='|'.join(map(str,H));r=pmap[key];prodsel=sorted(map(int,r.training_indices.split('|')))
  sumw=int(d['W'][np.ix_(sel,H)].sum());maxw=int(d['W'][np.ix_(sel,H)].max());sumd=int(d['D'][np.ix_(sel,H)].sum())
  obj=(maxw,sumw,-sumd);pobj=(int(r.objective_primary),int(r.objective_secondary),int(r.objective_tertiary));n+=1
  if sel!=prodsel or obj!=pobj:fail.append({'H':key,'expected_T':'|'.join(map(str,sel)),'production_T':r.training_indices,'expected_obj':str(obj),'production_obj':str(pobj)})
out=pd.DataFrame(fail);out.to_csv(ROOT/'tests/primary_leaf_exact_enumeration_independent_check_failures.csv',index=False)
pd.DataFrame([{'eligible_train':len(E),'eligible_heldout':len(L),'all_H':len(list(itertools.combinations(L,4))),'feasible_H_checked':n,'failures':len(fail),'pass':len(fail)==0}]).to_csv(ROOT/'tests/primary_leaf_exact_enumeration_independent_check_summary.csv',index=False)
print('checked',n,'failures',len(fail))

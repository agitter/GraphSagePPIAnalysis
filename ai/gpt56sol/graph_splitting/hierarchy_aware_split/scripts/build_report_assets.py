#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, json, struct
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

ROOT=Path('/mnt/data/GraphSAGE_PPI_hierarchy_aware_split_20260904');RES=ROOT/'results';PLOTS=ROOT/'plots';PLOTS.mkdir(exist_ok=True)
DT=np.dtype([('h0','<u2'),('h1','<u2'),('h2','<u2'),('h3','<u2'),('a0','<u2'),('a1','<u2'),('b0','<u2'),('b1','<u2'),('mean_wup','<f4'),('max_wup','<f4'),('mean_best_wup','<f4'),('mean_dist','<f4'),('min_dist','<f4'),('mean_nearest_dist','<f4'),('row_overlap','<f4'),('unique_overlap','<f4'),('pos_coverage','<f4'),('f1','<f4'),('roleA_mean_wup','<f4'),('roleB_mean_wup','<f4'),('role_worst_mean_wup','<f4'),('role_absdiff_wup','<f4'),('role_mean_f1','<f4'),('role_min_f1','<f4'),('role_max_f1','<f4'),('role_mean_row_overlap','<f4'),('role_min_row_overlap','<f4'),('role_max_row_overlap','<f4'),('train_nodes','<u4'),('train_edges','<u4'),('held_nodes','<u4'),('held_edges','<u4'),('train_branches','u1'),('held_branches','u1'),('train_internal','u1'),('held_internal','u1'),('related_pairs','u1')])
def load(name):
 p=ROOT/'work'/f'{name}.bin'
 with open(p,'rb') as f:f.read(8);N=struct.unpack('<Q',f.read(8))[0]
 return np.memmap(p,mode='r',dtype=DT,offset=32,shape=(N,))
leaf_u=load('leaf107__uniform');leaf_m=load('leaf107__matched_stratified');all_u=load('all144__uniform')
exact_leaf=pd.read_csv(RES/'leaf107__branch_distinct_node_stratified_minimax_wup.tsv.gz',sep='\t')
paired=pd.read_csv(RES/'paired_conditional_on_heldout_effects.csv');paired_leaf=paired[paired.universe=='leaf107']
opt=pd.read_csv(RES/'exact_global_optima_and_balanced_roles.csv')
actual=pd.read_csv(RES/'actual_split_metrics_recomputed.csv')

# Headline comparison table.
def meanrow(x,label,universe,method):
 return {'analysis_row':label,'universe':universe,'construction':method,'N_or_count':len(x),'mean_wup':float(np.mean(x['mean_wup'])),'max_wup_mean':float(np.mean(x['max_wup'])),'row_overlap':float(np.mean(x['row_overlap'])),'unique_overlap':float(np.mean(x['unique_overlap'])),'lookup_micro_f1':float(np.mean(x['f1'])),'training_nodes':float(np.mean(x['train_nodes'])),'training_edges':float(np.mean(x['train_edges']))}
rows=[meanrow(all_u,'All-144 manuscript-threshold random null','all144','uniform valid random'),meanrow(leaf_u,'Leaf-only manuscript-threshold random null','leaf107','uniform valid random'),meanrow(leaf_m,'Leaf-only matched random null','leaf107','branch-blocked + 4-per-size-stratum random')]
rows.append({'analysis_row':'Leaf-only hierarchy-optimized ensemble mean','universe':'leaf107','construction':'exact minimax hierarchy training for each of 1,062 feasible heldout quartets','N_or_count':len(exact_leaf),'mean_wup':exact_leaf.mean_wup.mean(),'max_wup_mean':exact_leaf.max_wup.mean(),'row_overlap':exact_leaf.row_overlap.mean(),'unique_overlap':exact_leaf.unique_overlap.mean(),'lookup_micro_f1':exact_leaf.lookup_micro_f1.mean(),'training_nodes':exact_leaf.training_node_total.mean(),'training_edges':exact_leaf.training_edge_total.mean()})
for method,label in [('branch_distinct_node_stratified_minimax_wup','Leaf-only global hierarchy optimum'),('branch_distinct_coverall_node_stratified_mean_wup','Leaf-only broad-coverage hierarchy optimum')]:
 d=opt[(opt.universe=='leaf107')&(opt.method==method)].iloc[0]
 rows.append({'analysis_row':label,'universe':'leaf107','construction':method,'N_or_count':1,'mean_wup':d.mean_wup,'max_wup_mean':d.max_wup,'row_overlap':d.row_overlap,'unique_overlap':d.unique_overlap,'lookup_micro_f1':d.lookup_micro_f1,'training_nodes':d.training_node_total,'training_edges':d.training_edge_total})
a=actual[actual.label=='actual_combined_heldout'].iloc[0]
rows.append({'analysis_row':'Released GraphSAGE split (external comparison)','universe':'released24','construction':'deposited split; does not satisfy current-release literal thresholds','N_or_count':1,'mean_wup':a.mean_wup,'max_wup_mean':a.max_wup,'row_overlap':a.row_overlap,'unique_overlap':a.unique_overlap,'lookup_micro_f1':a.lookup_micro_f1,'training_nodes':44906,'training_edges':633198})
pd.DataFrame(rows).to_csv(RES/'headline_comparison.csv',index=False)

# Method definitions.
methods=[
('Primary clean universe','leaf107','Only hierarchy leaves; raw edgelist thresholds >=15,000 training and >=35,000 heldout records.'),
('Primary hard separation','distinct heldout branches','Four heldout leaves must occupy four different root branches; all layers in those branches are excluded from training.'),
('Primary size control','five node-count strata','Twenty training layers, exactly four from each eligible-node-count quintile.'),
('Primary hierarchy objective','lexicographic minimax','Minimize the largest train-heldout Wu-Palmer similarity; then the total similarity; then maximize total path distance.'),
('Role assignment','balanced 2+2 partition','Partition four heldouts into two pairs by minimizing the worse train-to-role mean similarity, then the imbalance; no gene or label data used.'),
('Primary matched null','heldout-first conditional randomization','Same leaf, threshold, branch-blocking, and size-stratum constraints; random heldout quartet and random training selections within strata.'),
('Broad-coverage sensitivity','all remaining branches represented','Add at least one training leaf from every non-heldout root branch while retaining four-per-stratum size balance.'),
('Release-mirroring sensitivity','all144','Allow both hierarchy leaves and internal hierarchy layers, under the same raw thresholds.')]
pd.DataFrame(methods,columns=['component','name','definition']).to_csv(RES/'method_definitions.csv',index=False)

# Tissues for selected candidate splits.
spec=importlib.util.spec_from_file_location('enum',ROOT/'scripts/enumerate_hierarchy_aware_splits.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class A:pass
q=A();q.root=ROOT;q.networks_dir=Path('/mnt/data/_ohmnet/bio-tissue-networks');q.hierarchy_edges=Path('/mnt/data/_ohmnet/bio-tissue-hierarchy/tissue.edges');q.layer_classification=ROOT/'prepared/ohmnet_144_layer_classification.csv';dat=mod.load_data(q)
for method,outname in [('branch_distinct_node_stratified_minimax_wup','primary_leaf_split_tissues.tsv'),('branch_distinct_coverall_node_stratified_mean_wup','broad_coverage_leaf_split_tissues.tsv')]:
 r=opt[(opt.universe=='leaf107')&(opt.method==method)].iloc[0];roleA=set(str(r.role_pair_A).split('|'));roleB=set(str(r.role_pair_B).split('|'));T=set(str(r.training_tissues).split('|'));H=set(str(r.heldout_tissues).split('|'));out=[]
 for i,name in enumerate(dat['names']):
  if name not in T|H:continue
  role='train' if name in T else ('heldout_pair_A' if name in roleA else 'heldout_pair_B')
  x=dat['meta'].iloc[i]
  out.append({'role':role,'tissue':name,'hierarchy_node':x.hierarchy_node,'top_level_branch':dat['branches'][i],'node_count':int(x.node_count),'edge_count':int(x.edge_count),'meets_15000':int(x.meets_raw_15000),'meets_35000':int(x.meets_raw_35000),'layer_type':x.hierarchy_layer_type})
 pd.DataFrame(out).sort_values(['role','top_level_branch','tissue']).to_csv(RES/outname,sep='\t',index=False)

# Metric sensitivity: WUP and path distance agree for primary leaf solution.
sens=[]
for method in ['branch_distinct_node_stratified_minimax_wup','branch_distinct_node_stratified_maximin_distance','branch_distinct_node_stratified_mean_distance']:
 p=RES/f'leaf107__{method}.tsv.gz';d=pd.read_csv(p,sep='\t');keys=['objective_primary','objective_secondary','objective_tertiary'];b=min(map(tuple,d[keys].to_numpy()));r=d[(d[keys[0]]==b[0])&(d[keys[1]]==b[1])&(d[keys[2]]==b[2])].iloc[0]
 sens.append({'method':method,'heldout_tissues':r.heldout_tissues,'training_tissues':r.training_tissues,'mean_wup':r.mean_wup,'max_wup':r.max_wup,'mean_distance':r.mean_distance,'min_distance':r.min_distance,'row_overlap':r.row_overlap,'lookup_micro_f1':r.lookup_micro_f1})
pd.DataFrame(sens).to_csv(RES/'primary_metric_sensitivity.csv',index=False)

# Discrepancy register additions.
disc=[
{'id':'SPLIT-THRESHOLD-001','statement':'The deposited GraphSAGE PPI split is incompatible with the manuscript edge thresholds when applied to the current released OhmNet edgelists.','evidence':'astrocyte and basophil training layers have 9,635 and 4,193 raw edge records; midbrain test has 31,665, below the stated 15,000/35,000 cutoffs.','status':'confirmed for current release; historical cause unresolved','report_action':'state alongside other manuscript/data discrepancies'},
{'id':'SPLIT-HIERARCHY-002','statement':'The midbrain test layer is a descendant and exact subnetwork of the brain training layer.','evidence':'all 2,310 midbrain tissue-instance rows are present in training; midbrain nodes and edges are contained in brain.','status':'confirmed','report_action':'state as a limitation of the unseen-graph interpretation'}]
pd.DataFrame(disc).to_csv(RES/'discrepancy_register_additions.csv',index=False)

# Claims register.
claims=[
('HA-001','The literal raw-edge threshold universe contains 69 training-eligible and 15 heldout-eligible hierarchy leaves.','direct computation from supplied OhmNet files','high'),
('HA-002','The primary hierarchy-aware algorithm exactly evaluated all 1,062 feasible leaf heldout quartets.','exhaustive enumeration + independent reimplementation','high'),
('HA-003','Conditional on the same heldout quartet and hard constraints, hierarchy optimization lowers mean Wu-Palmer similarity for every quartet.','paired exact/Monte Carlo comparison','high'),
('HA-004','The corresponding reduction in node overlap is typical but not universal.','paired comparison; lower row overlap for 69.6% of quartets before and 66.8% after tight size matching','high'),
('HA-005','Hierarchy separation alone does not guarantee biological-entity separation.','broad-coverage hierarchy optimum has low WUP but 98.95% row overlap and lookup F1 0.99548','high'),
('HA-006','The primary leaf optimum is unchanged when optimizing Wu-Palmer similarity, mean path distance, or minimum path distance.','metric sensitivity enumeration','high'),
('HA-007','Lookup F1 outside the original 4,301-gene universe uses the recovered global GO transformation for 209 additional OhmNet genes.','derived-label sensitivity caveat','qualified')]
pd.DataFrame(claims,columns=['claim_id','claim','evidence_basis','confidence']).to_csv(RES/'CLAIMS.csv',index=False)

# Plots.
def ecdf(arr,max_points=25000):
 a=np.sort(np.asarray(arr,float));idx=np.linspace(0,len(a)-1,min(len(a),max_points),dtype=int);return a[idx],(idx+1)/len(a)
series=[('Hierarchy-blind random',leaf_u),('Matched random',leaf_m)]
for metric,xlabel,title,fn in [
 ('mean_wup','Mean train-heldout Wu-Palmer similarity','Hierarchy separation under leaf-only valid splits','leaf107_mean_wup_ecdf.png'),
 ('row_overlap','Heldout tissue-instance overlap with training','Biological node reuse under leaf-only valid splits','leaf107_row_overlap_ecdf.png'),
 ('f1','GeneID-lookup micro-F1','Lookup diagnostic under leaf-only valid splits','leaf107_lookup_f1_ecdf.png')]:
 plt.figure(figsize=(8,5.2))
 for lab,x in series:
  xx,yy=ecdf(x[metric]);plt.plot(xx,yy,label=lab)
 em={'mean_wup':'mean_wup','row_overlap':'row_overlap','f1':'lookup_micro_f1'}[metric];xx,yy=ecdf(exact_leaf[em]);plt.plot(xx,yy,label='Hierarchy-optimized training (all feasible heldouts)')
 plt.xlabel(xlabel);plt.ylabel('Cumulative fraction');plt.title(title);plt.legend();plt.grid(True,alpha=.25);plt.tight_layout();plt.savefig(PLOTS/fn,dpi=180);plt.close()

plt.figure(figsize=(8,5.5));plt.scatter(paired_leaf.difference_mean_wup,paired_leaf.difference_lookup_micro_f1,s=12,alpha=.45);plt.axhline(0,linewidth=1);plt.axvline(0,linewidth=1);plt.xlabel('Hierarchy-aware minus matched-random mean WUP');plt.ylabel('Hierarchy-aware minus matched-random lookup F1');plt.title('Paired effects conditional on the same heldout quartet');plt.grid(True,alpha=.25);plt.tight_layout();plt.savefig(PLOTS/'leaf107_paired_wup_vs_f1_effect.png',dpi=180);plt.close()

# Global-method tradeoff.
g=opt[opt.universe=='leaf107'].copy();g=g.sort_values(['method','lookup_micro_f1']).drop_duplicates('method')
plt.figure(figsize=(9,6));plt.scatter(g.mean_wup,g.lookup_micro_f1,s=55)
for _,r in g.iterrows():
 lab=str(r.method).replace('branch_distinct_','').replace('_wup','').replace('_',' ');plt.annotate(lab,(r.mean_wup,r.lookup_micro_f1),xytext=(4,4),textcoords='offset points',fontsize=8)
plt.xlabel('Mean train-heldout Wu-Palmer similarity');plt.ylabel('Combined heldout lookup micro-F1');plt.title('Hierarchy objective versus biological identity leakage');plt.grid(True,alpha=.25);plt.tight_layout();plt.savefig(PLOTS/'leaf107_global_method_tradeoff.png',dpi=180);plt.close()

print('wrote tables and plots')

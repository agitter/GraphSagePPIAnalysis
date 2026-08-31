#!/usr/bin/env python3
from __future__ import annotations
import csv, io, json, os, zipfile, gzip, hashlib, statistics
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, percentileofscore
STAMP=os.environ['B104A_STAMP']; ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/B104A_{STAMP}'); ANA=ROOT/'analysis'; DER=ROOT/'derived'; RI=ROOT/'retained_inputs'
ROW_MAP=RI/'graphsage_row_to_entrez_topology_features.csv'; TISSUES=RI/'tissue_partition.csv'; LABELS=RI/'collapsed_gene_labels_topology_features.csv'; MAP=RI/'B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz'; MISS=RI/'B102_GraphSAGE_genes_missing_from_historical_GPI_projection_20260827T162132Z.csv'; GZ=Path('/mnt/data/graphsage_ppi.zip')

def write_csv(p,rows):
 rows=list(rows); fields=list(rows[0]) if rows else []
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');
  if fields:w.writeheader();w.writerows(rows)
rm=pd.read_csv(ROW_MAP);td=pd.read_csv(TISSUES);ldf=pd.read_csv(LABELS);mdf=pd.read_csv(MAP)
with zipfile.ZipFile(GZ) as z:
 G=json.load(z.open('ppi/ppi-G.json')); feats=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')))
n=len(G['nodes']);degree=np.zeros(n,dtype=np.int64)
for e in G['links']:
 s=int(e['source']);t=int(e['target']);degree[s]+=1;degree[t]+=1
rt=np.empty(n,dtype=object);rs=np.empty(n,dtype=object);rg=np.zeros(n,dtype=np.int64)
for r in td.itertuples(index=False):
 a=int(r.row_start_inclusive);b=int(r.row_end_exclusive);rt[a:b]=str(r.ohmnet_tissue);rs[a:b]=str(r.split);rg[a:b]=int(r.graph_index)
rows=rm.graphsage_row.astype(int).to_numpy(); rm['degree']=degree[rows];rm['feature_nonzero_count']=(feats[rows]!=0).sum(axis=1);rm['tissue']=rt[rows];rm['split']=rs[rows];rm['graph_index']=rg[rows];rm=rm.dropna(subset=['entrez_gene_id']).copy();rm.entrez_gene_id=rm.entrez_gene_id.astype(int)
lcount={int(r.entrez_gene_id):sum(int(getattr(r,f'label_{j}')) for j in range(121)) for r in ldf.itertuples(index=False)}
gstats=[]
for g,s in rm.groupby('entrez_gene_id'):
 fc=sorted(set(map(int,s.feature_nonzero_count)))
 gstats.append({'GeneID':int(g),'mapped_graphsage_row_count':len(s),'graph_instance_count':len(s),'unique_tissue_count':s.tissue.nunique(),'tissues':'|'.join(sorted(s.tissue.unique())),'train_instance_count':int((s.split=='train').sum()),'valid_instance_count':int((s.split=='valid').sum()),'test_instance_count':int((s.split=='test').sum()),'mean_instance_degree':float(s.degree.mean()),'median_instance_degree':float(s.degree.median()),'max_instance_degree':int(s.degree.max()),'min_instance_degree':int(s.degree.min()),'total_instance_degree':int(s.degree.sum()),'feature_nonzero_counts':'|'.join(map(str,fc)),'feature_nonzero_count_consistent':int(len(fc)==1),'positive_label_count':lcount.get(int(g),0),'all_zero_label_vector':int(lcount.get(int(g),0)==0)})
write_csv(DER/f'B104A_all_resolved_gene_topology_feature_label_properties_{STAMP}.csv',gstats);gsdf=pd.DataFrame(gstats)
miss=pd.read_csv(MISS);mc='Entrez_GeneID' if 'Entrez_GeneID' in miss.columns else miss.columns[0];gap=set(miss[mc].astype(int)); mapped=set(mdf.GeneID.astype(int)); rem=pd.read_csv(ANA/f'B104A_remaining_gene_label_differences_by_policy_{STAMP}.csv');resgenes=set(rem[(rem.policy=='original_six_default_relations')&(rem.difference_type=='false_positive')].GeneID.astype(int));allzero=set(gsdf.loc[gsdf.all_zero_label_vector==1,'GeneID'].astype(int))
metrics=['graph_instance_count','unique_tissue_count','mean_instance_degree','max_instance_degree','total_instance_degree','positive_label_count']
props=[]
for g in sorted(gap|resgenes|({10159} if 10159 in set(gsdf.GeneID) else set())):
 h=gsdf[gsdf.GeneID==g]
 if h.empty:continue
 r=h.iloc[0].to_dict();groups=[]
 if g in gap:groups.append('missing_from_historical_direct_GPI_projection')
 if g not in mapped:groups.append('unmatched_after_component_aware_mapping')
 if g in resgenes:groups.append('gene_in_remaining_13_default_relation_false_positives')
 if g in allzero:groups.append('all_zero_label_gene')
 o={'GeneID':g,'analysis_groups':'|'.join(groups),**r}
 for m in metrics:o[f'{m}_percentile_among_all_resolved_genes']=float(percentileofscore(gsdf[m].astype(float),float(r[m]),kind='mean'))
 props.append(o)
write_csv(ANA/f'B104A_special_gene_topology_feature_label_properties_{STAMP}.csv',props)
comp=[]
for gn,gset in [('nine_missing_from_historical_direct_GPI_projection',gap),('one_unmatched_after_component_aware_mapping',set(gsdf.GeneID)-mapped),('genes_in_remaining_13_default_relation_false_positives',resgenes),('all_zero_label_genes',allzero)]:
 a=gsdf[gsdf.GeneID.isin(gset)];b=gsdf[~gsdf.GeneID.isin(gset)]
 for m in metrics:
  x=a[m].astype(float);y=b[m].astype(float)
  if len(x)==0 or len(y)==0:continue
  u=mannwhitneyu(x,y,alternative='two-sided');comp.append({'group':gn,'group_gene_count':len(a),'comparison_gene_count':len(b),'metric':m,'group_median':float(x.median()),'comparison_median':float(y.median()),'mann_whitney_u':float(u.statistic),'two_sided_p_value':float(u.pvalue),'caveat':'Small groups, especially one gene, do not support population inference.'})
write_csv(ANA/f'B104A_special_gene_group_comparisons_{STAMP}.csv',comp)
summary={'resolved_gene_count':len(gsdf),'all_zero_gene_count':len(allzero),'historical_direct_gap_genes':sorted(gap),'unmatched_after_component_aware_mapping':sorted(set(gsdf.GeneID)-mapped),'remaining_residual_genes':sorted(resgenes),'ATP6AP2_properties':next((r for r in props if r['GeneID']==10159),None),'FSBP_25788_properties':next((r for r in props if r['GeneID']==25788),None)}
(ROOT/f'B104A_gene_property_summary_{STAMP}.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

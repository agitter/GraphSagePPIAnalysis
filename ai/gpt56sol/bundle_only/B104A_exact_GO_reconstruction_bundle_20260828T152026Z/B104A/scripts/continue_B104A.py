#!/usr/bin/env python3
from __future__ import annotations
import collections,csv,functools,gzip,hashlib,io,itertools,json,math,os,statistics,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, percentileofscore

STAMP=os.environ['B104A_STAMP']
ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/B104A_{STAMP}')
ANA=ROOT/'analysis'; DER=ROOT/'derived'; RI=ROOT/'retained_inputs'
LABELMAP=RI/'B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv'
MAP_EDGES=RI/'B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz'
TERMS_FILE=RI/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES_FILE=RI/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
GAF_FILE=RI/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
LABELS_FILE=RI/'collapsed_gene_labels_topology_features.csv'
MISSING_HIST=RI/'B102_GraphSAGE_genes_missing_from_historical_GPI_projection_20260827T162132Z.csv'
ROW_MAP=RI/'graphsage_row_to_entrez_topology_features.csv'
TISSUES=RI/'tissue_partition.csv'
GRAPHSAGE=Path('/mnt/data/graphsage_ppi.zip')
ORIG=frozenset({'EXP','IDA','IMP','IGI','IEP','ISS'})
DEFAULT=frozenset({'involved_in','part_of','enables'})
ROOTS={'biological_process':'GO:0008150','cellular_component':'GO:0005575','molecular_function':'GO:0003674'}

def det_writer(p):
 raw=p.open('wb'); gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9); txt=io.TextIOWrapper(gz,encoding='utf-8',newline=''); return raw,gz,txt

def write_csv(p,rows,fields=None,delim=','):
 rows=list(rows); fields=fields or (list(rows[0]) if rows else [])
 with p.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter=delim,lineterminator='\n',extrasaction='ignore');
  if fields:w.writeheader(); w.writerows(rows)

def write_gz(p,rows,fields=None,delim=','):
 rows=list(rows); fields=fields or (list(rows[0]) if rows else []); raw,gz,txt=det_writer(p)
 try:
  w=csv.DictWriter(txt,fieldnames=fields,delimiter=delim,lineterminator='\n',extrasaction='ignore');
  if fields:w.writeheader(); w.writerows(rows)
 finally:
  txt.flush();txt.detach();gz.close();raw.close()

def bit_genes(b,genes):
 while b:
  l=b&-b;i=l.bit_length()-1;yield genes[i];b-=l

terms=pd.read_csv(TERMS_FILE,sep='\t',dtype=str).fillna(''); edges=pd.read_csv(EDGES_FILE,sep='\t',dtype=str)
name=dict(zip(terms.GO_ID,terms.name)); ns=dict(zip(terms.GO_ID,terms.namespace)); alt={}
for r in terms.itertuples(index=False):
 if r.alt_ids:
  for a in r.alt_ids.split('|'): alt[a]=r.GO_ID
parents=collections.defaultdict(set);children=collections.defaultdict(set)
for r in edges.itertuples(index=False):parents[r.child_GO_ID].add(r.parent_GO_ID);children[r.parent_GO_ID].add(r.child_GO_ID)
@functools.lru_cache(None)
def ancestors(g):
 o={g}
 for p in parents.get(g,()):o|=ancestors(p)
 return frozenset(o)
@functools.lru_cache(None)
def adist(g):
 d={g:0};q=collections.deque([g])
 while q:
  x=q.popleft()
  for p in parents.get(x,()):
   z=d[x]+1
   if p not in d or z<d[p]:d[p]=z;q.append(p)
 return d
@functools.lru_cache(None)
def mindepth(g):
 root=ROOTS.get(ns.get(g,''))
 if g==root:return 0
 ds=[mindepth(p) for p in parents.get(g,()) if ns.get(p)==ns.get(g)];ds=[x for x in ds if x is not None]
 return 1+min(ds) if ds else None
@functools.lru_cache(None)
def maxdepth(g):
 root=ROOTS.get(ns.get(g,''))
 if g==root:return 0
 ds=[maxdepth(p) for p in parents.get(g,()) if ns.get(p)==ns.get(g)];ds=[x for x in ds if x is not None]
 return 1+max(ds) if ds else None
@functools.lru_cache(None)
def descendants(g):
 o={g}
 for c in children.get(g,()):o|=descendants(c)
 return frozenset(o)

lm=pd.read_csv(LABELMAP); ldf=pd.read_csv(LABELS_FILE); genes=ldf.entrez_gene_id.astype(int).tolist(); idx={g:i for i,g in enumerate(genes)};N=len(genes);ALL=(1<<N)-1
obs=[]
for j in range(121):
 b=0
 for i,v in enumerate(ldf[f'label_{j}']):
  if int(v):b|=1<<i
 obs.append(b)
selected=lm.GO_ID.tolist(); uniq=sorted(set(selected)); cols=collections.defaultdict(list)
for j,t in enumerate(selected):cols[t].append(j)
mdf=pd.read_csv(MAP_EDGES); acc2g=collections.defaultdict(set);method={}
for r in mdf.itertuples(index=False):acc2g[str(r.UniProtKB_accession)].add(int(r.GeneID));method[(str(r.UniProtKB_accession),int(r.GeneID))]=str(r.mapping_method)

annotations=[]; dates=set();sources=set();evidences=set();
with gzip.open(GAF_FILE,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['Is_NOT']=='1' or 'NOT' in (r['Qualifier'] or '').split('|'):continue
  gs=acc2g.get(r['DB_Object_ID'],set()); b=0
  for g in gs:
   if g in idx:b|=1<<idx[g]
  if not b:continue
  go=alt.get(r['GO_ID'],r['GO_ID']); rec=(r['DB_Object_ID'],b,go,r['Evidence_Code'],r['Normalized_Relation'],r['Date'],r['Assigned_By'],r['DB_Reference'],r['Qualifier'])
  annotations.append(rec);dates.add(r['Date']);sources.add(r['Assigned_By']);evidences.add(r['Evidence_Code'])

# Aggregate by evidence/relation/term/distance, source, date.
bucket=collections.defaultdict(int); distbucket=collections.defaultdict(int); sourcebits=collections.defaultdict(int); datebits=collections.defaultdict(int)
for acc,b,go,e,rel,date,source,ref,qual in annotations:
 d=adist(go)
 for t in set(d)&set(uniq):
  bucket[(e,rel,t)]|=b;distbucket[(e,rel,t,d[t])]|=b
  if e in ORIG and rel in DEFAULT:
   sourcebits[(source,t)]|=b;datebits[(date,t)]|=b

def pred(codes,rels=DEFAULT,maxd=None):
 out={t:0 for t in uniq}
 if maxd is None:
  for e in codes:
   for rel in rels:
    for t in uniq:out[t]|=bucket.get((e,rel,t),0)
 else:
  for (e,rel,t,d),b in distbucket.items():
   if e in codes and rel in rels and d<=maxd.get(ns.get(t,''),99):out[t]|=b
 return out

def evaluate(P):
 fp=fn=exact=ge99=ge95=0;per=[]
 for j,t in enumerate(selected):
  p=P[t];o=obs[j];x=(p&~o&ALL).bit_count();y=(o&~p&ALL).bit_count();z=x+y
  fp+=x;fn+=y;exact+=z==0;ge99+=1-z/N>=.99;ge95+=1-z/N>=.95;per.append((j,t,x,y,z,p.bit_count(),o.bit_count()))
 return {'false_positives':fp,'false_negatives':fn,'total_mismatches':fp+fn,'exact_columns':exact,'at_least_99pct':ge99,'at_least_95pct':ge95,'per':per}
P=pred(ORIG); E=evaluate(P);assert E['total_mismatches']==13 and E['false_negatives']==0

# Read best evidence outputs from first-stage run.
global_top=pd.read_csv(ANA/f'B104A_default_relation_global_evidence_mask_search_top100_{STAMP}.csv');best_global=set(global_top.iloc[0].evidence_codes.split('|'))
aspect_top=pd.read_csv(ANA/f'B104A_default_relation_aspect_specific_evidence_mask_search_top50_each_{STAMP}.csv');best_aspect={}
for a,s in aspect_top.groupby('namespace',sort=False):best_aspect[a]=set(s.iloc[0].evidence_codes.split('|'))
Pglobal=pred(best_global);Eglobal=evaluate(Pglobal)
Pasp={t:0 for t in uniq}
for a,codes in best_aspect.items():
 q=pred(codes)
 for t in uniq:
  if ns.get(t)==a:Pasp[t]=q[t]
Easp=evaluate(Pasp)

# Date cumulative search (original six/default relation) and best-global if different.
date_rows=[]
for spec,codes in [('original_six',ORIG),('best_global_after_qualifier_exclusion',frozenset(best_global))]:
 bydate=collections.defaultdict(lambda:{t:0 for t in uniq})
 for acc,b,go,e,rel,date,source,ref,qual in annotations:
  if e not in codes or rel not in DEFAULT:continue
  for t in set(adist(go))&set(uniq):bydate[date][t]|=b
 cumulative={t:0 for t in uniq}
 for date in sorted(bydate):
  for t,b in bydate[date].items():cumulative[t]|=b
  ev=evaluate(cumulative)
  date_rows.append({'filter_spec':spec,'cutoff_date_inclusive':date,**{k:v for k,v in ev.items() if k!='per'}})
write_csv(ANA/f'B104A_default_relation_annotation_date_cutoff_search_{STAMP}.csv',date_rows)

# Source leave-one-out via source-aggregated predictions.
source_rows=[]
for spec,codes in [('original_six',ORIG),('best_global_after_qualifier_exclusion',frozenset(best_global))]:
 sb=collections.defaultdict(lambda:{t:0 for t in uniq})
 for acc,b,go,e,rel,date,source,ref,qual in annotations:
  if e not in codes or rel not in DEFAULT:continue
  for t in set(adist(go))&set(uniq):sb[source][t]|=b
 baseline=pred(codes);be=evaluate(baseline)
 for excluded in sorted(sb):
  out={t:0 for t in uniq}
  for source,dct in sb.items():
   if source==excluded:continue
   for t,b in dct.items():out[t]|=b
  ev=evaluate(out)
  source_rows.append({'filter_spec':spec,'excluded_assigned_by':excluded,**{k:v for k,v in ev.items() if k!='per'},'delta_total_mismatches_vs_no_exclusion':ev['total_mismatches']-be['total_mismatches']})
write_csv(ANA/f'B104A_default_relation_source_leave_one_out_{STAMP}.csv',source_rows)

# Distance limits.
distance_rows=[]
for spec,codes in [('original_six',ORIG),('best_global_after_qualifier_exclusion',frozenset(best_global))]:
 for kc in list(range(13))+[99]:
  for kf in list(range(13))+[99]:
   ev=evaluate(pred(codes,maxd={'biological_process':99,'cellular_component':kc,'molecular_function':kf}))
   distance_rows.append({'filter_spec':spec,'CC_max_is_a_distance':kc,'MF_max_is_a_distance':kf,**{k:v for k,v in ev.items() if k!='per'}})
write_csv(ANA/f'B104A_default_relation_propagation_distance_search_{STAMP}.csv',distance_rows)

# Remaining differences per policy.
policies={'original_six_default_relations':P,'best_global_evidence_default_relations':Pglobal,'best_aspect_specific_evidence_default_relations':Pasp}
rem=[]
for pn,Q in policies.items():
 for j,t in enumerate(selected):
  fp=Q[t]&~obs[j]&ALL;fn=obs[j]&~Q[t]&ALL
  for g in bit_genes(fp,genes):rem.append({'policy':pn,'difference_type':'false_positive','label_column':j,'GO_ID':t,'GO_name':name.get(t,''),'namespace':ns.get(t,''),'GeneID':g})
  for g in bit_genes(fn,genes):rem.append({'policy':pn,'difference_type':'false_negative','label_column':j,'GO_ID':t,'GO_name':name.get(t,''),'namespace':ns.get(t,''),'GeneID':g})
write_csv(ANA/f'B104A_remaining_gene_label_differences_by_policy_{STAMP}.csv',rem)
rem13=[r for r in rem if r['policy']=='original_six_default_relations' and r['difference_type']=='false_positive'];keys={(r['label_column'],r['GO_ID'],r['GeneID']) for r in rem13}

# Witnesses.
wits=[];support=collections.defaultdict(lambda:{'direct':set(),'dist':[],'rels':set(),'acc':set(),'methods':set(),'src':set(),'date':set(),'ev':set(),'refs':set()})
for acc,b,go,e,rel,date,source,ref,qual in annotations:
 if e not in ORIG or rel not in DEFAULT:continue
 d=adist(go)
 for t in set(d)&set(uniq):
  for j in cols[t]:
   for g in bit_genes(b,genes):
    k=(j,t,g)
    if k not in keys:continue
    m=method.get((acc,g),'')
    w={'label_column':j,'selected_GO_ID':t,'selected_GO_name':name.get(t,''),'selected_namespace':ns.get(t,''),'GeneID':g,'UniProtKB_accession':acc,'mapping_method':m,'direct_GO_ID':go,'direct_GO_name':name.get(go,''),'direct_GO_min_depth':mindepth(go),'selected_GO_min_depth':mindepth(t),'shortest_is_a_distance_to_selected':d[t],'evidence':e,'normalized_relation':rel,'qualifier':qual,'date':date,'assigned_by':source,'reference':ref,'direct_equals_selected':int(go==t)};wits.append(w)
    s=support[k];s['direct'].add(go);s['dist'].append(d[t]);s['rels'].add(rel);s['acc'].add(acc);s['methods'].add(m);s['src'].add(source);s['date'].add(date);s['ev'].add(e);s['refs'].add(ref)
write_gz(ANA/f'B104A_remaining_13_witness_annotation_rows_{STAMP}.csv.gz',wits)
cnt=collections.Counter((w['label_column'],w['selected_GO_ID'],w['GeneID']) for w in wits)
details=[]
for r in rem13:
 k=(r['label_column'],r['GO_ID'],r['GeneID']);s=support[k];t=r['GO_ID'];branches=set()
 for d in s['direct']:
  if d!=t:branches|=children.get(t,set())&ancestors(d)
 details.append({**r,'selected_min_depth':mindepth(t),'selected_max_depth':maxdepth(t),'selected_descendant_count':len(descendants(t)),'witness_annotation_rows':cnt[k],'distinct_direct_GO_terms':len(s['direct']),'direct_GO_terms':'|'.join(sorted(s['direct'])),'minimum_support_distance':min(s['dist']),'maximum_support_distance':max(s['dist']),'has_direct_selected_term_annotation':int(0 in s['dist']),'immediate_child_branches_to_selected':len(branches),'immediate_child_GO_IDs':'|'.join(sorted(branches)),'accessions':'|'.join(sorted(s['acc'])),'mapping_methods':'|'.join(sorted(s['methods'])),'evidence_codes':'|'.join(sorted(s['ev'])),'assigned_by':'|'.join(sorted(s['src'])),'annotation_dates':'|'.join(sorted(s['date'])),'references':'|'.join(sorted(s['refs'])),'can_ontology_is_a_edge_drift_alone_remove_pair':int(0 not in s['dist'])})
write_csv(ANA/f'B104A_remaining_13_pair_details_{STAMP}.csv',details)

# Immediate edge drift sensitivity.
gtd=collections.defaultdict(set)
for acc,b,go,e,rel,date,source,ref,qual in annotations:
 if e not in ORIG or rel not in DEFAULT:continue
 for t in set(adist(go))&set(uniq):
  for g in bit_genes(b,genes):gtd[(g,t)].add(go)
impact=[];opt=[]
for t in sorted({r['GO_ID'] for r in rem13}):
 j=cols[t][0];o=obs[j];p=P[t];ch=children.get(t,set());sup={}
 for g in bit_genes(p,genes):
  ds=gtd[(g,t)];direct=t in ds;br=set()
  for d in ds:
   if d!=t:br|=ch&ancestors(d)
  sup[g]=(direct,br)
 candidates=sorted(set().union(*(br for g,(direct,br) in sup.items() if not (o&(1<<idx[g])))) if sup else set())
 for c in candidates:
  fpr=tpl=0
  for g,(direct,br) in sup.items():
   after=direct or bool(br-{c})
   if not after:
    if o&(1<<idx[g]):tpl+=1
    else:fpr+=1
  impact.append({'selected_GO_ID':t,'selected_GO_name':name.get(t,''),'child_GO_ID':c,'child_GO_name':name.get(c,''),'false_positives_removed_if_only_this_edge_deleted':fpr,'true_positives_lost_if_only_this_edge_deleted':tpl,'net_change_in_total_mismatches':tpl-fpr})
 if len(candidates)<=20:
  best=None
  for mask in range(1<<len(candidates)):
   removed={candidates[i] for i in range(len(candidates)) if mask&(1<<i)};fp=fn=0
   for g,(direct,br) in sup.items():
    after=direct or bool(br-removed);ob=bool(o&(1<<idx[g]));fp+=after and not ob;fn+=ob and not after
   score=(fp+fn,fn,fp,len(removed),'|'.join(sorted(removed)))
   if best is None or score<best[0]:best=(score,removed)
  score,removed=best;opt.append({'selected_GO_ID':t,'selected_GO_name':name.get(t,''),'label_columns':'|'.join(map(str,cols[t])),'candidate_immediate_edges_considered':len(candidates),'search_method':'exhaustive','minimum_mismatches_after_edge_deletions':score[0],'false_negatives_at_optimum':score[1],'false_positives_at_optimum':score[2],'edges_deleted_at_optimum':score[3],'deleted_child_GO_IDs':'|'.join(sorted(removed))})
 else:opt.append({'selected_GO_ID':t,'selected_GO_name':name.get(t,''),'label_columns':'|'.join(map(str,cols[t])),'candidate_immediate_edges_considered':len(candidates),'search_method':'not_exhaustive_gt20','minimum_mismatches_after_edge_deletions':'','false_negatives_at_optimum':'','false_positives_at_optimum':'','edges_deleted_at_optimum':'','deleted_child_GO_IDs':''})
write_csv(ANA/f'B104A_immediate_is_a_edge_single_deletion_impacts_{STAMP}.csv',impact);write_csv(ANA/f'B104A_immediate_is_a_edge_subset_optimization_{STAMP}.csv',opt)

# Graph/topology properties.
rm=pd.read_csv(ROW_MAP);td=pd.read_csv(TISSUES)
with zipfile.ZipFile(GRAPHSAGE) as z:
 G=json.load(z.open('ppi/ppi-G.json'));feats=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')))
n=len(G['nodes']);degree=np.zeros(n,dtype=np.int64)
for e in G['links']:
 s=int(e['source']);t=int(e['target']);degree[s]+=1;degree[t]+=1
rt=['']*n;rs=['']*n;rg=[0]*n
for r in td.itertuples(index=False):
 for i in range(int(r.row_start_inclusive),int(r.row_end_exclusive)):rt[i]=str(r.ohmnet_tissue);rs[i]=str(r.split);rg[i]=int(r.graph_index)
rm['degree']=degree;rm['feature_nonzero_count']=(feats!=0).sum(axis=1);rm['tissue']=rt;rm['split']=rs;rm['graph_index']=rg;rm=rm.dropna(subset=['entrez_gene_id']).copy();rm.entrez_gene_id=rm.entrez_gene_id.astype(int)
lcount={int(r.entrez_gene_id):sum(int(getattr(r,f'label_{j}')) for j in range(121)) for r in ldf.itertuples(index=False)}
gstats=[]
for g,s in rm.groupby('entrez_gene_id'):
 fc=sorted(set(map(int,s.feature_nonzero_count)))
 gstats.append({'GeneID':int(g),'graph_instance_count':len(s),'unique_tissue_count':s.tissue.nunique(),'train_instance_count':int((s.split=='train').sum()),'valid_instance_count':int((s.split=='valid').sum()),'test_instance_count':int((s.split=='test').sum()),'mean_instance_degree':float(s.degree.mean()),'median_instance_degree':float(s.degree.median()),'max_instance_degree':int(s.degree.max()),'total_instance_degree':int(s.degree.sum()),'feature_nonzero_counts':'|'.join(map(str,fc)),'feature_nonzero_count_consistent':int(len(fc)==1),'positive_label_count':lcount.get(int(g),0),'all_zero_label_vector':int(lcount.get(int(g),0)==0)})
write_csv(DER/f'B104A_all_resolved_gene_topology_feature_label_properties_{STAMP}.csv',gstats);gsdf=pd.DataFrame(gstats)
miss=pd.read_csv(MISSING_HIST);mc='Entrez_GeneID' if 'Entrez_GeneID' in miss.columns else miss.columns[0];gap=set(miss[mc].astype(int));resgenes={int(r['GeneID']) for r in rem13};mapped=set(mdf.GeneID.astype(int));allzero=set(gsdf.loc[gsdf.all_zero_label_vector==1,'GeneID'].astype(int))
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

# Summary
summary={'default_relation_original_six':{k:v for k,v in E.items() if k!='per'},'best_global_evidence':{'codes':sorted(best_global),**{k:v for k,v in Eglobal.items() if k!='per'}},'best_aspect_specific_evidence':{'codes_by_namespace':{k:sorted(v) for k,v in best_aspect.items()},**{k:v for k,v in Easp.items() if k!='per'}},'remaining_13':{'pairs':len(rem13),'direct_pairs':sum(r['has_direct_selected_term_annotation'] for r in details),'ancestor_only_pairs':sum(1-r['has_direct_selected_term_annotation'] for r in details),'genes':sorted(resgenes),'terms':sorted({r['GO_ID'] for r in rem13})},'date_cutoff_best':sorted(date_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives']))[:20],'source_leave_one_out_best':sorted(source_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives']))[:20],'distance_rule_best':sorted(distance_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives'],-r['exact_columns']))[:20],'edge_optimization':opt,'unmatched_after_component_mapping':sorted(set(gsdf.GeneID)-mapped),'all_zero_gene_count':len(allzero)}
(ROOT/f'B104A_continuation_summary_{STAMP}.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

#!/usr/bin/env python3
from __future__ import annotations
import collections, csv, functools, gzip, hashlib, html, io, json, math, re, statistics, tarfile, zipfile
from pathlib import Path
import pandas as pd
import numpy as np

STAMP='20260828T160144Z'
ROOT=Path('/mnt/data/ppi_repro_corrected/batches/B104B_'+STAMP)
ANA=ROOT/'analysis'; DER=ROOT/'derived'
ANA.mkdir(parents=True,exist_ok=True); DER.mkdir(parents=True,exist_ok=True)
BASE=Path('/mnt/data/ppi_repro_corrected/batches/B104A_20260828T145842Z')
INP=BASE/'retained_inputs'
TERMS=INP/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES=INP/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
GAF159=INP/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
GPI159=INP/'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz'
HISTMAP=INP/'B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
ACCEPTED_MAP=INP/'B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz'
LABELS=INP/'collapsed_gene_labels_topology_features.csv'
EXACT_TERMS=BASE/'analysis/B104A_exact_GO_terms_for_each_label_column_20260828T145842Z.csv'
FINAL_MAP=BASE/'analysis/B104A_final_exact_121_column_mapping_20260828T145842Z.csv'
PARTITION=INP/'tissue_partition.csv'
NETWORKS=Path('/mnt/data/bio-tissue-networks.tar.gz')
MSIG={
 '5.1':Path('/mnt/data/msigdb_v5.1_files_to_download_locally.zip'),
 '5.2':Path('/mnt/data/msigdb_v5.2_files_to_download_locally.zip'),
 '6.0':Path('/mnt/data/msigdb_v6.0_files_to_download_locally.zip'),
}
ALLOWED={'EXP','IDA','IEP','IGI','IMP','ISS'}
REL={'involved_in','part_of','enables'}

def sha256_file(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def write_csv(path,rows,fields=None):
 rows=list(rows)
 if fields is None: fields=list(rows[0].keys()) if rows else []
 with open(path,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n',extrasaction='ignore')
  if fields:w.writeheader()
  w.writerows(rows)

# ontology
terms=pd.read_csv(TERMS,sep='\t',dtype=str).fillna('')
edges=pd.read_csv(EDGES,sep='\t',dtype=str).fillna('')
name_by=dict(zip(terms.GO_ID,terms.name)); ns_by=dict(zip(terms.GO_ID,terms.namespace))
alt={}
for r in terms.itertuples(index=False):
 for a in str(r.alt_ids).split('|') if r.alt_ids else []: alt[a]=r.GO_ID
parents=collections.defaultdict(set); children=collections.defaultdict(set)
for r in edges.itertuples(index=False):
 parents[r.child_GO_ID].add(r.parent_GO_ID); children[r.parent_GO_ID].add(r.child_GO_ID)
@functools.lru_cache(None)
def ancestors(x):
 out={x}
 for p in parents.get(x,()): out|=ancestors(p)
 return frozenset(out)
@functools.lru_cache(None)
def descendants(x):
 out={x}
 for c in children.get(x,()): out|=descendants(c)
 return frozenset(out)
roots={'biological_process':'GO:0008150','cellular_component':'GO:0005575','molecular_function':'GO:0003674'}
@functools.lru_cache(None)
def min_depth(x):
 ns=ns_by.get(x,''); root=roots.get(ns)
 if x==root:return 0
 ps=[p for p in parents.get(x,()) if ns_by.get(p)==ns]
 if not ps:return None
 ds=[min_depth(p) for p in ps]
 ds=[d for d in ds if d is not None]
 return 1+min(ds) if ds else None

# labels
lab=pd.read_csv(LABELS)
genes=lab.entrez_gene_id.astype(int).tolist(); gene_idx={g:i for i,g in enumerate(genes)}; N=len(genes); ALL=(1<<N)-1
obs=[]
for j in range(121):
 b=0
 for i in np.flatnonzero(lab[f'label_{j}'].to_numpy(dtype=np.uint8)): b|=1<<int(i)
 obs.append(b)
exact_df=pd.read_csv(EXACT_TERMS)
# The matrix identifies 121 distinct GO terms as a set: three duplicated vectors each have two exact GO terms.
target_ids=set()
term_to_columns=collections.defaultdict(set)
for r in exact_df.itertuples(index=False):
 for gid in str(r.exact_GO_IDs).split('|'):
  target_ids.add(gid); term_to_columns[gid].add(int(r.label_column))
assert len(target_ids)==121
selected_primary=set(pd.read_csv(FINAL_MAP).GO_ID)

# mapping, with O95073 -> 25788 excluded
accmap=collections.defaultdict(set)
for r in pd.read_csv(ACCEPTED_MAP).itertuples(index=False): accmap[str(r.UniProtKB_accession)].add(int(r.GeneID))
accmap['O95073'].discard(25788)

# Build all propagated GO term memberships on 4,268 resolved genes.
term_bits=collections.defaultdict(int); direct_bits=collections.defaultdict(int); accepted_rows=0
with gzip.open(GAF159,'rt',encoding='utf-8') as f:
 rd=csv.DictReader(f,delimiter='\t')
 for row in rd:
  if row['Is_NOT']=='1' or 'NOT' in row['Qualifier'].split('|'): continue
  if row['Evidence_Code'] not in ALLOWED or row['Normalized_Relation'] not in REL: continue
  bits=0
  for g in accmap.get(row['DB_Object_ID'],()):
   i=gene_idx.get(g)
   if i is not None: bits|=1<<i
  if not bits: continue
  go=alt.get(row['GO_ID'],row['GO_ID'])
  direct_bits[go]|=bits
  for a in ancestors(go): term_bits[a]|=bits
  accepted_rows+=1
# Verify all target IDs are exact to their relevant label vectors.
for gid in target_ids:
 if not any(term_bits[gid]==obs[c] for c in term_to_columns[gid]):
  raise RuntimeError(f'target term {gid} not exact to expected column set')

attr_re=re.compile(r'([A-Z_]+)="([^"]*)"')
go_re=re.compile(r'GO:\d{7}')

def parse_msig(version,path):
 rows=[]
 with zipfile.ZipFile(path) as zf:
  xml=[n for n in zf.namelist() if n.lower().endswith('.xml')]
  if len(xml)!=1: raise RuntimeError((version,xml))
  with zf.open(xml[0]) as fh:
   for order,raw in enumerate(fh):
    if b'<GENESET ' not in raw: continue
    attrs={k:html.unescape(v) for k,v in attr_re.findall(raw.decode('utf-8','replace'))}
    ez=[]
    for x in attrs.get('MEMBERS_EZID','').split(','):
     if x.isdigit(): ez.append(int(x))
    bits=0
    for g in ez:
     i=gene_idx.get(g)
     if i is not None: bits|=1<<i
    m=go_re.search(attrs.get('EXTERNAL_DETAILS_URL','')+' '+attrs.get('EXACT_SOURCE',''))
    gid=m.group(0) if m else ''
    rows.append({
     'version':version,'order':len(rows),'standard_name':attrs.get('STANDARD_NAME',''),
     'systematic_name':attrs.get('SYSTEMATIC_NAME',''),'category':attrs.get('CATEGORY_CODE',''),
     'subcategory':attrs.get('SUB_CATEGORY_CODE',''),'go_id':gid,'member_count_full':len(set(ez)),
     'member_count_graph':bits.bit_count(),'bits':bits,'chip':attrs.get('CHIP',''),
     'external_url':attrs.get('EXTERNAL_DETAILS_URL',''),'build_date':attrs.get('BUILD_DATE','')
    })
 return rows

msig_rows={}; msig_by_go={}
msig_summary=[]; direct_best_rows=[]; target_msig_rows=[]
for v,p in MSIG.items():
 rows=parse_msig(v,p); msig_rows[v]=rows
 bygo=collections.defaultdict(list)
 for r in rows:
  if r['go_id']: bygo[r['go_id']].append(r)
 msig_by_go[v]=bygo
 for scope in ['C5','ALL']:
  cand=[r for r in rows if scope=='ALL' or r['category']=='C5']
  exact_cols=set(); exact_pairs=0; best=[]
  for j,ob in enumerate(obs):
   bestmis=N+1; bestnames=[]
   for r in cand:
    mis=(r['bits']^ob).bit_count()
    if mis<bestmis: bestmis=mis; bestnames=[r['standard_name']]
    elif mis==bestmis: bestnames.append(r['standard_name'])
    if mis==0:
     exact_cols.add(j); exact_pairs+=1
   best.append(bestmis)
   direct_best_rows.append({'version':v,'scope':scope,'label_column':j,'best_mismatches':bestmis,
                            'best_agreement':1-bestmis/N,'best_set_names':'|'.join(bestnames[:20]),
                            'best_set_tie_count':len(bestnames)})
  msig_summary.append({'version':v,'scope':scope,'gene_sets_tested':len(cand),
                       'exact_label_columns':len(exact_cols),'exact_set_label_pairs':exact_pairs,
                       'columns_ge_99pct':sum(x<=math.floor(.01*N) for x in best),
                       'columns_ge_95pct':sum(x<=math.floor(.05*N) for x in best),
                       'closest_mismatch':min(best),'median_best_mismatch':statistics.median(best),
                       'mean_best_mismatch':sum(best)/len(best)})
 # target GO identity and direct membership
 present=0; direct_exact=0
 for gid in sorted(target_ids):
  candidates=bygo.get(gid,[])
  if candidates: present+=1
  goabits=term_bits[gid]
  if candidates:
   best=min((r['bits']^goabits).bit_count() for r in candidates)
   exact_any=any(r['bits']==goabits for r in candidates)
   direct_exact+=int(exact_any)
   names='|'.join(r['standard_name'] for r in candidates)
   full_sizes='|'.join(str(r['member_count_full']) for r in candidates)
   graph_sizes='|'.join(str(r['member_count_graph']) for r in candidates)
   orders='|'.join(str(r['order']) for r in candidates)
  else:
   best=''; exact_any=False; names=full_sizes=graph_sizes=orders=''
  target_msig_rows.append({'version':v,'GO_ID':gid,'GO_name':name_by.get(gid,''),'namespace':ns_by.get(gid,''),
                           'present_in_C5':int(bool(candidates)),'C5_set_names':names,'C5_orders':orders,
                           'C5_full_member_counts':full_sizes,'C5_graph_member_counts':graph_sizes,
                           'GOA_exact_graph_member_count':goabits.bit_count(),
                           'direct_C5_vs_exact_GOA_mismatches':best,'direct_C5_exact':int(exact_any)})
 msig_summary.append({'version':v,'scope':'TARGET_IDENTITY','gene_sets_tested':len([r for r in rows if r['category']=='C5']),
                      'exact_label_columns':'','exact_set_label_pairs':'','columns_ge_99pct':'','columns_ge_95pct':'',
                      'closest_mismatch':'','median_best_mismatch':'','mean_best_mismatch':'',
                      'target_GO_IDs_present':present,'target_GO_IDs_total':len(target_ids),
                      'target_GO_memberships_directly_exact':direct_exact})

write_csv(ANA/f'B104B_msigdb_direct_membership_summary_{STAMP}.csv',msig_summary)
write_csv(ANA/f'B104B_msigdb_direct_best_per_label_{STAMP}.csv',direct_best_rows)
write_csv(ANA/f'B104B_target_GO_ID_presence_and_membership_by_msigdb_version_{STAMP}.csv',target_msig_rows)

# Tissue gene occurrence weights.
part=pd.read_csv(PARTITION)
selected_tissues=set(part.ohmnet_tissue)
weights={'selected24':collections.Counter(),'all144':collections.Counter(),'train20':collections.Counter(),'valid2':collections.Counter(),'test2':collections.Counter()}
network_sets={}
with tarfile.open(NETWORKS,'r:gz') as tf:
 for m in tf.getmembers():
  if not m.isfile() or not m.name.endswith('.edgelist'): continue
  tissue=Path(m.name).stem
  gs=set()
  fh=tf.extractfile(m)
  for raw in fh:
   parts=raw.decode('utf-8','replace').strip().split()
   if len(parts)>=2:
    try: gs.add(int(parts[0])); gs.add(int(parts[1]))
    except ValueError: pass
  network_sets[tissue]=gs
  for g in gs: weights['all144'][g]+=1
for r in part.itertuples(index=False):
 gs=network_sets[str(r.ohmnet_tissue)]
 if len(gs)!=int(r.node_count):
  raise RuntimeError(f'node count mismatch {r.ohmnet_tissue}: {len(gs)} vs {r.node_count}')
 for g in gs:
  weights['selected24'][g]+=1
  weights[{'train':'train20','valid':'valid2','test':'test2'}[r.split]][g]+=1

# Candidate term feature table.
feature_rows=[]
for gid,bits in term_bits.items():
 if not bits: continue
 gs=[genes[i] for i in range(N) if (bits>>i)&1]
 row={'GO_ID':gid,'GO_name':name_by.get(gid,''),'namespace':ns_by.get(gid,''),
      'selected_target_121':int(gid in target_ids),'selected_primary_118':int(gid in selected_primary),
      'unique_resolved_genes':len(gs),'direct_unique_resolved_genes':direct_bits.get(gid,0).bit_count(),
      'min_is_a_depth':min_depth(gid),'descendant_terms_including_self':len(descendants(gid))}
 for key,w in weights.items(): row[key+'_mapped_node_instances']=sum(w[g] for g in gs)
 for v in MSIG:
  cs=msig_by_go[v].get(gid,[])
  row[f'msigdb_{v}_C5_present']=int(bool(cs))
  row[f'msigdb_{v}_C5_full_size_min']=min((x['member_count_full'] for x in cs),default='')
  row[f'msigdb_{v}_C5_graph_size_min']=min((x['member_count_graph'] for x in cs),default='')
  row[f'msigdb_{v}_C5_order_min']=min((x['order'] for x in cs),default='')
 feature_rows.append(row)
features=pd.DataFrame(feature_rows)
features.to_csv(DER/f'B104B_all_GOA_term_selection_features_{STAMP}.csv.gz',index=False,compression={'method':'gzip','mtime':0})

# Selection rule evaluation.
target=target_ids
rule_rows=[]
def eval_set(name,pred,metric='',notes=''):
 pred=set(pred); tp=len(pred&target); fp=len(pred-target); fn=len(target-pred)
 prec=tp/len(pred) if pred else 0; rec=tp/len(target); f1=2*prec*rec/(prec+rec) if prec+rec else 0
 rule_rows.append({'rule':name,'metric':metric,'predicted_term_count':len(pred),'target_count':len(target),
                   'true_positives':tp,'false_positives':fp,'false_negatives':fn,
                   'precision':prec,'recall':rec,'f1':f1,'jaccard':tp/len(pred|target) if pred|target else 1,
                   'notes':notes})

metric_cols=['unique_resolved_genes','selected24_mapped_node_instances','all144_mapped_node_instances',
             'train20_mapped_node_instances','valid2_mapped_node_instances','test2_mapped_node_instances',
             'direct_unique_resolved_genes']
for metric in metric_cols:
 d=features.sort_values([metric,'GO_ID'],ascending=[False,True])
 eval_set(f'top_121_all_terms_by_{metric}',d.head(121).GO_ID,metric)
 # best contiguous rank window of 121
 ids=d.GO_ID.tolist(); best=(-1,None)
 for i in range(0,max(1,len(ids)-120)):
  tp=len(set(ids[i:i+121])&target)
  if tp>best[0]: best=(tp,i)
 if best[1] is not None:
  eval_set(f'best_121_rank_window_all_terms_by_{metric}',ids[best[1]:best[1]+121],metric,f'ranks {best[1]+1}-{best[1]+121}')
 # top per observed aspect quota
 pred=set()
 for ns,k in [('biological_process',85),('cellular_component',26),('molecular_function',10)]:
  pred|=set(features[features.namespace==ns].sort_values([metric,'GO_ID'],ascending=[False,True]).head(k).GO_ID)
 eval_set(f'top_aspect_quotas_85_26_10_by_{metric}',pred,metric)

for v in MSIG:
 pres=f'msigdb_{v}_C5_present'
 pool=features[features[pres]==1]
 eval_set(f'all_C5_terms_msigdb_{v}',pool.GO_ID,notes='GO IDs represented in C5 and in the reconstructed GOA universe')
 for metric in metric_cols:
  d=pool.sort_values([metric,'GO_ID'],ascending=[False,True])
  eval_set(f'top_121_msigdb_{v}_C5_by_{metric}',d.head(121).GO_ID,metric)
  pred=set()
  for ns,k in [('biological_process',85),('cellular_component',26),('molecular_function',10)]:
   pred|=set(pool[pool.namespace==ns].sort_values([metric,'GO_ID'],ascending=[False,True]).head(k).GO_ID)
  eval_set(f'top_aspect_quotas_msigdb_{v}_C5_by_{metric}',pred,metric)

# Threshold rules: best single lower threshold, and exact-count thresholds if any.
threshold_rows=[]
for pool_name,pool in [('all',features)]+[(f'msigdb_{v}_C5',features[features[f'msigdb_{v}_C5_present']==1]) for v in MSIG]:
 for metric in metric_cols:
  vals=sorted(set(int(x) for x in pool[metric].dropna()),reverse=True)
  for t in vals:
   pred=set(pool.loc[pool[metric]>=t,'GO_ID'])
   tp=len(pred&target); fp=len(pred-target); fn=len(target-pred)
   f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0
   threshold_rows.append({'pool':pool_name,'metric':metric,'lower_threshold':t,'predicted_term_count':len(pred),
                          'true_positives':tp,'false_positives':fp,'false_negatives':fn,'f1':f1,
                          'exact_target_set':int(pred==target)})
th=pd.DataFrame(threshold_rows)
th.to_csv(ANA/f'B104B_single_threshold_rule_grid_{STAMP}.csv.gz',index=False,compression={'method':'gzip','mtime':0})
best_threshold=th.sort_values(['pool','metric','f1','true_positives'],ascending=[True,True,False,False]).groupby(['pool','metric'],as_index=False).first()
best_threshold.to_csv(ANA/f'B104B_best_single_threshold_rules_{STAMP}.csv',index=False)

# Column order comparison. Use first occurrence of each vector and all possible target GO terms.
# For duplicate vectors, report both candidates rather than forcing an assignment.
order_rows=[]
for v in MSIG:
 bygo=msig_by_go[v]
 for r in exact_df.itertuples(index=False):
  gids=str(r.exact_GO_IDs).split('|')
  for gid in gids:
   cs=bygo.get(gid,[])
   order_rows.append({'version':v,'label_column':int(r.label_column),'GO_ID_candidate':gid,
                      'GO_name':name_by.get(gid,''),'candidate_count_for_label_vector':len(gids),
                      'present_in_C5':int(bool(cs)),'C5_order_min':min((x['order'] for x in cs),default=''),
                      'C5_standard_names':'|'.join(x['standard_name'] for x in cs)})
write_csv(ANA/f'B104B_label_column_to_msigdb_order_candidates_{STAMP}.csv',order_rows)

# Explicit summaries of selected target identity by version and term-set ambiguity.
summary={
 'generated_utc':STAMP,
 'resolved_gene_count':N,
 'label_columns':121,
 'distinct_observed_label_vectors':len(set(obs)),
 'matrix_identifiable_primary_GO_IDs':len(selected_primary),
 'union_of_all_exact_GO_ID_candidates':len(target_ids),
 'extra_exact_GO_IDs_for_three_duplicate_vector_pairs':sorted(target_ids-selected_primary),
 'accepted_GAF_rows_used':accepted_rows,
 'candidate_propagated_GO_terms_with_at_least_one_resolved_gene':len(features),
 'network_gene_unions':{
  'selected24_union':len(set(weights['selected24'])),
  'selected24_union_covered_by_4268_resolved':len(set(weights['selected24'])&set(genes)),
  'all144_union':len(set(weights['all144'])),
  'all144_union_covered_by_4268_resolved':len(set(weights['all144'])&set(genes)),
  'selected24_node_instances':sum(weights['selected24'].values()),
  'selected24_node_instances_covered_by_4268_resolved':sum(w for g,w in weights['selected24'].items() if g in gene_idx),
 },
 'msigdb_archives':{v:{'path':str(p),'sha256':sha256_file(p),'parsed_gene_sets':len(msig_rows[v]),
                         'parsed_C5_gene_sets':sum(r['category']=='C5' for r in msig_rows[v]),
                         'target_GO_IDs_present_in_C5':sum(bool(msig_by_go[v].get(g)) for g in target_ids)} for v,p in MSIG.items()},
 'exact_target_rule_found_by_single_threshold':bool(th.exact_target_set.any()),
}
with open(ANA/f'B104B_analysis_summary_{STAMP}.json','w') as f: json.dump(summary,f,indent=2,sort_keys=True)

# Missing target IDs from each C5 version.
miss=[]
for v in MSIG:
 for gid in sorted(target_ids):
  if not msig_by_go[v].get(gid):
   miss.append({'version':v,'GO_ID':gid,'GO_name':name_by.get(gid,''),'namespace':ns_by.get(gid,''),
                'resolved_gene_count':term_bits[gid].bit_count()})
write_csv(ANA/f'B104B_target_GO_IDs_absent_from_msigdb_C5_{STAMP}.csv',miss)

print(json.dumps(summary,indent=2,sort_keys=True))
print('\nMSIG SUMMARY')
print(pd.DataFrame(msig_summary).to_string(index=False))
print('\nTOP RULES')
print(pd.DataFrame(rule_rows).sort_values(['f1','true_positives'],ascending=False).head(30).to_string(index=False))
print('\nBEST THRESHOLDS')
print(best_threshold.sort_values('f1',ascending=False).head(30).to_string(index=False))

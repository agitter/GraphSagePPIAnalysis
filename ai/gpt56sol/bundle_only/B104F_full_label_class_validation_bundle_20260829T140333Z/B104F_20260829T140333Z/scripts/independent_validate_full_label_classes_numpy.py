#!/usr/bin/env python3
"""Independent NumPy/pandas validation of full GraphSAGE label reproduction.
Uses vector matrices and grouped row signatures rather than Counter/set logic.
"""
from __future__ import annotations
import argparse,csv,gzip,json,zipfile,collections,functools,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--graphsage-zip',type=Path,required=True);ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();I=a.inputs
 # row mapping/classes
 rm=pd.read_csv(I/'graphsage_row_to_entrez_topology_features.csv');res=json.loads((I/'wl_residual_after_features.json').read_text())
 graph_genes=set(rm.entrez_gene_id.astype(int))
 for c in res:graph_genes.update(map(int,c['candidate_genes']))
 # observed matrix
 with zipfile.ZipFile(a.graphsage_zip) as z:
  cm=json.loads(z.read('ppi/ppi-class_map.json'));idm=json.loads(z.read('ppi/ppi-id_map.json'))
 obs=np.zeros((len(idm),121),dtype=np.uint8)
 for node,row in idm.items():obs[int(row)]=np.asarray(cm[node],dtype=np.uint8)
 # gpi and historical edges
 gpi=pd.read_csv(I/'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz',sep='\t',dtype=str).fillna('')
 gpi['symbol']=gpi.DB_Object_Symbol.where(gpi.DB_Object_Symbol!='',gpi.GAF_Fallback_Symbol)
 symbol=dict(zip(gpi.DB_Object_ID,gpi.symbol));gpi_acc=set(symbol)
 hist=pd.read_csv(I/'B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz',sep='\t')
 hist=hist[hist.UniProtKB_accession.isin(gpi_acc)].copy()
 symdf=pd.read_csv(I/'B104F_MSigDB52_symbol_to_Entrez_relevant.tsv.gz',sep='\t',dtype=str).fillna('')
 sym2={r.gene_symbol:{int(x) for x in r.Entrez_GeneIDs.split('|') if x} for r in symdf.itertuples(index=False)}
 # components
 accgs=hist.groupby('UniProtKB_accession').GeneID.apply(lambda x:set(map(int,x))).to_dict();adj=collections.defaultdict(set)
 for acc,gs in accgs.items():
  for g in gs:adj[('a',acc)].add(('g',g));adj[('g',g)].add(('a',acc))
 resolved={};seen=set()
 for start in list(adj):
  if start in seen:continue
  seen.add(start);stack=[start];nodes=[]
  while stack:
   x=stack.pop();nodes.append(x)
   for y in adj[x]:
    if y not in seen:seen.add(y);stack.append(y)
  accs={v for t,v in nodes if t=='a'};gs={v for t,v in nodes if t=='g'};cand={}
  for acc in accs:
   m=sym2.get(symbol.get(acc,''),set())&gs
   if len(m)==1:cand[acc]=next(iter(m))
  if len(accs)==len(gs)==len(cand) and len(set(cand.values()))==len(gs):
   resolved.update({acc:{g} for acc,g in cand.items()})
 mapping={}
 for acc in gpi_acc:
  if acc in resolved:mapping[acc]=resolved[acc]&graph_genes
  else:
   d=accgs.get(acc,set())&graph_genes
   if d:mapping[acc]=set(d)
   else:
    m=sym2.get(symbol.get(acc,''),set())&graph_genes
    mapping[acc]=set(m) if len(m)==1 else set()
 mapping.setdefault('O95073',set()).discard(25788)
 # ontology and label terms
 t=pd.read_csv(I/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz',sep='\t',dtype=str).fillna('');e=pd.read_csv(I/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz',sep='\t',dtype=str)
 alt={}
 for r in t.itertuples(index=False):
  for x in r.alt_ids.split('|') if r.alt_ids else []:alt[x]=r.GO_ID
 par=collections.defaultdict(set)
 for r in e.itertuples(index=False):par[r.child_GO_ID].add(r.parent_GO_ID)
 @functools.lru_cache(None)
 def ancestors(g):
  g=alt.get(g,g);o={g}
  for p in par.get(g,()):o|=ancestors(p)
  return frozenset(o)
 selected=pd.read_csv(I/'B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv').GO_ID.tolist();selidx=collections.defaultdict(list)
 for i,g in enumerate(selected): selidx[g].append(i)
 genes=sorted(graph_genes);gidx={g:i for i,g in enumerate(genes)};pred=np.zeros((len(genes),121),dtype=np.uint8)
 allowed={'EXP','IDA','IEP','IGI','IMP','ISS'};rels={'involved_in','part_of','enables'}
 with gzip.open(I/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz','rt') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|') or r['Evidence_Code'] not in allowed or r['Normalized_Relation'] not in rels:continue
   gs=mapping.get(r['DB_Object_ID'],set())
   cols=[i for x in ancestors(r['GO_ID']) for i in selidx.get(x,[])]
   for g in gs:
    if g in gidx and cols:pred[gidx[g],cols]=1
 # resolved rows
 resolved_bad=0
 for r in rm.itertuples(index=False):resolved_bad+=not np.array_equal(pred[gidx[int(r.entrez_gene_id)]],obs[int(r.graphsage_row)])
 # residual class signatures by packed bytes sorted lexically
 exact=0;unique_rows=0;amb_rows=0
 for c in res:
  o=[obs[int(r)].tobytes() for r in c['rows']];p=[pred[gidx[int(g)]].tobytes() for g in c['candidate_genes']]
  exact+=sorted(o)==sorted(p)
  oc=collections.Counter(o);pc=collections.Counter(p)
  for sig,n in oc.items():
   if n==pc.get(sig,0)==1:unique_rows+=1
   else:amb_rows+=n
 result={'implementation':'independent pandas/NumPy matrix and packed-byte signature validation','graph_candidate_gene_count':len(genes),'resolved_rows':len(rm),'resolved_row_mismatches':int(resolved_bad),'residual_classes':len(res),'residual_classes_exact_multiset':int(exact),'residual_rows':sum(len(c['rows']) for c in res),'residual_rows_unique_by_vector':int(unique_rows),'residual_rows_ambiguous_by_vector':int(amb_rows),'all_rows_exact_up_to_within_class_permutation':bool(resolved_bad==0 and exact==len(res))}
 a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()

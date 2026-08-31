import csv,gzip,collections,functools,json
from pathlib import Path
BASE=Path('/mnt/data/ppi_repro_corrected/batches/B104A_20260828T145842Z')
INP=BASE/'retained_inputs'
GAF=INP/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
MAP=INP/'B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
TERMS=INP/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES=INP/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
EXACT=BASE/'analysis/B104A_exact_GO_terms_for_each_label_column_20260828T145842Z.csv'
OUT=Path('/mnt/data/ppi_repro_corrected/batches/B104B_20260828T161826Z/analysis/B104B_term_selection_policy_robustness.csv')
# target
target=set()
with open(EXACT,newline='') as f:
 for r in csv.DictReader(f):target.update(r['exact_GO_IDs'].split('|'))
# ontology
alt={}
with gzip.open(TERMS,'rt',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  alt[r['GO_ID']]=r['GO_ID']
  for a in r['alt_ids'].split('|'):
   if a:alt[a]=r['GO_ID']
par=collections.defaultdict(set)
with gzip.open(EDGES,'rt',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['relation']=='is_a':par[r['child_GO_ID']].add(r['parent_GO_ID'])
@functools.lru_cache(None)
def anc(x):
 s={x};st=[x]
 while st:
  y=st.pop()
  for p in par.get(y,()):
   if p not in s:s.add(p);st.append(p)
 return frozenset(s)
# mapping
acc=collections.defaultdict(set)
with gzip.open(MAP,'rt',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['in_GPI159']=='1':acc[r['UniProtKB_accession']].add(int(r['GeneID']))
acc['O95073'].discard(25788)
# rows
rows=[]
with gzip.open(GAF,'rt',newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|'):continue
  gs=acc.get(r['DB_Object_ID'])
  if not gs:continue
  rows.append((r['Evidence_Code'],r['Normalized_Relation'],alt.get(r['GO_ID'],r['GO_ID']),frozenset(gs)))
all_ev=set(x[0] for x in rows)
EVSETS={
 'exact_six':{'EXP','IDA','IEP','IGI','IMP','ISS'},
 'classic_experimental':{'EXP','IDA','IPI','IMP','IGI','IEP'},
 'exact_plus_IPI':{'EXP','IDA','IEP','IGI','IMP','ISS','IPI'},
 'all_except_IEA':all_ev-{'IEA'},
 'all_evidence':all_ev,
 'exclude_IEA_ND':all_ev-{'IEA','ND'},
}
RELSETS={
 'default_only':{'involved_in','part_of','enables'},
 'all_positive':{'involved_in','part_of','enables','colocalizes_with','contributes_to'},
}
PROP={'is_a':True,'none':False}
out=[]
for ename,evs in EVSETS.items():
 for rname,rels in RELSETS.items():
  for pname,do_prop in PROP.items():
   mem=collections.defaultdict(set); accepted=0
   for ev,rel,go,gs in rows:
    if ev not in evs or rel not in rels:continue
    terms=anc(go) if do_prop else (go,)
    for t in terms:mem[t].update(gs)
    accepted+=1
   rank=sorted(((t,len(g)) for t,g in mem.items() if g),key=lambda x:(-x[1],x[0]))
   top=set(t for t,c in rank[:121]);th=set(t for t,c in rank if c>=1000)
   out.append({'evidence_policy':ename,'relation_policy':rname,'propagation':pname,'accepted_rows':accepted,'nonempty_terms':len(rank),'top121_overlap':len(top&target),'top121_extra':len(top-target),'top121_missing':len(target-top),'top121_exact':int(top==target),'threshold1000_size':len(th),'threshold1000_overlap':len(th&target),'threshold1000_extra':len(th-target),'threshold1000_missing':len(target-th),'threshold1000_exact':int(th==target),'rank121_count':rank[120][1] if len(rank)>120 else '', 'rank122_count':rank[121][1] if len(rank)>121 else ''})
with open(OUT,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),lineterminator='\n');w.writeheader();w.writerows(out)
print('wrote',OUT)
for r in out:
 if r['top121_exact'] or r['threshold1000_exact']:
  print(r)

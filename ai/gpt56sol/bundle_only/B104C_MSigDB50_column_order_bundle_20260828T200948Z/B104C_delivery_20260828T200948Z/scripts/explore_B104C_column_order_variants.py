#!/usr/bin/env python3
from __future__ import annotations
import bisect, collections, csv, functools, gzip, itertools, json, math
from pathlib import Path
import pandas as pd
from scipy.stats import kendalltau

STAMP='20260828T194921Z'
ROOT=Path('/mnt/data/ppi_repro_corrected/batches/B104C_'+STAMP)
ANA=ROOT/'analysis'; ANA.mkdir(parents=True,exist_ok=True)
BASE=Path('/mnt/data/ppi_repro_corrected/work_B104C/B104A')
INP=BASE/'retained_inputs'
GAF=INP/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
TERMS=INP/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES=INP/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
MAP=INP/'B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz'
LABELS=INP/'collapsed_gene_labels_topology_features.csv'
EXACT=BASE/'analysis/B104A_exact_GO_terms_for_each_label_column_20260828T145842Z.csv'
ALLOWED={'EXP','IDA','IEP','IGI','IMP','ISS'}
DEFAULT_REL={'involved_in','part_of','enables'}

# CPython 2.7 legacy hash/dict simulator, validated against ppi-class_map JSON ordering.
def py2_hash(s:str,bits:int=64)->int:
    b=s.encode('ascii'); mask=(1<<bits)-1
    if not b:return 0
    x=(b[0]<<7)&mask
    for c in b:x=((1000003*x)^c)&mask
    x^=len(b); x&=mask
    if x >= 1<<(bits-1):x-=1<<bits
    if x==-1:x=-2
    return x

class Py2Dict:
    def __init__(self,bits=64):
        self.bits=bits; self.mask=7; self.table=[None]*8; self.used=0; self.fill=0
    def lookup(self,key,h):
        mask=self.mask; i=h&mask; perturb=h & ((1<<self.bits)-1)
        first=None
        while True:
            e=self.table[i]
            if e is None:return first if first is not None else i
            if e[0]==key:return i
            i=(i*5+1+perturb)&mask; perturb>>=5
    def resize(self,minused):
        new=8
        while new<=minused:new<<=1
        old=[e for e in self.table if e is not None]
        self.mask=new-1; self.table=[None]*new; self.fill=self.used=0
        for key,h in old:self._insert_no_resize(key,h)
    def _insert_no_resize(self,key,h):
        i=self.lookup(key,h)
        if self.table[i] is None:
            self.table[i]=(key,h); self.used+=1; self.fill+=1
    def insert(self,key):
        h=py2_hash(key,self.bits); before=self.used; self._insert_no_resize(key,h)
        if self.used>before and self.fill*3 >= (self.mask+1)*2:
            self.resize((2 if self.used>50000 else 4)*self.used)
    def keys(self):return [e[0] for e in self.table if e is not None]
    def slots(self):return {e[0]:i for i,e in enumerate(self.table) if e is not None}

terms=pd.read_csv(TERMS,sep='\t',dtype=str).fillna('')
term_ids=set(terms.GO_ID)
alt={}
for r in terms.itertuples(index=False):
    for a in str(r.alt_ids).split('|') if r.alt_ids else []:alt[a]=r.GO_ID
edges=pd.read_csv(EDGES,sep='\t',dtype=str).fillna('')
parents_ordered=collections.defaultdict(list)
parents_set=collections.defaultdict(set)
for r in edges.itertuples(index=False):
    if r.parent_GO_ID not in parents_set[r.child_GO_ID]:
        parents_ordered[r.child_GO_ID].append(r.parent_GO_ID)
        parents_set[r.child_GO_ID].add(r.parent_GO_ID)
@functools.lru_cache(None)
def ancestors_set(g):
    out={g}
    for p in parents_set.get(g,()):out |= ancestors_set(p)
    return frozenset(out)
@functools.lru_cache(None)
def anc_dfs_pre(g):
    out=[]; seen=set()
    def rec(x):
        if x in seen:return
        seen.add(x); out.append(x)
        for p in parents_ordered.get(x,[]):rec(p)
    rec(g);return tuple(out)
@functools.lru_cache(None)
def anc_dfs_post(g):
    out=[]; seen=set()
    def rec(x):
        if x in seen:return
        seen.add(x)
        for p in parents_ordered.get(x,[]):rec(p)
        out.append(x)
    rec(g);return tuple(out)
@functools.lru_cache(None)
def anc_bfs(g):
    out=[]; seen={g}; q=collections.deque([g])
    while q:
        x=q.popleft(); out.append(x)
        for p in parents_ordered.get(x,[]):
            if p not in seen:seen.add(p);q.append(p)
    return tuple(out)

def py2_small_dict_order(keys):
    d=Py2Dict(64)
    for k in keys:d.insert(k)
    return tuple(d.keys())
@functools.lru_cache(None)
def anc_py2_set_proxy(g):
    # A plausible proxy for iteration of a Python2 hash container populated in sorted order.
    return py2_small_dict_order(sorted(ancestors_set(g)))

# Mapping coverage categories.
acc_to_gene=collections.defaultdict(set)
for r in pd.read_csv(MAP).itertuples(index=False):
    acc_to_gene[str(r.UniProtKB_accession)].add(int(r.GeneID))
acc_to_gene['O95073'].discard(25788)
graph_genes=set(pd.read_csv(LABELS).entrez_gene_id.astype(int))

# Exact target assignments.
exact=pd.read_csv(EXACT).sort_values('label_column')
allowed_by_col=[tuple(sorted(str(x).split('|'))) for x in exact.exact_GO_IDs]
target=set().union(*(set(x) for x in allowed_by_col))
dup_groups=[]; seen=[]
for ids in allowed_by_col:
    if len(ids)>1 and ids not in seen:
        dup_groups.append(([i for i,x in enumerate(allowed_by_col) if x==ids],ids));seen.append(ids)
def assignments():
    base=[ids[0] if len(ids)==1 else None for ids in allowed_by_col]
    for flips in itertools.product((0,1),repeat=len(dup_groups)):
        seq=base[:]
        for (cols,ids),flip in zip(dup_groups,flips):
            vals=list(ids)
            if flip:vals.reverse()
            for c,v in zip(cols,vals):seq[c]=v
        yield flips,seq

def lcs_len(a,b):
    rank={x:i for i,x in enumerate(a)}; tails=[]
    for x in b:
        if x not in rank:continue
        v=rank[x]; i=bisect.bisect_left(tails,v)
        if i==len(tails):tails.append(v)
        else:tails[i]=v
    return len(tails)
def score(pred):
    pred=[x for x in pred if x in target]
    if len(pred)!=121 or len(set(pred))!=121:return None
    best=None
    for flips,obs in assignments():
        rank={x:i for i,x in enumerate(pred)}; vals=[rank[x] for x in obs]
        tau=kendalltau(range(121),vals)
        exact=sum(a==b for a,b in zip(pred,obs)); prefix=0
        for a,b in zip(pred,obs):
            if a!=b:break
            prefix+=1
        lcs=lcs_len(pred,obs)
        key=(lcs,float(tau.statistic),exact,prefix)
        if best is None or key>best[0]:best=(key,flips,obs,float(tau.pvalue))
    (lcs,tau,exact,prefix),flips,obs,p=best
    concord=(tau+1)/2
    return dict(kendall_tau=tau,kendall_p=p,pairwise_concordance=concord,lcs=lcs,exact_positions=exact,exact_prefix=prefix,duplicate_flips=''.join(map(str,flips)),assigned_sequence=obs)

# Read all GAF rows minimally once.
rows=[]
with gzip.open(GAF,'rt',encoding='utf-8') as f:
    rd=csv.DictReader(f,delimiter='\t')
    for idx,r in enumerate(rd):
        go=alt.get(r['GO_ID'],r['GO_ID'])
        acc=r['DB_Object_ID']; genes=acc_to_gene.get(acc,set())
        rows.append(dict(idx=idx,go=go,acc=acc,symbol=r['DB_Object_Symbol'],ev=r['Evidence_Code'],rel=r['Normalized_Relation'],is_not=(r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|')),date=r['Date'],mapped=bool(genes),graph_mapped=bool(genes&graph_genes),min_gene=min(genes) if genes else 10**15))

filters={
 'exact_all_accessions':lambda r:(not r['is_not'] and r['ev'] in ALLOWED and r['rel'] in DEFAULT_REL),
 'exact_historical_mapped_accessions':lambda r:(not r['is_not'] and r['ev'] in ALLOWED and r['rel'] in DEFAULT_REL and r['mapped']),
 'all_evidence_default_rel':lambda r:(not r['is_not'] and r['rel'] in DEFAULT_REL),
}
orders={
 'source':lambda rs:rs,
 'GO_accession':lambda rs:sorted(rs,key=lambda r:(r['go'],r['acc'],r['idx'])),
 'GeneID_GO':lambda rs:sorted(rs,key=lambda r:(r['min_gene'],r['go'],r['acc'],r['idx'])),
 'date_accession_GO':lambda rs:sorted(rs,key=lambda r:(r['date'],r['acc'],r['go'],r['idx'])),
}
anc_modes={
 'sorted_asc':lambda g:tuple(sorted(ancestors_set(g))),
 'dfs_pre_obo':anc_dfs_pre,
 'py2_hash_container_proxy':anc_py2_set_proxy,
 'hash_low15':lambda g:tuple(sorted(ancestors_set(g),key=lambda x:(py2_hash(x)&32767,x))),
}

out=[]; best_detail=None
# Plausible per-row direct+ancestor insertion variants.
for fname,ff in filters.items():
    fr=[r for r in rows if ff(r)]
    for oname,of in orders.items():
        ordered=of(fr)
        # Only the first occurrence of each direct GO term can introduce new keys;
        # collapsing repeats preserves dictionary key-insertion history exactly.
        direct_first=[]; seen_direct=set()
        for r in ordered:
            if r['go'] not in seen_direct:
                seen_direct.add(r['go']); direct_first.append(r['go'])
        for aname,af in anc_modes.items():
            d=Py2Dict(64)
            for g in direct_first:
                d.insert(g)
                for a in af(g):d.insert(a)
            sc=score(d.keys())
            if sc:
                rec={k:v for k,v in sc.items() if k!='assigned_sequence'}
                rec.update(model='first_direct_occurrence_then_ancestors',filter=fname,row_order=oname,ancestor_order=aname,unique_keys=d.used,table_size=d.mask+1,direct_term_count=len(direct_first))
                out.append(rec)
                if best_detail is None or (rec['lcs'],rec['kendall_tau'],rec['exact_positions'],rec['exact_prefix'])>(best_detail[0],best_detail[1],best_detail[2],best_detail[3]):
                    best_detail=(rec['lcs'],rec['kendall_tau'],rec['exact_positions'],rec['exact_prefix'],rec,sc['assigned_sequence'],d.keys())

res=pd.DataFrame(out).sort_values(['lcs','kendall_tau','exact_positions','exact_prefix'],ascending=False)
res.to_csv(ANA/f'B104C_extended_python2_dictionary_order_simulation_grid_{STAMP}.csv',index=False)
# top summary
res.head(100).to_csv(ANA/f'B104C_extended_python2_dictionary_order_top100_{STAMP}.csv',index=False)

# table-size fingerprint for raw hash masks.
mask_rows=[]
for n in range(7,21):
    mask=(1<<n)-1
    pred=sorted(target,key=lambda x:(py2_hash(x)&mask,x))
    s=score(pred)
    mask_rows.append({k:v for k,v in s.items() if k!='assigned_sequence'}|{'mask_bits':n,'table_size':1<<n})
pd.DataFrame(mask_rows).to_csv(ANA/f'B104C_python2_hash_table_size_fingerprint_{STAMP}.csv',index=False)

# Best model per-column sequence and slot diagnostics.
if best_detail:
    _,_,_,_,rec,obs,pred_all=best_detail
    d=Py2Dict(64)
    # rebuild from stored best impossible without parsing rec, but score predicted target sequence from best simulation is in last run only.
    # Find target prediction by rerunning matching record through a helper not worth complexity; save best metadata and observed orientation.
    with open(ANA/f'B104C_extended_order_best_model_{STAMP}.json','w') as f:
        json.dump({**rec,'assigned_sequence':obs},f,indent=2,sort_keys=True)

summary={
 'simulations_scored':len(res),
 'best':res.iloc[0].to_dict() if len(res) else {},
 'models_with_table_size_32768':int((res.table_size==32768).sum()) if len(res) else 0,
 'all_models_duplicate_orientation_counts':res.duplicate_flips.value_counts().to_dict() if len(res) else {},
 'mask_fingerprint_best':pd.DataFrame(mask_rows).sort_values(['lcs','kendall_tau'],ascending=False).iloc[0].to_dict(),
}
with open(ANA/f'B104C_extended_order_simulation_summary_{STAMP}.json','w') as f:json.dump(summary,f,indent=2,sort_keys=True)
print(json.dumps(summary,indent=2,sort_keys=True))

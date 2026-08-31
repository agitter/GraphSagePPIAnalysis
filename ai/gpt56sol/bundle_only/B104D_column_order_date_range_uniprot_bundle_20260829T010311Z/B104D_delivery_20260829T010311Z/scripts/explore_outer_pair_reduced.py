#!/usr/bin/env python3
from __future__ import annotations
import bisect, collections, csv, functools, gzip, itertools, json, math, os, re, sys, hashlib
from pathlib import Path
import pandas as pd
from scipy.stats import kendalltau

ROOT=Path('/mnt/data/ppi_repro_corrected/work_B104D')
B104A=ROOT/'b104a_bundle/B104A'
B104C=ROOT/'b104c_bundle/B104C_delivery_20260828T200948Z'
OUT=ROOT/'analysis_outerpair'; OUT.mkdir(parents=True,exist_ok=True)
GAF=B104A/'retained_inputs/B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
TERMS=B104A/'retained_inputs/B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES=B104A/'retained_inputs/B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
MAP=B104A/'retained_inputs/B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz'
LABELS=B104A/'retained_inputs/collapsed_gene_labels_topology_features.csv'
EXACT=B104A/'analysis/B104A_exact_GO_terms_for_each_label_column_20260828T145842Z.csv'
ALLOWED={'EXP','IDA','IEP','IGI','IMP','ISS'}
DEFAULT_REL={'involved_in','part_of','enables'}
MASK64=(1<<64)-1

# ---------- exact CPython 2 style hashes / tables ----------
def py2_str_hash(s:str,bits:int=64)->int:
    b=s.encode('ascii','strict'); mask=(1<<bits)-1
    if not b:return 0
    x=(b[0]<<7)&mask
    for c in b:x=((1000003*x)^c)&mask
    x^=len(b); x&=mask
    if x >= 1<<(bits-1):x-=1<<bits
    if x==-1:x=-2
    return x

def py2_int_hash(n:int,bits:int=64)->int:
    # All GO numeric identifiers are small enough that Python2 int hash is identity.
    mask=(1<<bits)-1
    x=n & mask
    if x >= 1<<(bits-1): x-=1<<bits
    if x==-1:x=-2
    return x

def py2_tuple_hash(items, bits:int=64)->int:
    # CPython 2.7 tuplehash (historical algorithm).
    mask=(1<<bits)-1
    mult=1000003
    x=0x345678
    n=len(items)
    for i,item in enumerate(items):
        y=hash_key(item,bits)
        x=((x ^ (y & mask))*mult)&mask
        mult=(mult + 82520 + (n-i-1)*2)&mask
    x=(x+97531)&mask
    if x >= 1<<(bits-1):x-=1<<bits
    if x==-1:x=-2
    return x

def hash_key(k,bits=64):
    if isinstance(k,str):return py2_str_hash(k,bits)
    if isinstance(k,int):return py2_int_hash(k,bits)
    if isinstance(k,tuple):return py2_tuple_hash(k,bits)
    raise TypeError(type(k))

class Py2Dict:
    def __init__(self,bits=64):
        self.bits=bits; self.mask=7; self.table=[None]*8; self.used=0; self.fill=0
    def lookup(self,key,h):
        mask=self.mask; i=h&mask; perturb=h & ((1<<self.bits)-1)
        first_dummy=None
        while True:
            e=self.table[i]
            if e is None:return first_dummy if first_dummy is not None else i
            if e is DUMMY:
                if first_dummy is None:first_dummy=i
            elif e[0]==key:return i
            i=(i*5+1+perturb)&mask; perturb>>=5
    def resize(self,minused):
        new=8
        while new<=minused:new<<=1
        old=[e for e in self.table if e is not None and e is not DUMMY]
        self.mask=new-1; self.table=[None]*new; self.fill=self.used=0
        for key,h,val in old:self._insert_no_resize(key,h,val)
    def _insert_no_resize(self,key,h,val=None):
        i=self.lookup(key,h); old=self.table[i]
        if old is None:
            self.table[i]=(key,h,val); self.used+=1; self.fill+=1; return True
        if old is DUMMY:
            self.table[i]=(key,h,val); self.used+=1; return True
        self.table[i]=(key,h,val); return False
    def insert(self,key,val=None):
        h=hash_key(key,self.bits); added=self._insert_no_resize(key,h,val)
        if added and self.fill*3 >= (self.mask+1)*2:
            self.resize((2 if self.used>50000 else 4)*self.used)
        return added
    def insert_h(self,key,h,val=None):
        added=self._insert_no_resize(key,h,val)
        if added and self.fill*3 >= (self.mask+1)*2:
            self.resize((2 if self.used>50000 else 4)*self.used)
        return added
    def delete(self,key):
        h=hash_key(key,self.bits); i=self.lookup(key,h)
        if self.table[i] is None or self.table[i] is DUMMY:return False
        self.table[i]=DUMMY; self.used-=1; return True
    def keys(self):return [e[0] for e in self.table if e is not None and e is not DUMMY]
    def items(self):return [(e[0],e[2]) for e in self.table if e is not None and e is not DUMMY]
    def slots(self):return {e[0]:i for i,e in enumerate(self.table) if e is not None and e is not DUMMY}
    def __contains__(self,key):
        h=hash_key(key,self.bits); i=self.lookup(key,h); e=self.table[i]
        return e is not None and e is not DUMMY and e[0]==key
    def get(self,key,default=None):
        h=hash_key(key,self.bits); i=self.lookup(key,h); e=self.table[i]
        return default if e is None or e is DUMMY else e[2]
DUMMY=object()

class Py2Set:
    # Two load-policy variants are exposed because CPython set resize thresholds changed historically.
    def __init__(self,bits=64,policy='three_fifths'):
        self.bits=bits; self.policy=policy; self.mask=7; self.table=[None]*8; self.used=0; self.fill=0
    def lookup(self,key,h):
        mask=self.mask; i=h&mask; perturb=h & ((1<<self.bits)-1)
        first_dummy=None
        while True:
            e=self.table[i]
            if e is None:return first_dummy if first_dummy is not None else i
            if e is DUMMY:
                if first_dummy is None:first_dummy=i
            elif e[0]==key:return i
            i=(i*5+1+perturb)&mask; perturb>>=5
    def _insert_no_resize(self,key,h):
        i=self.lookup(key,h); old=self.table[i]
        if old is None:self.table[i]=(key,h);self.used+=1;self.fill+=1;return True
        if old is DUMMY:self.table[i]=(key,h);self.used+=1;return True
        return False
    def resize(self,minused):
        new=8
        while new<=minused:new<<=1
        old=[e for e in self.table if e is not None and e is not DUMMY]
        self.mask=new-1;self.table=[None]*new;self.used=self.fill=0
        for k,h in old:self._insert_no_resize(k,h)
    def add(self,key):
        h=hash_key(key,self.bits); added=self._insert_no_resize(key,h)
        if not added:return False
        n=self.mask+1
        trigger = self.fill*5 >= n*3 if self.policy=='three_fifths' else self.fill*3 >= n*2
        if trigger:self.resize(4*self.used if self.used<=50000 else 2*self.used)
        return True
    def keys(self):return [e[0] for e in self.table if e is not None and e is not DUMMY]

# ---------- load ontology / target ----------
print('stage terms',flush=True)
terms=pd.read_csv(TERMS,sep='\t',dtype=str).fillna('')
term_ids=list(terms.GO_ID)
term_set=set(term_ids)
name_by_id=dict(zip(terms.GO_ID,terms.name))
ns_by_id=dict(zip(terms.GO_ID,terms.namespace))
obo_index={g:i for i,g in enumerate(term_ids)}
alt={}
alt_list=[]; alt_by_primary=collections.defaultdict(list)
for r in terms.itertuples(index=False):
    aids=[a for a in str(r.alt_ids).split('|') if a]
    alt_list.extend((a,r.GO_ID) for a in aids); alt_by_primary[r.GO_ID].extend(aids)
    for a in aids:alt[a]=r.GO_ID
edges=pd.read_csv(EDGES,sep='\t',dtype=str).fillna('')
parents_ordered=collections.defaultdict(list)
parents_set=collections.defaultdict(set)
for r in edges.itertuples(index=False):
    c,p=r.child_GO_ID,r.parent_GO_ID
    if p not in parents_set[c]:parents_set[c].add(p);parents_ordered[c].append(p)
@functools.lru_cache(None)
def ancestors_set(g):
    out={g}
    for p in parents_set.get(g,()):out |= ancestors_set(p)
    return frozenset(out)
@functools.lru_cache(None)
def anc_dfs_pre(g):
    out=[];seen=set()
    def rec(x):
        if x in seen:return
        seen.add(x);out.append(x)
        for p in parents_ordered.get(x,[]):rec(p)
    rec(g);return tuple(out)
@functools.lru_cache(None)
def anc_dfs_post(g):
    out=[];seen=set()
    def rec(x):
        if x in seen:return
        seen.add(x)
        for p in parents_ordered.get(x,[]):rec(p)
        out.append(x)
    rec(g);return tuple(out)
@functools.lru_cache(None)
def anc_bfs_source(g):
    out=[];seen={g};q=collections.deque([g])
    while q:
        x=q.popleft();out.append(x)
        for p in parents_ordered.get(x,[]):
            if p not in seen:seen.add(p);q.append(p)
    return tuple(out)
@functools.lru_cache(None)
def anc_py2set_sorted_35(g):
    s=Py2Set(policy='three_fifths')
    for x in sorted(ancestors_set(g)):s.add(x)
    return tuple(s.keys())
@functools.lru_cache(None)
def anc_py2set_dfs_35(g):
    s=Py2Set(policy='three_fifths')
    for x in anc_dfs_pre(g):s.add(x)
    return tuple(s.keys())
@functools.lru_cache(None)
def anc_py2set_bfs_35(g):
    s=Py2Set(policy='three_fifths')
    for x in anc_bfs_source(g):s.add(x)
    return tuple(s.keys())
@functools.lru_cache(None)
def anc_py2set_sorted_23(g):
    s=Py2Set(policy='two_thirds')
    for x in sorted(ancestors_set(g)):s.add(x)
    return tuple(s.keys())

print('stage exact',flush=True)
exact=pd.read_csv(EXACT).sort_values('label_column')
allowed_by_col=[tuple(sorted(str(x).split('|'))) for x in exact.exact_GO_IDs]
target=set().union(*(set(x) for x in allowed_by_col))
dup_groups=[]; seen_dups=[]
for ids in allowed_by_col:
    if len(ids)>1 and ids not in seen_dups:
        dup_groups.append(([i for i,x in enumerate(allowed_by_col) if x==ids],ids));seen_dups.append(ids)
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
    rank={x:i for i,x in enumerate(a)};tails=[]
    for x in b:
        if x not in rank:continue
        v=rank[x];i=bisect.bisect_left(tails,v)
        if i==len(tails):tails.append(v)
        else:tails[i]=v
    return len(tails)

def score(pred_ids):
    pred=[x for x in pred_ids if x in target]
    # preserve first occurrence; key representations can collapse names.
    seen=set(); pred=[x for x in pred if not (x in seen or seen.add(x))]
    if len(pred)!=121 or len(set(pred))!=121:return None
    best=None
    for flips,obs in assignments():
        rank={x:i for i,x in enumerate(pred)}; vals=[rank[x] for x in obs]
        tau=kendalltau(range(121),vals)
        ex=sum(a==b for a,b in zip(pred,obs));pref=0
        for a,b in zip(pred,obs):
            if a!=b:break
            pref+=1
        l=lcs_len(pred,obs); key=(l,float(tau.statistic),ex,pref)
        if best is None or key>best[0]:best=(key,flips,obs,float(tau.pvalue))
    (l,t,ex,pref),flips,obs,p=best
    return dict(lcs=l,kendall_tau=t,kendall_p=p,pairwise_concordance=(t+1)/2,exact_positions=ex,exact_prefix=pref,duplicate_flips=''.join(map(str,flips)),assigned_observed='|'.join(obs),predicted='|'.join(pred))

# ---------- mapping / rows ----------
acc_to_gene=collections.defaultdict(set)
for r in pd.read_csv(MAP).itertuples(index=False):acc_to_gene[str(r.UniProtKB_accession)].add(int(r.GeneID))
acc_to_gene['O95073'].discard(25788)
graph_genes=set(pd.read_csv(LABELS).entrez_gene_id.astype(int))
print('stage rows',flush=True)
rows=[]
with gzip.open(GAF,'rt',encoding='utf-8') as f:
    rd=csv.DictReader(f,delimiter='\t')
    for idx,r in enumerate(rd):
        go=alt.get(r['GO_ID'],r['GO_ID']);acc=r['DB_Object_ID'];genes=acc_to_gene.get(acc,set())
        rows.append(dict(idx=idx,go=go,acc=acc,symbol=r['DB_Object_Symbol'],ev=r['Evidence_Code'],rel=r['Normalized_Relation'],is_not=(r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|')),date=r['Date'],assigned=r['Assigned_By'],mapped=bool(genes),graph_mapped=bool(genes&graph_genes),genes=tuple(sorted(genes)),min_gene=min(genes) if genes else 10**15,aspect=r['Aspect']))

print('rows loaded',len(rows),flush=True)
filters={
 'all_rows':lambda r:True,
 'positive_all_rel_all_ev':lambda r:not r['is_not'],
 'positive_default_rel_all_ev':lambda r:not r['is_not'] and r['rel'] in DEFAULT_REL,
 'positive_all_rel_allowed_ev':lambda r:not r['is_not'] and r['ev'] in ALLOWED,
 'exact_all_accessions':lambda r:not r['is_not'] and r['ev'] in ALLOWED and r['rel'] in DEFAULT_REL,
 'exact_historical_mapped':lambda r:not r['is_not'] and r['ev'] in ALLOWED and r['rel'] in DEFAULT_REL and r['mapped'],
 'exact_graph_mapped':lambda r:not r['is_not'] and r['ev'] in ALLOWED and r['rel'] in DEFAULT_REL and r['graph_mapped'],
 'all_ev_default_rel_mapped':lambda r:not r['is_not'] and r['rel'] in DEFAULT_REL and r['mapped'],
 'allowed_ev_all_rel_mapped':lambda r:not r['is_not'] and r['ev'] in ALLOWED and r['mapped'],
}
orders={
 'source':lambda rs:list(rs),
 'source_reverse':lambda rs:list(reversed(rs)),
 'GO_accession':lambda rs:sorted(rs,key=lambda r:(r['go'],r['acc'],r['idx'])),
 'accession_GO':lambda rs:sorted(rs,key=lambda r:(r['acc'],r['go'],r['idx'])),
 'GeneID_GO':lambda rs:sorted(rs,key=lambda r:(r['min_gene'],r['go'],r['acc'],r['idx'])),
 'GO_GeneID':lambda rs:sorted(rs,key=lambda r:(r['go'],r['min_gene'],r['acc'],r['idx'])),
 'date_accession_GO':lambda rs:sorted(rs,key=lambda r:(r['date'],r['acc'],r['go'],r['idx'])),
 'assigned_accession_GO':lambda rs:sorted(rs,key=lambda r:(r['assigned'],r['acc'],r['go'],r['idx'])),
 'evidence_accession_GO':lambda rs:sorted(rs,key=lambda r:(r['ev'],r['acc'],r['go'],r['idx'])),
 'symbol_GO':lambda rs:sorted(rs,key=lambda r:(r['symbol'],r['go'],r['acc'],r['idx'])),
}
anc_modes={
 'sorted_asc':lambda g:tuple(sorted(ancestors_set(g))),
 'sorted_desc':lambda g:tuple(sorted(ancestors_set(g),reverse=True)),
 'dfs_pre_obo':anc_dfs_pre,
 'dfs_post_obo':anc_dfs_post,
 'bfs_obo':anc_bfs_source,
 'py2set_sorted_35':anc_py2set_sorted_35,
 'py2set_dfs_35':anc_py2set_dfs_35,
 'py2set_bfs_35':anc_py2set_bfs_35,
 'py2set_sorted_23':anc_py2set_sorted_23,
 'hash_low15':lambda g:tuple(sorted(ancestors_set(g),key=lambda x:(py2_str_hash(x)&32767,x))),
}

# key representations. Return actual dict key; output term is carried separately.
def std_name(g):
    s=name_by_id.get(g,g).upper()
    s=re.sub(r'[^A-Z0-9]+','_',s).strip('_')
    return 'GO_'+s
key_reprs={
 'GO_ID':lambda g:g,
 'GO_ID_no_colon':lambda g:g.replace(':',''),
 'GO_ID_numeric_string':lambda g:str(int(g.split(':')[1])),
 'GO_numeric_int':lambda g:int(g.split(':')[1]),
 'GO_name':lambda g:name_by_id.get(g,g),
 'MSigDB_standard_name':std_name,
 'namespace_pipe_GO':lambda g:ns_by_id.get(g,'')+'|'+g,
 'GO_pipe_name':lambda g:g+'|'+name_by_id.get(g,''),
 'tuple_namespace_GO':lambda g:(ns_by_id.get(g,''),g),
}

_key_cache={}
_hash_cache={}
def dict_order_from_term_sequence(term_seq,key_mode):
    d=Py2Dict(64); keyfun=key_reprs[key_mode]; seen_fast=set()
    kc=_key_cache.setdefault(key_mode,{})
    hc=_hash_cache.setdefault(key_mode,{})
    for g in term_seq:
        k=kc.get(g)
        if k is None:
            k=keyfun(g);kc[g]=k;hc[g]=hash_key(k,64)
        if k in seen_fast:continue
        seen_fast.add(k); d.insert_h(k,hc[g],g)
    return [g for k,g in d.items()]

results=[]
def add_result(model, filter_name='', row_order='', ancestor_order='', key_mode='GO_ID', term_seq=None, extra=None):
    if term_seq is None:return
    pred=dict_order_from_term_sequence(term_seq,key_mode)
    sc=score(pred)
    if sc is None:return
    rec={k:v for k,v in sc.items() if k not in ('assigned_observed','predicted')}
    rec.update(model=model,filter=filter_name,row_order=row_order,ancestor_order=ancestor_order,key_mode=key_mode,input_terms=len(term_seq),unique_dictionary_terms=len(set(term_seq)))
    if extra:rec.update(extra)
    rec['assigned_observed']=sc['assigned_observed'];rec['predicted']=sc['predicted']
    results.append(rec)


# Optimized entity-first and annotation-pair models.

def build_outer_base_sequence(ordered, outer_kind, inner_kind, set_policy='three_fifths'):
    outer_first=[]; seen_outer=set(); term_lists=collections.defaultdict(list); term_seen=collections.defaultdict(set)
    for r in ordered:
        if outer_kind=='accession': oks=(r['acc'],)
        elif outer_kind=='symbol': oks=(r['symbol'],)
        elif outer_kind=='gene_int': oks=r['genes']
        elif outer_kind=='gene_str': oks=tuple(str(g) for g in r['genes'])
        else: raise ValueError(outer_kind)
        for ok in oks:
            if ok not in seen_outer:seen_outer.add(ok);outer_first.append(ok)
            if r['go'] not in term_seen[ok]:term_seen[ok].add(r['go']);term_lists[ok].append(r['go'])
    od=Py2Dict(64)
    for ok in outer_first:od.insert(ok,ok)
    seq=[]
    for ok,_ in od.items():
        if inner_kind=='dict':
            inn=Py2Dict(64)
            for g in term_lists[ok]:inn.insert(g,g)
            gos=[g for g,_ in inn.items()]
        else:
            inn=Py2Set(64,set_policy)
            for g in term_lists[ok]:inn.add(g)
            gos=inn.keys()
        seq.extend(gos)
    return seq,len(outer_first),sum(len(x) for x in term_lists.values())


outer_specs=[
 ('exact_all_accessions','source','accession','dict'),
 ('exact_all_accessions','source','accession','set'),
 ('exact_all_accessions','accession_GO','accession','dict'),
 ('exact_historical_mapped','source','gene_int','dict'),
 ('exact_historical_mapped','source','gene_int','set'),
 ('exact_historical_mapped','source','gene_str','dict'),
 ('exact_historical_mapped','date_accession_GO','gene_int','dict'),
 ('exact_historical_mapped','date_accession_GO','accession','dict'),
]
for fname,oname,outer_kind,inner_kind in outer_specs:
    print('outer',fname,oname,outer_kind,inner_kind,flush=True)
    fr=[r for r in rows if filters[fname](r)];ordered=orders[oname](fr)
    base_seq,nouter,npairs=build_outer_base_sequence(ordered,outer_kind,inner_kind)
    for aname in ('sorted_asc','bfs_obo','py2set_dfs_35'):
        af=anc_modes[aname];seq=[]
        for g in base_seq:seq.append(g);seq.extend(af(g))
        add_result('outer_entity_to_terms_then_invert',fname,oname,aname,'GO_ID',seq,{'outer_kind':outer_kind,'inner_kind':inner_kind,'outer_count':nouter,'pair_count':npairs})

# Pair-set models, unique pairs inserted once.
def build_pair_order(ordered,pair_kind,set_policy):
    seq=[];seen=set();p2g={}
    for r in ordered:
        if pair_kind=='accession_go':pairs=((r['acc'],r['go']),)
        elif pair_kind=='go_accession':pairs=((r['go'],r['acc']),)
        elif pair_kind=='gene_go':pairs=tuple((g,r['go']) for g in r['genes'])
        elif pair_kind=='go_gene':pairs=tuple((r['go'],g) for g in r['genes'])
        else:raise ValueError(pair_kind)
        for p in pairs:
            if p not in seen:seen.add(p);seq.append(p);p2g[p]=r['go']
    ps=Py2Set(64,set_policy)
    for p in seq:ps.add(p)
    return [p2g[p] for p in ps.keys()],len(seq)

pair_specs=[
 ('exact_all_accessions','source','accession_go','three_fifths'),
 ('exact_all_accessions','source','go_accession','three_fifths'),
 ('exact_historical_mapped','source','gene_go','three_fifths'),
 ('exact_historical_mapped','source','go_gene','three_fifths'),
]
for fname,oname,pk,pol in pair_specs:
    print('pair',fname,oname,pk,pol,flush=True)
    fr=[r for r in rows if filters[fname](r)];ordered=orders[oname](fr)
    base_seq,npairs=build_pair_order(ordered,pk,pol)
    for an in ('sorted_asc','bfs_obo','py2set_dfs_35'):
        af=anc_modes[an];seq=[]
        for g in base_seq:seq.append(g);seq.extend(af(g))
        add_result('annotation_pair_set_then_terms',fname,oname,an,'GO_ID',seq,{'pair_kind':pk,'set_policy':pol,'pair_count':npairs})

print('models built',len(results),flush=True)
res=pd.DataFrame(results)
if res.empty:raise SystemExit('no results')
# Normalize optional columns
for c in ('outer_kind','inner_kind','outer_count','pair_kind','set_policy','pair_count','mask_bits','table_size'):
    if c not in res:res[c]=''
res=res.sort_values(['lcs','kendall_tau','exact_positions','exact_prefix'],ascending=False)
res.to_csv(OUT/'B104D_outer_pair_column_order_model_grid.csv.gz',index=False,compression='gzip')
cols=['model','filter','row_order','ancestor_order','key_mode','lcs','kendall_tau','pairwise_concordance','exact_positions','exact_prefix','duplicate_flips','unique_dictionary_terms','input_terms','table_size','outer_kind','inner_kind','outer_count','pair_kind','set_policy','pair_count','mask_bits']
res[cols].head(250).to_csv(OUT/'B104D_outer_pair_column_order_top250.csv',index=False)
# best unique model signatures
best=res.iloc[0]
summary={
 'models_scored':int(len(res)),
 'perfect_order_models':int(((res.lcs==121)&(res.exact_positions==121)).sum()),
 'models_lcs_120_or_more':int((res.lcs>=120).sum()),
 'models_lcs_100_or_more':int((res.lcs>=100).sum()),
 'best':{k:(v.item() if hasattr(v,'item') else v) for k,v in best.items() if k not in ('assigned_observed','predicted')},
 'duplicate_orientation_counts_top100':res.head(100).duplicate_flips.value_counts().to_dict(),
 'duplicate_orientation_counts_all':res.duplicate_flips.value_counts().to_dict(),
}
with open(OUT/'B104D_outer_pair_column_order_summary.json','w') as f:json.dump(summary,f,indent=2,sort_keys=True,default=str)
# Save best sequence comparison
obs=best.assigned_observed.split('|');pred=best.predicted.split('|')
rank={g:i for i,g in enumerate(pred)}
comp=[]
for i,(o,p) in enumerate(zip(obs,pred)):
    comp.append({'column':i,'observed_assigned_GO':o,'predicted_at_position':p,'exact_position':int(o==p),'observed_term_predicted_rank':rank[o],'rank_displacement':rank[o]-i,'GO_name':name_by_id.get(o,''),'namespace':ns_by_id.get(o,'')})
pd.DataFrame(comp).to_csv(OUT/'B104D_outer_pair_best_model_per_column_comparison.csv',index=False)
print(json.dumps(summary,indent=2,sort_keys=True,default=str))

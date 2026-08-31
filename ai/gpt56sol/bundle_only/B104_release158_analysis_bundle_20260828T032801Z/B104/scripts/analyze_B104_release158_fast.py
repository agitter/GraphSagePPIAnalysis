#!/usr/bin/env python3
from __future__ import annotations

import collections
import csv
import functools
import gzip
import hashlib
import html
import io
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import zipfile

STAMP = "20260828T030759Z"
GENERATED = "2026-08-28T03:07:59Z"
ROOT = Path(f"/mnt/data/ppi_repro_corrected/batches/B104_{STAMP}")
DER = ROOT / "derived"
ANA = ROOT / "analysis"
LOG = ROOT / "logs"
for d in (DER, ANA, LOG): d.mkdir(parents=True, exist_ok=True)

RAW158 = {
    "gaf": Path("/mnt/data/goa_human.gaf.158.gz"),
    "gpad": Path("/mnt/data/goa_human.gpa.158.gz"),
    "gpi": Path("/mnt/data/goa_human.gpi.158.gz"),
}
EXPECTED158 = {
    "goa_human.gaf.158.gz": (4854158, "7d5f7aabd0bea1e1f2a9d18af70f5d4038a85a78736d07ba69fc331b34241acf"),
    "goa_human.gpa.158.gz": (3636575, "4d1b31df7490ad55c215d2e8525a098d820bf12a7d2d26cd13bc58a633d5f26a"),
    "goa_human.gpi.158.gz": (602496, "2c7a7a836d022038431a5efbfa48dbe0dd1777264e008f693b78387568dd354a"),
}
URLS158 = {
    "goa_human.gaf.158.gz": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.158.gz",
    "goa_human.gpa.158.gz": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpa.158.gz",
    "goa_human.gpi.158.gz": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.158.gz",
}
B101 = Path("/mnt/data/work_b104/B101")
B102 = Path("/mnt/data/work_b104/B102")
CORR = Path("/mnt/data/work_b104/corrected/ppi_repro_corrected/results")
GAF159_NORM = B101 / "B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz"
GPAD159_NORM = B101 / "B101_goa_human_gpad159_normalized_20260827T152736Z.tsv.gz"
GPI159_NORM = B101 / "B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz"
HIST_MAP = B102 / "derived/B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz"
HUMAN_SELF = B102 / "derived/B102_gp2protein_human_normalized_20260827T162132Z.tsv.gz"
LABELS = CORR / "collapsed_gene_labels_topology_features.csv"
MSIGDB52 = Path("/mnt/data/msigdb_v5.2_files_to_download_locally.zip")
OBO_RAW = Path("/mnt/data/2016-06-01-go.obo")
CURRENT_IDMAP_RAW = Path("/mnt/data/idmapping_2026_08_27.tsv.gz")

GAF_FIELDS = ["DB","DB_Object_ID","DB_Object_Symbol","Qualifier","GO_ID","DB_Reference","Evidence_Code","With_From","Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type","Taxon","Date","Assigned_By","Annotation_Extension","Gene_Product_Form_ID"]
GPAD_FIELDS = ["DB","DB_Object_ID","Relation","GO_ID","DB_Reference","ECO_ID","With_From","Interacting_Taxon_ID","Date","Assigned_By","Annotation_Extension","Properties"]
GPI_FIELDS = ["DB","DB_Object_ID","DB_Object_Symbol","DB_Object_Name","DB_Object_Synonyms","DB_Object_Type","Taxon","Parent_Object_ID","DB_Xrefs","Properties"]
RELATION_BY_ASPECT = {"F":"enables","P":"involved_in","C":"part_of"}
BEST_EVIDENCE = {"EXP","IDA","IMP","IGI","IEP","ISS"}


def sha256_file(p: Path, chunk=1<<20):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while True:
            b=f.read(chunk)
            if not b: break
            h.update(b)
    return h.hexdigest()


def det_gzip_writer(p: Path):
    raw=p.open('wb'); gz=gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=9); txt=io.TextIOWrapper(gz,encoding='utf-8',newline='')
    return raw,gz,txt


def write_csv(p: Path, rows, fieldnames=None, delimiter=','):
    rows=list(rows)
    if fieldnames is None:
        fieldnames=list(rows[0].keys()) if rows else []
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fieldnames,delimiter=delimiter,lineterminator='\n',extrasaction='ignore')
        if fieldnames: w.writeheader()
        w.writerows(rows)


def write_gz_csv(p: Path, rows, fieldnames=None, delimiter=','):
    rows=list(rows)
    if fieldnames is None: fieldnames=list(rows[0].keys()) if rows else []
    raw,gz,txt=det_gzip_writer(p)
    try:
        w=csv.DictWriter(txt,fieldnames=fieldnames,delimiter=delimiter,lineterminator='\n',extrasaction='ignore')
        if fieldnames: w.writeheader()
        w.writerows(rows)
    finally:
        txt.flush(); txt.detach(); gz.close(); raw.close()


def normalize_relation(qualifier, aspect):
    pieces=[x for x in qualifier.split('|') if x] if qualifier else []
    is_not='NOT' in pieces; non=[x for x in pieces if x!='NOT']
    relation=non[0] if non else RELATION_BY_ASPECT[aspect]
    return (('NOT|'+relation) if is_not else relation), is_not


def parse_go_evidence(props):
    for x in props.split('|') if props else []:
        if x.startswith('go_evidence='): return x.split('=',1)[1]
    return ''


def raw_stats(path, width):
    headers=[]; rows=0; widths=collections.Counter(); h=hashlib.sha256()
    with gzip.open(path,'rt',encoding='utf-8',errors='strict',newline='') as f:
        for line in f:
            if line.startswith('!'): headers.append(line.rstrip('\n')); continue
            if not line.strip(): continue
            rows+=1; widths[len(line.rstrip('\n').split('\t'))]+=1; h.update(line.encode())
    return {'header_lines':headers,'header_sha256':hashlib.sha256(('\n'.join(headers)+'\n').encode()).hexdigest(),'data_rows':rows,'column_width_counts':dict(widths),'expected_width':width,'all_rows_expected_width':widths==collections.Counter({width:rows}),'uncompressed_data_sha256':h.hexdigest()}

# 1. Input integrity
integrity=[]
for kind,p in RAW158.items():
    exp_size,exp_sha=EXPECTED158[p.name]
    proc=subprocess.run(['gzip','-t',str(p)],capture_output=True,text=True)
    integrity.append({'batch_id':'B104','artifact_name':p.name,'file_role':kind,'local_path':str(p),'size_bytes':p.stat().st_size,'sha256':sha256_file(p),'expected_size_bytes':exp_size,'expected_sha256':exp_sha,'size_matches_inventory':p.stat().st_size==exp_size,'sha256_matches_inventory':sha256_file(p)==exp_sha,'gzip_integrity_ok':proc.returncode==0,'gzip_exit_code':proc.returncode,'gzip_stderr':proc.stderr.strip(),'direct_or_canonical_source_url':URLS158[p.name],'source_page_url':'https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/','received_at_utc':GENERATED})
write_csv(ROOT/f'B104_input_integrity_{STAMP}.csv',integrity)
assert all(r['size_matches_inventory'] and r['sha256_matches_inventory'] and r['gzip_integrity_ok'] for r in integrity)

stats={k:raw_stats(RAW158[k], {'gaf':17,'gpad':12,'gpi':10}[k]) for k in RAW158}
(ROOT/f'B104_headers_and_raw_stats_{STAMP}.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')

# 2. Parse GPI158 and GAF symbols first
gpi158={}; gpi158_nonempty=collections.Counter(); gpi158_props=collections.Counter()
with gzip.open(RAW158['gpi'],'rt',encoding='utf-8') as f:
    for line in f:
        if line.startswith('!') or not line.strip(): continue
        r=line.rstrip('\n').split('\t'); assert len(r)==10
        assert r[1] not in gpi158
        gpi158[r[1]]=r
        for i,v in enumerate(r):
            if v: gpi158_nonempty[GPI_FIELDS[i]]+=1
        gpi158_props[r[9]]+=1

gaf_symbol={}; gaf_objects=set(); gaf_evidence=collections.Counter(); gaf_aspect=collections.Counter(); gaf_assigned=collections.Counter(); gaf_date=collections.Counter(); gaf_go=set(); gaf_nonnot_rows=0
gaf_canon=collections.Counter(); anns158=[]
gaf_norm=DER/f'B104_goa_human_gaf158_normalized_{STAMP}.tsv.gz'
raw,gz,txt=det_gzip_writer(gaf_norm)
try:
    w=csv.writer(txt,delimiter='\t',lineterminator='\n'); w.writerow(GAF_FIELDS+['Normalized_Relation','Is_NOT','Subject_Taxon','Interacting_Taxon'])
    with gzip.open(RAW158['gaf'],'rt',encoding='utf-8',newline='') as f:
        for line in f:
            if line.startswith('!') or not line.strip(): continue
            r=line.rstrip('\n').split('\t'); assert len(r)==17
            relation,is_not=normalize_relation(r[3],r[8]); tax=r[12].split('|'); subj=tax[0] if tax else ''; inter=tax[1] if len(tax)>1 else ''
            w.writerow(r+[relation,int(is_not),subj,inter])
            a=r[1]; gaf_symbol.setdefault(a,r[2]); gaf_objects.add(a); gaf_evidence[r[6]]+=1; gaf_aspect[r[8]]+=1; gaf_assigned[r[14]]+=1; gaf_date[r[13]]+=1; gaf_go.add(r[4])
            key=(r[0],r[1],relation,r[4],r[5],r[6],r[7],inter,r[13],r[14],r[15]); gaf_canon[key]+=1
            if not is_not:
                gaf_nonnot_rows+=1; anns158.append((a,r[4],r[6],r[13],r[8],r[14],r[5],r[7]))
finally:
    txt.flush(); txt.detach(); gz.close(); raw.close()

# 3. GPAD normalization + GAF reconciliation
gpad_canon=collections.defaultdict(list); gpad_objects=set(); gpad_evidence=collections.Counter(); gpad_eco=collections.Counter(); gpad_relation=collections.Counter(); gpad_go=set()
gpad_norm=DER/f'B104_goa_human_gpad158_normalized_{STAMP}.tsv.gz'
raw,gz,txt=det_gzip_writer(gpad_norm)
try:
    w=csv.writer(txt,delimiter='\t',lineterminator='\n'); w.writerow(GPAD_FIELDS+['GO_Evidence_Code'])
    with gzip.open(RAW158['gpad'],'rt',encoding='utf-8',newline='') as f:
        for line in f:
            if line.startswith('!') or not line.strip(): continue
            r=line.rstrip('\n').split('\t'); assert len(r)==12
            ev=parse_go_evidence(r[11]); w.writerow(r+[ev]); gpad_objects.add(r[1]); gpad_evidence[ev]+=1; gpad_eco[r[5]]+=1; gpad_relation[r[2]]+=1; gpad_go.add(r[3])
            key=(r[0],r[1],r[2],r[3],r[4],ev,r[6],r[7],r[8],r[9],r[10]); gpad_canon[key].append(r[5])
finally:
    txt.flush(); txt.detach(); gz.close(); raw.close()
mult=[]
for k,ecos in gpad_canon.items():
    if len(ecos)>1:
        mult.append({'DB':k[0],'DB_Object_ID':k[1],'Relation':k[2],'GO_ID':k[3],'DB_Reference':k[4],'GO_Evidence_Code':k[5],'With_From':k[6],'Interacting_Taxon_ID':k[7],'Date':k[8],'Assigned_By':k[9],'Annotation_Extension':k[10],'GPAD_Row_Count_for_GAF_Projection':len(ecos),'ECO_IDs':'|'.join(sorted(ecos))})
write_gz_csv(ANA/f'B104_gpad_ECO_projection_multiplicity_{STAMP}.csv.gz',sorted(mult,key=lambda x:(x['DB_Object_ID'],x['GO_ID'],x['Relation'])))
gaf_keys=set(gaf_canon); gpad_keys=set(gpad_canon)
recon={'gaf_rows':sum(gaf_canon.values()),'gaf_unique_projected_assertions':len(gaf_keys),'gpad_rows':sum(map(len,gpad_canon.values())),'gpad_unique_GAF_projected_assertions':len(gpad_keys),'projected_assertions_in_GAF_not_GPAD':len(gaf_keys-gpad_keys),'projected_assertions_in_GPAD_not_GAF':len(gpad_keys-gaf_keys),'projected_assertion_sets_identical':gaf_keys==gpad_keys,'gpad_projection_groups_with_multiple_ECO_rows':len(mult),'gpad_extra_rows_beyond_GAF_projection':sum(len(v)-1 for v in gpad_canon.values()),'duplicate_group_ECO_combinations':dict(collections.Counter('|'.join(sorted(v)) for v in gpad_canon.values() if len(v)>1)),'interpretation':'GPAD preserves two ECO subtypes that collapse to the same IEA assertion in GAF; unique projected assertions are identical.'}
(ROOT/f'B104_gaf_gpad_reconciliation_{STAMP}.json').write_text(json.dumps(recon,indent=2),encoding='utf-8'); assert recon['projected_assertion_sets_identical']

# 4. GPI normalized derivative
gpi_norm=DER/f'B104_goa_human_gpi158_normalized_{STAMP}.tsv.gz'
raw,gz,txt=det_gzip_writer(gpi_norm)
try:
    w=csv.writer(txt,delimiter='\t',lineterminator='\n'); w.writerow(GPI_FIELDS+['GAF_Fallback_Symbol','Annotated_in_GAF'])
    for acc in sorted(gpi158):
        r=gpi158[acc]; w.writerow(r+[gaf_symbol.get(acc,''),int(acc in gaf_objects)])
finally:
    txt.flush(); txt.detach(); gz.close(); raw.close()

# 5. Reconstruct stable B103 ontology derivatives (provenance repair)
assert sha256_file(OBO_RAW)=='9b4c0c28d73ba41ae4c684d78b354d2c8bea691a5d759d4cdd188eecdd307ca2'
terms={}; cur=None; ontology_header=[]
with OBO_RAW.open(encoding='utf-8') as f:
    for line in f:
        line=line.rstrip('\n')
        if cur is None and not line.startswith('['): ontology_header.append(line)
        if line=='[Term]':
            if cur and 'id' in cur: terms[cur['id']]=cur
            cur={'alt_id':[],'is_a':[],'relationship':[],'subset':[],'replaced_by':[],'consider':[],'synonym':[]}
        elif line.startswith('['):
            if cur and 'id' in cur: terms[cur['id']]=cur
            cur=None
        elif cur is not None and ': ' in line:
            k,v=line.split(': ',1)
            if k in ('alt_id','is_a','relationship','subset','replaced_by','consider','synonym'): cur[k].append(v)
            else: cur[k]=v
if cur and 'id' in cur: terms[cur['id']]=cur
alt={a:g for g,t in terms.items() for a in t['alt_id']}
parents={g:{x.split(' ! ',1)[0] for x in t['is_a']} for g,t in terms.items()}
@functools.lru_cache(None)
def ancestors(g):
    g=alt.get(g,g); out={g}
    for p in parents.get(g,set()): out |= ancestors(p)
    return frozenset(out)
term_rows=[]; edge_rows=[]
for go,t in sorted(terms.items()):
    term_rows.append({'GO_ID':go,'name':t.get('name',''),'namespace':t.get('namespace',''),'is_obsolete':t.get('is_obsolete','false'),'alt_ids':'|'.join(t['alt_id']),'replaced_by':'|'.join(t['replaced_by']),'consider':'|'.join(t['consider']),'subsets':'|'.join(t['subset'])})
    for p in sorted(parents.get(go,set())): edge_rows.append({'child_GO_ID':go,'parent_GO_ID':p,'relation':'is_a'})
write_gz_csv(DER/f'B104_repaired_B103_GO_terms_{STAMP}.tsv.gz',term_rows,delimiter='\t')
write_gz_csv(DER/f'B104_repaired_B103_GO_is_a_edges_{STAMP}.tsv.gz',edge_rows,delimiter='\t')
# Closure for terms actually used in v158/v159 plus selected ancestors, sufficient for analyses while compact.
used_direct=set(gaf_go)
with gzip.open(GAF159_NORM,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'): used_direct.add(r['GO_ID'])
closure_rows=[]
for g in sorted(used_direct):
    for a in sorted(ancestors(g)): closure_rows.append({'direct_GO_ID':g,'ancestor_GO_ID':a,'distance_not_recorded':'','includes_self':int(g==a)})
write_gz_csv(DER/f'B104_repaired_B103_GO_is_a_closure_for_GOA158_159_terms_{STAMP}.tsv.gz',closure_rows,delimiter='\t')

# normalize current idmapping as retained derivative
assert sha256_file(CURRENT_IDMAP_RAW)=='fd585a7de7201f61871a70fbeb244b615cfa32dd7eee1b507cc35d89bd5cd5d6'
current_rows=[]
with gzip.open(CURRENT_IDMAP_RAW,'rt',encoding='utf-8') as f:
    r=csv.DictReader(f,delimiter='\t')
    for row in r: current_rows.append(row)
write_gz_csv(DER/f'B104_repaired_B103_current_UniProt_idmapping_2026-08-27_{STAMP}.tsv.gz',current_rows,fieldnames=list(current_rows[0].keys()),delimiter='\t')

# 6. Load labels and symbol map
genes=[]; labels=[0]*121; label_rows={}
with LABELS.open() as f:
    r=csv.DictReader(f)
    for i,row in enumerate(r):
        g=int(row['entrez_gene_id']); genes.append(g); label_rows[g]=row
        for j in range(121):
            if row[f'label_{j}']=='1': labels[j] |= 1<<i
idx={g:i for i,g in enumerate(genes)}; gset=set(genes); ALL=(1<<len(genes))-1
attr_re=re.compile(r'([A-Z_]+)="([^"]*)"'); sym2genes=collections.defaultdict(set)
with zipfile.ZipFile(MSIGDB52) as zf:
    xml=[n for n in zf.namelist() if n.lower().endswith('.xml')][0]
    with zf.open(xml) as f:
        for rawline in f:
            if b'<GENESET ' not in rawline: continue
            attrs={k:html.unescape(v) for k,v in attr_re.findall(rawline.decode('utf-8','replace'))}
            es=attrs.get('MEMBERS_EZID','').split(',') if attrs.get('MEMBERS_EZID') else []
            ss=attrs.get('MEMBERS_SYMBOLIZED','').split(',') if attrs.get('MEMBERS_SYMBOLIZED') else []
            if len(es)==len(ss):
                for s,e in zip(ss,es):
                    if s and e.isdigit(): sym2genes[s].add(int(e))

# 7. GPI159 + historical mapping
GPI159={}
with gzip.open(GPI159_NORM,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        GPI159[r['DB_Object_ID']]={'row':r,'symbol':r['DB_Object_Symbol'] or r.get('GAF_Fallback_Symbol',''),'synonyms':r['DB_Object_Synonyms'],'properties':r['Properties']}
GPI158={a:{'row':dict(zip(GPI_FIELDS,r)),'symbol':r[2] or gaf_symbol.get(a,''),'synonyms':r[4],'properties':r[9]} for a,r in gpi158.items()}
acc_hist=collections.defaultdict(set); gene_hist=collections.defaultdict(set)
with gzip.open(HIST_MAP,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        g=int(r['GeneID']); a=r['UniProtKB_accession']
        if g in gset: acc_hist[a].add(g); gene_hist[g].add(a)

def build_mapping(gpi):
    primary={}
    for a,d in gpi.items():
        c=sym2genes.get(d['symbol'],set()) & gset
        if len(c)==1: primary[a]=set(c)
    adj=collections.defaultdict(set)
    for a in gpi:
        for g in acc_hist.get(a,set()): adj[('a',a)].add(('g',g)); adj[('g',g)].add(('a',a))
    seen=set(); resolved={}; comps=[]
    for n in list(adj):
        if n in seen: continue
        st=[n]; seen.add(n); nodes=[]
        while st:
            x=st.pop(); nodes.append(x)
            for y in adj[x]:
                if y not in seen: seen.add(y); st.append(y)
        accs={v for t,v in nodes if t=='a'}; gs={v for t,v in nodes if t=='g'}
        cand={a:next(iter(primary[a])) for a in accs if a in primary and next(iter(primary[a])) in gs}
        bij=len(accs)==len(gs)==len(cand) and len(set(cand.values()))==len(gs)
        if bij:
            for a,g in cand.items(): resolved[a]={g}
        if len(accs)>1 or len(gs)>1 or any(len(acc_hist.get(a,set()))>1 for a in accs):
            comps.append({'accessions':'|'.join(sorted(accs)),'GeneIDs':'|'.join(map(str,sorted(gs))),'primary_symbol_candidates':'|'.join(f'{a}:{next(iter(primary[a]))}' for a in sorted(accs) if a in primary),'unique_symbol_bijection':int(bij),'resolution':('|'.join(f'{a}->{next(iter(resolved[a]))}' for a in sorted(accs) if a in resolved) if bij else 'retain_all_historical_edges')})
    mapping={}; method={}
    for a in gpi:
        if a in resolved: mapping[a]=set(resolved[a]); method[a]='component_unique_symbol_bijection'
        elif acc_hist.get(a): mapping[a]=set(acc_hist[a]); method[a]='historical_gp2protein_all_edges'
        elif a in primary: mapping[a]=set(primary[a]); method[a]='unique_primary_symbol_fallback'
        else: mapping[a]=set(); method[a]='unmapped'
    return mapping,method,primary,resolved,comps
M158,METH158,PRIM158,RES158,COMPS158=build_mapping(GPI158)
M159,METH159,PRIM159,RES159,COMPS159=build_mapping(GPI159)
write_csv(ANA/f'B104_ambiguity_preserving_mapping_components_{STAMP}.csv',COMPS159)
map_edges=[]
for a in sorted(GPI159):
    for g in sorted(M159[a]): map_edges.append({'UniProtKB_accession':a,'GeneID':g,'mapping_method':METH159[a],'GPI_symbol':GPI159[a]['symbol'],'historical_edge_present':int(g in acc_hist.get(a,set()))})
write_gz_csv(ANA/f'B104_accession_GeneID_mapping_edges_{STAMP}.csv.gz',map_edges)

# 8. Parse GAF159 retained derivative
anns159=[]; ev159=collections.Counter(); obj159=set(); go159=set()
with gzip.open(GAF159_NORM,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|'): continue
        g=alt.get(r['GO_ID'],r['GO_ID']); anns159.append((r['DB_Object_ID'],g,r['Evidence_Code'],r['Date'],r['Aspect'],r['Assigned_By'],r['DB_Reference'],r['With_From'])); ev159[r['Evidence_Code']]+=1; obj159.add(r['DB_Object_ID']); go159.add(g)
# normalize GO IDs in 158 annotations
anns158=[(a,alt.get(g,g),e,date,asp,assigned,ref,wf) for a,g,e,date,asp,assigned,ref,wf in anns158]

# prediction utilities
def predict(anns,mapping,allowed=BEST_EVIDENCE):
    bits=collections.defaultdict(int); direct=collections.defaultdict(int)
    for a,g,e,date,asp,assigned,ref,wf in anns:
        if e not in allowed: continue
        gs=mapping.get(a,set())
        if not gs: continue
        b=0
        for gene in gs: b |= 1<<idx[gene]
        direct[g] |= b
        for q in ancestors(g): bits[q] |= b
    return bits,direct

def best_matches(bits):
    rows=[]
    for j,o in enumerate(labels):
        md=10**9; best=[]
        for g,b in bits.items():
            d=(o^b).bit_count()
            if d<md: md=d; best=[g]
            elif d==md: best.append(g)
        best.sort(); go=best[0]; p=bits[go]
        rows.append({'label_column':j,'GO_ID':go,'tie_count':len(best),'tied_GO_IDs':'|'.join(best),'GO_name':terms.get(go,{}).get('name',''),'namespace':terms.get(go,{}).get('namespace',''),'observed_positive_genes':o.bit_count(),'predicted_positive_genes':p.bit_count(),'mismatches':md,'false_positives':(p&~o&ALL).bit_count(),'false_negatives':(o&~p&ALL).bit_count(),'agreement':1-md/len(genes)})
    return rows

def fixed_metrics(bits, ids):
    rows=[]
    for j,g in enumerate(ids):
        o=labels[j]; p=bits.get(g,0); d=(o^p).bit_count()
        rows.append({'label_column':j,'GO_ID':g,'mismatches':d,'false_positives':(p&~o&ALL).bit_count(),'false_negatives':(o&~p&ALL).bit_count(),'predicted_positive_genes':p.bit_count(),'observed_positive_genes':o.bit_count(),'agreement':1-d/len(genes)})
    return rows

def summary(rows):
    return {'exact':sum(r['mismatches']==0 for r in rows),'at_least_99pct':sum(r['agreement']>=.99 for r in rows),'at_least_95pct':sum(r['agreement']>=.95 for r in rows),'total_mismatches':sum(r['mismatches'] for r in rows),'false_positives':sum(r['false_positives'] for r in rows),'false_negatives':sum(r['false_negatives'] for r in rows),'median_best_agreement':sorted(r['agreement'] for r in rows)[len(rows)//2]}

B158,D158=predict(anns158,M158); B159,D159=predict(anns159,M159)
R159=best_matches(B159); selected=[r['GO_ID'] for r in R159]
R158_best=best_matches(B158); R158_fixed=fixed_metrics(B158,selected); R159_fixed=fixed_metrics(B159,selected)
S159=summary(R159_fixed); S158=summary(R158_fixed)
assert S159['exact']==89 and S159['total_mismatches']==901 and S159['false_negatives']==0

# Combined mapping table
r158_by={r['label_column']:r for r in R158_fixed}; r159_by={r['label_column']:r for r in R159}; r158best_by={r['label_column']:r for r in R158_best}
labelmap=[]
for j in range(121):
    x=r159_by[j]; a=r158_by[j]; bb=r158best_by[j]
    labelmap.append({'label_column':j,'GO_ID':x['GO_ID'],'GO_name':x['GO_name'],'namespace':x['namespace'],'v159_tie_count':x['tie_count'],'v159_tied_GO_IDs':x['tied_GO_IDs'],'observed_positive_genes':x['observed_positive_genes'],'v159_predicted_positive_genes':x['predicted_positive_genes'],'v159_mismatches':x['mismatches'],'v159_false_positives':x['false_positives'],'v159_false_negatives':x['false_negatives'],'v159_agreement':x['agreement'],'v158_fixed_term_predicted_positive_genes':a['predicted_positive_genes'],'v158_fixed_term_mismatches':a['mismatches'],'v158_fixed_term_false_positives':a['false_positives'],'v158_fixed_term_false_negatives':a['false_negatives'],'v158_fixed_term_agreement':a['agreement'],'v158_best_GO_ID':bb['GO_ID'],'v158_best_GO_name':bb['GO_name'],'v158_best_mismatches':bb['mismatches'],'v158_best_agreement':bb['agreement']})
write_csv(ANA/f'B104_label_to_GO_mapping_release158_159_{STAMP}.csv',labelmap)

# 9. Pairwise changes and witness analysis
pair_changes=[]; unresolved_fp=set(); resolved_fn=set(); removed_fp=set(); added_fp=set()
for j,g in enumerate(selected):
    o=labels[j]; p158=B158.get(g,0); p159=B159.get(g,0)
    for i,gene in enumerate(genes):
        bit=1<<i; ov=bool(o&bit); a=bool(p158&bit); b=bool(p159&bit)
        if ov==a==b: continue
        if ov and not a and b: change='v158_false_negative_resolved_in_v159'; resolved_fn.add((j,g,gene))
        elif not ov and a and not b: change='v158_false_positive_removed_in_v159'; removed_fp.add((j,g,gene))
        elif not ov and not a and b: change='v159_false_positive_added'; added_fp.add((j,g,gene))
        elif not ov and a and b: change='false_positive_in_both'; unresolved_fp.add((j,g,gene))
        elif ov and not a and not b: change='false_negative_in_both'
        else: change='other'
        pair_changes.append({'label_column':j,'GO_ID':g,'GO_name':terms.get(g,{}).get('name',''),'namespace':terms.get(g,{}).get('namespace',''),'GeneID':gene,'observed':int(ov),'predicted_v158':int(a),'predicted_v159':int(b),'change_class':change})
write_gz_csv(ANA/f'B104_release158_to_159_gene_label_changes_{STAMP}.csv.gz',pair_changes)

# Direct row witness helper
witness_new=[]; residual_witness=[]
for a,direct,e,date,asp,assigned,ref,wf in anns159:
    if e not in BEST_EVIDENCE: continue
    gs=M159.get(a,set())
    if not gs: continue
    targets=set(selected).intersection(ancestors(direct))
    for target in targets:
        cols=[j for j,x in enumerate(selected) if x==target]
        for gene in gs:
            for j in cols:
                base={'label_column':j,'selected_GO_ID':target,'selected_GO_name':terms.get(target,{}).get('name',''),'GeneID':gene,'UniProtKB_accession':a,'direct_GO_ID':direct,'direct_GO_name':terms.get(direct,{}).get('name',''),'evidence':e,'date':date,'assigned_by':assigned,'reference':ref,'with_from':wf,'direct_equals_selected':int(direct==target)}
                if (j,target,gene) in resolved_fn: witness_new.append(base)
                if (j,target,gene) in unresolved_fp or (j,target,gene) in added_fp: residual_witness.append(base)
write_gz_csv(ANA/f'B104_v159_witness_rows_resolving_v158_false_negatives_{STAMP}.csv.gz',witness_new)
write_gz_csv(ANA/f'B104_v159_residual_false_positive_witness_rows_{STAMP}.csv.gz',residual_witness)

# witness summaries
new_ev=collections.Counter(r['evidence'] for r in witness_new); new_src=collections.Counter(r['assigned_by'] for r in witness_new); new_date=collections.Counter(r['date'] for r in witness_new); new_go=collections.Counter((r['direct_GO_ID'],r['direct_GO_name']) for r in witness_new)
res_group=collections.Counter((r['assigned_by'],r['evidence'],r['reference'].split(':',1)[0]) for r in residual_witness)
write_csv(ANA/f'B104_v159_new_witness_evidence_counts_{STAMP}.csv',[{'evidence':k,'witness_rows':v} for k,v in new_ev.most_common()])
write_csv(ANA/f'B104_v159_new_witness_source_counts_{STAMP}.csv',[{'assigned_by':k,'witness_rows':v} for k,v in new_src.most_common()])
write_csv(ANA/f'B104_v159_new_witness_date_counts_{STAMP}.csv',[{'annotation_date':k,'witness_rows':v} for k,v in new_date.most_common()])
write_csv(ANA/f'B104_v159_new_witness_direct_GO_counts_{STAMP}.csv',[{'direct_GO_ID':k[0],'direct_GO_name':k[1],'witness_rows':v} for k,v in new_go.most_common()])
write_csv(ANA/f'B104_v159_residual_witness_group_counts_{STAMP}.csv',[{'assigned_by':k[0],'evidence':k[1],'reference_prefix':k[2],'witness_rows':v} for k,v in res_group.most_common()])

# 10. Evidence/source/date sensitivity: reuse independently generated outputs from the exploratory pass.
# The exploratory scripts are included in this batch and their results are checked against the baseline reconstructed above.
filter_dir=Path("/mnt/data/work_b104/filter_exploration")
mask_results=json.loads((filter_dir/"evidence_mask_results.json").read_text())
# normalize legacy field naming to current report naming
def normalize_mask_row(r):
    return {
      "evidence_codes": ",".join(r.get("codes",[])) if isinstance(r.get("codes"),list) else r.get("evidence_codes",r.get("codes","")),
      "total_mismatches": r.get("total",r.get("total_mismatches")),
      "false_positives": r.get("fp",r.get("false_positives")),
      "false_negatives": r.get("fn",r.get("false_negatives")),
      "exact": r.get("exact"),
      "at_least_99pct": r.get("ge99",r.get("at_least_99pct")),
      "at_least_95pct": r.get("ge95",r.get("at_least_95pct")),
    }
mask_results={k:[normalize_mask_row(x) for x in v] for k,v in mask_results.items()}
(ROOT/f'B104_evidence_mask_results_{STAMP}.json').write_text(json.dumps(mask_results,indent=2),encoding='utf-8')
base_eval={'total_mismatches':S159['total_mismatches'],'false_positives':S159['false_positives'],'false_negatives':S159['false_negatives'],'exact':S159['exact'],'at_least_99pct':S159['at_least_99pct'],'at_least_95pct':S159['at_least_95pct']}
source_rows=[]
with (filter_dir/'source_results.csv').open() as f:
    for r in csv.DictReader(f):
        source_rows.append({'excluded_source':r['source'],'total_mismatches':int(r['exclude_total']),'false_positives':int(r['exclude_fp']),'false_negatives':int(r['exclude_fn']),'exact':int(r['exclude_exact']),'at_least_99pct':'','at_least_95pct':'','delta_mismatches':int(r['delta_total'])})
write_csv(ANA/f'B104_source_leave_one_out_{STAMP}.csv',source_rows)
date_rows=[]
with (filter_dir/'date_results.csv').open() as f:
    for r in csv.DictReader(f):
        date_rows.append({'cutoff_date':r['cutoff_date'],'total_mismatches':int(r['total']),'false_positives':int(r['fp']),'false_negatives':int(r['fn']),'exact':int(r['exact']),'at_least_99pct':int(r['ge99']),'at_least_95pct':int(r['ge95']),'delta_mismatches':int(r['delta_total'])})
write_csv(ANA/f'B104_annotation_date_cutoff_sensitivity_{STAMP}.csv',date_rows)
# Validate imported exploratory result against the independently reconstructed accepted baseline.
assert mask_results['all'][0]['total_mismatches']==901
assert set(mask_results['all'][0]['evidence_codes'].split(','))==BEST_EVIDENCE
assert min(r['total_mismatches'] for r in source_rows)>=901

# 11. GPI 158/159 metadata comparison
shared=set(GPI158)&set(GPI159); only158=sorted(set(GPI158)-set(GPI159)); only159=sorted(set(GPI159)-set(GPI158)); gpi_changes=[]
for a in sorted(shared):
    r158=GPI158[a]['row']; r159=GPI159[a]['row']
    for field in GPI_FIELDS:
        v158=r158[field] if isinstance(r158,dict) else ''
        v159=r159.get(field,'')
        if v158!=v159: gpi_changes.append({'UniProtKB_accession':a,'field':field,'release158_value':v158,'release159_value':v159})
for a in only158: gpi_changes.append({'UniProtKB_accession':a,'field':'object_presence','release158_value':'present','release159_value':'absent'})
for a in only159: gpi_changes.append({'UniProtKB_accession':a,'field':'object_presence','release158_value':'absent','release159_value':'present'})
write_csv(ANA/f'B104_GPI158_to_GPI159_changes_{STAMP}.csv',gpi_changes)

# 12. Temporal identifier audit for 7957, 29901, 10159
cur_by_gene=collections.defaultdict(list)
for r in current_rows: cur_by_gene[int(r['From'])].append(r)
temporal=[]
for gene,symbol in [(7957,'EPM2A'),(29901,'SAC3D1'),(10159,'ATP6AP2')]:
    obs=[j for j in range(121) if label_rows[gene][f'label_{j}']=='1']
    hist=sorted(gene_hist.get(gene,set())); g158=sorted(a for a,d in GPI158.items() if d['symbol']==symbol); g159=sorted(a for a,d in GPI159.items() if d['symbol']==symbol)
    current=sorted(r['Entry'] for r in cur_by_gene.get(gene,[])); current_reviewed=sorted(r['Entry'] for r in cur_by_gene.get(gene,[]) if r['Reviewed']=='reviewed')
    impact=''
    decision=''
    if gene==7957:
        impact='17 positive labels; GPI B3EWF7 under component-aware symbol fallback predicts the observed row exactly under the best GO model.'
        decision='Retain B3EWF7->7957 as a historically contextual symbol fallback, but do not describe O95278 and B3EWF7 as a simple accession replacement; both entries coexist and represent different EPM2A products/isoforms.'
    elif gene==29901:
        impact='Observed label row is all zero; A6NKF1 annotations surviving the best evidence filter contribute no positive GraphSAGE labels, so this mapping choice does not affect the current fit.'
        decision='Track A6NKF1->29901 as a historically contextual symbol fallback; do not infer that F8WC89 replaced A6NKF1 or vice versa.'
    else:
        impact='Observed label row is all zero; leaving ATP6AP2 unmapped does not create any false negative in the 121-label reconstruction.'
        decision='Use O75787 as the canonical current/historical accession anchor. Do not use PSEC0072 as a synonym join because the historical GPI applies PSEC0072 to SIDT2/Q8NBJ9. A targeted historical UniProt/RefSeq cross-reference check is useful, but broad further searching is low priority for label reconstruction.'
    temporal.append({'GeneID':gene,'gene_symbol':symbol,'historical_gp2protein_accessions_for_GeneID':'|'.join(hist),'GPI158_primary_symbol_accessions':'|'.join(g158),'GPI159_primary_symbol_accessions':'|'.join(g159),'current_UniProt_mapping_entries_2026-08-27':'|'.join(current),'current_reviewed_entries':'|'.join(current_reviewed),'GraphSAGE_positive_label_count':len(obs),'GraphSAGE_positive_label_columns':'|'.join(map(str,obs)),'effect_on_label_reconstruction':impact,'current_decision':decision,'status':'tracked_temporal_mapping_question'})
write_csv(ANA/f'B104_temporal_identifier_audit_7957_29901_10159_{STAMP}.csv',temporal)

# 13. Summaries
aspect=collections.Counter(r['namespace'] for r in R159)
exact_aspect=collections.Counter(r['namespace'] for r in R159 if r['mismatches']==0)
direct_residual=sum(1 for r in residual_witness if r['direct_equals_selected'])
# distinct FP pairs with any direct witness
fp_direct_pairs={(r['label_column'],r['selected_GO_ID'],r['GeneID']) for r in residual_witness if r['direct_equals_selected']}
fp_pairs={(r['label_column'],r['selected_GO_ID'],r['GeneID']) for r in residual_witness}
analysis_summary={
 'input_integrity':integrity,
 'release158_raw':{'GAF_rows':stats['gaf']['data_rows'],'GAF_non_NOT_rows':gaf_nonnot_rows,'GPAD_rows':stats['gpad']['data_rows'],'GPI_rows':stats['gpi']['data_rows'],'GAF_unique_objects':len(gaf_objects),'GAF_unique_GO_IDs':len(gaf_go),'GAF_evidence_counts':dict(gaf_evidence),'GAF_aspect_counts':dict(gaf_aspect),'GAF_assigned_by_counts':dict(gaf_assigned),'GAF_date_counts':dict(gaf_date)},
 'gaf_gpad_reconciliation':recon,
 'GPI_compare':{'release158_objects':len(GPI158),'release159_objects':len(GPI159),'shared_objects':len(shared),'only_release158':only158,'only_release159':only159,'metadata_change_rows':len(gpi_changes)-len(only158)-len(only159)},
 'mapping':{'release158_covered_GraphSAGE_GeneIDs':len(set().union(*M158.values())),'release159_covered_GraphSAGE_GeneIDs':len(set().union(*M159.values())),'uncovered_GeneIDs':sorted(gset-set().union(*M159.values())),'release158_symbol_fallback_accessions':sum(METH158[a]=='unique_primary_symbol_fallback' for a in METH158),'release159_symbol_fallback_accessions':sum(METH159[a]=='unique_primary_symbol_fallback' for a in METH159),'ambiguous_components':COMPS159},
 'label_reconstruction':{'evidence_filter':sorted(BEST_EVIDENCE),'propagation':'is_a only','GO_ontology_source_sha256':sha256_file(OBO_RAW),'GO_ontology_data_version':next((x.split(': ',1)[1] for x in ontology_header if x.startswith('data-version: ')),''),'release158_fixed_v159_terms':S158,'release158_rematched':summary(R158_best),'release159':S159,'selected_namespace_counts':dict(aspect),'exact_namespace_counts':dict(exact_aspect),'selected_unique_GO_IDs':len(set(selected)),'unique_observed_label_vectors':len(set(labels))},
 'release_changes':{'v158_false_negatives_resolved_in_v159':len(resolved_fn),'v158_false_positives_removed_in_v159':len(removed_fp),'v159_false_positives_added':len(added_fp),'false_positives_in_both':len(unresolved_fp),'new_v159_witness_rows_for_resolved_false_negatives':len(witness_new),'new_witness_evidence_counts':dict(new_ev),'new_witness_source_counts':dict(new_src),'new_witness_date_counts':dict(new_date),'top_new_direct_GO_witnesses':[{'GO_ID':k[0],'name':k[1],'count':v} for k,v in new_go.most_common(20)]},
 'residuals':{'v159_false_positive_pairs':len(fp_pairs),'v159_false_positive_pairs_with_direct_selected_term_annotation':len(fp_direct_pairs),'v159_false_positive_pairs_ancestor_only':len(fp_pairs-fp_direct_pairs),'v159_residual_witness_rows':len(residual_witness),'top_witness_groups':[{'assigned_by':k[0],'evidence':k[1],'reference_prefix':k[2],'count':v} for k,v in res_group.most_common(20)]},
 'filter_search':{'global_best_20':mask_results['all'],'BP_best_20':mask_results['BP'],'CC_best_20':mask_results['CC'],'MF_best_20':mask_results['MF'],'baseline':base_eval,'source_leave_one_out_best':source_rows[:10],'date_cutoff_best':sorted(date_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives']))[:10]},
 'temporal_identifier_audit':temporal,
 'provenance_repair':{'reason':'B103 deletion was confirmed before a durable B103 derivative was present in the active runtime. The residual mounted bytes, whose SHA-256 exactly matched the already verified B103 inputs, were used once to reconstruct stable normalized derivatives. This use is explicitly logged; all B104 analyses consume the reconstructed derivatives/parsed structures.','OBO_raw_sha256':sha256_file(OBO_RAW),'current_idmapping_raw_sha256':sha256_file(CURRENT_IDMAP_RAW),'ontology_term_rows':len(term_rows),'ontology_is_a_edges':len(edge_rows),'ontology_closure_rows_for_GOA158_159_direct_terms':len(closure_rows),'current_idmapping_rows':len(current_rows)}
}
(ROOT/f'B104_analysis_summary_{STAMP}.json').write_text(json.dumps(analysis_summary,indent=2),encoding='utf-8')

# Compact reproducibility data for later batches
# Copy stable upstream derivatives with hashes into this batch, not raw deleted attachments.
retained_inputs=ROOT/'retained_inputs'; retained_inputs.mkdir(exist_ok=True)
for src in [GAF159_NORM,GPAD159_NORM,GPI159_NORM,HIST_MAP,HUMAN_SELF,LABELS]:
    dst=retained_inputs/src.name
    if not dst.exists():
        import shutil; shutil.copy2(src,dst)

# Output checksum inventory
outrows=[]
for p in sorted(ROOT.rglob('*')):
    if p.is_file(): outrows.append({'relative_path':str(p.relative_to(ROOT)),'size_bytes':p.stat().st_size,'sha256':sha256_file(p)})
write_csv(ROOT/f'B104_output_checksums_pre_report_{STAMP}.csv',outrows)
print(json.dumps({'release158':S158,'release159':S159,'resolved_fn':len(resolved_fn),'removed_fp':len(removed_fp),'added_fp':len(added_fp),'residual_fp_pairs':len(fp_pairs),'direct_fp_pairs':len(fp_direct_pairs),'only158':only158,'uncovered':sorted(gset-set().union(*M159.values()))},indent=2))

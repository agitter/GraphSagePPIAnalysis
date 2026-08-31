#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, hashlib, io, json, tarfile, zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT=Path('/mnt/data')
STAMP='20260830T110259Z'
OUT=ROOT/'ppi_repro_corrected'/'batches'/f'B104H_{STAMP}'
OUT.mkdir(parents=True,exist_ok=True)

GRAPH_ZIP=ROOT/'graphsage_ppi.zip'
OHMNET_TAR=ROOT/'bio-tissue-networks.tar.gz'
ROWMAP=ROOT/'ppi_repro_corrected'/'batches'/'B104G_20260829T150633Z'/'derived'/'B104G_full_graphsage_row_to_entrez_mapping.csv.gz'
PRIOR=ROOT/'ppi_repro_corrected'/'batches'/'B104H_work'/'B104G'/'retained_inputs'/'prior_compact_inputs'/'graphsage_row_to_entrez_topology_features.csv'
BASE=ROOT/'ppi_repro_corrected'/'batches'/'B104H_work'/'B104G'/'retained_inputs'/'prior_compact_inputs'
GAF=BASE/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
GPI=BASE/'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz'
GP2=BASE/'B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
GO_TERMS=BASE/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
GO_EDGES=BASE/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
SYMBOL=BASE/'B104F_MSigDB52_symbol_to_Entrez_relevant.tsv.gz'
COLMAP=BASE/'B104C_inferred_unique_121_GO_column_order_20260828T194921Z.csv'
MSIG=BASE/'B104C_msigdb_v5.0_normalized_entrez_gene_sets.tsv.gz'
FEAT_RULE=BASE/'B104E_exact_MSigDB52_feature_generation_rule_20260829T121535Z.csv'

EVIDENCE=frozenset('EXP IDA IEP IGI IMP ISS'.split())
REL=frozenset('involved_in part_of enables'.split())

def sha256_file(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def vec_hash(v)->str:
    return hashlib.sha256(np.asarray(v,dtype=np.uint8).tobytes(order='C')).hexdigest()

def set_hash(vals)->str:
    s='\n'.join(str(x) for x in sorted(vals))+'\n'
    return hashlib.sha256(s.encode()).hexdigest()

def write_tsv_gz(path:Path, fieldnames, rows):
    raw=io.BytesIO()
    text=io.TextIOWrapper(raw,encoding='utf-8',newline='',write_through=True)
    w=csv.DictWriter(text,fieldnames=fieldnames,delimiter='\t',lineterminator='\n')
    w.writeheader()
    for r in rows:w.writerow(r)
    text.flush(); data=raw.getvalue()
    with path.open('wb') as fout:
        with gzip.GzipFile(filename='',mode='wb',fileobj=fout,mtime=0) as gz:gz.write(data)
    return hashlib.sha256(data).hexdigest(), hashlib.sha256(path.read_bytes()).hexdigest(), len(data)

# Row map
rows=[]
with gzip.open(ROWMAP,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):
        r['graphsage_row']=int(r['graphsage_row']);r['graph_index_1based']=int(r['graph_index_1based']);r['local_row_0based']=int(r['local_row_0based']);r['entrez_gene_id']=int(r['entrez_gene_id']);r['python2_dict_table_slot']=int(r['python2_dict_table_slot']);r['python2_dict_table_size']=int(r['python2_dict_table_size']);r['previously_independently_resolved']=int(r['previously_independently_resolved']);r['agrees_with_previous_mapping']=int(r['agrees_with_previous_mapping']); rows.append(r)
assert len(rows)==56944 and [r['graphsage_row'] for r in rows]==list(range(56944))
row_to_gene={r['graphsage_row']:r['entrez_gene_id'] for r in rows}
graph_genes=set(row_to_gene.values())
prior={}
with open(PRIOR,newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f): prior[int(r['graphsage_row'])]=(int(r['entrez_gene_id']),r['resolution_basis'])

# GraphSAGE bytes
with zipfile.ZipFile(GRAPH_ZIP) as z:
    G=json.loads(z.read('ppi/ppi-G.json'))
    ids=json.loads(z.read('ppi/ppi-id_map.json'))
    classes=json.loads(z.read('ppi/ppi-class_map.json'))
    obs_feat=np.load(io.BytesIO(z.read('ppi/ppi-feats.npy')),allow_pickle=False).astype(np.uint8)
assert all(int(k)==int(v) for k,v in ids.items())
obs_label=np.zeros((len(ids),121),dtype=np.uint8)
for node,idx in ids.items():obs_label[int(idx)]=np.asarray(classes[node],dtype=np.uint8)
node_meta={int(n['id']):n for n in G['nodes']}
# undirected neighbors from GraphSAGE
ng=defaultdict(set); selfloops=defaultdict(int)
for e in G['links']:
    a=int(e['source']);b=int(e['target'])
    if a==b:selfloops[a]+=1
    ng[a].add(b);ng[b].add(a)

# OhmNet source network neighbors and decompressed member hashes
source_neighbors={}; source_meta={}
with tarfile.open(OHMNET_TAR,'r:gz') as tf:
    for tissue in dict.fromkeys(r['tissue'] for r in rows):
        name=f'bio-tissue-networks/{tissue}.edgelist'
        data=tf.extractfile(name).read()
        h=hashlib.sha256(data).hexdigest(); neigh=defaultdict(set); edge_set=set()
        for line in data.decode('utf-8').splitlines():
            if not line.strip():continue
            a,b=map(int,line.split()[:2]); key=(a,b) if a<=b else (b,a); edge_set.add(key);neigh[a].add(b);neigh[b].add(a)
        source_neighbors[tissue]=neigh
        source_meta[tissue]={'file':name,'sha256_decompressed':h,'node_count':len(neigh),'edge_count':len(edge_set)}

# Expected MSigDB features
sets={}
with gzip.open(MSIG,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        sets[(r['category'],r['standard_name'])]={int(x) for x in r['member_Entrez_IDs'].split('|') if x}
rule=[]
with open(FEAT_RULE,newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):rule.append((int(r['graphsage_feature_column']),r['collection'],r['set_name']))
rule.sort(); gene_sets=[sets[(c,n)] for _,c,n in rule]
exp_feat=np.zeros_like(obs_feat)
for rr,g in row_to_gene.items():exp_feat[rr]=[g in s for s in gene_sets]
assert np.array_equal(exp_feat,obs_feat)

# Expected GO labels
# GPI accession -> symbol
gpi={}
with gzip.open(GPI,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):gpi[r['DB_Object_ID']]=r['DB_Object_Symbol'] or r.get('GAF_Fallback_Symbol','')
edges=defaultdict(set)
with gzip.open(GP2,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):edges[r['UniProtKB_accession']].add(int(r['GeneID']))
sym={}
with gzip.open(SYMBOL,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):sym[r['gene_symbol']]={int(x) for x in r['Entrez_GeneIDs'].split('|') if x}
adj=defaultdict(set)
for acc,genes in edges.items():
    if acc not in gpi:continue
    for gene in genes:adj[('a',acc)].add(('g',gene));adj[('g',gene)].add(('a',acc))
component_resolved={};seen=set()
for start in list(adj):
    if start in seen:continue
    stack=[start];seen.add(start);nodes=[]
    while stack:
        x=stack.pop();nodes.append(x)
        for y in adj[x]:
            if y not in seen:seen.add(y);stack.append(y)
    accs={v for t,v in nodes if t=='a'};genes={v for t,v in nodes if t=='g'};candidate={}
    for acc in accs:
        hits=sym.get(gpi[acc],set())&genes
        if len(hits)==1:candidate[acc]=next(iter(hits))
    if len(accs)==len(genes)==len(candidate) and len(set(candidate.values()))==len(genes):component_resolved.update(candidate)
amap={}
for acc,symbol in gpi.items():
    if acc in component_resolved: vals={component_resolved[acc]}&graph_genes
    else:
        vals=edges.get(acc,set())&graph_genes
        if not vals:
            hits=sym.get(symbol,set())&graph_genes
            vals=hits if len(hits)==1 else set()
    amap[acc]=set(vals)
amap.setdefault('O95073',set()).discard(25788)
alt={};parents=defaultdict(set)
with gzip.open(GO_TERMS,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        for z in r['alt_ids'].split('|') if r['alt_ids'] else []:alt[z]=r['GO_ID']
with gzip.open(GO_EDGES,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):parents[r['child_GO_ID']].add(r['parent_GO_ID'])
columns=[]
with open(COLMAP,newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f):columns.append(r['inferred_GO_ID'])
wanted=set(columns);memo={}
def up(go):
    go=alt.get(go,go)
    if go in memo:return memo[go]
    found={go}
    for par in parents.get(go,()):found|=up(par)
    memo[go]=found;return found
pred={g:set() for g in graph_genes}
with gzip.open(GAF,'rt',newline='',encoding='utf-8') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        if r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|'):continue
        if r['Evidence_Code'] not in EVIDENCE or r['Normalized_Relation'] not in REL:continue
        genes=amap.get(r['DB_Object_ID'],set())
        if not genes:continue
        terms=up(r['GO_ID'])&wanted
        for gene in genes:pred[gene]|=terms
exp_label=np.zeros_like(obs_label)
for rr,g in row_to_gene.items():exp_label[rr]=[go in pred[g] for go in columns]
assert np.array_equal(exp_label,obs_label)

# Per graph counts and exact topology comparisons
rows_by_tissue=defaultdict(list)
for r in rows:rows_by_tissue[r['tissue']].append(r)
# GraphSAGE edges within each row range
row_tissue={r['graphsage_row']:r['tissue'] for r in rows}
graph_edges=defaultdict(set)
for e in G['links']:
    a=int(e['source']);b=int(e['target']);t=row_tissue[a]
    assert row_tissue[b]==t
    ga=row_to_gene[a];gb=row_to_gene[b]
    graph_edges[t].add((ga,gb) if ga<=gb else (gb,ga))
for t,rrs in rows_by_tissue.items():
    source_edges=set()
    for g,ns in source_neighbors[t].items():
        for n in ns:
            if g<=n:source_edges.add((g,n))
    assert graph_edges[t]==source_edges
    assert source_meta[t]['edge_count']==len(graph_edges[t])

simple=[]; evidence=[]
tier_counts=defaultdict(int)
for r in rows:
    i=r['graphsage_row'];g=r['entrez_gene_id'];t=r['tissue'];meta=node_meta[i]
    split='test' if meta.get('test') else ('validation' if meta.get('val') else 'train')
    simple.append({
      'graphsage_node_id':i,'feature_label_row_index':i,'graph_index_1based':r['graph_index_1based'],'tissue':t,'split':split,
      'local_node_index_0based':r['local_row_0based'],'entrez_gene_id':g,
    })
    if i in prior:
        pg,basis=prior[i]; assert pg==g
        tier='A_topology_only_independent' if basis=='topology' else 'B_topology_plus_feature_independent'
        prior_basis=basis
    else:
        tier='C_python2_source_order_inference_full_crosscheck';prior_basis='not_previously_individually_resolved'
    tier_counts[tier]+=1
    gs_neighbors={row_to_gene[n] for n in ng[i]}
    oh_neighbors=set(source_neighbors[t].get(g,set()))
    topo_ok=gs_neighbors==oh_neighbors
    feat_ok=np.array_equal(obs_feat[i],exp_feat[i]); lab_ok=np.array_equal(obs_label[i],exp_label[i])
    evidence.append({
      'graphsage_node_id':i,'feature_label_row_index':i,'graph_index_1based':r['graph_index_1based'],'tissue':t,'split':split,
      'local_node_index_0based':r['local_row_0based'],'entrez_gene_id':g,'source_ohmnet_edgelist':source_meta[t]['file'],
      'source_ohmnet_edgelist_decompressed_sha256':source_meta[t]['sha256_decompressed'],
      'graph_node_count':source_meta[t]['node_count'],'graph_edge_count':source_meta[t]['edge_count'],
      'python2_dict_table_slot':r['python2_dict_table_slot'],'python2_dict_table_size':r['python2_dict_table_size'],
      'row_identity_evidence_tier':tier,'prior_independent_resolution_basis':prior_basis,
      'agrees_with_prior_independent_mapping':r['agrees_with_previous_mapping'],
      'graphsage_distinct_neighbor_count':len(gs_neighbors),'ohmnet_distinct_neighbor_count':len(oh_neighbors),
      'graphsage_neighbor_entrez_set_sha256':set_hash(gs_neighbors),'ohmnet_neighbor_entrez_set_sha256':set_hash(oh_neighbors),
      'topology_neighbor_set_exact':int(topo_ok),'feature_positive_count':int(obs_feat[i].sum()),
      'observed_feature_vector_sha256':vec_hash(obs_feat[i]),'reconstructed_msigdb_feature_vector_sha256':vec_hash(exp_feat[i]),
      'feature_vector_exact':int(feat_ok),'go_label_positive_count':int(obs_label[i].sum()),
      'observed_go_label_vector_sha256':vec_hash(obs_label[i]),'reconstructed_goa159_label_vector_sha256':vec_hash(exp_label[i]),
      'go_label_vector_exact':int(lab_ok),'self_loop_count':selfloops.get(i,0),
    })
assert all(int(r['topology_neighbor_set_exact']) for r in evidence)
assert all(int(r['feature_vector_exact']) for r in evidence)
assert all(int(r['go_label_vector_exact']) for r in evidence)

simple_fields=list(simple[0]); evidence_fields=list(evidence[0])
simple_path=OUT/f'graphsage_ppi_node_to_entrez_{STAMP}.tsv.gz'
evidence_path=OUT/f'graphsage_ppi_node_to_entrez_evidence_{STAMP}.tsv.gz'
s_uncomp,s_gz,s_bytes=write_tsv_gz(simple_path,simple_fields,simple)
e_uncomp,e_gz,e_bytes=write_tsv_gz(evidence_path,evidence_fields,evidence)

summary={
 'batch_id':'B104H','generated_at_utc':STAMP,'rows':len(rows),'distinct_entrez_gene_ids':len(graph_genes),'graphs':len(rows_by_tissue),
 'graphsage_edges':len(G['links']),'feature_cells':int(obs_feat.size),'label_cells':int(obs_label.size),
 'evidence_tier_counts':dict(tier_counts),'topology_exact_rows':sum(int(r['topology_neighbor_set_exact']) for r in evidence),
 'feature_exact_rows':sum(int(r['feature_vector_exact']) for r in evidence),'label_exact_rows':sum(int(r['go_label_vector_exact']) for r in evidence),
 'simple_mapping_file':str(simple_path),'simple_mapping_gzip_sha256':s_gz,'simple_mapping_uncompressed_tsv_sha256':s_uncomp,'simple_mapping_uncompressed_bytes':s_bytes,
 'evidence_mapping_file':str(evidence_path),'evidence_mapping_gzip_sha256':e_gz,'evidence_mapping_uncompressed_tsv_sha256':e_uncomp,'evidence_mapping_uncompressed_bytes':e_bytes,
 'graphsage_ppi_zip_sha256':sha256_file(GRAPH_ZIP),'ohmnet_networks_tar_sha256':sha256_file(OHMNET_TAR),'b104g_row_map_sha256':sha256_file(ROWMAP),
 'observed_feature_matrix_sha256':hashlib.sha256(obs_feat.tobytes(order='C')).hexdigest(),
 'reconstructed_feature_matrix_sha256':hashlib.sha256(exp_feat.tobytes(order='C')).hexdigest(),
 'observed_label_matrix_sha256':hashlib.sha256(obs_label.tobytes(order='C')).hexdigest(),
 'reconstructed_label_matrix_sha256':hashlib.sha256(exp_label.tobytes(order='C')).hexdigest(),
 'all_checks_pass':True,
}
summary_path=OUT/f'B104H_MAPPING_VALIDATION_{STAMP}.json'
summary_path.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2,sort_keys=True))

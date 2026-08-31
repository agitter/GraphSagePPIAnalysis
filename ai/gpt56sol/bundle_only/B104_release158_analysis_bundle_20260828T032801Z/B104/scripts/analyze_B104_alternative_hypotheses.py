#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, collections, json
from pathlib import Path

STAMP='20260828T030759Z'
ROOT=Path('/mnt/data/ppi_repro_corrected/batches/B104_20260828T030759Z')
ANA=ROOT/'analysis'
LABELS=ROOT/'retained_inputs/collapsed_gene_labels_topology_features.csv'
GAF159=ROOT/'retained_inputs/B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
GPI159=ROOT/'retained_inputs/B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz'
HIST=ROOT/'retained_inputs/B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
MAP_EDGES=ANA/f'B104_accession_GeneID_mapping_edges_{STAMP}.csv.gz'
TERM=ROOT/f'derived/B104_repaired_B103_GO_terms_{STAMP}.tsv.gz'
CLOSURE=ROOT/f'derived/B104_repaired_B103_GO_is_a_closure_for_GOA158_159_terms_{STAMP}.tsv.gz'
LABELMAP=ANA/f'B104_label_to_GO_mapping_release158_159_{STAMP}.csv'
WIT=ANA/f'B104_v159_residual_false_positive_witness_rows_{STAMP}.csv.gz'
BEST={'EXP','IDA','IEP','IGI','IMP','ISS'}

# Labels
rows=[]; genes=[]; observed=[]; zero_genes=set()
with LABELS.open() as f:
    for r in csv.DictReader(f):
        g=int(r['entrez_gene_id']); genes.append(g)
        vals=[int(r[f'label_{j}']) for j in range(121)]
        rows.append(vals)
        if sum(vals)==0: zero_genes.add(g)
idx={g:i for i,g in enumerate(genes)}; gset=set(genes); ALL=(1<<len(genes))-1
for j in range(121):
    b=0
    for i,r in enumerate(rows):
        if r[j]: b|=1<<i
    observed.append(b)

# GO alt IDs and closure
alt={}; name={}; ns={}
with gzip.open(TERM,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        name[r['GO_ID']]=r['name']; ns[r['GO_ID']]=r['namespace']
        for a in (r['alt_ids'] or '').split('|'):
            if a: alt[a]=r['GO_ID']
closure=collections.defaultdict(set)
with gzip.open(CLOSURE,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        closure[r['direct_GO_ID']].add(r['ancestor_GO_ID'])

def canonical(go): return alt.get(go,go)

# GPI accessions and mappings
accs=set()
with gzip.open(GPI159,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'): accs.add(r['DB_Object_ID'])
baseline=collections.defaultdict(set)
with gzip.open(MAP_EDGES,'rt') as f:
    for r in csv.DictReader(f): baseline[r['UniProtKB_accession']].add(int(r['GeneID']))
hist=collections.defaultdict(set); gene_to_acc=collections.defaultdict(set)
with gzip.open(HIST,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        g=int(r['GeneID']); a=r['UniProtKB_accession']
        if a in accs and g in gset:
            hist[a].add(g); gene_to_acc[g].add(a)

# GAF rows. Keep NOT separately.
pos=[]; not_rows=[]
with gzip.open(GAF159,'rt') as f:
    for r in csv.DictReader(f,delimiter='\t'):
        go=canonical(r['GO_ID']); isnot=(r.get('Is_NOT')=='1' or 'NOT' in (r.get('Qualifier') or '').split('|'))
        t=(r['DB_Object_ID'],go,r['Evidence_Code'],r['Date'],r['Assigned_By'],r['DB_Reference'])
        if isnot: not_rows.append(t)
        else: pos.append(t)

selected=[]
with LABELMAP.open() as f:
    for r in csv.DictReader(f): selected.append(canonical(r['GO_ID']))

# Predictor for mapping policies

def predict(mapping):
    bits=collections.defaultdict(int)
    for a,go,ev,date,assigned,ref in pos:
        if ev not in BEST: continue
        gs=mapping.get(a,set())
        if not gs: continue
        b=0
        for g in gs: b|=1<<idx[g]
        ancs=closure.get(go,{go})
        for q in ancs: bits[q]|=b
    return bits

def evaluate(mapping):
    bits=predict(mapping); fp=fn=exact=0; per=[]
    for j,go in enumerate(selected):
        o=observed[j]; p=bits.get(go,0); fpi=(p&~o&ALL).bit_count(); fni=(o&~p&ALL).bit_count(); d=fpi+fni
        fp+=fpi; fn+=fni; exact+=d==0
        per.append((j,d,fpi,fni))
    return {'exact_columns':exact,'at_least_99pct':sum((d/len(genes))<=.01 for _,d,_,_ in per),'at_least_95pct':sum((d/len(genes))<=.05 for _,d,_,_ in per),'total_mismatches':fp+fn,'false_positives':fp,'false_negatives':fn}

# Policies
policies={}
policies['component_aware_hybrid']=baseline
policies['historical_all_edges_no_symbol_fallback']=hist
policies['historical_accession_degree_1_only']={a:set(gs) for a,gs in hist.items() if len(gs)==1}
policies['historical_global_bijective_only']={a:set(gs) for a,gs in hist.items() if len(gs)==1 and len(gene_to_acc[next(iter(gs))])==1}
map_results=[]
for pname,m in policies.items():
    e=evaluate(m); e['mapping_policy']=pname; e['mapped_accessions']=sum(bool(v) for v in m.values()); e['covered_GraphSAGE_GeneIDs']=len(set().union(*[v for v in m.values() if v]) if m else set()); map_results.append(e)

# Residual pairs
fp_pairs=set(); fp_pair_witness_direct=set()
with gzip.open(WIT,'rt') as f:
    for r in csv.DictReader(f):
        p=(int(r['GeneID']),canonical(r['selected_GO_ID']),int(r['label_column']))
        fp_pairs.add(p)
        if str(r['direct_equals_selected']) in {'1','True','true'}: fp_pair_witness_direct.add(p)

# NOT mappings and overlaps
not_direct=set(); not_up=set(); mapped_not_rows=0
for a,go,ev,date,assigned,ref in not_rows:
    gs=baseline.get(a,set())
    if not gs: continue
    mapped_not_rows+=1
    for g in gs:
        not_direct.add((g,go))
        for q in closure.get(go,{go}): not_up.add((g,q))
fp_direct_not={p for p in fp_pairs if (p[0],p[1]) in not_direct}
fp_up_not={p for p in fp_pairs if (p[0],p[1]) in not_up}

zero_fp={p for p in fp_pairs if p[0] in zero_genes}
nonzero_fp=fp_pairs-zero_fp

# Direct vs ancestor residual counts by namespace/label
residual_by_ns=collections.Counter(); residual_by_label=collections.Counter(); residual_genes=set()
for g,go,j in fp_pairs:
    residual_by_ns[ns.get(go,'')]+=1; residual_by_label[j]+=1; residual_genes.add(g)

summary={
 'mapping_policy_sensitivity':map_results,
 'all_zero_gene_posthoc_mask':{
   'all_zero_GraphSAGE_genes':len(zero_genes),
   'residual_false_positive_pairs_on_all_zero_genes':len(zero_fp),
   'all_zero_genes_in_residuals':len({p[0] for p in zero_fp}),
   'residual_false_positive_pairs_remaining_after_invalid_posthoc_mask':len(nonzero_fp),
   'interpretation':'Masking genes using their observed all-zero label rows is target leakage and removes only a small fraction of residuals.'
 },
 'NOT_annotation_check':{
   'GAF159_NOT_rows':len(not_rows),
   'NOT_rows_with_component_aware_GeneID_mapping':mapped_not_rows,
   'unique_mapped_direct_GeneID_GO_NOT_pairs':len(not_direct),
   'residual_false_positive_pairs_with_direct_NOT_to_selected_term':len(fp_direct_not),
   'residual_false_positive_pairs_with_NOT_reaching_selected_term_by_is_a':len(fp_up_not),
   'interpretation':'Direct NOT annotations explain no residual false positives. Propagating NOT upward through is_a would remove only a few pairs and is not a generally valid inference.'
 },
 'residual_structure':{
   'residual_false_positive_pairs':len(fp_pairs),
   'unique_residual_genes':len(residual_genes),
   'pairs_with_direct_selected_term_witness':len(fp_pair_witness_direct),
   'pairs_ancestor_only':len(fp_pairs-fp_pair_witness_direct),
   'namespace_counts':dict(residual_by_ns),
   'labels_with_residuals':len(residual_by_label)
 }
}

out=ANA/f'B104_alternative_hypothesis_checks_{STAMP}.json'
out.write_text(json.dumps(summary,indent=2),encoding='utf-8')
with (ANA/f'B104_mapping_policy_sensitivity_{STAMP}.csv').open('w',newline='') as f:
    fn=['mapping_policy','mapped_accessions','covered_GraphSAGE_GeneIDs','exact_columns','at_least_99pct','at_least_95pct','total_mismatches','false_positives','false_negatives']
    w=csv.DictWriter(f,fieldnames=fn); w.writeheader(); w.writerows(map_results)
with (ANA/f'B104_residual_false_positives_on_all_zero_genes_{STAMP}.csv').open('w',newline='') as f:
    fn=['GeneID','selected_GO_ID','label_column']; w=csv.DictWriter(f,fieldnames=fn); w.writeheader()
    for g,go,j in sorted(zero_fp): w.writerow({'GeneID':g,'selected_GO_ID':go,'label_column':j})
print(json.dumps(summary,indent=2))

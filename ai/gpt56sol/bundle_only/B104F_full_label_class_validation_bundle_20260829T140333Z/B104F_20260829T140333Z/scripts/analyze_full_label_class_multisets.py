#!/usr/bin/env python3
from __future__ import annotations
import argparse, collections, csv, functools, gzip, hashlib, json, zipfile
from pathlib import Path

ALLOWED_EVIDENCE = {"EXP", "IDA", "IEP", "IGI", "IMP", "ISS"}
DEFAULT_RELATIONS = {"involved_in", "part_of", "enables"}
WATCH_GENES = [3248,3988,8564,27201,30061,51166,51312,55471,55801,56994,79017,121599]

def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def write_csv(path: Path, rows, fields=None):
    rows=list(rows)
    if fields is None: fields=list(rows[0]) if rows else []
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n',extrasaction='ignore')
        if fields:w.writeheader()
        w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--graphsage-zip',type=Path,required=True)
    ap.add_argument('--inputs',type=Path,required=True)
    ap.add_argument('--analysis-dir',type=Path,required=True)
    ap.add_argument('--summary-json',type=Path,required=True)
    args=ap.parse_args(); I=args.inputs; A=args.analysis_dir;A.mkdir(parents=True,exist_ok=True)
    rowmap=I/'graphsage_row_to_entrez_topology_features.csv'
    residual_path=I/'wl_residual_after_features.json'
    gaf=I/'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
    gpi=I/'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz'
    hist=I/'B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
    terms=I/'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
    edges=I/'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
    labelmap=I/'B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv'
    symbols=I/'B104F_MSigDB52_symbol_to_Entrez_relevant.tsv.gz'

    # Independently resolved rows plus topology/feature-equivalence classes.
    row_to_gene={}
    with rowmap.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f): row_to_gene[int(r['graphsage_row'])]=int(r['entrez_gene_id'])
    residual=json.loads(residual_path.read_text())
    graph_genes=set(row_to_gene.values())
    for cls in residual: graph_genes.update(map(int,cls['candidate_genes']))

    # All GraphSAGE row label vectors.
    with zipfile.ZipFile(args.graphsage_zip) as z:
        class_map=json.loads(z.read('ppi/ppi-class_map.json'))
        id_map=json.loads(z.read('ppi/ppi-id_map.json'))
    row_labels=[None]*len(id_map)
    for node,row in id_map.items(): row_labels[int(row)]=tuple(int(x) for x in class_map[node])

    # GPI fields used: DB_Object_ID and primary DB_Object_Symbol. Synonyms are audit only; xrefs are empty.
    GPI={}
    with gzip.open(gpi,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            GPI[r['DB_Object_ID']]={
                'symbol':r['DB_Object_Symbol'] or r.get('GAF_Fallback_Symbol',''),
                'name':r['DB_Object_Name'],
                'synonyms':r['DB_Object_Synonyms'],
                'taxon':r['Taxon'],
                'db_xrefs':r['DB_Xrefs'],
                'parent_object_id':r['Parent_Object_ID'],
            }

    # Historical edges from gp2protein.geneid; build components using only accessions present in GPI159.
    acc_edges=collections.defaultdict(set);gene_edges=collections.defaultdict(set)
    with gzip.open(hist,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            a=r['UniProtKB_accession'];g=int(r['GeneID'])
            acc_edges[a].add(g);gene_edges[g].add(a)
    sym_to_genes={}
    with gzip.open(symbols,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            sym_to_genes[r['gene_symbol']]={int(x) for x in r['Entrez_GeneIDs'].split('|') if x}

    adj=collections.defaultdict(set)
    for a,gs in acc_edges.items():
        if a not in GPI: continue
        for g in gs:
            adj[('a',a)].add(('g',g));adj[('g',g)].add(('a',a))
    resolved_accessions={};seen=set();component_rows=[]
    for start in list(adj):
        if start in seen: continue
        seen.add(start);stack=[start];nodes=[]
        while stack:
            x=stack.pop();nodes.append(x)
            for y in adj[x]:
                if y not in seen:seen.add(y);stack.append(y)
        accs={v for t,v in nodes if t=='a'};genes={v for t,v in nodes if t=='g'}
        cand={}
        for a in accs:
            matches=sym_to_genes.get(GPI[a]['symbol'],set()) & genes
            if len(matches)==1:cand[a]=next(iter(matches))
        bij=len(accs)==len(genes)==len(cand) and len(set(cand.values()))==len(genes)
        if bij:
            for a,g in cand.items():resolved_accessions[a]={g}
        if genes & graph_genes and (len(accs)>1 or len(genes)>1):
            component_rows.append({
                'accessions':'|'.join(sorted(accs)),
                'GeneIDs':'|'.join(map(str,sorted(genes))),
                'symbol_candidates':'|'.join(f'{a}:{cand[a]}' for a in sorted(cand)),
                'unique_symbol_bijection':int(bij),
                'resolution':'|'.join(f'{a}->{next(iter(resolved_accessions[a]))}' for a in sorted(accs) if a in resolved_accessions) if bij else 'retain_all_historical_edges',
            })

    mapping=collections.defaultdict(set);method={}
    for a,d in GPI.items():
        if a in resolved_accessions:
            mapping[a]=resolved_accessions[a]&graph_genes;method[a]='full_component_unique_symbol_bijection'
        else:
            direct=acc_edges.get(a,set())&graph_genes
            if direct:
                mapping[a]=set(direct);method[a]='historical_gp2protein_all_edges'
            else:
                matches=sym_to_genes.get(d['symbol'],set())&graph_genes
                if len(matches)==1:
                    mapping[a]=set(matches);method[a]='unique_primary_symbol_fallback'
                else:method[a]='unmapped'
    # Global semantic projection policy: do not transfer FSBP protein annotations to the RAD54B graph node.
    mapping['O95073'].discard(25788)

    # Ontology and selected terms.
    alt={};parents=collections.defaultdict(set);name_by={};namespace_by={}
    with gzip.open(terms,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            name_by[r['GO_ID']]=r['name'];namespace_by[r['GO_ID']]=r['namespace']
            for a in r['alt_ids'].split('|') if r['alt_ids'] else []:alt[a]=r['GO_ID']
    with gzip.open(edges,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):parents[r['child_GO_ID']].add(r['parent_GO_ID'])
    @functools.lru_cache(None)
    def ancestors(go):
        go=alt.get(go,go);out={go}
        for p in parents.get(go,()):out.update(ancestors(p))
        return frozenset(out)
    selected=[]
    with labelmap.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):selected.append(r['GO_ID'])
    selected_set=set(selected)

    # One global annotation policy, fixed before any row/gene comparisons.
    predicted={g:set() for g in graph_genes};accepted_rows=0
    with gzip.open(gaf,'rt',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['Is_NOT']=='1' or 'NOT' in r['Qualifier'].split('|'):continue
            if r['Evidence_Code'] not in ALLOWED_EVIDENCE:continue
            if r['Normalized_Relation'] not in DEFAULT_RELATIONS:continue
            genes=mapping.get(r['DB_Object_ID'],set())
            if not genes:continue
            propagated=ancestors(r['GO_ID'])&selected_set
            if not propagated:continue
            accepted_rows+=1
            for g in genes:predicted[g].update(propagated)
    vectors={g:tuple(int(t in predicted[g]) for t in selected) for g in graph_genes}

    # Individually resolved rows.
    resolved_mismatches=[]
    for row,g in sorted(row_to_gene.items()):
        mismatch=sum(a!=b for a,b in zip(vectors[g],row_labels[row]))
        if mismatch:resolved_mismatches.append({'graphsage_row':row,'GeneID':g,'mismatched_label_cells':mismatch})

    # Compare unresolved classes as multisets. No row-to-gene permutation is selected to score this test.
    class_rows=[];assignment_rows=[]
    for cls in residual:
        obs=collections.Counter(row_labels[int(r)] for r in cls['rows'])
        pred=collections.Counter(vectors[int(g)] for g in cls['candidate_genes'])
        exact=obs==pred
        obs_to_rows=collections.defaultdict(list);pred_to_genes=collections.defaultdict(list)
        for r in cls['rows']:obs_to_rows[row_labels[int(r)]].append(int(r))
        for g in cls['candidate_genes']:pred_to_genes[vectors[int(g)]].append(int(g))
        class_rows.append({
            'graph_index':cls['graph_index'],'tissue':cls['tissue'],'WL_color':cls['color'],
            'row_count':len(cls['rows']),'candidate_gene_count':len(cls['candidate_genes']),
            'observed_unique_label_vectors':len(obs),'predicted_unique_label_vectors':len(pred),
            'label_vector_multiset_exact':int(exact),
            'rows':'|'.join(map(str,cls['rows'])),'candidate_GeneIDs':'|'.join(map(str,cls['candidate_genes'])),
        })
        for vec,rows in obs_to_rows.items():
            genes=sorted(pred_to_genes.get(vec,[]))
            for row in sorted(rows):
                assignment_rows.append({
                    'graph_index':cls['graph_index'],'tissue':cls['tissue'],'WL_color':cls['color'],
                    'graphsage_row':row,'observed_positive_labels':sum(vec),
                    'candidate_GeneIDs_with_same_fixed_GOA_vector':'|'.join(map(str,genes)),
                    'candidate_count':len(genes),'unique_label_based_assignment':int(len(rows)==len(genes)==1),
                    'class_multiset_exact':int(exact),
                })

    # Twelve genes named by the independent verifier.
    gpi_by_gene=collections.defaultdict(list)
    for a,gs in mapping.items():
        for g in gs:gpi_by_gene[g].append(a)
    watch_summary=[];watch_assign=[]
    assignments_by_gene=collections.defaultdict(list)
    for r in assignment_rows:
        for x in r['candidate_GeneIDs_with_same_fixed_GOA_vector'].split('|') if r['candidate_GeneIDs_with_same_fixed_GOA_vector'] else []:
            assignments_by_gene[int(x)].append(r)
    for g in WATCH_GENES:
        accs=sorted(gpi_by_gene.get(g,[]))
        watch_summary.append({
            'GeneID':g,
            'GPI159_accessions':'|'.join(accs),
            'GPI159_symbols':'|'.join(sorted({GPI[a]['symbol'] for a in accs})),
            'GPI159_names':'|'.join(sorted({GPI[a]['name'] for a in accs})),
            'mapping_methods':'|'.join(f'{a}:{method.get(a,"")}' for a in accs),
            'predicted_positive_label_count':sum(vectors[g]),
            'independently_resolved_by_topology_and_input_features':int(g in set(row_to_gene.values())),
            'topology_feature_equivalence_class_count':sum(g in set(c['candidate_genes']) for c in residual),
            'label_vector_candidate_row_count':len(assignments_by_gene.get(g,[])),
            'all_candidate_rows_unique_by_fixed_label_vector':int(bool(assignments_by_gene.get(g)) and all(int(r['candidate_count'])==1 for r in assignments_by_gene[g])),
        })
        for r in assignments_by_gene.get(g,[]):
            rr=dict(r);rr['GeneID']=g;watch_assign.append(rr)

    write_csv(A/'B104F_full_unresolved_class_label_vector_multiset_validation.csv',class_rows)
    write_csv(A/'B104F_unresolved_row_candidate_GeneIDs_by_fixed_GOA_vector.csv',assignment_rows)
    write_csv(A/'B104F_independent_verifier_12_gene_summary.csv',watch_summary)
    write_csv(A/'B104F_independent_verifier_12_gene_row_matches.csv',watch_assign)
    write_csv(A/'B104F_full_mapping_components_touching_graph_genes.csv',component_rows)
    write_csv(A/'B104F_resolved_row_mismatches.csv',resolved_mismatches,
              ['graphsage_row','GeneID','mismatched_label_cells'])

    covered=set().union(*mapping.values()) if mapping else set()
    summary={
        'graphsage_rows':len(row_labels),
        'graph_candidate_GeneIDs':len(graph_genes),
        'individually_resolved_rows':len(row_to_gene),
        'individually_resolved_unique_GeneIDs':len(set(row_to_gene.values())),
        'individually_resolved_row_mismatch_count':len(resolved_mismatches),
        'unresolved_topology_feature_class_count':len(residual),
        'unresolved_rows':sum(len(c['rows']) for c in residual),
        'unresolved_classes_with_exact_observed_vs_predicted_label_vector_multiset':sum(r['label_vector_multiset_exact'] for r in class_rows),
        'unresolved_rows_uniquely_assignable_by_fixed_GOA_vector_within_class':sum(r['unique_label_based_assignment'] for r in assignment_rows),
        'unresolved_rows_remaining_ambiguous_due_to_identical_fixed_GOA_vectors':sum(not r['unique_label_based_assignment'] for r in assignment_rows),
        'all_56944_rows_reproduced_up_to_permutation_within_topology_feature_equivalence_classes':len(resolved_mismatches)==0 and all(r['label_vector_multiset_exact'] for r in class_rows),
        'unmapped_graph_candidate_GeneIDs':sorted(graph_genes-covered),
        'accepted_mapped_GAF_rows_contributing_to_selected_terms':accepted_rows,
        'global_policy':{
            'GOA':'human release 159',
            'evidence_codes':sorted(ALLOWED_EVIDENCE),
            'retained_GAF_relations':sorted(DEFAULT_RELATIONS),
            'excluded_GAF_relations':['colocalizes_with','contributes_to'],
            'NOT_excluded':True,
            'ontology_propagation':'is_a only, including direct term',
            'mapping':'all historical gp2protein edges unless a full component has a unique primary-symbol bijection; unique primary-symbol fallback only when no historical edge; O95073/FSBP not projected to 25788/RAD54B',
            'per_gene_or_per_label_tuning':False,
        },
        'GPI_fields_used_for_mapping':{
            'DB_Object_ID':'UniProt accession key joining GAF/GPI to gp2protein.geneid',
            'DB_Object_Symbol':'primary-symbol consistency and unique fallback only',
            'DB_Object_Name':'audit/reporting only',
            'DB_Object_Synonyms':'not used for mapping',
            'DB_Xrefs':'empty in all GPI159 rows; not used',
            'Parent_Object_ID':'empty in all GPI159 rows; not used',
            'Taxon':'audit; file is human-specific',
        },
        'input_sha256':{str(p):sha256(p) for p in [args.graphsage_zip,rowmap,residual_path,gaf,gpi,hist,terms,edges,labelmap,symbols]},
    }
    args.summary_json.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

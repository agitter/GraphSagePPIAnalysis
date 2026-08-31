#!/usr/bin/env python3
"""Reconstruct GraphSAGE PPI features from MSigDB v5.2 Entrez GMT source order.

Inputs are immutable user-supplied archives plus an independently recovered
GraphSAGE-row-to-Entrez map. Outputs record a complete cell-level comparison.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np

GMT_MEMBERS = [
    ("C1", "msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c1.all.v5.2.entrez.gmt"),
    ("C3", "msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c3.all.v5.2.entrez.gmt"),
    ("C7", "msigdb_v5.2_files_to_download_locally/msigdb_v5.2_GMTs/c7.all.v5.2.entrez.gmt"),
]

def sha256_path(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()

def read_mapping(path: Path):
    row_to_gene={}
    with path.open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            row_to_gene[int(r['graphsage_row'])]=int(r['entrez_gene_id'])
    return row_to_gene

def read_features(zip_path: Path):
    with zipfile.ZipFile(zip_path) as zf:
        data=zf.read('ppi/ppi-feats.npy')
    return np.load(io.BytesIO(data), allow_pickle=False)

def read_gmt_sets(zip_path: Path):
    all_sets=[]
    per_collection={}
    with zipfile.ZipFile(zip_path) as zf:
        for collection, member in GMT_MEMBERS:
            rows=[]
            text=zf.read(member).decode('utf-8')
            for idx,line in enumerate(text.splitlines()):
                if not line.strip():
                    continue
                parts=line.rstrip('\n\r').split('\t')
                name=parts[0]
                description=parts[1] if len(parts)>1 else ''
                genes=[]
                for x in parts[2:]:
                    try: genes.append(int(x))
                    except ValueError: pass
                gene_set=set(genes)
                rec={
                    'collection':collection,
                    'source_member':member,
                    'source_row_index_0based':idx,
                    'source_row_index_1based':idx+1,
                    'set_name':name,
                    'description':description,
                    'entrez_ids':gene_set,
                    'unique_entrez_count':len(gene_set),
                    'raw_entrez_field_count':len(genes),
                }
                rows.append(rec); all_sets.append(rec)
            per_collection[collection]=rows
    return all_sets, per_collection

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--graphsage-zip', required=True, type=Path)
    ap.add_argument('--msigdb-zip', required=True, type=Path)
    ap.add_argument('--row-map', required=True, type=Path)
    ap.add_argument('--output-csv', required=True, type=Path)
    ap.add_argument('--summary-json', required=True, type=Path)
    ap.add_argument('--min-size', type=int, default=200)
    ap.add_argument('--max-features', type=int, default=50)
    args=ap.parse_args()
    row_to_gene=read_mapping(args.row_map)
    feats=read_features(args.graphsage_zip)
    if feats.shape[1] != args.max_features:
        raise SystemExit(f'unexpected feature columns: {feats.shape}')
    if not np.all((feats==0)|(feats==1)):
        raise SystemExit('GraphSAGE features are not binary')
    all_sets, per_collection=read_gmt_sets(args.msigdb_zip)
    selected=[]
    qualifying_counts={}
    for collection,_ in GMT_MEMBERS:
        qualifying=[r for r in per_collection[collection] if r['unique_entrez_count']>=args.min_size]
        qualifying_counts[collection]=len(qualifying)
        for r in qualifying:
            if len(selected) >= args.max_features: break
            selected.append(r)
        if len(selected) >= args.max_features: break
    if len(selected) != args.max_features:
        raise SystemExit(f'only selected {len(selected)} features')

    resolved_rows=np.array(sorted(row_to_gene), dtype=int)
    genes=np.array([row_to_gene[int(i)] for i in resolved_rows], dtype=int)
    observed=feats[resolved_rows,:]
    expected=np.zeros_like(observed)
    for j,r in enumerate(selected):
        expected[:,j]=np.array([g in r['entrez_ids'] for g in genes], dtype=observed.dtype)
    diff=(observed!=expected)

    # Check repeated appearances of the same Entrez gene have one identical vector.
    gene_to_first={}; repeated_conflicts=[]
    for row,g in zip(resolved_rows,genes):
        vec=feats[row,:]
        if g in gene_to_first:
            if not np.array_equal(vec, gene_to_first[g][1]):
                repeated_conflicts.append((int(g), int(gene_to_first[g][0]), int(row)))
        else:
            gene_to_first[int(g)]=(int(row), vec.copy())

    out=[]
    for j,r in enumerate(selected):
        col_diff=diff[:,j]
        resolved_unique_genes=set(int(x) for x in genes)
        record={k:v for k,v in r.items() if k!='entrez_ids'}
        record.update({
            'graphsage_feature_column':j,
            'selection_threshold_min_unique_entrez':args.min_size,
            'selection_global_cap':args.max_features,
            'selected_rank_1based':j+1,
            'resolved_unique_graph_genes_in_set':len(resolved_unique_genes & r['entrez_ids']),
            'observed_positive_resolved_rows':int(observed[:,j].sum()),
            'expected_positive_resolved_rows':int(expected[:,j].sum()),
            'mismatch_resolved_rows':int(col_diff.sum()),
            'exact_on_all_resolved_rows':int(not col_diff.any()),
            'all_zero_on_all_graphsage_rows':int(feats[:,j].sum()==0),
            'all_zero_on_resolved_rows':int(observed[:,j].sum()==0),
        })
        out.append(record)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader(); w.writerows(out)

    summary={
        'input_sha256':{
            'graphsage_ppi_zip':sha256_path(args.graphsage_zip),
            'msigdb_v5_2_zip':sha256_path(args.msigdb_zip),
            'row_map':sha256_path(args.row_map),
        },
        'feature_matrix_shape':[int(x) for x in feats.shape],
        'feature_matrix_dtype':str(feats.dtype),
        'resolved_rows':int(len(resolved_rows)),
        'resolved_unique_entrez_genes':int(len(set(int(x) for x in genes))),
        'resolved_feature_cells_compared':int(diff.size),
        'total_resolved_feature_mismatches':int(diff.sum()),
        'exact_feature_columns':int(sum(not diff[:,j].any() for j in range(diff.shape[1]))),
        'all_zero_feature_columns_all_rows':[int(j) for j in range(feats.shape[1]) if feats[:,j].sum()==0],
        'repeated_gene_feature_conflicts':repeated_conflicts,
        'selection_rule':{
            'collection_order':[x[0] for x in GMT_MEMBERS],
            'minimum_unique_entrez_count_in_full_msigdb_set':args.min_size,
            'global_feature_cap':args.max_features,
            'preserve_gmt_source_order':True,
            'stop_when_global_cap_reached':True,
        },
        'qualifying_set_count_by_collection':qualifying_counts,
        'selected_set_count_by_collection':dict((c,sum(r['collection']==c for r in selected)) for c,_ in GMT_MEMBERS),
        'first_unselected_qualifying_set':None,
        'selected_feature_csv':str(args.output_csv),
    }
    # Identify the next qualifying set after the first 50 under the same collection stream.
    stream=[]
    for c,_ in GMT_MEMBERS:
        stream += [r for r in per_collection[c] if r['unique_entrez_count']>=args.min_size]
    if len(stream)>args.max_features:
        r=stream[args.max_features]
        summary['first_unselected_qualifying_set']={k:v for k,v in r.items() if k!='entrez_ids'}
    args.summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    if summary['total_resolved_feature_mismatches']!=0 or summary['exact_feature_columns']!=50 or repeated_conflicts:
        raise SystemExit('feature reconstruction validation failed')
    print(json.dumps(summary, indent=2, sort_keys=True))
if __name__=='__main__': main()

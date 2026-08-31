#!/usr/bin/env python3
from __future__ import annotations

import collections
import csv
import functools
import gzip
import hashlib
import io
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import sys
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, percentileofscore, spearmanr

STAMP = os.environ.get('B104A_STAMP', '20260828T145842Z')
ROOT = Path(f'/mnt/data/ppi_repro_corrected/batches/B104A_{STAMP}')
ANA = ROOT / 'analysis'
DER = ROOT / 'derived'
LOG = ROOT / 'logs'
RI = ROOT / 'retained_inputs'
for d in (ANA, DER, LOG, RI):
    d.mkdir(parents=True, exist_ok=True)

B104_BUNDLE = Path('/mnt/data/B104_release158_analysis_bundle_20260828T032801Z.zip')
CORE_BUNDLE = Path('/mnt/data/ppi_reproduction_corrected_bundle.zip')
GRAPHSAGE_ZIP = Path('/mnt/data/graphsage_ppi.zip')
EXPECTED_BUNDLES = {
    str(B104_BUNDLE): 'c131fb1f65d4b8c242b60d3db3ea815dc873674b18162d8dfd164f32f5609570',
    str(CORE_BUNDLE): 'ffa1997c4311ad859f109de01a0fd5580f6fd301bee0e60023c1ebe817a338e4',
    str(GRAPHSAGE_ZIP): '53aeb76e54fd41b645e7edb48b62929240b89839495396b048086fd212503fbd',
}

LABELMAP = RI / 'B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv'
MAP_EDGES = RI / 'B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz'
TERMS_FILE = RI / 'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES_FILE = RI / 'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
GAF_FILE = RI / 'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
LABELS_FILE = RI / 'collapsed_gene_labels_topology_features.csv'
MISSING_HIST = RI / 'B102_GraphSAGE_genes_missing_from_historical_GPI_projection_20260827T162132Z.csv'
ROW_MAP = RI / 'graphsage_row_to_entrez_topology_features.csv'
TISSUES = RI / 'tissue_partition.csv'

ORIGINAL_EVIDENCE = frozenset({'EXP','IDA','IMP','IGI','IEP','ISS'})
DEFAULT_RELATIONS = frozenset({'involved_in','part_of','enables'})
QUALIFIED_RELATIONS = frozenset({'colocalizes_with','contributes_to'})
ALL_POSITIVE_RELATIONS = DEFAULT_RELATIONS | QUALIFIED_RELATIONS
ROOTS = {
    'biological_process': 'GO:0008150',
    'cellular_component': 'GO:0005575',
    'molecular_function': 'GO:0003674',
}


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def det_gzip_text_writer(path: Path):
    raw = path.open('wb')
    gz = gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9)
    txt = io.TextIOWrapper(gz, encoding='utf-8', newline='')
    return raw, gz, txt


def write_csv(path: Path, rows, fieldnames=None, delimiter=','):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter=delimiter, lineterminator='\n', extrasaction='ignore')
        if fieldnames:
            w.writeheader()
        w.writerows(rows)


def write_gz_csv(path: Path, rows, fieldnames=None, delimiter=','):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    raw, gz, txt = det_gzip_text_writer(path)
    try:
        w = csv.DictWriter(txt, fieldnames=fieldnames, delimiter=delimiter, lineterminator='\n', extrasaction='ignore')
        if fieldnames:
            w.writeheader()
        w.writerows(rows)
    finally:
        txt.flush()
        txt.detach()
        gz.close()
        raw.close()


def bits_to_genes(bits: int, genes: list[int]):
    while bits:
        lsb = bits & -bits
        i = lsb.bit_length() - 1
        yield genes[i]
        bits -= lsb


# ---------- input verification ----------
input_rows = []
for p, expected in EXPECTED_BUNDLES.items():
    path = Path(p)
    actual = sha256_file(path)
    input_rows.append({
        'artifact_name': path.name,
        'artifact_role': 'frozen_upstream_bundle' if path != GRAPHSAGE_ZIP else 'raw_core_dataset_retained_from_initial_upload',
        'local_path': str(path),
        'size_bytes': path.stat().st_size,
        'sha256': actual,
        'expected_sha256': expected,
        'hash_matches_expected': actual == expected,
    })
    assert actual == expected, (path, actual, expected)
for p in [LABELMAP, MAP_EDGES, TERMS_FILE, EDGES_FILE, GAF_FILE, LABELS_FILE, MISSING_HIST, ROW_MAP, TISSUES]:
    input_rows.append({
        'artifact_name': p.name,
        'artifact_role': 'retained_derivative_extracted_from_frozen_prior_bundle',
        'local_path': str(p),
        'size_bytes': p.stat().st_size,
        'sha256': sha256_file(p),
        'expected_sha256': '',
        'hash_matches_expected': '',
    })
write_csv(ROOT / f'B104A_input_integrity_{STAMP}.csv', input_rows)

# ---------- ontology ----------
terms_df = pd.read_csv(TERMS_FILE, sep='\t', dtype=str).fillna('')
edges_df = pd.read_csv(EDGES_FILE, sep='\t', dtype=str).fillna('')
name_by = dict(zip(terms_df.GO_ID, terms_df.name))
namespace_by = dict(zip(terms_df.GO_ID, terms_df.namespace))
alt_to_primary: dict[str, str] = {}
for r in terms_df.itertuples(index=False):
    if r.alt_ids:
        for a in r.alt_ids.split('|'):
            if a:
                alt_to_primary[a] = r.GO_ID
parents: dict[str, set[str]] = collections.defaultdict(set)
children: dict[str, set[str]] = collections.defaultdict(set)
for r in edges_df.itertuples(index=False):
    parents[r.child_GO_ID].add(r.parent_GO_ID)
    children[r.parent_GO_ID].add(r.child_GO_ID)

@functools.lru_cache(None)
def ancestors(go: str) -> frozenset[str]:
    out = {go}
    for p in parents.get(go, ()):
        out.update(ancestors(p))
    return frozenset(out)

@functools.lru_cache(None)
def descendants(go: str) -> frozenset[str]:
    out = {go}
    for c in children.get(go, ()):
        out.update(descendants(c))
    return frozenset(out)

@functools.lru_cache(None)
def min_depth(go: str):
    root = ROOTS.get(namespace_by.get(go, ''))
    if go == root:
        return 0
    ds = [min_depth(p) for p in parents.get(go, ()) if namespace_by.get(p) == namespace_by.get(go)]
    ds = [d for d in ds if d is not None]
    return 1 + min(ds) if ds else None

@functools.lru_cache(None)
def max_depth(go: str):
    root = ROOTS.get(namespace_by.get(go, ''))
    if go == root:
        return 0
    ds = [max_depth(p) for p in parents.get(go, ()) if namespace_by.get(p) == namespace_by.get(go)]
    ds = [d for d in ds if d is not None]
    return 1 + max(ds) if ds else None

@functools.lru_cache(None)
def ancestor_distances(go: str) -> dict[str, int]:
    d = {go: 0}
    q = collections.deque([go])
    while q:
        x = q.popleft()
        for p in parents.get(x, ()):
            nd = d[x] + 1
            if p not in d or nd < d[p]:
                d[p] = nd
                q.append(p)
    return d

# ---------- labels and selected terms ----------
labelmap = pd.read_csv(LABELMAP)
labels_df = pd.read_csv(LABELS_FILE)
genes = labels_df.entrez_gene_id.astype(int).tolist()
gene_index = {g: i for i, g in enumerate(genes)}
N = len(genes)
ALL_BITS = (1 << N) - 1
observed_bits: list[int] = []
for j in range(121):
    b = 0
    for i, v in enumerate(labels_df[f'label_{j}'].values):
        if int(v):
            b |= 1 << i
    observed_bits.append(b)
selected_terms = labelmap.GO_ID.tolist()
unique_selected = sorted(set(selected_terms))
columns_by_term: dict[str, list[int]] = collections.defaultdict(list)
for j, t in enumerate(selected_terms):
    columns_by_term[t].append(j)

# ---------- accepted mapping ----------
map_df = pd.read_csv(MAP_EDGES)
accession_to_genes: dict[str, set[int]] = collections.defaultdict(set)
accession_gene_method: dict[tuple[str,int], str] = {}
for r in map_df.itertuples(index=False):
    g = int(r.GeneID)
    accession_to_genes[str(r.UniProtKB_accession)].add(g)
    accession_gene_method[(str(r.UniProtKB_accession), g)] = str(r.mapping_method)

# ---------- parse GAF and aggregate selected-term support ----------
# annotation tuple fields: accession, gene_bits, direct_go, evidence, relation, date, source, reference, qualifier
annotations = []
all_evidence_codes = set()
all_sources = set()
all_dates = set()
with gzip.open(GAF_FILE, 'rt', encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        if r['Is_NOT'] == '1' or 'NOT' in (r['Qualifier'] or '').split('|'):
            continue
        gs = accession_to_genes.get(r['DB_Object_ID'], set())
        if not gs:
            continue
        bits = 0
        for g in gs:
            if g in gene_index:
                bits |= 1 << gene_index[g]
        if not bits:
            continue
        go = alt_to_primary.get(r['GO_ID'], r['GO_ID'])
        evidence = r['Evidence_Code']
        relation = r['Normalized_Relation']
        date = r['Date']
        source = r['Assigned_By']
        annotations.append((r['DB_Object_ID'], bits, go, evidence, relation, date, source, r['DB_Reference'], r['Qualifier']))
        all_evidence_codes.add(evidence)
        all_sources.add(source)
        all_dates.add(date)
all_evidence_codes = sorted(all_evidence_codes)
all_sources = sorted(all_sources)
all_dates = sorted(all_dates)

# support[evidence][relation][selected term] => bitset; also distance and direct-term support
bucket: dict[tuple[str,str,str], int] = collections.defaultdict(int)
distance_bucket: dict[tuple[str,str,str,int], int] = collections.defaultdict(int)
direct_term_gene_bits: dict[tuple[str,str,str], int] = collections.defaultdict(int)  # evidence, relation, direct go
for accession, bits, direct_go, evidence, relation, date, source, ref, qualifier in annotations:
    dmap = ancestor_distances(direct_go)
    targets = set(dmap).intersection(unique_selected)
    for t in targets:
        bucket[(evidence, relation, t)] |= bits
        distance_bucket[(evidence, relation, t, dmap[t])] |= bits
    direct_term_gene_bits[(evidence, relation, direct_go)] |= bits


def predictions(evidence_codes: set[str] | frozenset[str], relations: set[str] | frozenset[str], max_distance_by_namespace=None, qualified_direct_only=False) -> dict[str, int]:
    out = {t: 0 for t in unique_selected}
    if max_distance_by_namespace is None and not qualified_direct_only:
        for e in evidence_codes:
            for rel in relations:
                for t in unique_selected:
                    out[t] |= bucket.get((e, rel, t), 0)
        return out
    for (e, rel, t, d), bits in distance_bucket.items():
        if e not in evidence_codes or rel not in relations:
            continue
        if qualified_direct_only and rel in QUALIFIED_RELATIONS and d != 0:
            continue
        maxd = 10**9 if max_distance_by_namespace is None else max_distance_by_namespace.get(namespace_by.get(t,''), 10**9)
        if d <= maxd:
            out[t] |= bits
    return out


def eval_predictions(pred: dict[str, int]):
    rows = []
    fp_total = fn_total = exact = ge99 = ge95 = 0
    for j, t in enumerate(selected_terms):
        p = pred.get(t, 0)
        o = observed_bits[j]
        fp = (p & ~o & ALL_BITS).bit_count()
        fn = (o & ~p & ALL_BITS).bit_count()
        d = fp + fn
        agreement = 1 - d / N
        rows.append({'label_column': j, 'GO_ID': t, 'false_positives': fp, 'false_negatives': fn, 'mismatches': d, 'agreement': agreement, 'predicted_positive_genes': p.bit_count(), 'observed_positive_genes': o.bit_count()})
        fp_total += fp
        fn_total += fn
        exact += d == 0
        ge99 += agreement >= .99
        ge95 += agreement >= .95
    return {'false_positives': fp_total, 'false_negatives': fn_total, 'total_mismatches': fp_total + fn_total, 'exact_columns': exact, 'at_least_99pct': ge99, 'at_least_95pct': ge95, 'per_column': rows}

# Validate prior accepted baseline.
pred_all_rel = predictions(ORIGINAL_EVIDENCE, ALL_POSITIVE_RELATIONS)
eval_all_rel = eval_predictions(pred_all_rel)
assert eval_all_rel['false_positives'] == 901 and eval_all_rel['false_negatives'] == 0 and eval_all_rel['exact_columns'] == 89, eval_all_rel

# ---------- relation/qualifier policies ----------
relation_policies = {
    'all_positive_relations': ALL_POSITIVE_RELATIONS,
    'exclude_colocalizes_with': ALL_POSITIVE_RELATIONS - {'colocalizes_with'},
    'exclude_contributes_to': ALL_POSITIVE_RELATIONS - {'contributes_to'},
    'default_relations_only': DEFAULT_RELATIONS,
}
relation_summary = []
relation_evals = {}
for name, rels in relation_policies.items():
    ev = eval_predictions(predictions(ORIGINAL_EVIDENCE, rels))
    relation_evals[name] = ev
    relation_summary.append({'policy': name, 'included_relations': '|'.join(sorted(rels)), **{k:v for k,v in ev.items() if k != 'per_column'}})
# qualified rows direct-only, but default rows fully propagated
pred_qualified_direct_only = predictions(ORIGINAL_EVIDENCE, ALL_POSITIVE_RELATIONS, qualified_direct_only=True)
ev_qdirect = eval_predictions(pred_qualified_direct_only)
relation_summary.append({'policy':'default_propagation_plus_qualified_direct_only','included_relations':'all; qualified relations only at distance 0', **{k:v for k,v in ev_qdirect.items() if k!='per_column'}})
write_csv(ANA / f'B104A_relation_qualifier_policy_summary_{STAMP}.csv', relation_summary)

pred_default = predictions(ORIGINAL_EVIDENCE, DEFAULT_RELATIONS)
eval_default = eval_predictions(pred_default)
assert eval_default['false_positives'] == 13 and eval_default['false_negatives'] == 0 and eval_default['exact_columns'] == 108, eval_default

# Baseline FP support decomposition.
fp_decomp = []
for j, t in enumerate(selected_terms):
    o = observed_bits[j]
    p_all = pred_all_rel[t]
    p_def = pred_default[t]
    p_qual = predictions(ORIGINAL_EVIDENCE, QUALIFIED_RELATIONS)[t]
    fp_all = p_all & ~o & ALL_BITS
    for g in bits_to_genes(fp_all, genes):
        bit = 1 << gene_index[g]
        fp_decomp.append({
            'label_column': j,
            'GO_ID': t,
            'GO_name': name_by.get(t,''),
            'namespace': namespace_by.get(t,''),
            'GeneID': g,
            'supported_by_default_relation': int(bool(p_def & bit)),
            'supported_by_qualified_relation': int(bool(p_qual & bit)),
            'eliminated_when_qualifier_rows_excluded': int(not bool(p_def & bit)),
        })
write_gz_csv(ANA / f'B104A_original_901_false_positive_decomposition_{STAMP}.csv.gz', fp_decomp)

# ---------- term depth/breadth metrics ----------
# direct selected bits under all and default policies.
def direct_bits_for(term: str, evidence: set[str] | frozenset[str], relations: set[str] | frozenset[str]):
    b = 0
    for e in evidence:
        for rel in relations:
            b |= distance_bucket.get((e, rel, term, 0), 0)
    return b

term_metric_rows = []
for row in labelmap.itertuples(index=False):
    j = int(row.label_column)
    t = row.GO_ID
    p_all = pred_all_rel[t]
    p_def = pred_default[t]
    o = observed_bits[j]
    fp_all = (p_all & ~o & ALL_BITS).bit_count()
    fp_def = (p_def & ~o & ALL_BITS).bit_count()
    md = min_depth(t)
    xd = max_depth(t)
    term_metric_rows.append({
        'label_column': j,
        'GO_ID': t,
        'GO_name': name_by.get(t,''),
        'namespace': namespace_by.get(t,''),
        'min_is_a_depth_from_namespace_root': md,
        'max_is_a_depth_from_namespace_root': xd,
        'depth_span': (xd-md if md is not None and xd is not None else ''),
        'ancestor_count_including_self': len(ancestors(t)),
        'descendant_count_including_self': len(descendants(t)),
        'observed_positive_genes': o.bit_count(),
        'all_relation_predicted_positive_genes': p_all.bit_count(),
        'all_relation_false_positives': fp_all,
        'default_relation_predicted_positive_genes': p_def.bit_count(),
        'default_relation_false_positives': fp_def,
        'false_positives_eliminated_by_excluding_qualified_relations': fp_all-fp_def,
        'all_relation_exact': int((p_all ^ o).bit_count()==0),
        'default_relation_exact': int((p_def ^ o).bit_count()==0),
        'direct_positive_genes_all_relations': direct_bits_for(t, ORIGINAL_EVIDENCE, ALL_POSITIVE_RELATIONS).bit_count(),
        'direct_positive_genes_default_relations': direct_bits_for(t, ORIGINAL_EVIDENCE, DEFAULT_RELATIONS).bit_count(),
        'predicted_information_content_bits_default': -math.log2(p_def.bit_count()/N) if p_def else '',
    })
write_csv(ANA / f'B104A_selected_term_ontology_depth_and_residual_metrics_{STAMP}.csv', term_metric_rows)
term_metrics = pd.DataFrame(term_metric_rows)

comparison_rows = []
def cliffs_delta(a, b):
    a = list(a); b = list(b)
    if not a or not b:
        return None
    gt = sum(x > y for x in a for y in b)
    lt = sum(x < y for x in a for y in b)
    return (gt - lt) / (len(a)*len(b))
for scope, sub in [('all_aspects', term_metrics), ('cellular_component_and_molecular_function_only', term_metrics[term_metrics.namespace != 'biological_process'])]:
    nonexact = sub[sub.all_relation_exact == 0]
    exact = sub[sub.all_relation_exact == 1]
    for metric in ['min_is_a_depth_from_namespace_root','max_is_a_depth_from_namespace_root','descendant_count_including_self','all_relation_predicted_positive_genes','predicted_information_content_bits_default']:
        x = pd.to_numeric(nonexact[metric], errors='coerce').dropna().tolist()
        y = pd.to_numeric(exact[metric], errors='coerce').dropna().tolist()
        if x and y:
            test = mannwhitneyu(x, y, alternative='two-sided')
            comparison_rows.append({
                'scope': scope, 'group_definition': 'nonexact_under_all_relations_vs_exact_under_all_relations', 'metric': metric,
                'nonexact_n': len(x), 'exact_n': len(y), 'nonexact_median': statistics.median(x), 'exact_median': statistics.median(y),
                'mann_whitney_u': float(test.statistic), 'two_sided_p_value': float(test.pvalue), 'cliffs_delta_nonexact_minus_exact': cliffs_delta(x,y),
                'caveat': 'Aspect composition is confounded in all_aspects; the CC/MF comparison has only four exact columns.'
            })
write_csv(ANA / f'B104A_depth_breadth_group_comparisons_{STAMP}.csv', comparison_rows)

# ---------- evidence mask search with default relations ----------
ev_codes = all_evidence_codes
# bitsets by evidence and selected term, default relations only
ev_term_bits = {e: {t: 0 for t in unique_selected} for e in ev_codes}
for e in ev_codes:
    for rel in DEFAULT_RELATIONS:
        for t in unique_selected:
            ev_term_bits[e][t] |= bucket.get((e, rel, t), 0)

label_indices_by_aspect = collections.defaultdict(list)
for j, t in enumerate(selected_terms):
    label_indices_by_aspect[namespace_by.get(t,'')].append(j)

def evaluate_mask(mask: int, indices=None):
    included = [ev_codes[i] for i in range(len(ev_codes)) if mask & (1 << i)]
    if indices is None:
        indices = range(121)
    fp = fn = exact = ge99 = ge95 = 0
    for j in indices:
        t = selected_terms[j]
        p = 0
        for e in included:
            p |= ev_term_bits[e][t]
        o = observed_bits[j]
        x = (p & ~o & ALL_BITS).bit_count()
        y = (o & ~p & ALL_BITS).bit_count()
        d = x+y
        fp += x; fn += y; exact += d==0; ge99 += 1-d/N >= .99; ge95 += 1-d/N >= .95
    return included, fp, fn, exact, ge99, ge95

mask_rows = []
for mask in range(1, 1 << len(ev_codes)):
    included, fp, fn, exact, ge99, ge95 = evaluate_mask(mask)
    mask_rows.append({'mask_integer':mask,'evidence_codes':'|'.join(included),'evidence_count':len(included),'false_positives':fp,'false_negatives':fn,'total_mismatches':fp+fn,'exact_columns':exact,'at_least_99pct':ge99,'at_least_95pct':ge95})
mask_rows.sort(key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives'],-r['exact_columns'],r['evidence_count'],r['evidence_codes']))
write_gz_csv(ANA / f'B104A_default_relation_global_evidence_mask_search_all_{STAMP}.csv.gz', mask_rows)
write_csv(ANA / f'B104A_default_relation_global_evidence_mask_search_top100_{STAMP}.csv', mask_rows[:100])
best_global = mask_rows[0]
best_global_evidence = frozenset(best_global['evidence_codes'].split('|'))

aspect_mask_rows = []
best_aspect = {}
for aspect, indices in label_indices_by_aspect.items():
    rows = []
    for mask in range(1, 1 << len(ev_codes)):
        included, fp, fn, exact, ge99, ge95 = evaluate_mask(mask, indices)
        rows.append({'namespace':aspect,'mask_integer':mask,'evidence_codes':'|'.join(included),'evidence_count':len(included),'false_positives':fp,'false_negatives':fn,'total_mismatches':fp+fn,'exact_columns':exact,'namespace_column_count':len(indices),'at_least_99pct':ge99,'at_least_95pct':ge95})
    rows.sort(key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives'],-r['exact_columns'],r['evidence_count'],r['evidence_codes']))
    best_aspect[aspect] = rows[0]
    aspect_mask_rows.extend(rows)
write_gz_csv(ANA / f'B104A_default_relation_aspect_specific_evidence_mask_search_all_{STAMP}.csv.gz', aspect_mask_rows)
write_csv(ANA / f'B104A_default_relation_aspect_specific_evidence_mask_search_top50_each_{STAMP}.csv', [r for aspect in sorted(best_aspect) for r in sorted((x for x in aspect_mask_rows if x['namespace']==aspect), key=lambda x:(x['total_mismatches'],x['false_negatives'],x['false_positives'],-x['exact_columns'],x['evidence_count'],x['evidence_codes']))[:50]])

# Aspect-specific combined prediction.
pred_aspect_best = {t: 0 for t in unique_selected}
for aspect, best in best_aspect.items():
    codes = set(best['evidence_codes'].split('|'))
    for t in unique_selected:
        if namespace_by.get(t) != aspect:
            continue
        for e in codes:
            pred_aspect_best[t] |= ev_term_bits[e][t]
eval_aspect_best = eval_predictions(pred_aspect_best)

pred_best_global = predictions(best_global_evidence, DEFAULT_RELATIONS)
eval_best_global = eval_predictions(pred_best_global)

# ---------- date cutoff and source leave-one-out under best global and original filters ----------
def build_predictions_from_annotations(evidence_codes, excluded_source=None, cutoff_date=None):
    out = {t:0 for t in unique_selected}
    for accession, bits, direct_go, evidence, relation, date, source, ref, qualifier in annotations:
        if evidence not in evidence_codes or relation not in DEFAULT_RELATIONS:
            continue
        if excluded_source is not None and source == excluded_source:
            continue
        if cutoff_date is not None and date > cutoff_date:
            continue
        for t in set(ancestor_distances(direct_go)).intersection(unique_selected):
            out[t] |= bits
    return out

filter_specs = [('original_six', ORIGINAL_EVIDENCE), ('best_global_after_qualifier_exclusion', best_global_evidence)]
date_rows = []
for spec_name, codes in filter_specs:
    for cutoff in all_dates:
        ev = eval_predictions(build_predictions_from_annotations(codes, cutoff_date=cutoff))
        date_rows.append({'filter_spec':spec_name,'cutoff_date_inclusive':cutoff, **{k:v for k,v in ev.items() if k!='per_column'}})
write_csv(ANA / f'B104A_default_relation_annotation_date_cutoff_search_{STAMP}.csv', date_rows)

source_rows = []
for spec_name, codes in filter_specs:
    baseline = eval_predictions(build_predictions_from_annotations(codes))
    for source in all_sources:
        ev = eval_predictions(build_predictions_from_annotations(codes, excluded_source=source))
        source_rows.append({'filter_spec':spec_name,'excluded_assigned_by':source, **{k:v for k,v in ev.items() if k!='per_column'}, 'delta_total_mismatches_vs_no_source_exclusion':ev['total_mismatches']-baseline['total_mismatches']})
write_csv(ANA / f'B104A_default_relation_source_leave_one_out_{STAMP}.csv', source_rows)

# ---------- propagation-distance search, default relations ----------
distance_rows = []
for spec_name, codes in filter_specs:
    for cc_max in list(range(0,13)) + [99]:
        for mf_max in list(range(0,13)) + [99]:
            pred = predictions(codes, DEFAULT_RELATIONS, {'biological_process':99,'cellular_component':cc_max,'molecular_function':mf_max})
            ev = eval_predictions(pred)
            distance_rows.append({'filter_spec':spec_name,'CC_max_is_a_distance':cc_max,'MF_max_is_a_distance':mf_max, **{k:v for k,v in ev.items() if k!='per_column'}})
write_csv(ANA / f'B104A_default_relation_propagation_distance_search_{STAMP}.csv', distance_rows)

# ---------- remaining pairs under three candidate policies ----------
policies = {
    'original_six_default_relations': pred_default,
    'best_global_evidence_default_relations': pred_best_global,
    'best_aspect_specific_evidence_default_relations': pred_aspect_best,
}
remaining_pair_rows = []
for policy_name, pred in policies.items():
    for j, t in enumerate(selected_terms):
        p = pred[t]; o = observed_bits[j]
        fp_bits = p & ~o & ALL_BITS
        fn_bits = o & ~p & ALL_BITS
        for g in bits_to_genes(fp_bits, genes):
            remaining_pair_rows.append({'policy':policy_name,'difference_type':'false_positive','label_column':j,'GO_ID':t,'GO_name':name_by.get(t,''),'namespace':namespace_by.get(t,''),'GeneID':g})
        for g in bits_to_genes(fn_bits, genes):
            remaining_pair_rows.append({'policy':policy_name,'difference_type':'false_negative','label_column':j,'GO_ID':t,'GO_name':name_by.get(t,''),'namespace':namespace_by.get(t,''),'GeneID':g})
write_csv(ANA / f'B104A_remaining_gene_label_differences_by_policy_{STAMP}.csv', remaining_pair_rows)

# Witness rows and path support for original-six/default relation 13 FPs.
remaining13 = [r for r in remaining_pair_rows if r['policy']=='original_six_default_relations' and r['difference_type']=='false_positive']
remaining_keys = {(r['label_column'],r['GO_ID'],r['GeneID']) for r in remaining13}
witness_rows = []
pair_support = collections.defaultdict(lambda: {'direct_terms':set(),'distances':[],'relations':set(),'accessions':set(),'mapping_methods':set(),'sources':set(),'dates':set(),'evidence':set(),'references':set()})
for accession, bits, direct_go, evidence, relation, date, source, ref, qualifier in annotations:
    if evidence not in ORIGINAL_EVIDENCE or relation not in DEFAULT_RELATIONS:
        continue
    dmap = ancestor_distances(direct_go)
    for t in set(dmap).intersection(unique_selected):
        for j in columns_by_term[t]:
            for g in bits_to_genes(bits, genes):
                key = (j,t,g)
                if key not in remaining_keys:
                    continue
                method = accession_gene_method.get((accession,g),'')
                witness_rows.append({
                    'label_column':j,'selected_GO_ID':t,'selected_GO_name':name_by.get(t,''),'selected_namespace':namespace_by.get(t,''),'GeneID':g,
                    'UniProtKB_accession':accession,'mapping_method':method,'direct_GO_ID':direct_go,'direct_GO_name':name_by.get(direct_go,''),
                    'direct_GO_min_depth':min_depth(direct_go),'selected_GO_min_depth':min_depth(t),'shortest_is_a_distance_to_selected':dmap[t],
                    'evidence':evidence,'normalized_relation':relation,'qualifier':qualifier,'date':date,'assigned_by':source,'reference':ref,
                    'direct_equals_selected':int(direct_go==t),
                })
                s=pair_support[key]
                s['direct_terms'].add(direct_go); s['distances'].append(dmap[t]); s['relations'].add(relation); s['accessions'].add(accession); s['mapping_methods'].add(method); s['sources'].add(source); s['dates'].add(date); s['evidence'].add(evidence); s['references'].add(ref)
write_gz_csv(ANA / f'B104A_remaining_13_witness_annotation_rows_{STAMP}.csv.gz', witness_rows)

pair_detail_rows = []
for r in remaining13:
    key=(r['label_column'],r['GO_ID'],r['GeneID']); s=pair_support[key]
    t=r['GO_ID']; g=r['GeneID']
    immediate_children=set()
    for d in s['direct_terms']:
        if d==t: continue
        immediate_children |= children.get(t,set()).intersection(ancestors(d))
    pair_detail_rows.append({
        **r,
        'selected_min_depth':min_depth(t),'selected_max_depth':max_depth(t),'selected_descendant_count':len(descendants(t)),
        'witness_annotation_rows':sum(1 for x in witness_rows if x['label_column']==r['label_column'] and x['GO_ID'] if False),
        'distinct_direct_GO_terms':len(s['direct_terms']),'direct_GO_terms':'|'.join(sorted(s['direct_terms'])),
        'minimum_support_distance':min(s['distances']) if s['distances'] else '', 'maximum_support_distance':max(s['distances']) if s['distances'] else '',
        'has_direct_selected_term_annotation':int(0 in s['distances']),
        'immediate_child_branches_to_selected':len(immediate_children),'immediate_child_GO_IDs':'|'.join(sorted(immediate_children)),
        'accessions':'|'.join(sorted(s['accessions'])),'mapping_methods':'|'.join(sorted(s['mapping_methods'])),'evidence_codes':'|'.join(sorted(s['evidence'])),'assigned_by':'|'.join(sorted(s['sources'])),'annotation_dates':'|'.join(sorted(s['dates'])),'references':'|'.join(sorted(s['references'])),
        'can_ontology_is_a_edge_drift_alone_remove_pair':int(0 not in s['distances']),
    })
# Correct witness counts separately.
counts=collections.Counter((x['label_column'],x['selected_GO_ID'],x['GeneID']) for x in witness_rows)
for x in pair_detail_rows:
    x['witness_annotation_rows']=counts[(x['label_column'],x['GO_ID'],x['GeneID'])]
write_csv(ANA / f'B104A_remaining_13_pair_details_{STAMP}.csv', pair_detail_rows)

# ---------- immediate-parent-edge drift sensitivity for the remaining 13 ----------
# Build direct term support under original evidence/default relations for every gene and selected term.
gene_term_directs = collections.defaultdict(set)
for accession, bits, direct_go, evidence, relation, date, source, ref, qualifier in annotations:
    if evidence not in ORIGINAL_EVIDENCE or relation not in DEFAULT_RELATIONS:
        continue
    dmap=ancestor_distances(direct_go)
    for t in set(dmap).intersection(unique_selected):
        for g in bits_to_genes(bits,genes):
            gene_term_directs[(g,t)].add(direct_go)

edge_impact_rows=[]
edge_optimization_rows=[]
for t in sorted({r['GO_ID'] for r in remaining13}):
    js=columns_by_term[t]
    # duplicate columns have identical observed vectors; use first for optimization and record all columns.
    j=js[0]; o=observed_bits[j]; p=pred_default[t]
    child_set=children.get(t,set())
    support_by_gene={}
    for g in bits_to_genes(p,genes):
        directs=gene_term_directs[(g,t)]
        direct_t=t in directs
        branches=set()
        for d in directs:
            if d!=t:
                branches |= child_set.intersection(ancestors(d))
        support_by_gene[g]=(direct_t,branches)
    candidate_edges=sorted(set().union(*(branches for g,(direct_t,branches) in support_by_gene.items() if not (o & (1<<gene_index[g])))) if support_by_gene else set())
    for c in candidate_edges:
        fp_removed=tp_lost=0
        for g,(direct_t,branches) in support_by_gene.items():
            before=True
            after=direct_t or bool(branches-{c})
            if before and not after:
                if o & (1<<gene_index[g]): tp_lost+=1
                else: fp_removed+=1
        edge_impact_rows.append({'selected_GO_ID':t,'selected_GO_name':name_by.get(t,''),'child_GO_ID':c,'child_GO_name':name_by.get(c,''),'edge':'%s is_a %s'%(c,t),'false_positives_removed_if_only_this_edge_deleted':fp_removed,'true_positives_lost_if_only_this_edge_deleted':tp_lost,'net_mismatch_change':tp_lost-fp_removed})
    # Exact subset search over candidate edges when tractable; candidate set concerns only branches used by remaining FPs.
    if len(candidate_edges)<=20:
        best=None
        for mask in range(1<<len(candidate_edges)):
            removed={candidate_edges[i] for i in range(len(candidate_edges)) if mask&(1<<i)}
            fp=fn=0
            for g,(direct_t,branches) in support_by_gene.items():
                after=direct_t or bool(branches-removed)
                observed=bool(o & (1<<gene_index[g]))
                fp += after and not observed
                fn += observed and not after
            score=(fp+fn,fn,fp,len(removed),'|'.join(sorted(removed)))
            if best is None or score<best[0]: best=(score,removed)
        score,removed=best
        edge_optimization_rows.append({'selected_GO_ID':t,'selected_GO_name':name_by.get(t,''),'label_columns':'|'.join(map(str,js)),'candidate_immediate_edges_considered':len(candidate_edges),'search_method':'exhaustive_all_subsets','minimum_mismatches_after_immediate_edge_deletions':score[0],'false_negatives_at_optimum':score[1],'false_positives_at_optimum':score[2],'edges_deleted_at_optimum':score[3],'deleted_child_GO_IDs':'|'.join(sorted(removed))})
    else:
        edge_optimization_rows.append({'selected_GO_ID':t,'selected_GO_name':name_by.get(t,''),'label_columns':'|'.join(map(str,js)),'candidate_immediate_edges_considered':len(candidate_edges),'search_method':'not_exhaustive_candidate_count_gt_20','minimum_mismatches_after_immediate_edge_deletions':'','false_negatives_at_optimum':'','false_positives_at_optimum':'','edges_deleted_at_optimum':'','deleted_child_GO_IDs':''})
write_csv(ANA / f'B104A_immediate_is_a_edge_single_deletion_impacts_{STAMP}.csv', edge_impact_rows)
write_csv(ANA / f'B104A_immediate_is_a_edge_subset_optimization_{STAMP}.csv', edge_optimization_rows)

# ---------- topology / feature properties for historically unmapped and residual genes ----------
row_map_df = pd.read_csv(ROW_MAP)
tissue_df = pd.read_csv(TISSUES)
# GraphSAGE row IDs are contiguous and correspond to feature rows.
with zipfile.ZipFile(GRAPHSAGE_ZIP) as zf:
    G = json.load(zf.open('ppi/ppi-G.json'))
    feat_bytes = io.BytesIO(zf.read('ppi/ppi-feats.npy'))
    feats = np.load(feat_bytes)
row_count = len(G['nodes'])
assert row_count == len(row_map_df) == feats.shape[0]
degrees = np.zeros(row_count, dtype=np.int64)
for e in G['links']:
    s=int(e['source']); t=int(e['target'])
    if s==t:
        degrees[s]+=2
    else:
        degrees[s]+=1; degrees[t]+=1
row_tissue=['']*row_count; row_split=['']*row_count; row_graph=[0]*row_count
for r in tissue_df.itertuples(index=False):
    for i in range(int(r.row_start_inclusive),int(r.row_end_exclusive)):
        row_tissue[i]=str(r.ohmnet_tissue); row_split[i]=str(r.split); row_graph[i]=int(r.graph_index)
row_map_df['degree']=degrees
row_map_df['feature_nonzero_count']=(feats!=0).sum(axis=1)
row_map_df['tissue']=row_tissue
row_map_df['split']=row_split
row_map_df['graph_index']=row_graph
row_map_df=row_map_df.dropna(subset=['entrez_gene_id']).copy(); row_map_df['entrez_gene_id']=row_map_df.entrez_gene_id.astype(int)
label_count_by_gene={int(r.entrez_gene_id):sum(int(getattr(r,f'label_{j}')) for j in range(121)) for r in labels_df.itertuples(index=False)}

gene_stats=[]
for g,sub in row_map_df.groupby('entrez_gene_id'):
    fcounts=sorted(set(int(x) for x in sub.feature_nonzero_count))
    gene_stats.append({'GeneID':int(g),'graph_instance_count':len(sub),'unique_tissue_count':sub.tissue.nunique(),'train_instance_count':int((sub.split=='train').sum()),'valid_instance_count':int((sub.split=='valid').sum()),'test_instance_count':int((sub.split=='test').sum()),'mean_instance_degree':float(sub.degree.mean()),'median_instance_degree':float(sub.degree.median()),'max_instance_degree':int(sub.degree.max()),'total_instance_degree':int(sub.degree.sum()),'feature_nonzero_counts':'|'.join(map(str,fcounts)),'feature_nonzero_count_consistent':int(len(fcounts)==1),'positive_label_count':label_count_by_gene.get(int(g),0),'all_zero_label_vector':int(label_count_by_gene.get(int(g),0)==0)})
gene_stats_df=pd.DataFrame(gene_stats)
write_csv(DER / f'B104A_all_resolved_gene_topology_feature_label_properties_{STAMP}.csv', gene_stats)

missing_df=pd.read_csv(MISSING_HIST)
missing_col='Entrez_GeneID' if 'Entrez_GeneID' in missing_df.columns else 'GeneID'
historical_gap_genes=set(missing_df[missing_col].astype(int))
remaining_residual_genes={int(r['GeneID']) for r in remaining13}
accepted_mapped_genes=set(map_df.GeneID.astype(int))
all_zero_genes=set(gene_stats_df.loc[gene_stats_df.all_zero_label_vector==1,'GeneID'].astype(int))
metrics_for_percentile=['graph_instance_count','unique_tissue_count','mean_instance_degree','max_instance_degree','total_instance_degree','positive_label_count']
property_rows=[]
for g in sorted(historical_gap_genes | remaining_residual_genes | ({10159} if 10159 in set(gene_stats_df.GeneID) else set())):
    hit=gene_stats_df[gene_stats_df.GeneID==g]
    if hit.empty: continue
    r=hit.iloc[0].to_dict()
    groups=[]
    if g in historical_gap_genes: groups.append('missing_from_historical_direct_GPI_projection')
    if g not in accepted_mapped_genes: groups.append('unmatched_after_component_aware_mapping')
    if g in remaining_residual_genes: groups.append('gene_in_remaining_13_default_relation_false_positives')
    if g in all_zero_genes: groups.append('all_zero_label_gene')
    out={'GeneID':g,'analysis_groups':'|'.join(groups)}
    out.update(r)
    for metric in metrics_for_percentile:
        vals=gene_stats_df[metric].astype(float).values
        out[f'{metric}_percentile_among_all_resolved_genes']=float(percentileofscore(vals,float(r[metric]),kind='mean'))
    property_rows.append(out)
write_csv(ANA / f'B104A_special_gene_topology_feature_label_properties_{STAMP}.csv', property_rows)

# Group comparisons for unmatched/historical gap and residual genes.
group_compare=[]
for group_name, gset in [
    ('nine_missing_from_historical_direct_GPI_projection',historical_gap_genes),
    ('one_unmatched_after_component_aware_mapping',set(gene_stats_df.GeneID)-accepted_mapped_genes),
    ('genes_in_remaining_13_default_relation_false_positives',remaining_residual_genes),
    ('all_zero_label_genes',all_zero_genes),
]:
    sub=gene_stats_df[gene_stats_df.GeneID.isin(gset)]
    rest=gene_stats_df[~gene_stats_df.GeneID.isin(gset)]
    for metric in metrics_for_percentile:
        x=sub[metric].astype(float).tolist(); y=rest[metric].astype(float).tolist()
        if not x or not y: continue
        test=mannwhitneyu(x,y,alternative='two-sided')
        group_compare.append({'group':group_name,'group_gene_count':len(sub),'comparison_gene_count':len(rest),'metric':metric,'group_median':float(np.median(x)),'comparison_median':float(np.median(y)),'mann_whitney_u':float(test.statistic),'two_sided_p_value':float(test.pvalue),'interpretation_caveat':'Small groups, especially the single unmatched gene, do not support population-level inference.'})
write_csv(ANA / f'B104A_special_gene_group_comparisons_{STAMP}.csv', group_compare)

# ---------- concise machine-readable summary ----------
# Summaries of key depth metrics.
def med(sub, col):
    vals=pd.to_numeric(sub[col],errors='coerce').dropna()
    return float(vals.median()) if len(vals) else None
all_nonexact=term_metrics[term_metrics.all_relation_exact==0]
all_exact=term_metrics[term_metrics.all_relation_exact==1]
ccmf=term_metrics[term_metrics.namespace!='biological_process']
ccmf_nonexact=ccmf[ccmf.all_relation_exact==0]
ccmf_exact=ccmf[ccmf.all_relation_exact==1]

summary={
    'definition':{
        'residual_gene_label_difference':'A gene-label pair for which the reconstructed membership differs from the deposited GraphSAGE binary label. Under the prior release-159 model all 901 were false positives: GOA plus mapping and propagation predicted 1 while GraphSAGE stored 0.',
        'largest_residual_terms':'The selected GO label columns with the greatest number of those false-positive gene-label pairs, not the GO terms with the largest annotation sets in general.',
    },
    'validated_prior_baseline':{k:v for k,v in eval_all_rel.items() if k!='per_column'},
    'qualifier_relation_result':{
        'default_relations_only':{k:v for k,v in eval_default.items() if k!='per_column'},
        'original_false_positive_pairs_eliminated_by_excluding_qualified_relations':sum(int(r['eliminated_when_qualifier_rows_excluded']) for r in fp_decomp),
        'original_false_positive_pairs_remaining':sum(int(r['supported_by_default_relation']) for r in fp_decomp),
        'interpretation':'Every deposited positive remains supported after excluding contributes_to and colocalizes_with rows. Those qualified rows account for 888 of the 901 prior false positives.'
    },
    'ontology_depth_and_breadth':{
        'all_relation_nonexact_term_median_min_depth':med(all_nonexact,'min_is_a_depth_from_namespace_root'),
        'all_relation_exact_term_median_min_depth':med(all_exact,'min_is_a_depth_from_namespace_root'),
        'CC_MF_nonexact_term_median_descendant_count':med(ccmf_nonexact,'descendant_count_including_self'),
        'CC_MF_exact_term_median_descendant_count':med(ccmf_exact,'descendant_count_including_self'),
        'caveat':'Broad/shallow CC and MF ancestors amplify qualified annotations, but broadness is not sufficient: broad Biological Process labels are exact because these two qualifiers are not BP relations.'
    },
    'evidence_search':{
        'all_evidence_codes_considered':ev_codes,
        'best_global':best_global,
        'best_by_namespace':best_aspect,
        'combined_best_by_namespace':{k:v for k,v in eval_aspect_best.items() if k!='per_column'},
    },
    'remaining_original_six_default_relation_pairs':{
        'count':len(remaining13),
        'direct_selected_term_pairs':sum(int(r['has_direct_selected_term_annotation']) for r in pair_detail_rows),
        'ancestor_only_pairs':sum(1-int(r['has_direct_selected_term_annotation']) for r in pair_detail_rows),
        'ontology_drift_can_in_principle_affect_only_ancestor_only_pairs':True,
    },
    'propagation_distance_best':sorted(distance_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives'],-r['exact_columns']))[:20],
    'source_leave_one_out_best':sorted(source_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives']))[:20],
    'date_cutoff_best':sorted(date_rows,key=lambda r:(r['total_mismatches'],r['false_negatives'],r['false_positives']))[:20],
    'unmatched_genes':{
        'historical_direct_projection_gap_gene_count':len(historical_gap_genes),
        'unmatched_after_component_aware_mapping':sorted(set(gene_stats_df.GeneID)-accepted_mapped_genes),
        'all_zero_gene_count':len(all_zero_genes),
    },
}
(ROOT / f'B104A_analysis_summary_{STAMP}.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')

# output hashes before report
output_rows=[]
for p in sorted(ROOT.rglob('*')):
    if p.is_file():
        output_rows.append({'relative_path':str(p.relative_to(ROOT)),'size_bytes':p.stat().st_size,'sha256':sha256_file(p)})
write_csv(ROOT / f'B104A_output_checksums_pre_report_{STAMP}.csv', output_rows)

print(json.dumps({
    'baseline_all_relations':{k:v for k,v in eval_all_rel.items() if k!='per_column'},
    'default_relations':{k:v for k,v in eval_default.items() if k!='per_column'},
    'qualified_direct_only':{k:v for k,v in ev_qdirect.items() if k!='per_column'},
    'best_global_evidence':best_global,
    'best_aspect_evidence':best_aspect,
    'combined_aspect_best':{k:v for k,v in eval_aspect_best.items() if k!='per_column'},
    'remaining13_count':len(remaining13),
    'remaining13_direct':sum(int(r['has_direct_selected_term_annotation']) for r in pair_detail_rows),
    'unmatched_after_component_mapping':sorted(set(gene_stats_df.GeneID)-accepted_mapped_genes),
}, indent=2))

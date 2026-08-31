#!/usr/bin/env python3
"""Independent, set-based validation of the full-human GO term selection result.

This intentionally avoids pandas, NumPy, and the integer-bitset implementation
used by the first analysis. It reads the normalized historical inputs with only
Python's standard library and represents memberships as ordinary Python sets.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, Set

ROOT = Path('/mnt/data/ppi_repro_corrected/batches/B104A_20260828T145842Z')
INP = ROOT / 'retained_inputs'
GAF = INP / 'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz'
MAP = INP / 'B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz'
TERMS = INP / 'B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz'
EDGES = INP / 'B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz'
EXACT = ROOT / 'analysis/B104A_exact_GO_terms_for_each_label_column_20260828T145842Z.csv'
OUT = Path('/mnt/data/ppi_repro_corrected/batches/B104B_20260828T160144Z/analysis/B104B_independent_global_term_selection_validation.json')

ACCEPTED_EVIDENCE = {'EXP', 'IDA', 'IEP', 'IGI', 'IMP', 'ISS'}
ACCEPTED_RELATIONS = {'involved_in', 'part_of', 'enables'}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def split_pipe(value: str) -> Iterable[str]:
    for item in value.split('|'):
        item = item.strip()
        if item:
            yield item


# Parse ontology term metadata and alternate IDs.
term_name: Dict[str, str] = {}
namespace: Dict[str, str] = {}
canonical: Dict[str, str] = {}
with gzip.open(TERMS, 'rt', newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        go = row['GO_ID']
        term_name[go] = row['name']
        namespace[go] = row['namespace']
        canonical[go] = go
        for alt in split_pipe(row.get('alt_ids', '')):
            canonical[alt] = go

# Parse is_a graph.
parents: Dict[str, Set[str]] = defaultdict(set)
with gzip.open(EDGES, 'rt', newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row.get('relation') == 'is_a':
            parents[row['child_GO_ID']].add(row['parent_GO_ID'])

# Iterative memoized closure to avoid sharing code with the recursive first implementation.
ancestor_cache: Dict[str, FrozenSet[str]] = {}

def ancestors_including_self(term: str) -> FrozenSet[str]:
    cached = ancestor_cache.get(term)
    if cached is not None:
        return cached
    seen = {term}
    stack = [term]
    while stack:
        child = stack.pop()
        for parent in parents.get(child, ()):
            if parent not in seen:
                seen.add(parent)
                stack.append(parent)
    result = frozenset(seen)
    ancestor_cache[term] = result
    return result

# Parse all historical GPI159-linked GeneID edges. Preserve many-to-many mappings.
accession_to_geneids: Dict[str, Set[int]] = defaultdict(set)
with gzip.open(MAP, 'rt', newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['in_GPI159'] != '1':
            continue
        accession_to_geneids[row['UniProtKB_accession']].add(int(row['GeneID']))

# Correct the historically inconsistent cross-symbol edge before restricting any universe.
accession_to_geneids['O95073'].discard(25788)
if not accession_to_geneids['O95073']:
    del accession_to_geneids['O95073']

all_geneids = set()
for genes in accession_to_geneids.values():
    all_geneids.update(genes)

# Build direct term memberships from the final accepted relation/evidence policy.
direct_members: Dict[str, Set[int]] = defaultdict(set)
accepted_rows = 0
with gzip.open(GAF, 'rt', newline='') as fh:
    for row in csv.DictReader(fh, delimiter='\t'):
        if row['Is_NOT'] == '1' or 'NOT' in set(split_pipe(row.get('Qualifier', ''))):
            continue
        if row['Evidence_Code'] not in ACCEPTED_EVIDENCE:
            continue
        if row['Normalized_Relation'] not in ACCEPTED_RELATIONS:
            continue
        geneids = accession_to_geneids.get(row['DB_Object_ID'])
        if not geneids:
            continue
        term = canonical.get(row['GO_ID'], row['GO_ID'])
        direct_members[term].update(geneids)
        accepted_rows += 1

# Propagate each direct membership set through is_a only.
propagated_members: Dict[str, Set[int]] = defaultdict(set)
for direct_term, genes in direct_members.items():
    for ancestor in ancestors_including_self(direct_term):
        propagated_members[ancestor].update(genes)

# Candidate IDs are the union of every exact GO ID for every deposited label column.
target_ids: Set[str] = set()
with EXACT.open(newline='') as fh:
    for row in csv.DictReader(fh):
        target_ids.update(split_pipe(row['exact_GO_IDs']))

# Deterministic ranking: decreasing full-human GeneID prevalence, then GO ID.
ranked = sorted(
    ((go, len(genes)) for go, genes in propagated_members.items() if genes),
    key=lambda item: (-item[1], item[0]),
)
top121 = {go for go, _ in ranked[:121]}
threshold1007 = {go for go, count in ranked if count >= 1007}

# Aspect quota sensitivity.
aspect_quota = set()
for aspect, quota in (
    ('biological_process', 85),
    ('cellular_component', 26),
    ('molecular_function', 10),
):
    aspect_ranked = [(go, count) for go, count in ranked if namespace.get(go) == aspect]
    aspect_quota.update(go for go, _ in aspect_ranked[:quota])

# Boundary rows and direct-count control.
direct_ranked = sorted(
    ((go, len(genes)) for go, genes in direct_members.items() if genes),
    key=lambda item: (-item[1], item[0]),
)
direct_top121 = {go for go, _ in direct_ranked[:121]}

report = {
    'implementation': 'independent_standard_library_set_based',
    'input_sha256': {str(p): sha256(p) for p in (GAF, MAP, TERMS, EDGES, EXACT)},
    'accepted_evidence_codes': sorted(ACCEPTED_EVIDENCE),
    'accepted_relations': sorted(ACCEPTED_RELATIONS),
    'mapping_correction': 'removed O95073 -> GeneID 25788 while retaining other many-to-many historical edges',
    'historical_GPI159_linked_geneids': len(all_geneids),
    'accepted_GAF_rows': accepted_rows,
    'direct_terms_with_members': len(direct_members),
    'propagated_terms_with_members': len(propagated_members),
    'target_exact_GO_ID_union_size': len(target_ids),
    'top121_size': len(top121),
    'top121_intersection_size': len(top121 & target_ids),
    'top121_extra_ids': sorted(top121 - target_ids),
    'top121_missing_ids': sorted(target_ids - top121),
    'top121_exact_set_match': top121 == target_ids,
    'threshold_1007_size': len(threshold1007),
    'threshold_1007_exact_set_match': threshold1007 == target_ids,
    'aspect_quota_size': len(aspect_quota),
    'aspect_quota_exact_set_match': aspect_quota == target_ids,
    'direct_top121_intersection_size': len(direct_top121 & target_ids),
    'direct_top121_exact_set_match': direct_top121 == target_ids,
    'rank_118_to_124': [
        {
            'rank': i + 1,
            'GO_ID': go,
            'GO_name': term_name.get(go, ''),
            'namespace': namespace.get(go, ''),
            'global_gene_count': count,
            'is_target': go in target_ids,
        }
        for i, (go, count) in enumerate(ranked)
        if 117 <= i <= 123
    ],
}

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open('w') as fh:
    json.dump(report, fh, indent=2, sort_keys=True)
    fh.write('\n')

print(json.dumps(report, indent=2, sort_keys=True))
if not report['top121_exact_set_match']:
    raise SystemExit('Independent top-121 validation did not match the candidate GO term set')

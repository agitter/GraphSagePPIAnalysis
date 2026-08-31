#!/usr/bin/env python3
from __future__ import annotations
import csv, gzip, hashlib, json, os
from pathlib import Path
import pandas as pd

BATCH_STAMP='20260828T145842Z'; FINAL_STAMP='20260828T152026Z'
ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/B104A_{BATCH_STAMP}')
RESULTS=Path('/mnt/data/ppi_repro_corrected/results')
checks=[]
def ck(name,passed,details=''): checks.append({'check':name,'passed':bool(passed),'details':details})

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while True:
   b=f.read(1<<20)
   if not b:break
   h.update(b)
 return h.hexdigest()

exact=json.loads((ROOT/f'B104A_EXACT_RECONSTRUCTION_VALIDATION_{BATCH_STAMP}.json').read_text())
ind=json.loads((ROOT/f'B104A_INDEPENDENT_SET_VALIDATION_{BATCH_STAMP}.json').read_text())
ck('primary_exact_121', exact['cell_mismatch_count']==0 and exact['columns_with_at_least_one_exact_GO_term']==121, json.dumps({k:exact[k] for k in ['cell_mismatch_count','columns_with_at_least_one_exact_GO_term']}))
ck('matrix_hashes_identical', exact['matrix_hashes_identical'], exact['deposited_matrix_uint8_C_order_sha256'])
ck('independent_set_validation', ind['all_121_columns_exact'] and ind['total_false_positives']==0 and ind['total_false_negatives']==0 and ind['matrix_hashes_identical'])

factor=pd.read_csv(ROOT/f'analysis/B104A_two_factor_mapping_qualifier_reconstruction_{BATCH_STAMP}.csv')
final=factor[(factor.mapping_policy=='symbol_resolved_full_component_projection')&(factor.relation_policy=='default_aspect_relations_only')]
ck('factorial_final_zero',len(final)==1 and int(final.iloc[0].total_mismatches)==0 and int(final.iloc[0].exact_columns)==121)
ck('factorial_baseline_901',int(factor[(factor.mapping_policy=='original_component_aware_projection')&(factor.relation_policy=='all_positive_GAF_relations')].iloc[0].total_mismatches)==901)

release=pd.read_csv(ROOT/f'analysis/B104A_release158_159_final_policy_comparison_{BATCH_STAMP}.csv')
ck('release159_temporal_specificity',int(release[release.GOA_release==159].iloc[0].total_mismatches)==0 and int(release[release.GOA_release==158].iloc[0].total_mismatches)==846)

component=pd.read_csv(ROOT/f'analysis/B104A_O95073_Q9Y620_historical_mapping_component_{BATCH_STAMP}.csv')
needed={(25788,'O95073'),(25788,'Q9Y620'),(100861412,'O95073')}
seen={(int(r.GeneID),r.UniProtKB_accession) for r in component.itertuples(index=False)}
ck('full_mapping_component_present',needed.issubset(seen),str(sorted(seen)))

man_path=RESULTS/f'actual_input_file_manifest_through_B104A_{FINAL_STAMP}.csv'
ledger_path=RESULTS/f'source_ledger_through_B104A_{FINAL_STAMP}.csv'
events_path=RESULTS/f'provenance_events_through_B104A_{FINAL_STAMP}.csv'
for p in [man_path,ledger_path,events_path]: ck(f'exists_{p.name}',p.exists(),str(p))
man=pd.read_csv(man_path).fillna('')
ret=man[man.record_type=='retained_input']
missing=[];bad=[]
for r in ret.itertuples(index=False):
 p=Path(r.local_path)
 if not p.exists(): missing.append(str(p))
 elif sha(p)!=r.sha256: bad.append(str(p))
ck('retained_input_paths_and_hashes',not missing and not bad,f'missing={missing}; bad_hash={bad}')

events=pd.read_csv(events_path).fillna('')
for n in ['goa_human.gaf.158.gz','goa_human.gpa.158.gz','goa_human.gpi.158.gz']:
 q=events[(events.batch_id=='B104')&(events.event_type=='user_deletion_confirmed')&(events.artifact_name==n)]
 ck(f'B104_deletion_confirmed_{n}',len(q)>=1)

# Logical deletion discipline: no B104A script reads deleted raw paths.
forbidden=['/mnt/data/goa_human.gaf.158.gz','/mnt/data/goa_human.gpa.158.gz','/mnt/data/goa_human.gpi.158.gz','/mnt/data/2016-06-01-go.obo','/mnt/data/2016-06-01-gp2protein.geneid.gz']
hits=[]
for p in (ROOT/'scripts').glob('*.py'):
 if p.name == 'validate_B104A_delivery.py':
  continue
 txt=p.read_text()
 for x in forbidden:
  if x in txt:hits.append(f'{p.name}:{x}')
ck('no_deleted_raw_paths_in_B104A_scripts',not hits,str(hits))

report=(ROOT/f'B104A_REPORT_{BATCH_STAMP}.md').read_text()
for phrase in ['What “residual” meant','Ontology depth and breadth','Uniform qualifier handling','121 / 121','Revised assessment of ontology drift','Remaining unmatched']:
 ck(f'report_contains_{phrase}',phrase in report)

# Matrix derivative re-read and 4,268 rows.
matrix_path=ROOT/f'derived/B104A_reconstructed_exact_label_matrix_{BATCH_STAMP}.csv.gz'
with gzip.open(matrix_path,'rt',encoding='utf-8',newline='') as fh:
 r=csv.reader(fh); header=next(r); count=sum(1 for _ in r)
ck('matrix_derivative_shape',count==4268 and len(header)==122,f'rows={count}; columns={len(header)}')

out={'batch_id':'B104A','generated_at_utc':'2026-08-28T15:20:26Z','all_checks_passed':all(c['passed'] for c in checks),'checks':checks}
(ROOT/f'B104A_DELIVERY_VALIDATION_{FINAL_STAMP}.json').write_text(json.dumps(out,indent=2))
if not out['all_checks_passed']:
 print(json.dumps(out,indent=2)); raise SystemExit(1)
print(json.dumps({'all_checks_passed':True,'check_count':len(checks)},indent=2))

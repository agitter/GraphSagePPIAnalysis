#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, gzip, hashlib, json, zipfile
from pathlib import Path
import pandas as pd

EXPECTED_RAW_SHA='d372fc23f229cbb79656d824e0519587db6110963d22d1f4c95e5154963a32d2'
EXPECTED_RAW_SIZE=115484475
EXPECTED_DER_SHA='6775fb96be44080c768bd5789a0dbb0c802a1a0faa45927aa2a07d70af9f7c1f'
EXPECTED_DER_SIZE=2982724

def sha(p:Path):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--stage',type=Path,required=True)
 ap.add_argument('--raw',type=Path,default=Path('/mnt/data/msigdb_v5.0_files_to_download_locally.zip'))
 ap.add_argument('--bundle',type=Path)
 ap.add_argument('--output',type=Path,required=True)
 args=ap.parse_args(); S=args.stage
 checks={}
 def check(name,value):
  checks[name]=bool(value)
  if not value: raise AssertionError(name)
 check('stage_exists',S.is_dir())
 for name in [
  'B104C_MSIGDB50_COLUMN_ORDER_REPORT_20260828T200948Z.md',
  'B104C_EXECUTION_DIAGNOSTICS_20260828T200948Z.md',
  'B104C_DELETION_CLEARANCE_20260828T200948Z.md',
  'B104C_UNIPROT_AUDIT_INSTRUCTIONS_20260828T200948Z.md']:
  check('exists_'+name,(S/name).is_file())
 check('raw_size',args.raw.stat().st_size==EXPECTED_RAW_SIZE)
 check('raw_sha256',sha(args.raw)==EXPECTED_RAW_SHA)
 der=S/'retained_inputs/B104C_msigdb_v5.0_normalized_entrez_gene_sets_20260828T194921Z.tsv.gz'
 check('derivative_size',der.stat().st_size==EXPECTED_DER_SIZE)
 check('derivative_sha256',sha(der)==EXPECTED_DER_SHA)
 with gzip.open(der,'rt',encoding='utf-8') as f:
  header=f.readline().rstrip('\n').split('\t'); rows=0;c5=0;msum=0
  for line in f:
   rows+=1; p=line.rstrip('\n').split('\t');
   if p[4]=='C5':c5+=1
   if p[9]:msum+=len(p[9].split('|'))
 check('normalized_header',header==['order_all','order_C5','standard_name','systematic_name','category','subcategory','GO_ID','chip','external_details_url','member_Entrez_IDs'])
 check('normalized_rows',rows==10348)
 check('normalized_C5_rows',c5==1454)
 check('normalized_membership_sum',msum==1315074)
 direct=pd.read_csv(S/'analysis/B104C_msigdb_direct_membership_summary_20260828T194921Z.csv')
 v50=direct[direct.version==5.0]
 check('v50_two_scopes',set(v50.scope)=={'C5','ALL'})
 check('v50_no_exact',(v50.exact_label_columns==0).all())
 check('v50_closest_307',(v50.closest_mismatch==307).all())
 check('v50_median_823',(v50.median_best_mismatch==823).all())
 cross=pd.read_csv(S/'analysis/B104C_msigdb_cross_version_C5_target_presence_20260828T194921Z.csv')
 got=dict(zip(cross.version.astype(str),cross.target_GO_IDs_present_in_C5.astype(int)))
 check('corrected_C5_presence',got=={'5.0':57,'5.1':57,'5.2':6,'6.0':6})
 cm=pd.read_csv(S/'analysis/B104C_graphsage_class_map_python2_dict_validation_20260828T194921Z.csv')
 r64=cm[cm['bits']==64].iloc[0]
 check('class_map_56944_exact',int(r64.exact_positions)==56944 and int(r64.node_keys)==56944 and int(r64.complete_sequence_exact)==1)
 ext=json.loads((S/'analysis/B104C_extended_order_simulation_summary_20260828T194921Z.json').read_text())
 check('extended_models_48',ext['simulations_scored']==48)
 check('extended_orientation_unanimous',ext['all_models_duplicate_orientation_counts']=={'001':48})
 check('best_lcs_94',ext['best']['lcs']==94)
 check('best_table_32768',ext['best']['table_size']==32768)
 order=pd.read_csv(S/'analysis/B104C_inferred_unique_121_GO_column_order_20260828T194921Z.csv')
 check('unique_order_rows',len(order)==121)
 check('unique_order_distinct_GO',order.inferred_GO_ID.nunique()==121)
 expected={24:'GO:0043228',71:'GO:0043232',39:'GO:0006464',63:'GO:0036211',48:'GO:1903561',70:'GO:0043230'}
 check('duplicate_assignments',all(order.loc[order.label_column==c,'inferred_GO_ID'].iloc[0]==g for c,g in expected.items()))
 # Script test logs.
 check('uniprot_self_test_passed','"self_test": "passed"' in (S/'logs/uniprot_audit_self_test.stdout').read_text())
 check('uniprot_dry_run_three_releases',sum(1 for x in (S/'logs/uniprot_audit_dry_run.stdout').read_text().splitlines() if x.startswith('[2016_'))==3)
 # Raw archives must not be in staging.
 bad=[p for p in S.rglob('*') if p.is_file() and (p.name=='msigdb_v5.0_files_to_download_locally.zip' or p.name.startswith('uniprot_sprot-only'))]
 check('stage_excludes_large_raw_archives',not bad)
 bundle_info={}
 if args.bundle:
  check('bundle_exists',args.bundle.is_file())
  with zipfile.ZipFile(args.bundle) as z:
   check('bundle_zip_test',z.testzip() is None)
   names=z.namelist()
   check('bundle_excludes_msigdb50_raw',not any(n.endswith('/msigdb_v5.0_files_to_download_locally.zip') or n=='msigdb_v5.0_files_to_download_locally.zip' for n in names))
   check('bundle_excludes_uniprot_raw',not any('uniprot_sprot-only' in n and n.endswith('.tar.gz') for n in names))
   check('bundle_contains_report',any(n.endswith('B104C_MSIGDB50_COLUMN_ORDER_REPORT_20260828T200948Z.md') for n in names))
  bundle_info={'path':str(args.bundle),'size_bytes':args.bundle.stat().st_size,'sha256':sha(args.bundle)}
 out={'validation':'passed','checks':checks,'bundle':bundle_info,'stage':str(S),'raw_input_sha256':sha(args.raw),'normalized_derivative_sha256':sha(der)}
 args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
 print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()

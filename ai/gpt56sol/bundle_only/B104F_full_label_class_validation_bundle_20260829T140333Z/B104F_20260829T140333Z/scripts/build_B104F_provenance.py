#!/usr/bin/env python3
from __future__ import annotations
import csv,hashlib,datetime,json
from pathlib import Path
import pandas as pd
STAMP='20260829T140333Z';NOW='2026-08-29T14:03:33Z';BATCH='B104F'
BASE=Path('/mnt/data');ROOT=BASE/f'ppi_repro_corrected/batches/{BATCH}_{STAMP}';RES=BASE/'ppi_repro_corrected/results'
ACT_IN=RES/'actual_input_file_manifest_through_B104E_20260829T122734Z.csv'
LED_IN=RES/'source_ledger_through_B104E_FINAL_20260829T122734Z.csv'
EV_IN=RES/'provenance_events_through_B104E_FINAL_20260829T122734Z.csv'
ACT_OUT=RES/f'actual_input_file_manifest_through_{BATCH}_{STAMP}.csv'
LED_OUT=RES/f'source_ledger_through_{BATCH}_PREBUNDLE_{STAMP}.csv'
EV_OUT=RES/f'provenance_events_through_{BATCH}_PREBUNDLE_{STAMP}.csv'

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def blank(cols):return {c:'' for c in cols}
# Current-state actual-input manifest
act=pd.read_csv(ACT_IN,dtype=str).fillna('')
for name in ['uniprot_2016_mapping_audit_ledger.zip','goa_date_screen_results.zip']:
 m=act.artifact_name.eq(name)
 act.loc[m,'user_deletion_confirmed_at_utc']=NOW
 act.loc[m,'raw_retention_status']='conversation_attachment_deleted_by_user; compact retained derivatives preserved'
 act.loc[m,'raw_available_in_conversation']='no_user_confirmed_deleted'
 act.loc[m,'deletion_notes']='User reported Deleted B104E on 2026-08-29; exact compact members and provenance remain retained.'
audit=pd.read_csv(ROOT/'analysis/B104F_uploaded_dhimmel_file_audit.csv',dtype=str).fillna('')
for r in audit.to_dict('records'):
 row=blank(act.columns)
 row.update({
  'record_type':'actual_input_rejected', 'artifact_name':r['artifact_name'],'local_path':r['local_path'],
  'local_status':'present_verified_wrong_media_type','origin_in_this_run':'user_upload',
  'analysis_role':'intended Entrez-native GO annotation control; rejected because bytes are GitHub HTML, not TSV',
  'used_by':'B104F input audit only','direct_or_canonical_source_url':r['commit_pinned_raw_url'],
  'source_page_url':r['commit_pinned_raw_url'].replace('raw.githubusercontent.com','github.com').replace(f'/962a5e12f8590400c2891cde93fd6a783b26e02e/','/blob/962a5e12f8590400c2891cde93fd6a783b26e02e/'),
  'url_status':'commit-pinned canonical raw URL recorded; uploaded bytes do not match it',
  'retrieval_status':'uploaded_by_user_as_saved_GitHub_HTML_page','retrieved_at_utc':NOW,
  'size_bytes':r['size_bytes'],'sha256':r['sha256'],'parent_or_derivation':'saved GitHub blob HTML page',
  'notes':r['notes'],'batch_id':BATCH,'received_at_utc':NOW,'integrity_verified_at_utc':NOW,
  'analysis_completed_at_utc':NOW,'deletion_clearance_issued_at_utc':NOW,
  'raw_retention_status':'rejected_wrong_media_type; audit retained; safe to delete conversation copy',
  'raw_available_in_conversation':'yes_pending_user_deletion','retained_derivative_paths':str(ROOT/'analysis/B104F_uploaded_dhimmel_file_audit.csv'),
  'derivative_sha256s':sha(ROOT/'analysis/B104F_uploaded_dhimmel_file_audit.csv'),
  'parser_version_or_script_sha256':sha(ROOT/'scripts/audit_uploaded_dhimmel_files.py'),
  'raw_to_derived_reconciliation':'HTML signature and title recorded; zero annotation rows parsed',
  'reacquisition_url':r['commit_pinned_raw_url'],'deletion_notes':'Delete uploaded conversation copy; replace local file with validated raw download before analysis.'
 })
 act=pd.concat([act,pd.DataFrame([row])],ignore_index=True)
act.to_csv(ACT_OUT,index=False)
# Append-only source ledger
led=pd.read_csv(LED_IN,dtype=str).fillna('');new=[]
for name in ['uniprot_2016_mapping_audit_ledger.zip','goa_date_screen_results.zip']:
 row=blank(led.columns);row.update({'record_type':'user_deletion_confirmation','artifact_name':name,'local_status':'conversation_attachment_deleted_by_user','origin_in_this_run':'user_statement','analysis_role':'append-only retention-state transition','used_by':'future provenance','notes':'User reported Deleted B104E. Compact retained derivatives remain available.','batch_id':BATCH,'deletion_state':'deletion_confirmed_by_user','event_recorded_at_utc':NOW,'hash_authority':'prior verified input record','runtime_verification_status':'raw conversation copy treated as unavailable','user_deletion_confirmed_at_utc':NOW});new.append(row)
for r in audit.to_dict('records'):
 row=blank(led.columns);row.update({'record_type':'actual_input_rejected','artifact_name':r['artifact_name'],'local_path':r['local_path'],'local_status':'present_wrong_media_type','origin_in_this_run':'user_upload','analysis_role':'intended historical Entrez-native GO annotation; rejected','used_by':'B104F audit only','direct_or_canonical_source_url':r['commit_pinned_raw_url'],'source_page_url':r['commit_pinned_raw_url'].replace('raw.githubusercontent.com','github.com').replace('/962a5e12f8590400c2891cde93fd6a783b26e02e/','/blob/962a5e12f8590400c2891cde93fd6a783b26e02e/'),'url_status':'canonical raw URL recorded','retrieval_status':'user uploaded HTML page','retrieved_at_utc':NOW,'size_bytes':r['size_bytes'],'sha256':r['sha256'],'parent_or_derivation':'GitHub blob HTML page','notes':r['notes'],'batch_id':BATCH,'deletion_state':'clearance_issued','event_recorded_at_utc':NOW,'hash_authority':'runtime SHA-256','runtime_verification_status':'HTML detected; TSV validation failed','planned_batch':'future raw re-upload optional','retained_derivative_paths':str(ROOT/'analysis/B104F_uploaded_dhimmel_file_audit.csv'),'parser_script_sha256':sha(ROOT/'scripts/audit_uploaded_dhimmel_files.py'),'raw_to_derived_reconciliation':'zero annotation rows; exact HTML hash and source title retained'});new.append(row)
# Generated/retained B104F artifacts except logs
for p in sorted(ROOT.rglob('*')):
 if not p.is_file() or '/logs/' in str(p):continue
 rel=str(p.relative_to(ROOT))
 if rel.startswith('retained_inputs/'):
  rt='retained_input';role='compact prior derivative or B104F symbol crosswalk required for reproduction'
 elif rel.startswith('scripts/'):
  rt='analysis_script';role='B104F reproducibility script'
 elif rel.startswith('analysis/'):
  rt='analysis_output';role='B104F scientific analysis output'
 else:
  rt='batch_report_or_validation';role='B104F report, response, diagnostics, or validation'
 row=blank(led.columns);row.update({'record_type':rt,'artifact_name':p.name,'local_path':str(p),'local_status':'present_verified','origin_in_this_run':'generated_or_retained_by_B104F','analysis_role':role,'used_by':'B104F delivery and future reproduction','retrieval_status':'generated_locally' if rt!='retained_input' else 'copied_from_verified_retained_derivative','retrieved_at_utc':NOW,'size_bytes':str(p.stat().st_size),'sha256':sha(p),'parent_or_derivation':'B104F fixed-policy full-class analysis and prior frozen derivatives','batch_id':BATCH,'deletion_state':'retained','event_recorded_at_utc':NOW,'hash_authority':'runtime SHA-256','runtime_verification_status':'present','retained_derivative_paths':str(p),'parser_script_sha256':sha(ROOT/'scripts/analyze_full_label_class_multisets.py') if rt=='analysis_output' else '','raw_to_derived_reconciliation':'see B104F report and diagnostics'});new.append(row)
led=pd.concat([led,pd.DataFrame(new)],ignore_index=True);led.to_csv(LED_OUT,index=False)
# Events
ev=pd.read_csv(EV_IN,dtype=str).fillna('');events=[]
def add(batch,event,artifact,status,details):events.append({'event_time_utc':NOW,'batch_id':batch,'event_type':event,'artifact_name':artifact,'status':status,'details':details})
for n in ['uniprot_2016_mapping_audit_ledger.zip','goa_date_screen_results.zip']:add('B104E','user_deletion_confirmed',n,'deleted_from_conversation','User reported Deleted B104E; compact retained members remain.')
for r in audit.to_dict('records'):
 add(BATCH,'input_received',r['artifact_name'],'received',f"bytes={r['size_bytes']}; sha256={r['sha256']}")
 add(BATCH,'input_validation',r['artifact_name'],'rejected_wrong_media_type','HTML GitHub blob page detected; no TSV rows parsed.')
 add(BATCH,'deletion_clearance_issued',r['artifact_name'],'safe_to_delete_conversation_attachment','Audit and exact SHA-256 retained; local file should be replaced by commit-pinned raw TSV.')
add(BATCH,'analysis_completed','full GraphSAGE label-class validation','passed','56,411 resolved rows exact; 183/183 unresolved classes exact as vector multisets; 438 unique and 95 ambiguous unresolved rows.')
add(BATCH,'independent_validation','NumPy packed-signature implementation','passed','Independent implementation agrees with primary full-class result.')
ev=pd.concat([ev,pd.DataFrame(events)],ignore_index=True);ev.to_csv(EV_OUT,index=False)
print(ACT_OUT);print(LED_OUT);print(EV_OUT)

#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path

STAMP='20260829T010311Z'
ISO='2026-08-29T01:03:11Z'
BATCH='B104D'
BASE=Path('/mnt/data/ppi_repro_corrected')
ROOT=BASE/f'batches/B104D_{STAMP}'
RESULTS=BASE/'results'
PREV_ACT=RESULTS/'actual_input_file_manifest_through_B104C_20260828T200948Z.csv'
PREV_LED=RESULTS/'source_ledger_through_B104C_FINAL_20260828T200948Z.csv'
PREV_EVT=RESULTS/'provenance_events_through_B104C_FINAL_20260828T200948Z.csv'
UNIPROT=Path('/mnt/data/uniprot_2016_mapping.zip')
SCREENSHOT=Path('/mnt/data/fbe92bb6-4a87-4f70-834e-2a920376c88b.png')
LEDGER=ROOT/'provenance/uniprot_2016_mapping_audit_ledger.csv'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def readcsv(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def writecsv(p,rows,fields):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)

def blank(fields):return {k:'' for k in fields}

# Current-state actual input manifest.
act=readcsv(PREV_ACT); af=list(act[0])
for r in act:
    if r.get('artifact_name')=='msigdb_v5.0_files_to_download_locally.zip' and r.get('record_type')=='actual_input':
        r['local_status']='logically_deleted_by_user_from_conversation_after_B104C_clearance'
        r['user_deletion_confirmed_at_utc']=ISO
        r['raw_retention_status']='normalized_complete_derivative_retained; user_deletion_confirmed'
        r['raw_available_in_conversation']='no_logically_deleted_by_user'
        r['deletion_notes']=(r.get('deletion_notes','')+' User confirmed `Deleted B104C`; residual runtime mounts, if any, are not treated as available inputs.').strip()

u=blank(af);u.update({
'record_type':'actual_input','artifact_name':UNIPROT.name,'local_path':str(UNIPROT),'local_status':'present_verified_ledger_only_package',
'origin_in_this_run':'user_upload','analysis_role':'user-side UniProt 2016 audit ledger package for O95073 and Q9Y620','used_by':'B104D package review',
'direct_or_canonical_source_url':'multiple official UniProt release URLs recorded in retained ledger','source_page_url':'https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/',
'url_status':'source URLs recorded per release in ledger','retrieval_status':'uploaded_by_user; parent archives downloaded and deleted on user system',
'retrieved_at_utc':'','size_bytes':str(UNIPROT.stat().st_size),'sha256':sha(UNIPROT),'parent_or_derivation':'output of download_extract_uniprot_2016_mapping_audit.py',
'notes':'ZIP contains only uniprot_2016_mapping_audit_ledger.csv; extracted DAT/TSV/provenance files are absent and remain requested.','batch_id':BATCH,
'received_at_utc':ISO,'integrity_verified_at_utc':ISO,'analysis_completed_at_utc':ISO,
'raw_retention_status':'ledger retained; uploaded wrapper ZIP not required after deletion clearance','raw_available_in_conversation':'yes_pending_user_deletion',
'retained_derivative_paths':str(LEDGER),'derivative_sha256s':sha(LEDGER),
'parser_version_or_script_sha256':'B104D package review; package_uniprot_audit_outputs.py:'+sha(ROOT/'scripts/package_uniprot_audit_outputs.py'),
'raw_to_derived_reconciliation':'ZIP has one member; retained ledger bytes and SHA-256 match the member exactly',
'reacquisition_url':'not applicable; regenerate from retained user-side audit output directory','deletion_notes':'Safe to delete conversation ZIP after B104D bundle validation; extracted output package still needed.'})
act.append(u)

s=blank(af);s.update({
'record_type':'diagnostic_input','artifact_name':SCREENSHOT.name,'local_path':str(SCREENSHOT),'local_status':'present_verified',
'origin_in_this_run':'user_upload','analysis_role':'screenshot of ChatGPT Library storage display','used_by':'B104D storage-limit explanation',
'retrieval_status':'uploaded_by_user','size_bytes':str(SCREENSHOT.stat().st_size),'sha256':sha(SCREENSHOT),'batch_id':BATCH,
'received_at_utc':ISO,'integrity_verified_at_utc':ISO,'analysis_completed_at_utc':ISO,
'raw_retention_status':'diagnostic only; no scientific dependency','raw_available_in_conversation':'yes',
'raw_to_derived_reconciliation':'visually inspected; no data extraction required','deletion_notes':'May be deleted at any time; not needed for reproduction.'})
act.append(s)

l=blank(af);l.update({
'record_type':'retained_input','artifact_name':LEDGER.name,'local_path':str(LEDGER),'local_status':'present_verified',
'origin_in_this_run':'extracted byte-for-byte from user-uploaded ZIP','analysis_role':'append-only user-side UniProt archive verification ledger','used_by':'B104D package review and future extracted-record audit',
'direct_or_canonical_source_url':'per-release URLs inside ledger','source_page_url':'https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/',
'url_status':'official URLs recorded in ledger','retrieval_status':'retained_from_verified_user_package','retrieved_at_utc':ISO,
'size_bytes':str(LEDGER.stat().st_size),'sha256':sha(LEDGER),'parent_or_derivation':UNIPROT.name+':'+sha(UNIPROT),
'notes':'Three successful releases; expected/observed byte sizes and MD5 values match; extracted record payloads not present in uploaded ZIP.','batch_id':BATCH,
'received_at_utc':ISO,'integrity_verified_at_utc':ISO,'analysis_completed_at_utc':ISO,
'raw_retention_status':'required compact provenance record','raw_available_in_conversation':'not_applicable_retained_derivative',
'retained_derivative_paths':str(LEDGER),'derivative_sha256s':sha(LEDGER),
'raw_to_derived_reconciliation':'exact ZIP member extraction; member count and bytes verified'})
act.append(l)
ACT_OUT=RESULTS/f'actual_input_file_manifest_through_B104D_{STAMP}.csv'
writecsv(ACT_OUT,act,af)

# Append-only source ledger.
led=readcsv(PREV_LED); lf=list(led[0])
def add_ledger(**kwargs):
    r=blank(lf);r.update({k:str(v) for k,v in kwargs.items()});led.append(r)
add_ledger(record_type='user_deletion_confirmation',artifact_name='msigdb_v5.0_files_to_download_locally.zip',local_status='logically_deleted_by_user',origin_in_this_run='user statement',analysis_role='state correction',batch_id=BATCH,user_local_status='deleted_from_user_local_storage_and_conversation_as_reported',deletion_state='user_deletion_confirmed',event_recorded_at_utc=ISO,supersedes_prior_record_type='B104C deletion clearance pending confirmation',hash_authority='prior verified SHA-256',runtime_verification_status='residual mount ignored',user_deletion_confirmed_at_utc=ISO,notes='User reported Deleted B104C. Earlier message also reported older MSigDB files deleted locally for space.')
add_ledger(record_type='user_uploaded_audit_package',artifact_name=UNIPROT.name,local_path=str(UNIPROT),local_status='present_verified_ledger_only',origin_in_this_run='user upload',analysis_role='UniProt 2016 parent-archive audit ledger',used_by='B104D',source_page_url='https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/',url_status='official URLs inside ledger',retrieval_status='uploaded_by_user',size_bytes=UNIPROT.stat().st_size,sha256=sha(UNIPROT),parent_or_derivation='download_extract_uniprot_2016_mapping_audit.py outputs',notes='Only ledger was included; extracted DAT/TSV/JSON outputs absent.',batch_id=BATCH,deletion_state='safe_to_delete_after_B104D_bundle',event_recorded_at_utc=ISO,hash_authority='container_sha256',runtime_verification_status='ZIP integrity and sole member verified',retained_derivative_paths=str(LEDGER),raw_to_derived_reconciliation='one ZIP member extracted byte-for-byte')
add_ledger(record_type='retained_input',artifact_name=LEDGER.name,local_path=str(LEDGER),local_status='present_verified',origin_in_this_run='extracted from user package',analysis_role='UniProt archive verification ledger',used_by='B104D and future mapping audit',direct_or_canonical_source_url='per-release official URLs in CSV',source_page_url='https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/',retrieval_status='retained_derivative_from_user_package',size_bytes=LEDGER.stat().st_size,sha256=sha(LEDGER),parent_or_derivation=UNIPROT.name+':'+sha(UNIPROT),notes='Parent archives verified by user-side script; record content cannot be independently inspected until small outputs are repackaged.',batch_id=BATCH,deletion_state='retain',event_recorded_at_utc=ISO,hash_authority='container_sha256',runtime_verification_status='present_and_verified')
add_ledger(record_type='diagnostic_input',artifact_name=SCREENSHOT.name,local_path=str(SCREENSHOT),local_status='present_verified',origin_in_this_run='user upload',analysis_role='Library storage screenshot',used_by='B104D storage explanation',retrieval_status='uploaded_by_user',size_bytes=SCREENSHOT.stat().st_size,sha256=sha(SCREENSHOT),batch_id=BATCH,deletion_state='not_required_for_reproduction',event_recorded_at_utc=ISO,hash_authority='container_sha256',runtime_verification_status='present_and_verified')
# Generated/retained B104D artifacts; bundle row is appended later after bundle creation.
for p,role in [
(ROOT/f'B104D_COLUMN_ORDER_DATE_RANGE_AND_UNIPROT_REPORT_{STAMP}.md','scientific report'),
(ROOT/f'B104D_GOA_DATE_SCREEN_INSTRUCTIONS_{STAMP}.md','date-screen instructions'),
(ROOT/f'B104D_UNIPROT_REPACK_INSTRUCTIONS_{STAMP}.md','UniProt repack instructions'),
(ROOT/f'B104D_GOA_date_screen_reference_pack_{STAMP}.zip','date-screen compact reference pack'),
(ROOT/'scripts/screen_goa_release_date_range.py','sequential GOA date-screen script'),
(ROOT/'scripts/package_uniprot_audit_outputs.py','UniProt small-output packager'),
(ROOT/'analysis/B104D_column_order_extended_model_summary.json','extended model summary'),
(ROOT/'analysis/B104D_column_order_model_family_summary.csv','model-family summary'),
(ROOT/'analysis/B104D_graphsage_zip_member_timestamps.csv','GraphSAGE ZIP timestamp audit'),
(ROOT/'analysis/B104D_uniprot_uploaded_package_review.json','UniProt package completeness review'),
]:
    add_ledger(record_type='generated_output',artifact_name=p.name,local_path=str(p),local_status='present_and_hash_verified',origin_in_this_run=BATCH,analysis_role=role,used_by='user delivery and future reproduction',retrieval_status='generated_in_runtime',retrieved_at_utc=ISO,size_bytes=p.stat().st_size,sha256=sha(p),parent_or_derivation='verified prior derivatives plus B104D user inputs',notes='',batch_id=BATCH,deletion_state='retain_generated_output',event_recorded_at_utc=ISO,hash_authority='container_sha256',runtime_verification_status='present_and_verified')
LED_OUT=RESULTS/f'source_ledger_through_B104D_PREBUNDLE_{STAMP}.csv'
writecsv(LED_OUT,led,lf)

# Append-only event stream.
events=readcsv(PREV_EVT); ef=list(events[0])
def ev(event_type,artifact,status,details,batch=BATCH):
    events.append({'event_time_utc':ISO,'batch_id':batch,'event_type':event_type,'artifact_name':artifact,'status':status,'details':details})
ev('user_deletion_confirmed','msigdb_v5.0_files_to_download_locally.zip','confirmed','User reported `Deleted B104C`; compact normalized derivative and frozen B104C bundle retained.')
ev('input_received',UNIPROT.name,'accepted',f'{UNIPROT.stat().st_size} bytes; SHA-256 {sha(UNIPROT)}.')
ev('input_integrity',UNIPROT.name,'passed','ZIP integrity passed; exactly one member: uniprot_2016_mapping_audit_ledger.csv.')
ev('audit_ledger_review',LEDGER.name,'partial_success','Three rows report successful size/MD5 verification, both targets found, and source archives deleted. Extracted DAT/TSV/JSON payloads absent from upload.')
ev('column_order_model_expansion','B104D extended grids','completed','492 new configurations plus 48 prior configurations; no exact order and no LCS >=100; best LCS 94 unchanged.')
ev('duplicate_orientation_stability','three duplicate-vector pairs','strongly_supported_provisional','All 161 configurations with LCS >=80 selected orientation 001; original ordered source still absent.')
ev('graphsage_timestamp_audit','graphsage_ppi.zip','completed','ppi-class_map.json extended UTC timestamp 2017-05-10T19:25:23Z; practical source-date upper bound, not cryptographic proof.')
ev('date_screen_script_validation','screen_goa_release_date_range.py','passed','Offline release-158/159 self-test reproduced B104A: v158 846 mismatches/8 exact columns; v159 0 mismatches/121 exact; best order LCS 94 and orientation 001.')
ev('uniprot_repack_script_validation','package_uniprot_audit_outputs.py','passed','Compiled; synthetic Windows-path/hash/integrity/deletion-safe packaging test passed.')
EVT_OUT=RESULTS/f'provenance_events_through_B104D_PREBUNDLE_{STAMP}.csv'
writecsv(EVT_OUT,events,ef)
print(json.dumps({'actual_manifest':str(ACT_OUT),'source_ledger':str(LED_OUT),'events':str(EVT_OUT),'rows':{'actual':len(act),'ledger':len(led),'events':len(events)}},indent=2))

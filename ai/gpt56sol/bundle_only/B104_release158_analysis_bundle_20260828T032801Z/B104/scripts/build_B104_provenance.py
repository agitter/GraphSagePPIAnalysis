#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json, os, shutil
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd

BATCH_STAMP='20260828T030759Z'
DELIVERY_STAMP='20260828T032243Z'
NOW='2026-08-28T03:22:43Z'
B102_DELETE_OBSERVED='2026-08-28T02:52:31Z'
B103_DELETE_OBSERVED='2026-08-28T02:52:31Z'
ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/B104_{BATCH_STAMP}')
RESULTS=Path('/mnt/data/ppi_repro_corrected/results')
RESULTS.mkdir(parents=True,exist_ok=True)
ANA=ROOT/'analysis'; DER=ROOT/'derived'; SCRIPTS=ROOT/'scripts'
BASE_MAN=RESULTS/'actual_input_file_manifest_through_B102_20260827T163101Z.csv'
BASE_LEDGER=RESULTS/'source_ledger_through_B102_20260827T163101Z.csv'
BASE_EVENTS=RESULTS/'provenance_events_through_B102_20260827T163101Z.csv'
BASE_INV=RESULTS/'user_local_inventory_full_enriched_20260827T163101Z.csv'


def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def hash_paths(paths):
 return '|'.join(sha(Path(p)) for p in paths)

def strpaths(paths): return '|'.join(str(Path(p)) for p in paths)

# Stable retained input copies/derivatives already exist under ROOT.
script_main=SCRIPTS/'analyze_B104_release158_fast.py'
script_alt=SCRIPTS/'analyze_B104_alternative_hypotheses.py'
parser_hash=sha(script_main)

# ---------- Actual input manifest ----------
man=pd.read_csv(BASE_MAN,dtype=str).fillna('')
# B102 deletion confirmation
mask=man['batch_id'].eq('B102')
man.loc[mask,'local_status']='conversation_copy_deletion_confirmed_by_user; residual_runtime_mount_ignored'
man.loc[mask,'user_deletion_confirmed_at_utc']=B102_DELETE_OBSERVED
man.loc[mask,'raw_retention_status']='user_local_master_retained; conversation_attachment_deleted_by_user'
man.loc[mask,'raw_available_in_conversation']='false'
man.loc[mask,'deletion_notes']=man.loc[mask,'deletion_notes'].astype(str).str.cat(pd.Series([' User explicitly reported Deleted B102; event time approximated from receipt of the immediately following B103 upload because an exact UI timestamp was unavailable.']*mask.sum(),index=man.index[mask]),sep='')

cols=list(man.columns)
new=[]
def add_man(**kw):
 row={c:'' for c in cols}; row.update({k:str(v) for k,v in kw.items()}); new.append(row)

# B103 rows, now deletion confirmed and represented by durable repaired derivatives.
obo_der=[DER/f'B104_repaired_B103_GO_terms_{BATCH_STAMP}.tsv.gz',DER/f'B104_repaired_B103_GO_is_a_edges_{BATCH_STAMP}.tsv.gz',DER/f'B104_repaired_B103_GO_is_a_closure_for_GOA158_159_terms_{BATCH_STAMP}.tsv.gz']
id_der=[DER/f'B104_repaired_B103_current_UniProt_idmapping_2026-08-27_{BATCH_STAMP}.tsv.gz']
add_man(record_type='actual_input',artifact_name='2016-06-01-go.obo',local_path='/mnt/data/2016-06-01-go.obo',local_status='conversation_copy_deletion_confirmed_by_user; residual_runtime_mount_used_once_for_explicit_provenance_repair_then_ignored',origin_in_this_run='user_upload_B103',analysis_role='historical GO ontology for term metadata and is_a propagation',used_by='analyze_B104_release158_fast.py',direct_or_canonical_source_url='https://release.geneontology.org/2016-06-01/ontology/go.obo',source_page_url='https://release.geneontology.org/2016-06-01/',url_status='official_GO_archive_path_recorded; uploaded_bytes_not_remote_compared',retrieval_status='supplied_by_user; SHA-256 matched full user-local inventory and prior B103 verification',retrieved_at_utc='',size_bytes='34762536',sha256='9b4c0c28d73ba41ae4c684d78b354d2c8bea691a5d759d4cdd188eecdd307ca2',parent_or_derivation='local_upload_inventory_full_20260827T160408Z.csv',notes='B103 deletion confirmation recorded. Residual mount was used once to reconstruct durable derivatives because the previously claimed B103 retained artifacts were absent; exact raw hash matched the already verified B103 hash.',batch_id='B103',inventory_sha256='4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4',received_at_utc='exact receipt timestamp unavailable; B103 completed 2026-08-27T21:06:59Z',integrity_verified_at_utc='2026-08-28T03:07:59Z (reverified for repair)',analysis_completed_at_utc='2026-08-27T21:06:59Z; derivative repair completed 2026-08-28T03:07:59Z',deletion_clearance_issued_at_utc='2026-08-27T21:06:59Z',user_deletion_confirmed_at_utc=B103_DELETE_OBSERVED,raw_retention_status='user_local_master_retained; conversation_attachment_deleted_by_user; residual_mount_logically_ignored after repair',raw_available_in_conversation='false',retained_derivative_paths=strpaths(obo_der),derivative_sha256s=hash_paths(obo_der),parser_version_or_script_sha256=parser_hash,raw_to_derived_reconciliation='44,797 ontology terms; 73,691 is_a edges; 238,711 closure rows for GOA158/159 direct terms; source raw SHA-256 reverified',reacquisition_url='https://release.geneontology.org/2016-06-01/ontology/go.obo',deletion_notes='User explicitly reported Deleted B103. Event time approximated from receipt of B104 uploads because exact UI timestamp was unavailable.')
add_man(record_type='actual_input',artifact_name='idmapping_2026_08_27.tsv.gz',local_path='/mnt/data/idmapping_2026_08_27.tsv.gz',local_status='conversation_copy_deletion_confirmed_by_user; residual_runtime_mount_used_once_for_explicit_provenance_repair_then_ignored',origin_in_this_run='user_upload_B103; user-generated_current_UniProt_ID_mapping_output',analysis_role='current temporal sensitivity mapping for nine GeneIDs',used_by='analyze_B104_release158_fast.py; build_B104_identifier_watchlist.py',direct_or_canonical_source_url='',source_page_url='https://www.uniprot.org/id-mapping',url_status='no_stable_job_result_URL_supplied; source_tool_page_recorded',retrieval_status='generated_by_user_with_current_UniProt_ID_mapping_tool_on_2026-08-27',retrieved_at_utc='2026-08-27',size_bytes='1219',sha256='fd585a7de7201f61871a70fbeb244b615cfa32dd7eee1b507cc35d89bd5cd5d6',parent_or_derivation='current UniProt ID Mapping query for GeneIDs 176,337,3108,4018,7957,10159,29901,55125,84919',notes='Used only for temporal sensitivity, never as a substitute for 2016 mapping. B103 deletion confirmation recorded; residual mount used once to reconstruct a durable normalized derivative.',batch_id='B103',inventory_sha256='',received_at_utc='exact receipt timestamp unavailable; B103 completed 2026-08-27T21:06:59Z',integrity_verified_at_utc='2026-08-28T03:07:59Z (reverified for repair)',analysis_completed_at_utc='2026-08-27T21:06:59Z; derivative repair completed 2026-08-28T03:07:59Z',deletion_clearance_issued_at_utc='2026-08-27T21:06:59Z',user_deletion_confirmed_at_utc=B103_DELETE_OBSERVED,raw_retention_status='user_local_master_retained; conversation_attachment_deleted_by_user; residual_mount_logically_ignored after repair',raw_available_in_conversation='false',retained_derivative_paths=strpaths(id_der),derivative_sha256s=hash_paths(id_der),parser_version_or_script_sha256=parser_hash,raw_to_derived_reconciliation='15 rows retained with original columns; derivative re-read; raw SHA-256 reverified',reacquisition_url='not_reacquirable_as_same_job_output_without_query_parameters_and_UniProt_version',deletion_notes='User explicitly reported Deleted B103. Event time approximated from receipt of B104 uploads because exact UI timestamp was unavailable.')

# B104 raw rows
summary=json.load(open(ROOT/f'B104_analysis_summary_{BATCH_STAMP}.json'))
role={'gaf':'historical GOA release-158 GAF annotations','gpad':'historical GOA release-158 GPAD annotations','gpi':'historical GOA release-158 gene-product information'}
der_for={
 'goa_human.gaf.158.gz':[DER/f'B104_goa_human_gaf158_normalized_{BATCH_STAMP}.tsv.gz'],
 'goa_human.gpa.158.gz':[DER/f'B104_goa_human_gpad158_normalized_{BATCH_STAMP}.tsv.gz'],
 'goa_human.gpi.158.gz':[DER/f'B104_goa_human_gpi158_normalized_{BATCH_STAMP}.tsv.gz'],
}
recon_for={
 'goa_human.gaf.158.gz':'388,218 rows retained; 17-column width check passed; GAF unique projection reconciled exactly to GPAD',
 'goa_human.gpa.158.gz':'389,235 rows retained; 12-column width check passed; 388,218 unique GAF-projected assertions plus 1,017 ECO subtype duplicate projections',
 'goa_human.gpi.158.gz':'21,005 rows retained; 10-column width check passed; all rows re-read; compared against 21,002 GPI159 objects',
}
for r in summary['input_integrity']:
 ds=der_for[r['artifact_name']]
 add_man(record_type='actual_input',artifact_name=r['artifact_name'],local_path=r['local_path'],local_status='present_as_B104_conversation_attachment; deletion_clearance_issued',origin_in_this_run='user_upload_B104',analysis_role=role[r['file_role']],used_by='analyze_B104_release158_fast.py',direct_or_canonical_source_url=r['direct_or_canonical_source_url'],source_page_url=r['source_page_url'],url_status='official_EBI_GOA_archive_exact_filename; uploaded_bytes_not_remote_compared',retrieval_status='supplied_by_user; size and SHA-256 matched authoritative full local inventory; gzip integrity passed',retrieved_at_utc='2026-08-28T02:52:31Z (runtime mount time; exact UI receipt time unavailable)',size_bytes=r['size_bytes'],sha256=r['sha256'],parent_or_derivation='local_upload_inventory_full_20260827T160408Z.csv',notes='B104 fully analyzed. Clearance applies only to conversation attachment; retain user-local master.',batch_id='B104',inventory_sha256='4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4',received_at_utc='2026-08-28T02:52:31Z (runtime mount time)',integrity_verified_at_utc='2026-08-28T03:07:59Z',analysis_completed_at_utc=NOW,deletion_clearance_issued_at_utc=NOW,user_deletion_confirmed_at_utc='',raw_retention_status='user_local_master_retained; conversation_attachment_pending_user_deletion',raw_available_in_conversation='true_until_user_confirms_Deleted_B104',retained_derivative_paths=strpaths(ds),derivative_sha256s=hash_paths(ds),parser_version_or_script_sha256=parser_hash,raw_to_derived_reconciliation=recon_for[r['artifact_name']],reacquisition_url=r['direct_or_canonical_source_url'],deletion_notes='SAFE TO DELETE B104 issued; awaiting user confirmation.')

man2=pd.concat([man,pd.DataFrame(new,columns=cols)],ignore_index=True)
man_out=RESULTS/f'actual_input_file_manifest_through_B104_{DELIVERY_STAMP}.csv'
man2.to_csv(man_out,index=False)
# readable markdown, one row per actual input with key fields
md_out=RESULTS/f'actual_input_file_manifest_through_B104_{DELIVERY_STAMP}.md'
with md_out.open('w',encoding='utf-8') as f:
 f.write('# Actual input manifest through B104\n\n')
 f.write(f'Generated: `{NOW}`. Rows: {len(man2)}. This manifest records actual bytes used or explicitly retained derivatives; historical candidates not materialized remain in the source ledger instead.\n\n')
 f.write('| Batch | Artifact | Bytes | SHA-256 | Current raw state | Source URL | Retained derivative(s) |\n|---|---|---:|---|---|---|---|\n')
 for _,r in man2.iterrows():
  f.write(f"| {r['batch_id']} | `{r['artifact_name']}` | {r['size_bytes']} | `{r['sha256']}` | {r['raw_retention_status']} | {r['direct_or_canonical_source_url'] or r['source_page_url']} | `{r['retained_derivative_paths']}` |\n")

# ---------- Local inventory augmented ----------
inv=pd.read_csv(BASE_INV,dtype=str).fillna('')
# update runtime states for B103/B104
for name,h,status,batch in [
 ('2016-06-01-go.obo','9b4c0c28d73ba41ae4c684d78b354d2c8bea691a5d759d4cdd188eecdd307ca2','B103_uploaded_hash_reverified; deletion_confirmed','B103'),
 ('goa_human.gaf.158.gz','7d5f7aabd0bea1e1f2a9d18af70f5d4038a85a78736d07ba69fc331b34241acf','B104_uploaded_hash_matches_inventory; gzip_passed','B104'),
 ('goa_human.gpa.158.gz','4d1b31df7490ad55c215d2e8525a098d820bf12a7d2d26cd13bc58a633d5f26a','B104_uploaded_hash_matches_inventory; gzip_passed','B104'),
 ('goa_human.gpi.158.gz','2c7a7a836d022038431a5efbfa48dbe0dd1777264e008f693b78387568dd354a','B104_uploaded_hash_matches_inventory; gzip_passed','B104')]:
 m=inv['artifact_name'].eq(name)
 inv.loc[m,'runtime_verification_status']=status
 inv.loc[m,'runtime_recomputed_sha256']=h
 inv.loc[m,'planned_batch']=batch
 inv.loc[m,'last_updated_at_utc']=NOW
# append current idmapping generated after initial inventory
if not inv['artifact_name'].eq('idmapping_2026_08_27.tsv.gz').any():
 row={c:'' for c in inv.columns}
 row.update({'artifact_name':'idmapping_2026_08_27.tsv.gz','relative_path':'idmapping_2026_08_27.tsv.gz','size_bytes':'1219','sha256':'fd585a7de7201f61871a70fbeb244b615cfa32dd7eee1b507cc35d89bd5cd5d6','mtime_utc':'2026-08-28T02:52:31.629755+00:00','file_suffixes':'.tsv.gz','direct_or_canonical_source_url':'','source_page_url':'https://www.uniprot.org/id-mapping','obtained_from':'user-generated current UniProt ID Mapping output','acquired_at_utc':'2026-08-27','provenance_notes':'Query for nine GeneIDs; exact job URL and UniProt release identifier not supplied. Used only for temporal sensitivity.','url_status':'source_tool_page_recorded; no_stable_result_URL','source_notes':'post-inventory file added from verified B103 upload','inventory_semantics':'post_inventory_user_local_file_declaration_verified_on_upload','inventory_file':'local_upload_inventory_full_20260827T160408Z.csv plus B103 upload','inventory_file_sha256':'4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4','runtime_verification_status':'B103_uploaded_hash_verified; deletion_confirmed','runtime_recomputed_sha256':'fd585a7de7201f61871a70fbeb244b615cfa32dd7eee1b507cc35d89bd5cd5d6','planned_batch':'B103','last_updated_at_utc':NOW})
 inv=pd.concat([inv,pd.DataFrame([row])],ignore_index=True)
inv_out=RESULTS/f'user_local_inventory_augmented_through_B104_{DELIVERY_STAMP}.csv'
inv.to_csv(inv_out,index=False)

# ---------- Source ledger append-only ----------
ledger=pd.read_csv(BASE_LEDGER,dtype=str).fillna('')
lcols=list(ledger.columns); lnew=[]
def add_led(**kw):
 r={c:'' for c in lcols}; r.update({k:str(v) for k,v in kw.items()}); lnew.append(r)
# Deletion confirmations
for name in ['2016-06-01-annotations-README','2016-06-01-gp2protein.geneid.gz','2016-06-01-gp2protein.human.gz','local_upload_inventory_full_20260827T160408Z.csv']:
 add_led(record_type='deletion_confirmation',artifact_name=name,local_status='conversation_attachment_deleted_by_user',origin_in_this_run='user_confirmation',analysis_role='raw lifecycle tracking',used_by='provenance protocol',notes='User reported Deleted B102; exact event time unavailable and approximated from next upload receipt.',batch_id='B102',deletion_state='deletion_confirmed_by_user',event_recorded_at_utc=B102_DELETE_OBSERVED,hash_authority='prior B102 verified hash',runtime_verification_status='raw_mount_ignored_after_confirmation',user_deletion_confirmed_at_utc=B102_DELETE_OBSERVED)
# B103 actual inputs and repair
for _,r in pd.DataFrame(new[:2]).iterrows():
 add_led(record_type='actual_input',artifact_name=r['artifact_name'],local_path=r['local_path'],local_status=r['local_status'],origin_in_this_run=r['origin_in_this_run'],analysis_role=r['analysis_role'],used_by=r['used_by'],direct_or_canonical_source_url=r['direct_or_canonical_source_url'],source_page_url=r['source_page_url'],url_status=r['url_status'],retrieval_status=r['retrieval_status'],size_bytes=r['size_bytes'],sha256=r['sha256'],parent_or_derivation=r['parent_or_derivation'],notes=r['notes'],batch_id='B103',user_local_relative_path=r['artifact_name'],user_local_status='retained_by_user_locally',user_local_size_bytes=r['size_bytes'],user_local_sha256=r['sha256'],inventory_file='local_upload_inventory_full_20260827T160408Z.csv' if r['artifact_name'].endswith('.obo') else 'post_inventory_B103_upload',inventory_sha256=r['inventory_sha256'],deletion_state='deletion_confirmed_by_user',event_recorded_at_utc=NOW,hash_authority='runtime_recomputed_and_prior_B103_verified',runtime_verification_status='hash_reverified_during_B104_provenance_repair',planned_batch='B103',retained_derivative_paths=r['retained_derivative_paths'],parser_script_sha256=parser_hash,raw_to_derived_reconciliation=r['raw_to_derived_reconciliation'],user_deletion_confirmed_at_utc=B103_DELETE_OBSERVED)
 add_led(record_type='provenance_repair',artifact_name=r['artifact_name'],local_status='residual_mount_used_once_then_logically_ignored',origin_in_this_run='B104_explicit_repair',analysis_role='reconstruct missing durable B103 derivative',used_by='analyze_B104_release158_fast.py',notes='This does not reverse user deletion confirmation. Source bytes matched previously verified SHA-256. Durable normalized derivative created and all later analyses use derivative.',batch_id='B104',deletion_state='B103_deletion_still_confirmed',event_recorded_at_utc='2026-08-28T03:07:59Z',hash_authority='runtime_recomputed_equals_prior_verified',runtime_verification_status='complete',retained_derivative_paths=r['retained_derivative_paths'],parser_script_sha256=parser_hash,raw_to_derived_reconciliation=r['raw_to_derived_reconciliation'],user_deletion_confirmed_at_utc=B103_DELETE_OBSERVED)
# B104 raw
for r in new[2:]:
 add_led(record_type='actual_input',artifact_name=r['artifact_name'],local_path=r['local_path'],local_status=r['local_status'],origin_in_this_run=r['origin_in_this_run'],analysis_role=r['analysis_role'],used_by=r['used_by'],direct_or_canonical_source_url=r['direct_or_canonical_source_url'],source_page_url=r['source_page_url'],url_status=r['url_status'],retrieval_status=r['retrieval_status'],size_bytes=r['size_bytes'],sha256=r['sha256'],parent_or_derivation=r['parent_or_derivation'],notes=r['notes'],batch_id='B104',user_local_relative_path=r['artifact_name'],user_local_status='present_on_user_machine_per_full_inventory',user_local_size_bytes=r['size_bytes'],user_local_sha256=r['sha256'],inventory_file='local_upload_inventory_full_20260827T160408Z.csv',inventory_sha256=r['inventory_sha256'],deletion_state='deletion_clearance_issued_pending_user_confirmation',event_recorded_at_utc=NOW,hash_authority='runtime_recomputed_matches_user_full_inventory',runtime_verification_status='gzip_integrity_and_complete_parse_passed',planned_batch='B104',retained_derivative_paths=r['retained_derivative_paths'],parser_script_sha256=parser_hash,raw_to_derived_reconciliation=r['raw_to_derived_reconciliation'])
# Important outputs
output_records=[
 (ROOT/f'B104_REPORT_{BATCH_STAMP}.md','accepted scientific report'),
 (ROOT/f'B104_EXECUTION_DIAGNOSTICS_{BATCH_STAMP}.md','complete execution diagnostics'),
 (ROOT/f'B104_analysis_summary_{BATCH_STAMP}.json','machine-readable accepted summary'),
 (ANA/f'B104_label_to_GO_mapping_release158_159_{BATCH_STAMP}.csv','121-column GO label map'),
 (ANA/f'B104_identifier_mapping_watchlist_{BATCH_STAMP}.csv','temporal identifier mapping watchlist'),
 (ANA/f'B104_alternative_hypothesis_checks_{BATCH_STAMP}.json','mapping, NOT, and all-zero sensitivity checks'),
]
for p,role_ in output_records:
 add_led(record_type='generated_output',artifact_name=p.name,local_path=str(p),local_status='retained_and_hashed',origin_in_this_run='B104_analysis',analysis_role=role_,used_by='B104 report and future batches',size_bytes=p.stat().st_size,sha256=sha(p),parent_or_derivation='B104 verified inputs plus retained B101/B102/B103 derivatives',notes='Accepted output; raw B104 files excluded from final bundle.',batch_id='B104',deletion_state='retain_generated_output',event_recorded_at_utc=NOW,hash_authority='runtime_recomputed',runtime_verification_status='present_and_hashed',retained_derivative_paths=str(p),parser_script_sha256=parser_hash)
# Mapping decisions and web references
for gene,artifact,note,url in [
 ('7957','GeneID 7957 / EPM2A temporal mapping decision','B3EWF7->7957 retained as historically contextual symbol fallback; exact 17-label row; no replacement claim between B3EWF7 and O95278.','https://www.uniprot.org/uniprotkb/B3EWF7'),
 ('29901','GeneID 29901 / SAC3D1 temporal mapping decision','A6NKF1->29901 tracked as contextual fallback; observed row all zero; F8WC89 and A6NKF1 treated as coexisting products, not inferred replacements.','https://www.uniprot.org/uniprotkb/A6NKF1'),
 ('10159','GeneID 10159 / ATP6AP2 temporal mapping decision','O75787 is current/historical anchor; secondary accessions not independent primary records; PSEC0072 synonym join prohibited due SIDT2 collision.','https://www.uniprot.org/uniprotkb/O75787'),
]:
 add_led(record_type='mapping_decision',artifact_name=artifact,local_status='decision_record',origin_in_this_run='B104_analysis_plus_official_current_sources',analysis_role='temporal identifier audit',used_by='future mapping reconstruction',direct_or_canonical_source_url=url,source_page_url=url,url_status='official_current_UniProt_page',retrieval_status='web_inspected_2026-08-28',notes=note,batch_id='B104',event_recorded_at_utc=NOW,hash_authority='not_applicable_decision_record',runtime_verification_status='local historical data cross-checked with official current page',planned_batch='watchlist',retained_derivative_paths=str(ANA/f'B104_identifier_mapping_watchlist_{BATCH_STAMP}.csv'),parser_script_sha256=sha(SCRIPTS/'build_B104_identifier_watchlist.py'))
# Exact ontology candidate
add_led(record_type='historical_candidate',artifact_name='2016-06-29 GO ontology matching GOA159 header',local_status='not_materialized_in_runtime',origin_in_this_run='identified_from_GAF159_header',analysis_role='highest-priority ontology-version test for 878 ancestor-only residuals',used_by='planned B105 exact ontology test',direct_or_canonical_source_url='http://purl.obolibrary.org/obo/go/releases/2016-06-29/go.obo',source_page_url='https://release.geneontology.org/2016-07-01/',url_status='exact PURL inferred from GAF header OWL version; monthly archive fallback recorded; bytes not obtained',retrieval_status='not_materialized',notes='Prefer exact June-29 serialization; otherwise upload July-01 archive go.obo and verify internal data-version before analysis.',batch_id='B105_planned',deletion_state='not_applicable',event_recorded_at_utc=NOW,hash_authority='not_available',runtime_verification_status='not_downloaded',planned_batch='B105')
ledger2=pd.concat([ledger,pd.DataFrame(lnew,columns=lcols)],ignore_index=True)
ledger_out=RESULTS/f'source_ledger_through_B104_{DELIVERY_STAMP}.csv'
ledger2.to_csv(ledger_out,index=False)
ledger_md=RESULTS/f'source_ledger_through_B104_{DELIVERY_STAMP}.md'
with ledger_md.open('w',encoding='utf-8') as f:
 f.write('# Source ledger through B104\n\n')
 f.write(f'Generated: `{NOW}`. Append-only rows: {len(ledger2)}. Earlier rows are preserved even when superseded by later deletion confirmations or corrections.\n\n')
 f.write('| Record type | Batch | Artifact | Status | SHA-256 | Source | Notes |\n|---|---|---|---|---|---|---|\n')
 for _,r in ledger2.tail(80).iterrows():
  f.write(f"| {r['record_type']} | {r['batch_id']} | `{r['artifact_name']}` | {r['local_status'] or r['deletion_state']} | `{r['sha256']}` | {r['direct_or_canonical_source_url'] or r['source_page_url']} | {r['notes']} |\n")

# ---------- Events ----------
events=pd.read_csv(BASE_EVENTS,dtype=str).fillna('')
evnew=[]
def ev(t,b,e,a,s,d): evnew.append({'event_time_utc':t,'batch_id':b,'event_type':e,'artifact_name':a,'status':s,'details':d})
for name in ['2016-06-01-annotations-README','2016-06-01-gp2protein.geneid.gz','2016-06-01-gp2protein.human.gz','local_upload_inventory_full_20260827T160408Z.csv']:
 ev(B102_DELETE_OBSERVED,'B102','user_deletion_confirmed',name,'conversation_copy_deleted','User reported Deleted B102; time approximated from immediate next upload receipt because exact UI timestamp unavailable.')
for name,h in [('2016-06-01-go.obo','9b4c0c28d73ba41ae4c684d78b354d2c8bea691a5d759d4cdd188eecdd307ca2'),('idmapping_2026_08_27.tsv.gz','fd585a7de7201f61871a70fbeb244b615cfa32dd7eee1b507cc35d89bd5cd5d6')]:
 ev(B103_DELETE_OBSERVED,'B103','user_deletion_confirmed',name,'conversation_copy_deleted','User reported Deleted B103; time approximated from B104 upload receipt because exact UI timestamp unavailable.')
 ev('2026-08-28T03:07:59Z','B104','prior_batch_derivative_repair',name,'completed','Residual mounted bytes matched the prior verified hash and were used once to reconstruct a durable derivative; logical deletion remained in force.')
for r in summary['input_integrity']:
 ev('2026-08-28T02:52:31Z','B104','input_received',r['artifact_name'],'received','Runtime mount time recorded; exact UI receipt time unavailable.')
 ev('2026-08-28T03:07:59Z','B104','integrity_verified',r['artifact_name'],'passed',f"Size and SHA-256 matched full local inventory; gzip -t passed; SHA-256 {r['sha256']}.")
ev('2026-08-28T03:07:59Z','B104','analysis_attempt', 'analyze_B104_release158.py','timed_out','Exhaustive initial implementation exceeded 120 seconds; no partial outputs accepted.')
ev('2026-08-28T03:07:59Z','B104','analysis_completed','analyze_B104_release158_fast.py','accepted','Exit 0; release158 1,733 differences, release159 901; normalized derivatives and reconciliation outputs written.')
ev(NOW,'B104','alternative_hypotheses_completed','analyze_B104_alternative_hypotheses.py','accepted','Exit 0; mapping-policy, NOT, all-zero-gene, and residual-structure tests written.')
ev(NOW,'B104','identifier_watchlist_created','B104_identifier_mapping_watchlist','accepted','GeneIDs 7957, 29901, and 10159 tracked without replacement assumptions.')
ev(NOW,'B104','deletion_clearance_issued','B104 raw conversation attachments','safe_to_delete','All three B104 uploads passed integrity, parsing, analysis, derivative retention, and final validation; awaiting user confirmation Deleted B104.')
events2=pd.concat([events,pd.DataFrame(evnew)],ignore_index=True)
ev_out=RESULTS/f'provenance_events_through_B104_{DELIVERY_STAMP}.csv'
events2.to_csv(ev_out,index=False)

print(man_out); print(md_out); print(inv_out); print(ledger_out); print(ledger_md); print(ev_out)

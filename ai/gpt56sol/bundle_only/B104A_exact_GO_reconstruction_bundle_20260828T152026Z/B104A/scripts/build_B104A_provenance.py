#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import pandas as pd

BATCH='B104A'
BATCH_STAMP='20260828T145842Z'
FINAL_STAMP='20260828T152026Z'
EVENT_TIME='2026-08-28T15:20:26Z'
ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/{BATCH}_{BATCH_STAMP}')
RESULTS=Path('/mnt/data/ppi_repro_corrected/results')
BASE_MAN=RESULTS/'actual_input_file_manifest_through_B104_20260828T032243Z.csv'
BASE_LEDGER=RESULTS/'source_ledger_through_B104_FINAL_20260828T032801Z.csv'
BASE_EVENTS=RESULTS/'provenance_events_through_B104_FINAL_20260828T032801Z.csv'

OUT_MAN=RESULTS/f'actual_input_file_manifest_through_B104A_{FINAL_STAMP}.csv'
OUT_LEDGER=RESULTS/f'source_ledger_through_B104A_{FINAL_STAMP}.csv'
OUT_EVENTS=RESULTS/f'provenance_events_through_B104A_{FINAL_STAMP}.csv'


def sha(path: Path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while True:
            b=f.read(1<<20)
            if not b: break
            h.update(b)
    return h.hexdigest()

script_hashes={p.name:sha(p) for p in sorted((ROOT/'scripts').glob('*.py'))}
parser_hash='|'.join(f'{k}:{v}' for k,v in script_hashes.items())

# -------- actual input manifest --------
man=pd.read_csv(BASE_MAN,dtype=str).fillna('')
# Current-state update for B104 raw conversation copies.
for name in ['goa_human.gaf.158.gz','goa_human.gpa.158.gz','goa_human.gpi.158.gz']:
    mask=man.artifact_name.eq(name)
    man.loc[mask,'local_status']='conversation_copy_deleted_by_user; user_local_master_retained'
    man.loc[mask,'user_deletion_confirmed_at_utc']=EVENT_TIME
    man.loc[mask,'raw_retention_status']='user_local_master_retained; conversation_attachment_deleted_by_user'
    man.loc[mask,'raw_available_in_conversation']='false'
    man.loc[mask,'deletion_notes']='User reported Deleted B104. Logical deletion remains in force; B104A used only frozen bundles and retained derivatives.'

# Append B104A script usage to existing source rows.
usage_updates={
    'graphsage_ppi.zip':'run_core_verification.py; finish_B104A_graph_summary.py',
    'msigdb_v5.2_files_to_download_locally.zip':'run_core_verification.py; run_local_label_source_screen.py; finalize_B104A_exact_reconstruction.py',
    'goa_human.gaf.159.gz':'analyze_B101_goa159.py; analyze_B104A_residuals.py; continue_B104A.py; finalize_B104A_exact_reconstruction.py; validate_B104A_exact_reconstruction_set_based.py',
    '2016-06-01-gp2protein.geneid.gz':'analyze_B102_gp2protein_final.py; finalize_B104A_exact_reconstruction.py',
    '2016-06-01-go.obo':'analyze_B104_release158_fast.py; analyze_B104A_residuals.py; continue_B104A.py; finalize_B104A_exact_reconstruction.py; validate_B104A_exact_reconstruction_set_based.py',
    'collapsed_gene_labels_topology_features.csv':'analyze_B102_gp2protein_final.py; all B104A reconstruction/validation scripts',
}
for name,used in usage_updates.items():
    man.loc[man.artifact_name.eq(name),'used_by']=used

cols=list(man.columns)
def row_template(): return {c:'' for c in cols}
retained_specs=[
 ('B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz','release-159 normalized GAF','B101 raw GAF159','analyze_B104A_residuals.py; continue_B104A.py; finalize_B104A_exact_reconstruction.py; validate_B104A_exact_reconstruction_set_based.py','https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.159.gz','https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/'),
 ('B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz','release-159 normalized GPI and primary symbols','B101 raw GPI159','finalize_B104A_exact_reconstruction.py','https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gpi.159.gz','https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/'),
 ('B102_gp2protein_geneid_relevant_subset_20260827T162132Z.tsv.gz','full relevant May-2016 GeneID-UniProt component, including GeneIDs outside GraphSAGE','B102 raw gp2protein.geneid','finalize_B104A_exact_reconstruction.py','https://release.geneontology.org/2016-06-01/annotations/gp2protein/gp2protein.geneid.gz','https://release.geneontology.org/2016-06-01/annotations/gp2protein/'),
 ('B104_accession_GeneID_mapping_edges_20260828T030759Z.csv.gz','prior accepted accession-GeneID projection for factorial comparison','B104 mapping analysis','all B104A reconstruction/validation scripts','',''),
 ('B104_label_to_GO_mapping_release158_159_20260828T030759Z.csv','previously recovered GO IDs and column order','B104 label analysis','all B104A reconstruction/validation scripts','',''),
 ('B104_repaired_B103_GO_terms_20260828T030759Z.tsv.gz','normalized May-31 GO term metadata','B103 ontology through B104 repair','all B104A ontology/reconstruction scripts','https://release.geneontology.org/2016-06-01/ontology/go.obo','https://release.geneontology.org/2016-06-01/ontology/'),
 ('B104_repaired_B103_GO_is_a_edges_20260828T030759Z.tsv.gz','normalized May-31 GO is_a edges','B103 ontology through B104 repair','all B104A ontology/reconstruction scripts','https://release.geneontology.org/2016-06-01/ontology/go.obo','https://release.geneontology.org/2016-06-01/ontology/'),
 ('B104_goa_human_gaf158_normalized_20260828T030759Z.tsv.gz','release-158 normalized GAF for temporal control','B104 raw GAF158','finalize_B104A_exact_reconstruction.py','https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/goa_human.gaf.158.gz','https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/'),
 ('collapsed_gene_labels_topology_features.csv','collapsed deposited labels by independently recovered GeneID','corrected core analysis','all B104A reconstruction/validation scripts','',''),
 ('graphsage_row_to_entrez_topology_features.csv','row-level topology-derived GraphSAGE-to-GeneID map','corrected core analysis','finish_B104A_graph_summary.py','',''),
 ('tissue_partition.csv','GraphSAGE graph/tissue row ranges','corrected core analysis','finish_B104A_graph_summary.py','',''),
]
add=[]
for filename,role,parent,used,url,page in retained_specs:
    p=ROOT/'retained_inputs'/filename
    r=row_template()
    r.update({
      'record_type':'retained_input','artifact_name':filename,'local_path':str(p),'local_status':'present_verified',
      'origin_in_this_run':'copied_from_frozen_prior_bundle_or_stable_prior_derivative','analysis_role':role,'used_by':used,
      'direct_or_canonical_source_url':url,'source_page_url':page,'url_status':'prior_source_url_recorded' if url else 'not_applicable_generated_derivative',
      'retrieval_status':'retained_derivative_from_verified_prior_batch','retrieved_at_utc':'','size_bytes':str(p.stat().st_size),'sha256':sha(p),
      'parent_or_derivation':parent,'notes':'B104A consumed this retained copy; no logically deleted raw attachment was used.','batch_id':BATCH,
      'received_at_utc':EVENT_TIME,'integrity_verified_at_utc':EVENT_TIME,'analysis_completed_at_utc':EVENT_TIME,
      'raw_retention_status':'retained_derivative_required_for_reproduction','raw_available_in_conversation':'not_applicable_retained_derivative',
      'retained_derivative_paths':str(p),'derivative_sha256s':sha(p),'parser_version_or_script_sha256':parser_hash,
      'raw_to_derived_reconciliation':'Input hash recorded and file successfully re-read by accepted B104A scripts.','reacquisition_url':url,
    })
    add.append(r)
man=pd.concat([man,pd.DataFrame(add,columns=cols)],ignore_index=True)
man.to_csv(OUT_MAN,index=False)

# -------- append-only provenance events --------
events=pd.read_csv(BASE_EVENTS,dtype=str).fillna('')
new_events=[]
def ev(event_type,artifact,status,details,time=EVENT_TIME,batch=BATCH):
    new_events.append({'event_time_utc':time,'batch_id':batch,'event_type':event_type,'artifact_name':artifact,'status':status,'details':details})
for name in ['goa_human.gaf.158.gz','goa_human.gpa.158.gz','goa_human.gpi.158.gz']:
    ev('user_deletion_confirmed',name,'conversation_copy_deleted','User reported Deleted B104; user-local master retained. B104A used frozen bundle/derivatives only.',batch='B104')
ev('analysis_attempt','analyze_B104A_residuals.py','timed_out','Timed out after writing completed relation/depth/evidence outputs; process itself not accepted as a complete run.')
ev('analysis_attempt','continue_B104A.py','partial_then_failed','Completed residual/date/source/edge outputs, then failed on row-array length mismatch; no outputs after the failing statement accepted.')
ev('analysis_completed','finish_B104A_graph_summary.py','accepted','Corrected graph row indexing via explicit graphsage_row; gene-property summaries completed.')
ev('qualifier_hypothesis_tested','GAF contributes_to and colocalizes_with','accepted','Uniform exclusion from ordinary binary membership reduced release159 differences from 901 to 13 without false negatives.')
ev('mapping_component_resolved','O95073/Q9Y620/GeneID25788/100861412','accepted','Full historical component and independent primary symbols identify O95073 as FSBP and Q9Y620 as RAD54B; O95073->25788 excluded before GraphSAGE projection. No label values used for the mapping decision.')
ev('exact_reconstruction_completed','finalize_B104A_exact_reconstruction.py','accepted','GOA159 final model reproduced all 516,428 binary cells exactly; release158 control retained 846 differences.')
ev('independent_validation_completed','validate_B104A_exact_reconstruction_set_based.py','accepted','Separate set-based implementation reproduced 121 exact columns, zero FP/FN, and identical matrix hash.')
ev('report_generated',f'B104A_REPORT_{BATCH_STAMP}.md','accepted','Detailed residual definition, ontology depth analysis, qualifier decomposition, mapping audit, exact reconstruction, and B105 implications documented.')
events=pd.concat([events,pd.DataFrame(new_events,columns=events.columns)],ignore_index=True)
events.to_csv(OUT_EVENTS,index=False)

# -------- append-only source ledger --------
ledger=pd.read_csv(BASE_LEDGER,dtype=str).fillna('')
ledger_cols=list(ledger.columns)
def led_template(): return {c:'' for c in ledger_cols}
new_ledger=[]
# deletion confirmations
for name in ['goa_human.gaf.158.gz','goa_human.gpa.158.gz','goa_human.gpi.158.gz']:
    r=led_template(); r.update({'record_type':'deletion_confirmation','artifact_name':name,'local_status':'conversation_copy_deleted_by_user','origin_in_this_run':'user_confirmation','analysis_role':'provenance state transition','used_by':'B104A provenance only','notes':'Logical deletion confirmed; no B104A raw use.','batch_id':'B104','deletion_state':'deleted_by_user','event_recorded_at_utc':EVENT_TIME,'supersedes_prior_record_type':'B104 pending deletion state','hash_authority':'prior verified B104 manifest','runtime_verification_status':'not_rechecked_after_logical_deletion','user_deletion_confirmed_at_utc':EVENT_TIME})
    new_ledger.append(r)
# retained inputs
for r0 in add:
    r=led_template()
    for c in ledger_cols:
        if c in r0: r[c]=r0[c]
    r['record_type']='retained_input'
    r['hash_authority']='container_sha256_of_retained_copy'
    r['runtime_verification_status']='present_and_re_read_by_accepted_B104A_scripts'
    r['event_recorded_at_utc']=EVENT_TIME
    r['retained_derivative_paths']=r0['local_path']
    r['parser_script_sha256']=parser_hash
    new_ledger.append(r)
# web references
web_refs=[
 ('NCBI Gene 25788 RAD54B','https://www.ncbi.nlm.nih.gov/gene/25788','current independent identity check for GeneID 25788/RAD54B'),
 ('NCBI Gene 100861412 FSBP','https://www.ncbi.nlm.nih.gov/gene/100861412','current independent identity check for GeneID 100861412/FSBP and O95073'),
 ('UniProt O95073 FSBP','https://www.uniprot.org/uniprotkb/O95073','current independent identity check for O95073/FSBP'),
 ('UniProt Q9Y620 RAD54B','https://www.uniprot.org/uniprotkb/Q9Y620/entry','current independent identity check for Q9Y620/RAD54B'),
 ('GO GAF 2.1 format','https://geneontology.org/docs/go-annotation-file-gaf-format-2.1/','official qualifier semantics for NOT, contributes_to, and colocalizes_with'),
]
for name,url,role in web_refs:
    r=led_template(); r.update({'record_type':'web_reference','artifact_name':name,'local_status':'not_materialized','origin_in_this_run':'web_verification','analysis_role':role,'used_by':f'B104A_REPORT_{BATCH_STAMP}.md','direct_or_canonical_source_url':url,'source_page_url':url,'url_status':'official_primary_source_opened','retrieval_status':'inspected_online_not_downloaded','notes':'Current records are sensitivity/identity evidence; they do not by themselves prove 2016 cross-reference state.','batch_id':BATCH,'event_recorded_at_utc':EVENT_TIME,'hash_authority':'not_applicable_web_reference','runtime_verification_status':'official page inspected'})
    new_ledger.append(r)
# generated outputs existing at this point
for p in sorted(ROOT.rglob('*')):
    if not p.is_file() or 'retained_inputs' in p.parts: continue
    r=led_template(); r.update({'record_type':'generated_output','artifact_name':p.name,'local_path':str(p),'local_status':'present','origin_in_this_run':'B104A_internal_analysis','analysis_role':'B104A output, script, report, log, or validation artifact','used_by':'B104A delivery and future B105 comparisons','size_bytes':str(p.stat().st_size),'sha256':sha(p),'parent_or_derivation':'B104A retained inputs and frozen prior bundles','notes':'Generated without new raw user uploads.','batch_id':BATCH,'event_recorded_at_utc':EVENT_TIME,'hash_authority':'container_sha256','runtime_verification_status':'present','parser_script_sha256':parser_hash})
    new_ledger.append(r)
ledger=pd.concat([ledger,pd.DataFrame(new_ledger,columns=ledger_cols)],ignore_index=True)
ledger.to_csv(OUT_LEDGER,index=False)

print(json.dumps({'manifest_rows':len(man),'ledger_rows':len(ledger),'event_rows':len(events),'outputs':[str(OUT_MAN),str(OUT_LEDGER),str(OUT_EVENTS)]},indent=2))

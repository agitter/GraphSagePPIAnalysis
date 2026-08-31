#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

BATCH = 'B104C'
START = '20260828T194921Z'
DELIVERY = '20260828T200948Z'
EVENT_TIME = '2026-08-28T20:09:48Z'
RECEIVED = '2026-08-28T17:52:44Z'
ROOT = Path(f'/mnt/data/ppi_repro_corrected/batches/{BATCH}_{START}')
RESULTS = Path('/mnt/data/ppi_repro_corrected/results')
BASE_MANIFEST = RESULTS/'actual_input_file_manifest_through_B104B_20260828T161826Z.csv'
BASE_LEDGER = RESULTS/'source_ledger_through_B104B_FINAL2_20260828T161826Z.csv'
BASE_EVENTS = RESULTS/'provenance_events_through_B104B_FINAL2_20260828T161826Z.csv'
RAW = Path('/mnt/data/msigdb_v5.0_files_to_download_locally.zip')
DER = ROOT/'retained_inputs/B104C_msigdb_v5.0_normalized_entrez_gene_sets_20260828T194921Z.tsv.gz'
RECON = ROOT/'analysis/B104C_msigdb_v5.0_normalization_reconciliation_20260828T194921Z.json'
MAIN = ROOT/'scripts/analyze_B104C_msigdb50_column_order.py'
EXT = ROOT/'scripts/explore_B104C_column_order_variants.py'
NORM = ROOT/'scripts/normalize_msigdb50.py'
UNI = ROOT/'scripts/download_extract_uniprot_2016_mapping_audit.py'
REPORT = ROOT/f'B104C_MSIGDB50_COLUMN_ORDER_REPORT_{DELIVERY}.md'
DIAG = ROOT/f'B104C_EXECUTION_DIAGNOSTICS_{DELIVERY}.md'
UNI_INST = ROOT/f'B104C_UNIPROT_AUDIT_INSTRUCTIONS_{DELIVERY}.md'
UNIQUE_ORDER = ROOT/'analysis/B104C_inferred_unique_121_GO_column_order_20260828T194921Z.csv'


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def md_table(df: pd.DataFrame) -> str:
    vals = df.fillna('').astype(str)
    headers = vals.columns.tolist()
    out = ['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |']
    for row in vals.itertuples(index=False, name=None):
        out.append('| ' + ' | '.join(str(x).replace('|', '\\|').replace('\n', '<br>') for x in row) + ' |')
    return '\n'.join(out) + '\n'

# Input receipt.
receipt = ROOT/f'B104C_INPUT_RECEIPT_{DELIVERY}.md'
receipt.write_text(f'''# B104C input receipt\n\nReceived and verified:\n\n- `msigdb_v5.0_files_to_download_locally.zip`\n- Bytes: {RAW.stat().st_size:,}\n- SHA-256: `{sha(RAW)}`\n- User inventory match: yes\n- ZIP/XML parse: passed\n- Internal XML build date: Apr 27, 2015\n- MSigDB v5.0 public release-notes label: Mar 2015\n\nRetained normalized derivative:\n\n- `{DER}`\n- Bytes: {DER.stat().st_size:,}\n- SHA-256: `{sha(DER)}`\n- Reconciliation: 10,348 rows; 1,454 C5 rows; 1,315,074 summed unique memberships.\n''', encoding='utf-8')

# Actual-input manifest.
manifest = pd.read_csv(BASE_MANIFEST, dtype=str).fillna('')
cols = manifest.columns.tolist()
script_hashes = '|'.join(f'{p.name}:{sha(p)}' for p in [MAIN, EXT, NORM])
raw_row = {c: '' for c in cols}
raw_row.update({
    'record_type': 'actual_input',
    'artifact_name': RAW.name,
    'local_path': str(RAW),
    'local_status': 'present_and_hash_verified_at_analysis_time',
    'origin_in_this_run': 'user_upload; predeclared in full local inventory',
    'analysis_role': 'MSigDB v5.0 direct label-membership test, cross-version C5 identity comparison, and candidate source-order analysis',
    'used_by': 'analyze_B104C_msigdb50_column_order.py; normalize_msigdb50.py',
    'direct_or_canonical_source_url': 'https://www.gsea-msigdb.org/gsea/downloads.jsp',
    'source_page_url': 'https://docs.gsea-msigdb.org/MSigDB/Release_Notes/MSigDB_5.0/',
    'url_status': 'official_authenticated_download_page; direct archive URL unavailable without login',
    'retrieval_status': 'supplied_by_user; not downloaded by analysis runtime',
    'size_bytes': str(RAW.stat().st_size),
    'sha256': sha(RAW),
    'parent_or_derivation': 'MSigDB v5.0 authenticated local-download bundle',
    'notes': 'Archive matches the B000B full-inventory size and SHA-256. XML internal BUILD_DATE is Apr 27, 2015; public release notes label v5.0 as Mar 2015.',
    'batch_id': BATCH,
    'inventory_sha256': '4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4',
    'received_at_utc': RECEIVED,
    'integrity_verified_at_utc': EVENT_TIME,
    'analysis_completed_at_utc': EVENT_TIME,
    'deletion_clearance_issued_at_utc': EVENT_TIME,
    'raw_retention_status': 'normalized_complete_derivative_retained; deletion_clearance_issued',
    'raw_available_in_conversation': 'yes_at_analysis_time; pending_user_deletion_confirmation',
    'retained_derivative_paths': str(DER),
    'derivative_sha256s': sha(DER),
    'parser_version_or_script_sha256': script_hashes,
    'raw_to_derived_reconciliation': '10,348 XML gene-set rows preserved; 1,454 C5 rows; 1,315,074 summed unique Entrez memberships; source and C5 row orders retained',
    'reacquisition_url': 'https://www.gsea-msigdb.org/gsea/downloads.jsp',
    'deletion_notes': 'Clearance applies to the conversation attachment after frozen bundle validation; user may keep a local/backed-up master because the source requires registration.',
})
der_row = {c: '' for c in cols}
der_row.update({
    'record_type': 'retained_input',
    'artifact_name': DER.name,
    'local_path': str(DER),
    'local_status': 'present_and_hash_verified',
    'origin_in_this_run': BATCH,
    'analysis_role': 'complete compact normalized representation of MSigDB v5.0 Entrez gene-set memberships and source order',
    'used_by': 'future MSigDB comparisons and provenance reproduction',
    'direct_or_canonical_source_url': 'https://www.gsea-msigdb.org/gsea/downloads.jsp',
    'source_page_url': 'https://docs.gsea-msigdb.org/MSigDB/Release_Notes/MSigDB_5.0/',
    'url_status': 'derived_from_verified_user_archive',
    'retrieval_status': 'generated_and_reconciled',
    'retrieved_at_utc': EVENT_TIME,
    'size_bytes': str(DER.stat().st_size),
    'sha256': sha(DER),
    'parent_or_derivation': f'{RAW.name}:{sha(RAW)}',
    'notes': 'Preserves names, systematic IDs, collection/subcollection, GO ID, chip namespace, URL, XML order, C5 order, and sorted unique Entrez memberships.',
    'batch_id': BATCH,
    'received_at_utc': EVENT_TIME,
    'integrity_verified_at_utc': EVENT_TIME,
    'analysis_completed_at_utc': EVENT_TIME,
    'raw_retention_status': 'required_compact_derivative',
    'raw_available_in_conversation': 'not_applicable_retained_derivative',
    'retained_derivative_paths': str(DER),
    'derivative_sha256s': sha(DER),
    'parser_version_or_script_sha256': f'{NORM.name}:{sha(NORM)}',
    'raw_to_derived_reconciliation': 'gzip valid; 10,348 rows + header; C5 1,454; membership sum 1,315,074',
    'reacquisition_url': 'https://www.gsea-msigdb.org/gsea/downloads.jsp',
})
manifest = pd.concat([manifest, pd.DataFrame([raw_row, der_row], columns=cols)], ignore_index=True)
man_path = RESULTS/f'actual_input_file_manifest_through_B104C_{DELIVERY}.csv'
man_md = RESULTS/f'actual_input_file_manifest_through_B104C_{DELIVERY}.md'
manifest.to_csv(man_path, index=False)
man_md.write_text('# Actual input manifest through B104C\n\n' + md_table(manifest), encoding='utf-8')

# Collection-level storage status update, intentionally not file-specific.
storage = pd.DataFrame([{
    'event_time_utc': EVENT_TIME,
    'scope': 'user-local storage',
    'artifact_group': 'previously held MSigDB archives predating the new v5.0 upload',
    'reported_status': 'deleted locally to free space',
    'file_level_status_changed': 'no',
    'reason': 'The user did not enumerate exact deleted filenames; preserve prior per-file records and add this collection-level event.',
    'conversation_runtime_status': 'previously uploaded copies were still available at B104C analysis time',
    'reupload_needed_now': 'no; compact results and cross-version metadata retained',
}])
storage_path = RESULTS/f'user_local_storage_status_update_B104C_{DELIVERY}.csv'
storage.to_csv(storage_path, index=False)

# Source ledger prebundle.
ledger = pd.read_csv(BASE_LEDGER, dtype=str).fillna('')
lcols = ledger.columns.tolist()

def lrow(**kwargs):
    r = {c: '' for c in lcols}
    r.update(kwargs)
    return r

new_ledger = []
new_ledger.append(lrow(
    record_type='actual_input', artifact_name=RAW.name, local_path=str(RAW),
    local_status='present_and_hash_verified_at_analysis_time', origin_in_this_run='user_upload_B104C',
    analysis_role='MSigDB v5.0 direct membership, cross-version C5 identity, and source-order investigation',
    used_by=f'{MAIN.name}; {NORM.name}', direct_or_canonical_source_url='https://www.gsea-msigdb.org/gsea/downloads.jsp',
    source_page_url='https://docs.gsea-msigdb.org/MSigDB/Release_Notes/MSigDB_5.0/',
    url_status='official_authenticated_download_page; direct_archive_URL_not_recorded',
    retrieval_status='supplied_by_user', size_bytes=str(RAW.stat().st_size), sha256=sha(RAW),
    parent_or_derivation='B000B full inventory declaration plus B104C received bytes',
    notes='Runtime hash and size match full user inventory. XML build date Apr 27, 2015.', batch_id=BATCH,
    user_local_relative_path=RAW.name, user_local_status='present_when_uploaded_B104C; later retention not assumed',
    user_local_size_bytes=str(RAW.stat().st_size), user_local_sha256=sha(RAW),
    inventory_file='local_upload_inventory_full_20260827T160408Z.csv',
    inventory_sha256='4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4',
    deletion_state='conversation_deletion_clearance_issued', event_recorded_at_utc=EVENT_TIME,
    supersedes_prior_record_type='historical_candidate_or_user_inventory_declaration_for_same_filename',
    hash_authority='container_sha256_matching_user_inventory', runtime_verification_status='size_hash_zip_xml_verified',
    planned_batch=BATCH, retained_derivative_paths=str(DER), parser_script_sha256=script_hashes,
    raw_to_derived_reconciliation='10,348 rows; 1,454 C5; membership sum 1,315,074',
))
new_ledger.append(lrow(
    record_type='retained_input', artifact_name=DER.name, local_path=str(DER), local_status='present_and_hash_verified',
    origin_in_this_run=BATCH, analysis_role='complete compact MSigDB v5.0 normalized Entrez gene-set table',
    used_by='future source and membership analyses', direct_or_canonical_source_url='https://www.gsea-msigdb.org/gsea/downloads.jsp',
    source_page_url='https://docs.gsea-msigdb.org/MSigDB/Release_Notes/MSigDB_5.0/', url_status='derived_from_verified_archive',
    retrieval_status='generated_and_reconciled', retrieved_at_utc=EVENT_TIME, size_bytes=str(DER.stat().st_size), sha256=sha(DER),
    parent_or_derivation=f'{RAW.name}:{sha(RAW)}', notes='All 10,348 XML rows retained with source order and normalized memberships.',
    batch_id=BATCH, deletion_state='retain_required_derivative', event_recorded_at_utc=EVENT_TIME,
    hash_authority='container_sha256', runtime_verification_status='gzip_and_reconciliation_passed',
    retained_derivative_paths=str(DER), parser_script_sha256=sha(NORM),
    raw_to_derived_reconciliation='10,348 rows; 1,454 C5; membership sum 1,315,074',
))
new_ledger.append(lrow(
    record_type='user_local_storage_status_update', artifact_name='previously held MSigDB archives (scope not itemized)',
    local_path='user_local_only', local_status='user_reported_deleted_locally_to_free_space', origin_in_this_run='user_statement_B104C',
    analysis_role='storage/provenance status only', used_by='batch planning and reacquisition decisions',
    source_page_url='https://www.gsea-msigdb.org/gsea/downloads.jsp', url_status='official_authenticated_download_page',
    retrieval_status='not_applicable', parent_or_derivation='user statement in B104C request',
    notes='Exact deleted filenames were not enumerated. Prior per-file records are not silently overwritten; this collection-level event records the statement faithfully.',
    batch_id=BATCH, user_local_status='collection_level_deletion_reported', deletion_state='user_local_deletion_reported_scope_not_itemized',
    event_recorded_at_utc=EVENT_TIME, hash_authority='user_statement', runtime_verification_status='not_applicable',
))
new_ledger.append(lrow(
    record_type='planned_batch_status', artifact_name='B105 late-June ontology robustness test', local_path='not_materialized_in_B104C',
    local_status='deferred_by_user', origin_in_this_run='B104C user instruction', analysis_role='future ontology/order robustness test',
    used_by='future B105', direct_or_canonical_source_url='https://release.geneontology.org/2016-07-01/ontology/go.obo',
    source_page_url='https://release.geneontology.org/2016-07-01/', url_status='official historical GO release',
    retrieval_status='not_downloaded; deferred', notes='B105 explicitly deferred because of storage constraints. B104C did not use a new ontology.',
    batch_id=BATCH, deletion_state='not_applicable', event_recorded_at_utc=EVENT_TIME, hash_authority='not_applicable',
    runtime_verification_status='not_materialized', planned_batch='B105',
))

# Metalink source records, adding explicit provenance endpoints used by the audit script.
for rel, date, size, md5 in [
    ('2016_04','2016-04-13','1516525310','e607b83de1ac87e6f63b13715c049a3f'),
    ('2016_05','2016-05-11','1504161063','fe9525832026b03ab34f0971b43c0c81'),
    ('2016_06','2016-06-08','1504963399','e3a5ac5a166efc95e9ad06465d5bd2c4'),
]:
    new_ledger.append(lrow(
        record_type='historical_candidate_verification_source', artifact_name=f'UniProt {rel} RELEASE.metalink',
        local_path='not_materialized_in_runtime', local_status='web_verified_source', origin_in_this_run=BATCH,
        analysis_role=f'official release, archive-size, and MD5 authority for sequential O95073/Q9Y620 audit ({date})',
        used_by=UNI.name,
        direct_or_canonical_source_url=f'https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-{rel}/knowledgebase/RELEASE.metalink',
        source_page_url=f'https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-{rel}/knowledgebase/',
        url_status='official_UniProt_archive_metadata', retrieval_status='inspected_via_web; script_will_download_locally',
        size_bytes='', sha256='', parent_or_derivation=f'uniprot_sprot-only{rel}.tar.gz expected bytes {size}; MD5 {md5}',
        notes='The audit script re-downloads and parses this official metalink before downloading the large archive.', batch_id=BATCH,
        deletion_state='not_materialized', event_recorded_at_utc=EVENT_TIME, hash_authority='official_metalink_MD5_and_size',
        runtime_verification_status='web_source_verified', planned_batch='future UniProt mapping audit',
    ))

# Generated outputs and scripts included in the prebundle ledger.
output_paths = [
    REPORT, DIAG, UNI_INST, receipt, UNIQUE_ORDER,
    ROOT/'analysis/B104C_analysis_summary_20260828T194921Z.json',
    ROOT/'analysis/B104C_msigdb_direct_membership_summary_20260828T194921Z.csv',
    ROOT/'analysis/B104C_msigdb_direct_best_per_label_20260828T194921Z.csv',
    ROOT/'analysis/B104C_graphsage_class_map_python2_dict_validation_20260828T194921Z.csv',
    ROOT/'analysis/B104C_column_order_model_comparison_20260828T194921Z.csv',
    ROOT/'analysis/B104C_duplicate_vector_GO_disambiguation_20260828T194921Z.csv',
    ROOT/'analysis/B104C_extended_python2_dictionary_order_simulation_grid_20260828T194921Z.csv',
    ROOT/'analysis/B104C_extended_order_best_model_20260828T194921Z.json',
    ROOT/'analysis/B104C_msigdb_cross_version_C5_target_presence_20260828T194921Z.csv',
    ROOT/'analysis/B104C_msigdb_v5.0_normalization_reconciliation_20260828T194921Z.json',
    MAIN, EXT, NORM, UNI,
]
for p in output_paths:
    new_ledger.append(lrow(
        record_type='generated_output', artifact_name=p.name, local_path=str(p), local_status='present_and_hash_verified',
        origin_in_this_run=BATCH, analysis_role='B104C analysis, reproducibility, provenance, or user delivery artifact',
        used_by='user delivery and future reproduction', retrieval_status='generated_in_runtime', retrieved_at_utc=EVENT_TIME,
        size_bytes=str(p.stat().st_size), sha256=sha(p), parent_or_derivation='B104A/B104B retained derivatives plus verified MSigDB v5.0 archive',
        notes='Included in B104C prebundle file inventory unless explicitly excluded.', batch_id=BATCH,
        deletion_state='retain_generated_output', event_recorded_at_utc=EVENT_TIME, hash_authority='container_sha256',
        runtime_verification_status='present_and_verified',
    ))
ledger = pd.concat([ledger, pd.DataFrame(new_ledger, columns=lcols)], ignore_index=True)
ledger_path = RESULTS/f'source_ledger_through_B104C_PREBUNDLE_{DELIVERY}.csv'
ledger_md = RESULTS/f'source_ledger_through_B104C_PREBUNDLE_{DELIVERY}.md'
ledger.to_csv(ledger_path, index=False)
ledger_md.write_text('# Source ledger through B104C — prebundle\n\n' + md_table(ledger), encoding='utf-8')

# Append-only provenance events prebundle.
events = pd.read_csv(BASE_EVENTS, dtype=str).fillna('')
evcols = events.columns.tolist()

def ev(event_type, artifact, status, details):
    return {'event_time_utc': EVENT_TIME, 'batch_id': BATCH, 'event_type': event_type,
            'artifact_name': artifact, 'status': status, 'details': details}
new_events = [
    ev('raw_input_received', RAW.name, 'accepted', f'{RAW.stat().st_size} bytes; SHA-256 {sha(RAW)}; matches B000B full inventory.'),
    ev('raw_input_integrity_verified', RAW.name, 'accepted', 'ZIP parsed; one XML file found; internal XML BUILD_DATE Apr 27, 2015.'),
    ev('msigdb_v5_direct_membership_test', RAW.name, 'accepted', '10,348 all sets and 1,454 C5 sets tested; zero exact GraphSAGE label columns; closest mismatch 307 genes.'),
    ev('msigdb_C5_scope_correction', 'cross-version MSigDB target identity counts', 'accepted_superseding_prior_count', 'C5-scoped target GO presence is v5.0=57, v5.1=57, v5.2=6, v6.0=6. Prior v5.2/v6.0 count of 23 counted non-C5 XML rows.'),
    ev('graphsage_python2_serialization_validation', 'ppi-class_map.json', 'accepted', '64-bit unrandomized CPython2-style dictionary simulation reproduces all 56,944 JSON keys exactly; 32-bit model matches only 5,410 positions.'),
    ev('column_order_conventional_models_tested', 'GraphSAGE 121-label GO order', 'accepted_negative_result', 'GO ID, name, prevalence, OBO order, GAF first occurrence, and MSigDB XML orders do not explain the sequence.'),
    ev('column_order_legacy_dictionary_model', 'GraphSAGE 121-label GO order', 'accepted_strong_inference', 'Large accepted-GAF Python2 dictionary model: tau 0.772176, p 3.54e-36, 88.61% pairwise concordance, LCS 92/121, first five exact.'),
    ev('column_order_extended_grid', '48 GAF-derived Python2 simulations', 'accepted', 'Best LCS 94/121; all simulations had a 32,768-slot table and selected duplicate orientation 001.'),
    ev('duplicate_vector_order_disambiguation', UNIQUE_ORDER.name, 'strongly_supported_not_source_code_proven', 'Columns 24/71 = GO:0043228/GO:0043232; 39/63 = GO:0006464/GO:0036211; 48/70 = GO:1903561/GO:0043230.'),
    ev('msigdb_v5_normalized_derivative_created', DER.name, 'accepted', f'10,348 rows; 1,454 C5; SHA-256 {sha(DER)}.'),
    ev('uniprot_sequential_audit_script_created', UNI.name, 'accepted', 'Processes 2016_04, 2016_05, and 2016_06 sequentially; verifies metalink size/MD5, records SHA-256, extracts O95073/Q9Y620, deletes archive only after validation.'),
    ev('uniprot_audit_script_self_test', UNI.name, 'passed', 'Synthetic archive with 100,004 records scanned; both target records extracted; post-validation archive deletion gate passed.'),
    ev('uniprot_audit_script_dry_run', UNI.name, 'passed', 'All three official release URLs, sizes, MD5 values, and scratch paths printed; no large download initiated.'),
    ev('user_local_storage_deletion_report', 'previously held MSigDB archives', 'recorded_collection_level', 'User reported deleting older MSigDB files to free space; exact filenames not enumerated, so per-file local statuses were not guessed.'),
    ev('B105_status', 'late-June ontology robustness test', 'deferred_by_user', 'No B105 input requested or analyzed in B104C.'),
    ev('deletion_clearance_issued', RAW.name, 'safe_after_bundle_validation', f'Complete normalized derivative retained; raw SHA-256 {sha(RAW)}.'),
]
events = pd.concat([events, pd.DataFrame(new_events, columns=evcols)], ignore_index=True)
events_path = RESULTS/f'provenance_events_through_B104C_PREBUNDLE_{DELIVERY}.csv'
events.to_csv(events_path, index=False)

# Copy key manifests into batch for frozen bundle.
for src in [man_path, man_md, ledger_path, ledger_md, events_path, storage_path]:
    target = ROOT/src.name
    target.write_bytes(src.read_bytes())

print(json.dumps({
    'manifest': str(man_path), 'ledger_prebundle': str(ledger_path), 'events_prebundle': str(events_path),
    'receipt': str(receipt), 'storage_update': str(storage_path),
    'rows': {'manifest': len(manifest), 'ledger': len(ledger), 'events': len(events)},
}, indent=2))

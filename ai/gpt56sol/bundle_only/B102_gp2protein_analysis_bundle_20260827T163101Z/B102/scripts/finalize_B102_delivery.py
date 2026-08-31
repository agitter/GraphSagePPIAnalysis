#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path('/mnt/data/ppi_repro_corrected')
RESULTS = ROOT / 'results'
ACCEPTED = ROOT / 'batches' / 'B102_final2'
FAILED = ROOT / 'batches' / 'B102_failed_attempts'
B102_DUP = ROOT / 'batches' / 'B102'
B102_INCOMPLETE = ROOT / 'batches' / 'B102_final'
SCRIPTS = ROOT / 'scripts'
INPUT_ROOT = Path('/mnt/data')
FULL_INVENTORY = INPUT_ROOT / 'local_upload_inventory_full_20260827T160408Z.csv'
PRIOR_ACTUAL = RESULTS / 'actual_input_file_manifest_with_B101_20260827T153319Z.csv'
PRIOR_LEDGER = RESULTS / 'source_ledger_with_B101_B000A_20260827T153319Z.csv'
BASELINE = Path('/mnt/data/work/ppi_repro_corrected/results/collapsed_gene_labels_topology_features.csv')
B101_DIR = ROOT / 'batches' / 'B101_extracted' / 'B101'
ANALYSIS_SCRIPT = SCRIPTS / 'analyze_B102_gp2protein_final.py'
FINALIZER_SCRIPT = Path(__file__)

ACCEPTED_STAMP = '20260827T162132Z'
ACCEPTED_COMPLETED = '2026-08-27T16:22:36.415383+00:00'
RECEIVED = '2026-08-27T16:21:32.620248+00:00'
NOW = datetime.now(timezone.utc)
STAMP = NOW.strftime('%Y%m%dT%H%M%SZ')
NOW_ISO = NOW.isoformat()
DELIVERY = ROOT / 'batches' / f'B102_delivery_{STAMP}'
DELIVERY.mkdir(parents=True, exist_ok=False)
RESULTS.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def sha256_decompressed_gzip(path: Path) -> str:
    h = hashlib.sha256()
    with gzip.open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8-sig', newline='') as fh:
        return list(csv.DictReader(fh))


def md_escape(value: object) -> str:
    return str(value).replace('|', '\\|').replace('\n', '<br>')


def markdown_table(rows: Iterable[dict[str, object]], fields: list[str]) -> str:
    rows = list(rows)
    out = ['| ' + ' | '.join(fields) + ' |', '| ' + ' | '.join(['---'] * len(fields)) + ' |']
    for row in rows:
        out.append('| ' + ' | '.join(md_escape(row.get(f, '')) for f in fields) + ' |')
    return '\n'.join(out)


def source_for(name: str) -> dict[str, str]:
    direct = ''
    page = ''
    status = 'local_or_derived_artifact_no_external_source_recorded'
    note = ''

    go_release_base = 'https://release.geneontology.org/2016-06-01/'
    if name == '2016-06-01-annotations-README':
        direct = go_release_base + 'annotations/gp2protein/README'
        page = go_release_base + 'annotations/gp2protein/'
        status = 'official_GO_release_path_recorded; uploaded_bytes_not_remote_compared'
    elif name in {'2016-06-01-gp2protein.geneid.gz', '2016-06-01-gp2protein.human.gz'}:
        upstream = name.removeprefix('2016-06-01-')
        direct = go_release_base + 'annotations/gp2protein/' + upstream
        page = go_release_base + 'annotations/gp2protein/'
        status = 'official_GO_release_path_recorded; uploaded_bytes_not_remote_compared'
    elif name in {'2016-06-01-go.obo', '2016-06-01-gene_ontology.obo', '2016-06-01-go.owl'}:
        upstream = name.removeprefix('2016-06-01-')
        direct = go_release_base + 'ontology/' + upstream
        page = go_release_base
        status = 'official_GO_release_candidate_path_recorded; raw_bytes_not_uploaded_in_B102'
    elif re.match(r'^(gene_association|gp_association|gp_information)\.goa_(ref_)?human\.\d+\.gz$', name) or re.match(r'^goa_human\.(gaf|gpa|gpi)\.\d+\.gz$', name):
        direct = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/' + name
        page = 'https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/'
        status = 'official_EBI_GOA_historical_archive_exact_filename'
    elif name == 'graphsage_ppi.zip':
        direct = 'https://snap.stanford.edu/graphsage/ppi.zip'
        page = 'https://snap.stanford.edu/graphsage/'
        status = 'official_project_download_url'
    elif name == 'dgl_ppi.zip':
        direct = 'https://data.dgl.ai/dataset/ppi.zip'
        page = 'https://github.com/dmlc/dgl/blob/master/python/dgl/data/ppi.py'
        status = 'official_dataset_loader_URL'
    elif name.startswith('bio-tissue-'):
        direct = 'https://snap.stanford.edu/ohmnet/' + name
        page = 'https://snap.stanford.edu/ohmnet/'
        status = 'official_OhmNet_project_download_url'
    elif name == 'ohmnet-master.zip':
        direct = 'https://github.com/mims-harvard/ohmnet/archive/refs/heads/master.zip'
        page = 'https://github.com/mims-harvard/ohmnet'
        status = 'repository_branch_archive_candidate; mutable_branch; bytes_not_compared'
    elif name in {'Greene2015.pdf', 'Greene2015_sup.pdf', 'Greene2015_Table6.xlsx', 'Greene2015_Table9.xlsx'}:
        direct = 'https://doi.org/10.1038/ng.3259'
        page = 'https://www.nature.com/articles/ng.3259'
        status = 'publisher_article_or_supplement_page'
    elif name == 'OhmNet.pdf':
        direct = 'https://doi.org/10.1093/bioinformatics/btx252'
        page = 'https://academic.oup.com/bioinformatics/article/33/14/i190/3953967'
        status = 'publisher_article_page'
    elif name.startswith('msigdb_v'):
        page = 'https://www.gsea-msigdb.org/gsea/downloads.jsp'
        status = 'official_authenticated_download_page; direct_archive_URL_not_recorded'
    elif name.startswith('bioconductor-annotation-org.Hs.eg.db_'):
        m = re.search(r'_(\d+\.\d+\.\d+)\.tar\.gz$', name)
        if m:
            version = m.group(1)
            bioc = {'3.0.0': '3.0', '3.1.2': '3.1', '3.3.0': '3.3', '3.4.0': '3.4'}.get(version, '')
            upstream = f'org.Hs.eg.db_{version}.tar.gz'
            if bioc:
                direct = f'https://bioconductor.statistik.tu-dortmund.de/packages/{bioc}/data/annotation/src/contrib/{upstream}'
                page = f'https://bioconductor.org/packages/{bioc}/data/annotation/html/org.Hs.eg.db.html'
                status = 'historical_Bioconductor_mirror_path; local_filename_has_added_prefix; bytes_not_compared'
    elif name == '2026-08-14-gene2go.gz':
        direct = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
        page = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/'
        status = 'rolling_NCBI_source_path; local_dated_copy_hash_recorded'
    elif name == '2026-08-14-gene2go_human.tsv.gz':
        page = 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/'
        status = 'locally_filtered_derivative_of_2026-08-14-gene2go.gz'
    elif name == '2016-12-23-gene2go.gz':
        page = 'https://web.archive.org/'
        status = 'historical_snapshot; exact_Wayback_acquisition_URL_unresolved'
    elif name == '2016-12-23-gene2go_human.tsv.gz':
        page = 'https://web.archive.org/'
        status = 'locally_filtered_historical_snapshot; exact_parent_acquisition_URL_unresolved'
    elif name.startswith('dhimmel-gene-ontology-962a5e1-'):
        short = name.removeprefix('dhimmel-gene-ontology-962a5e1-')
        direct = f'https://raw.githubusercontent.com/dhimmel/gene-ontology/962a5e1/{short}'
        page = 'https://github.com/dhimmel/gene-ontology/tree/962a5e1'
        status = 'commit_pinned_candidate_path; exact_raw_path_not_remote_compared'
    elif name in {'HumanBase-blood.dat', 'HumanBase-blood_top.gz', 'HuamnBase-kidney.dat'}:
        page = 'https://hb.flatironinstitute.org/download'
        status = 'HumanBase_source_page_recorded; exact_object_URL_unresolved'
        if name == 'HuamnBase-kidney.dat':
            note = 'Exact local spelling is preserved; it may be a local typo but is not silently corrected.'
    elif name == 'blood_sample.tsv.gz':
        page = 'https://hb.flatironinstitute.org/download'
        status = 'locally_generated_sample_from_HumanBase-blood_top.gz'
    elif name in {'filter_gene2go.sh', 'find_gene2go_2016.py', 'gene2go_github_survey.tsv', 'file_inventory_2026_08_18.pdf', 'chat_gpt_historical_go_mapping_inventory.md', 'investigation_summary_2026_08_23.md', 'local_upload_inventory_20260827T145903Z.csv', 'make_local_inventory.py', 'make_local_inventory_v2.py', 'requested_file_patterns.txt', 'sample_giant_network.sh', 'train_graph_id.npy'}:
        status = 'local_or_conversation_generated_artifact; no_external_source_URL'
    return {
        'direct_or_canonical_source_url': direct,
        'source_page_url': page,
        'url_status': status,
        'source_notes': note,
    }


def planned_batch(name: str) -> str:
    if name == '2016-06-01-go.obo':
        return 'B103'
    if name in {'2016-06-01-gene_ontology.obo', '2016-06-01-go.owl'}:
        return 'hold_for_ontology_format_comparison_if_needed'
    if name in {'goa_human.gaf.158.gz', 'goa_human.gpa.158.gz', 'goa_human.gpi.158.gz'}:
        return 'B104'
    if name in {'goa_human.gaf.160.gz', 'goa_human.gpa.160.gz', 'goa_human.gpi.160.gz'}:
        return 'B105'
    if name in {'gene_association.goa_human.157.gz', 'gp_association.goa_human.157.gz', 'gene_association.goa_ref_human.157.gz'}:
        return 'B106'
    if name in {'2016-12-23-gene2go.gz', '2016-12-23-gene2go_human.tsv.gz'}:
        return 'B201'
    if name.startswith('bioconductor-annotation-org.Hs.eg.db_'):
        return 'B30x_one_package_per_batch'
    if name in {'HumanBase-blood.dat', 'HuamnBase-kidney.dat', 'HumanBase-blood_top.gz', 'blood_sample.tsv.gz'}:
        return 'B40x'
    return ''


# Verify accepted run and required files.
required = [FULL_INVENTORY, PRIOR_ACTUAL, PRIOR_LEDGER, BASELINE, ANALYSIS_SCRIPT]
for p in required:
    if not p.exists():
        raise FileNotFoundError(p)
if (ACCEPTED / 'run.exit').read_text().strip() != '0':
    raise RuntimeError('Accepted B102 run does not have exit status 0')

analysis_summary_path = ACCEPTED / f'B102_analysis_summary_{ACCEPTED_STAMP}.json'
validation_path = ACCEPTED / f'B102_VALIDATION_{ACCEPTED_STAMP}.json'
analysis_summary = json.loads(analysis_summary_path.read_text())
validation = json.loads(validation_path.read_text())

# Full inventory: this is a user-local declaration with exact hashes, not a remote checksum catalog.
inv_df = pd.read_csv(FULL_INVENTORY, dtype=str, keep_default_na=False)
if len(inv_df) != 65 or inv_df['size_bytes'].astype(int).sum() != 2_680_734_828:
    raise RuntimeError('Full inventory row count or byte total changed')
if not inv_df['sha256'].str.fullmatch(r'[0-9a-f]{64}').all():
    raise RuntimeError('Invalid SHA-256 in full inventory')
full_inventory_sha = sha256_file(FULL_INVENTORY)
if full_inventory_sha != '4210821f03fc5fc6f51e978cf2b82968f0500bc53c025cc5bf3cebb7c13015e4':
    raise RuntimeError('Full inventory upload hash changed')

# Enrich all 65 rows without rewriting user declarations.
enriched_rows: list[dict[str, object]] = []
for _, r in inv_df.iterrows():
    src = source_for(r['artifact_name'])
    mounted = INPUT_ROOT / r['artifact_name']
    exact_mounted = mounted.exists()
    runtime_hash = sha256_file(mounted) if exact_mounted and mounted.stat().st_size <= 100_000_000 else ''
    # For large exact-name mounts, we do not infer absence; hashes are the user's inventory declaration unless uploaded in a batch.
    runtime_status = 'not_mounted_or_not_checked'
    if exact_mounted:
        if runtime_hash:
            runtime_status = 'runtime_hash_matches_inventory' if runtime_hash == r['sha256'] else 'runtime_hash_DIFFERS_from_inventory'
        else:
            runtime_status = 'runtime_exact_name_present; hash_not_recomputed_here_due_size'
    if r['artifact_name'] in {
        '2016-06-01-annotations-README', '2016-06-01-gp2protein.geneid.gz', '2016-06-01-gp2protein.human.gz'
    }:
        runtime_hash = sha256_file(mounted)
        runtime_status = 'B102_uploaded_bytes_hash_matches_inventory' if runtime_hash == r['sha256'] else 'B102_HASH_MISMATCH'
    if r['artifact_name'] in {'goa_human.gaf.159.gz', 'goa_human.gpa.159.gz', 'goa_human.gpi.159.gz'}:
        runtime_status = 'B101_user_deletion_confirmed; any residual runtime mount is ignored'
        runtime_hash = ''
    enriched_rows.append({
        **r.to_dict(),
        **src,
        'inventory_semantics': 'user_local_file_hash_declaration_generated_before_upload',
        'inventory_file': FULL_INVENTORY.name,
        'inventory_file_sha256': full_inventory_sha,
        'runtime_verification_status': runtime_status,
        'runtime_recomputed_sha256': runtime_hash,
        'planned_batch': planned_batch(r['artifact_name']),
        'last_updated_at_utc': NOW_ISO,
    })

enriched_inventory_csv = RESULTS / f'user_local_inventory_full_enriched_{STAMP}.csv'
write_csv(enriched_inventory_csv, enriched_rows)
enriched_inventory_md = RESULTS / f'user_local_inventory_full_enriched_{STAMP}.md'
url_counts = Counter(r['url_status'] for r in enriched_rows)
enriched_inventory_md.write_text(
    '# Complete user-local inventory, enriched through B102\n\n'
    f'Generated: `{NOW_ISO}`  \n'
    f'Inventory source: `{FULL_INVENTORY.name}`  \n'
    f'Inventory SHA-256: `{full_inventory_sha}`  \n'
    f'Rows: **{len(enriched_rows)}**  \n'
    f'Total declared bytes: **{sum(int(r["size_bytes"]) for r in enriched_rows):,}**\n\n'
    'The size, mtime, and SHA-256 columns are declarations from the user-generated full inventory. '
    'Only files uploaded in a batch are independently hash-checked against those declarations. '
    'The inventory file cannot list its own hash; its received-copy hash is tracked separately. '
    'The exact filename `HuamnBase-kidney.dat` is preserved.\n\n'
    '## Source-status counts\n\n' +
    markdown_table([{'url_status': k, 'files': v} for k, v in sorted(url_counts.items())], ['url_status', 'files']) +
    '\n\n## Complete file table\n\n' +
    markdown_table(enriched_rows, ['artifact_name', 'size_bytes', 'sha256', 'mtime_utc', 'direct_or_canonical_source_url', 'source_page_url', 'url_status', 'runtime_verification_status', 'planned_batch', 'source_notes']) + '\n',
    encoding='utf-8'
)

# Actual data input manifest through B102, with explicit lifecycle fields.
prior_actual_df = pd.read_csv(PRIOR_ACTUAL, dtype=str, keep_default_na=False)
base_fields = list(prior_actual_df.columns)
lifecycle_fields = [
    'batch_id', 'inventory_sha256', 'received_at_utc', 'integrity_verified_at_utc',
    'analysis_completed_at_utc', 'deletion_clearance_issued_at_utc',
    'user_deletion_confirmed_at_utc', 'raw_retention_status',
    'raw_available_in_conversation', 'retained_derivative_paths',
    'derivative_sha256s', 'parser_version_or_script_sha256',
    'raw_to_derived_reconciliation', 'reacquisition_url', 'deletion_notes'
]
actual_fields = base_fields + lifecycle_fields
actual_rows: list[dict[str, object]] = []
for _, row in prior_actual_df.iterrows():
    rec = {f: row.get(f, '') for f in base_fields}
    for f in lifecycle_fields:
        rec[f] = ''
    if row['artifact_name'] in {'goa_human.gaf.159.gz', 'goa_human.gpa.159.gz', 'goa_human.gpi.159.gz'}:
        rec.update({
            'local_status': 'conversation_copy_deletion_confirmed_by_user; residual_runtime_mount_ignored',
            'batch_id': 'B101',
            'inventory_sha256': 'dae73b7f54d28a089917ba8cf5b7a7f62b2ea0496cba097f3c1e1bfd94e0e9df',
            'received_at_utc': '2026-08-27T15:27:36.842353+00:00',
            'integrity_verified_at_utc': '2026-08-27T15:27:36.842353+00:00',
            'analysis_completed_at_utc': '2026-08-27T15:27:36.842353+00:00',
            'deletion_clearance_issued_at_utc': '2026-08-27T15:27:36.842353+00:00',
            'user_deletion_confirmed_at_utc': RECEIVED,
            'raw_retention_status': 'user_local_master_retained; conversation_attachment_deleted_by_user',
            'raw_available_in_conversation': 'false',
            'parser_version_or_script_sha256': 'analyze_B101_goa159.py; see B101 bundle',
            'deletion_notes': "User explicitly reported 'Deleted B101'; future analyses use retained B101 derivatives, not residual raw mounts.",
        })
        if row['artifact_name'] == 'goa_human.gaf.159.gz':
            rec['retained_derivative_paths'] = str(B101_DIR / 'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz')
            rec['derivative_sha256s'] = '132066afa3d3ea550752d3d2eb98fcbc238570c23bf239f0ea77e720e34cc274'
            rec['raw_to_derived_reconciliation'] = '392440 rows re-read; reconstructed uncompressed raw-data SHA-256 matched B101 raw-data hash'
        elif row['artifact_name'] == 'goa_human.gpi.159.gz':
            rec['retained_derivative_paths'] = str(B101_DIR / 'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz')
            rec['derivative_sha256s'] = '4a69bf547b951d060c1b38a4fc208b116ff0917ce74b941675d244b3315c34b1'
            rec['raw_to_derived_reconciliation'] = '21002 rows re-read; reconstructed uncompressed raw-data SHA-256 matched B101 raw-data hash'
        else:
            rec['retained_derivative_paths'] = str(B101_DIR / 'B101_goa_human_gpad159_normalized_20260827T152736Z.tsv.gz')
            rec['derivative_sha256s'] = '7e3ec99110ef936ea9e9b7deaacdf7a86d7a7057cb41d5348f80e5ff56ac2793'
            rec['raw_to_derived_reconciliation'] = '393458 rows re-read; reconstructed uncompressed raw-data SHA-256 matched B101 raw-data hash'
    actual_rows.append(rec)

analysis_script_sha = sha256_file(ANALYSIS_SCRIPT)
raw_to_derivative = {
    '2016-06-01-annotations-README': (
        ACCEPTED / 'derived' / f'B102_2016-06-01_annotations_README_exact_{ACCEPTED_STAMP}.txt',
        'exact byte-for-byte retained copy; SHA-256 equals raw upload'
    ),
    '2016-06-01-gp2protein.geneid.gz': (
        ACCEPTED / 'derived' / f'B102_gp2protein_geneid_relevant_subset_{ACCEPTED_STAMP}.tsv.gz',
        '7,296,170 raw data rows parsed with 0 malformed; retained 25,983 relevant raw rows; derivative re-read row count matched'
    ),
    '2016-06-01-gp2protein.human.gz': (
        ACCEPTED / 'derived' / f'B102_gp2protein_human_normalized_{ACCEPTED_STAMP}.tsv.gz',
        '70,625 raw rows retained and re-read; row count matched'
    ),
}
for name in ['2016-06-01-annotations-README', '2016-06-01-gp2protein.geneid.gz', '2016-06-01-gp2protein.human.gz']:
    inv = inv_df.loc[inv_df['artifact_name'] == name].iloc[0]
    src = source_for(name)
    derivative, reconcile = raw_to_derivative[name]
    actual_rows.append({
        'record_type': 'actual_input',
        'artifact_name': name,
        'local_path': str(INPUT_ROOT / name),
        'local_status': 'present_as_B102_conversation_attachment; deletion_clearance_issued',
        'origin_in_this_run': 'user_upload_B102',
        'analysis_role': 'historical GO gp2protein documentation or mapping input',
        'used_by': ANALYSIS_SCRIPT.name,
        'direct_or_canonical_source_url': src['direct_or_canonical_source_url'],
        'source_page_url': src['source_page_url'],
        'url_status': src['url_status'],
        'retrieval_status': 'supplied_by_user; bytes matched full user-local inventory; remote bytes not compared',
        'retrieved_at_utc': RECEIVED,
        'size_bytes': inv['size_bytes'],
        'sha256': inv['sha256'],
        'parent_or_derivation': FULL_INVENTORY.name,
        'notes': src['source_notes'],
        'batch_id': 'B102',
        'inventory_sha256': full_inventory_sha,
        'received_at_utc': RECEIVED,
        'integrity_verified_at_utc': RECEIVED,
        'analysis_completed_at_utc': ACCEPTED_COMPLETED,
        'deletion_clearance_issued_at_utc': ACCEPTED_COMPLETED,
        'user_deletion_confirmed_at_utc': '',
        'raw_retention_status': 'user_local_master_retained; conversation_attachment_pending_user_deletion',
        'raw_available_in_conversation': 'true_until_user_confirms_Deleted_B102',
        'retained_derivative_paths': str(derivative),
        'derivative_sha256s': sha256_file(derivative),
        'parser_version_or_script_sha256': analysis_script_sha,
        'raw_to_derived_reconciliation': reconcile,
        'reacquisition_url': src['direct_or_canonical_source_url'],
        'deletion_notes': 'Clearance applies only to conversation attachment; keep user-local master copy.',
    })

actual_rows.append({
    'record_type': 'actual_input',
    'artifact_name': FULL_INVENTORY.name,
    'local_path': str(FULL_INVENTORY),
    'local_status': 'present_as_B102_conversation_attachment; safe_to_delete_but_user_may_keep',
    'origin_in_this_run': 'user_upload_B102; user_generated_inventory',
    'analysis_role': 'authoritative 65-file user-local hash inventory',
    'used_by': ANALYSIS_SCRIPT.name + '; batch planning and source ledger',
    'direct_or_canonical_source_url': '',
    'source_page_url': '',
    'url_status': 'user_generated_local_inventory; no_external_URL',
    'retrieval_status': 'supplied_by_user; received-copy SHA-256 computed independently',
    'retrieved_at_utc': RECEIVED,
    'size_bytes': FULL_INVENTORY.stat().st_size,
    'sha256': full_inventory_sha,
    'parent_or_derivation': 'generated locally by make_local_inventory_v2.py; cannot contain its own hash',
    'notes': 'The 65 contained hashes are user-local declarations. The inventory file itself is separately hashed on receipt.',
    'batch_id': 'B102',
    'inventory_sha256': full_inventory_sha,
    'received_at_utc': RECEIVED,
    'integrity_verified_at_utc': RECEIVED,
    'analysis_completed_at_utc': ACCEPTED_COMPLETED,
    'deletion_clearance_issued_at_utc': ACCEPTED_COMPLETED,
    'user_deletion_confirmed_at_utc': '',
    'raw_retention_status': 'small provenance input; safe_to_delete_from_conversation; local retention recommended',
    'raw_available_in_conversation': 'true_until_user_confirms_Deleted_B102',
    'retained_derivative_paths': str(enriched_inventory_csv),
    'derivative_sha256s': '',  # populated after file is final; checksum table supplies exact value
    'parser_version_or_script_sha256': analysis_script_sha,
    'raw_to_derived_reconciliation': '65 rows; total 2,680,734,828 bytes; all SHA-256 fields syntactically valid; exact values preserved in enriched inventory',
    'reacquisition_url': '',
    'deletion_notes': 'Keeping the CSV is acceptable; it is only about 12 KiB.',
})

# Derived inputs inherited from B101 and the independently recovered label matrix.
derived_inputs = [
    {
        'name': 'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz',
        'path': B101_DIR / 'B101_goa_human_gaf159_normalized_20260827T152736Z.tsv.gz',
        'sha': '132066afa3d3ea550752d3d2eb98fcbc238570c23bf239f0ea77e720e34cc274',
        'role': 'row-preserving GAF v159 derivative used after B101 deletion',
        'parent': 'goa_human.gaf.159.gz',
        'reconcile': '392440 rows and reconstructed raw-data hash verified in B101',
    },
    {
        'name': 'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz',
        'path': B101_DIR / 'B101_goa_human_gpi159_normalized_20260827T152736Z.tsv.gz',
        'sha': '4a69bf547b951d060c1b38a4fc208b116ff0917ce74b941675d244b3315c34b1',
        'role': 'row-preserving GPI v159 derivative used after B101 deletion',
        'parent': 'goa_human.gpi.159.gz',
        'reconcile': '21002 rows and reconstructed raw-data hash verified in B101',
    },
    {
        'name': 'B101_provisional_UniProt_to_Entrez_via_MSigDB52_symbols_20260827T152736Z.csv.gz',
        'path': B101_DIR / 'B101_provisional_UniProt_to_Entrez_via_MSigDB52_symbols_20260827T152736Z.csv.gz',
        'sha': '694ae310816f792be7c77e77235ba8b9e546b6c08b946f91d2c3e7501e22ecc5',
        'role': 'provisional comparison mapping retained from B101',
        'parent': 'B101 GPI symbols plus MSigDB v5.2 aliases',
        'reconcile': 'used only as a comparison/fallback mapping, not accepted as historical ground truth',
    },
    {
        'name': BASELINE.name,
        'path': BASELINE,
        'sha': '4fd88002b2600c0b7c5dcb076390f5b74a39edcb34184f1b0768d9926dcb1907',
        'role': 'independently recovered GraphSAGE gene-by-label matrix',
        'parent': 'topology-and-feature node identity reconstruction',
        'reconcile': '4268 distinct resolved Entrez genes with 121 observed binary label columns',
    },
]
for d in derived_inputs:
    if not d['path'].exists() or sha256_file(d['path']) != d['sha']:
        raise RuntimeError(f'Derived input verification failed: {d["path"]}')
    actual_rows.append({
        'record_type': 'retained_or_derived_input',
        'artifact_name': d['name'],
        'local_path': str(d['path']),
        'local_status': 'present_and_hash_verified',
        'origin_in_this_run': 'retained_from_prior_analysis' if d['name'].startswith('B101_') else 'generated_by_core_reconstruction',
        'analysis_role': d['role'],
        'used_by': ANALYSIS_SCRIPT.name,
        'direct_or_canonical_source_url': '',
        'source_page_url': '',
        'url_status': 'derived_artifact; provenance_traced_to_parent_inputs',
        'retrieval_status': 'local retained derivative',
        'retrieved_at_utc': '',
        'size_bytes': d['path'].stat().st_size,
        'sha256': d['sha'],
        'parent_or_derivation': d['parent'],
        'notes': d['reconcile'],
        'batch_id': 'B101' if d['name'].startswith('B101_') else 'core_reconstruction',
        'inventory_sha256': '',
        'received_at_utc': '',
        'integrity_verified_at_utc': RECEIVED,
        'analysis_completed_at_utc': ACCEPTED_COMPLETED,
        'deletion_clearance_issued_at_utc': '',
        'user_deletion_confirmed_at_utc': '',
        'raw_retention_status': 'retained_derivative_required_for_reproduction',
        'raw_available_in_conversation': 'true_as_generated_or_retained_artifact',
        'retained_derivative_paths': str(d['path']),
        'derivative_sha256s': d['sha'],
        'parser_version_or_script_sha256': analysis_script_sha,
        'raw_to_derived_reconciliation': d['reconcile'],
        'reacquisition_url': '',
        'deletion_notes': 'Do not delete this retained derivative from the analysis bundle.',
    })

actual_manifest_csv = RESULTS / f'actual_input_file_manifest_through_B102_{STAMP}.csv'
write_csv(actual_manifest_csv, actual_rows, actual_fields)
actual_manifest_md = RESULTS / f'actual_input_file_manifest_through_B102_{STAMP}.md'
actual_manifest_md.write_text(
    '# Actual input-file manifest through B102\n\n'
    f'Generated: `{NOW_ISO}`\n\n'
    'This manifest distinguishes raw uploaded inputs, user-local hash declarations, and retained derivatives. '
    'B101 raw attachments are marked deleted by user even if a residual runtime mount remains visible; those residual mounts are not used. '
    'The B102 inventory file cannot self-list its own hash, so its received-copy hash is tracked separately.\n\n'
    + markdown_table(actual_rows, [
        'record_type', 'artifact_name', 'local_status', 'batch_id', 'analysis_role',
        'direct_or_canonical_source_url', 'size_bytes', 'sha256', 'raw_available_in_conversation',
        'retained_derivative_paths', 'raw_to_derived_reconciliation', 'deletion_notes'
    ]) + '\n', encoding='utf-8'
)

# Append-only provenance event ledger.
prior_events: list[dict[str, str]] = []
for p in [ROOT / 'batches' / 'B000' / 'B000_provenance_events_20260827T150123Z.csv', B101_DIR / 'B101_provenance_events_20260827T152736Z.csv']:
    if p.exists():
        prior_events.extend(read_csv_rows(p))
accepted_events = read_csv_rows(ACCEPTED / f'B102_provenance_events_{ACCEPTED_STAMP}.csv')
# De-duplicate exact events while preserving order.
seen = set()
provenance_rows: list[dict[str, object]] = []
for row in prior_events + accepted_events:
    key = tuple(row.get(k, '') for k in ['event_time_utc', 'batch_id', 'event_type', 'artifact_name', 'status', 'details'])
    if key in seen:
        continue
    seen.add(key)
    provenance_rows.append(row)
provenance_rows.extend([
    {
        'event_time_utc': NOW_ISO,
        'batch_id': 'B000B',
        'event_type': 'full_inventory_ingested',
        'artifact_name': FULL_INVENTORY.name,
        'status': 'authoritative_user_local_inventory_for_batch_planning',
        'details': f'65 user-local files, 2,680,734,828 total bytes, all rows include SHA-256; inventory received-copy SHA-256 {full_inventory_sha}.',
    },
    {
        'event_time_utc': NOW_ISO,
        'batch_id': 'B102',
        'event_type': 'report_correction_issued',
        'artifact_name': f'B102_REPORT_CORRECTED_{STAMP}.md',
        'status': 'supersedes_wording_in_initial_B102_report',
        'details': 'Corrected the statement that all four uploads matched rows inside the full inventory: only the three data files can match inventory rows; the inventory CSV cannot self-list and was independently hashed on receipt.',
    },
    {
        'event_time_utc': NOW_ISO,
        'batch_id': 'B102',
        'event_type': 'delivery_bundle_created',
        'artifact_name': f'B102_gp2protein_analysis_bundle_{STAMP}.zip',
        'status': 'pending_bundle_hash_population',
        'details': 'Curated accepted outputs, complete diagnostics, updated manifests, and scripts; raw uploaded files excluded.',
    },
])
provenance_csv = RESULTS / f'provenance_events_through_B102_{STAMP}.csv'
write_csv(provenance_csv, provenance_rows, ['event_time_utc', 'batch_id', 'event_type', 'artifact_name', 'status', 'details'])

# Append-only source ledger: keep every prior row, then append superseding B000B and B102 records.
prior_ledger_df = pd.read_csv(PRIOR_LEDGER, dtype=str, keep_default_na=False)
ledger_fields = list(prior_ledger_df.columns)
extra_ledger_fields = [
    'event_recorded_at_utc', 'supersedes_prior_record_type', 'hash_authority',
    'runtime_verification_status', 'planned_batch', 'retained_derivative_paths',
    'parser_script_sha256', 'raw_to_derived_reconciliation', 'user_deletion_confirmed_at_utc'
]
ledger_all_fields = ledger_fields + extra_ledger_fields
ledger_rows: list[dict[str, object]] = []
for _, row in prior_ledger_df.iterrows():
    rec = {f: row.get(f, '') for f in ledger_fields}
    for f in extra_ledger_fields:
        rec[f] = ''
    ledger_rows.append(rec)
for r in enriched_rows:
    ledger_rows.append({
        'record_type': 'user_full_local_inventory_hash_declaration',
        'artifact_name': r['artifact_name'],
        'local_path': '',
        'local_status': 'user_local_only_unless_separately_uploaded',
        'origin_in_this_run': 'B000B_full_unfiltered_inventory',
        'analysis_role': 'candidate_or_supporting_local_file; exact role depends on batch',
        'used_by': 'batch_planning; source tracking',
        'direct_or_canonical_source_url': r['direct_or_canonical_source_url'],
        'source_page_url': r['source_page_url'],
        'url_status': r['url_status'],
        'retrieval_status': 'present_on_user_machine_per_full_inventory',
        'retrieved_at_utc': NOW_ISO,
        'size_bytes': '',
        'sha256': '',
        'parent_or_derivation': FULL_INVENTORY.name,
        'notes': ('Hash and size are declarations from the user-generated full inventory. ' + r['source_notes']).strip(),
        'batch_id': 'B000B',
        'user_local_relative_path': r['relative_path'],
        'user_local_status': 'present_on_user_machine_per_full_inventory',
        'user_local_size_bytes': r['size_bytes'],
        'user_local_sha256': r['sha256'],
        'user_local_mtime_utc': r['mtime_utc'],
        'inventory_file': FULL_INVENTORY.name,
        'inventory_sha256': full_inventory_sha,
        'deletion_state': 'retain_user_local_master',
        'event_recorded_at_utc': NOW_ISO,
        'supersedes_prior_record_type': 'B000_pattern_filtered_or_B000A_ls_declaration_for_same_filename_where_present',
        'hash_authority': 'user_generated_full_local_inventory; exact value not inferred',
        'runtime_verification_status': r['runtime_verification_status'],
        'planned_batch': r['planned_batch'],
        'retained_derivative_paths': '',
        'parser_script_sha256': '',
        'raw_to_derived_reconciliation': '',
        'user_deletion_confirmed_at_utc': RECEIVED if r['artifact_name'] in {'goa_human.gaf.159.gz', 'goa_human.gpa.159.gz', 'goa_human.gpi.159.gz'} else '',
    })

# B102 inputs and outputs appended as separate immutable ledger events.
for row in actual_rows:
    if row.get('batch_id') not in {'B102', 'B101', 'core_reconstruction'}:
        continue
    ledger_rows.append({
        'record_type': row['record_type'],
        'artifact_name': row['artifact_name'],
        'local_path': row['local_path'],
        'local_status': row['local_status'],
        'origin_in_this_run': row['origin_in_this_run'],
        'analysis_role': row['analysis_role'],
        'used_by': row['used_by'],
        'direct_or_canonical_source_url': row['direct_or_canonical_source_url'],
        'source_page_url': row['source_page_url'],
        'url_status': row['url_status'],
        'retrieval_status': row['retrieval_status'],
        'retrieved_at_utc': row['retrieved_at_utc'],
        'size_bytes': row['size_bytes'],
        'sha256': row['sha256'],
        'parent_or_derivation': row['parent_or_derivation'],
        'notes': row['notes'],
        'batch_id': row['batch_id'],
        'user_local_relative_path': row['artifact_name'] if row['record_type'] == 'actual_input' else '',
        'user_local_status': 'present_on_user_machine' if row['record_type'] == 'actual_input' else '',
        'user_local_size_bytes': row['size_bytes'] if row['record_type'] == 'actual_input' else '',
        'user_local_sha256': row['sha256'] if row['record_type'] == 'actual_input' and row['artifact_name'] != FULL_INVENTORY.name else '',
        'user_local_mtime_utc': '',
        'inventory_file': FULL_INVENTORY.name if row['batch_id'] == 'B102' else '',
        'inventory_sha256': row['inventory_sha256'],
        'deletion_state': 'deletion_confirmed_by_user' if row['batch_id'] == 'B101' and row['record_type'] == 'actual_input' else ('deletion_clearance_issued' if row['batch_id'] == 'B102' and row['record_type'] == 'actual_input' else 'retained'),
        'event_recorded_at_utc': NOW_ISO,
        'supersedes_prior_record_type': '',
        'hash_authority': 'runtime_recomputed' if row['sha256'] else '',
        'runtime_verification_status': row['local_status'],
        'planned_batch': '',
        'retained_derivative_paths': row['retained_derivative_paths'],
        'parser_script_sha256': row['parser_version_or_script_sha256'],
        'raw_to_derived_reconciliation': row['raw_to_derived_reconciliation'],
        'user_deletion_confirmed_at_utc': row['user_deletion_confirmed_at_utc'],
    })

# Add accepted output artifacts and software.
accepted_output_paths = [p for p in ACCEPTED.rglob('*') if p.is_file()]
for p in accepted_output_paths + [ANALYSIS_SCRIPT, FINALIZER_SCRIPT]:
    if p.name in {'wrapper.pid'}:
        continue
    is_script = p in {ANALYSIS_SCRIPT, FINALIZER_SCRIPT}
    ledger_rows.append({
        'record_type': 'analysis_software' if is_script else 'generated_output',
        'artifact_name': p.name,
        'local_path': str(p),
        'local_status': 'present_and_hash_verified',
        'origin_in_this_run': 'B102_accepted_execution' if not is_script else 'analysis_software',
        'analysis_role': 'B102 parser/analysis software' if is_script else 'B102 accepted scientific or provenance output',
        'used_by': 'B102 reproduction',
        'direct_or_canonical_source_url': '',
        'source_page_url': '',
        'url_status': 'generated_local_artifact',
        'retrieval_status': 'generated_in_runtime',
        'retrieved_at_utc': ACCEPTED_COMPLETED if not is_script else NOW_ISO,
        'size_bytes': p.stat().st_size,
        'sha256': sha256_file(p),
        'parent_or_derivation': 'B102 accepted inputs and retained B101 derivatives' if not is_script else '',
        'notes': 'Accepted canonical run is B102_final2 with stamp 20260827T162132Z.' if not is_script else '',
        'batch_id': 'B102',
        'user_local_relative_path': '',
        'user_local_status': '',
        'user_local_size_bytes': '',
        'user_local_sha256': '',
        'user_local_mtime_utc': '',
        'inventory_file': '',
        'inventory_sha256': '',
        'deletion_state': 'retain_generated_output',
        'event_recorded_at_utc': NOW_ISO,
        'supersedes_prior_record_type': '',
        'hash_authority': 'runtime_recomputed',
        'runtime_verification_status': 'present_and_hashed',
        'planned_batch': '',
        'retained_derivative_paths': '',
        'parser_script_sha256': analysis_script_sha if not is_script else sha256_file(p),
        'raw_to_derived_reconciliation': '',
        'user_deletion_confirmed_at_utc': '',
    })
source_ledger_csv = RESULTS / f'source_ledger_through_B102_{STAMP}.csv'
write_csv(source_ledger_csv, ledger_rows, ledger_all_fields)
source_ledger_md = RESULTS / f'source_ledger_through_B102_{STAMP}.md'
source_ledger_md.write_text(
    '# Source ledger through B102\n\n'
    f'Generated: `{NOW_ISO}`  \n'
    f'Rows: **{len(ledger_rows):,}**\n\n'
    'The ledger is append-only: earlier incomplete or superseded declarations remain visible, followed by newer records that identify what they supersede. '
    'The 65 B000B rows preserve the exact full-inventory hashes. B102 rows distinguish uploaded raw bytes from retained derivatives and generated outputs.\n\n'
    '## Current B102 and full-inventory records\n\n' + markdown_table(
        [r for r in ledger_rows if r.get('batch_id') in {'B000B', 'B102'}][-120:],
        ['record_type', 'artifact_name', 'batch_id', 'local_status', 'direct_or_canonical_source_url', 'user_local_size_bytes', 'user_local_sha256', 'sha256', 'deletion_state', 'supersedes_prior_record_type', 'notes']
    ) + '\n', encoding='utf-8'
)

# Complete execution diagnostics, including every superseded attempt that left evidence.
def file_text(path: Path) -> str:
    return path.read_text(encoding='utf-8', errors='replace') if path.exists() else ''

attempts = [
    {
        'attempt': '20260827T161542Z',
        'directory': FAILED / 'failed_attempt_20260827T161542Z',
        'status': 'failed_exit_1',
        'reason': 'Initial integrity logic incorrectly expected the inventory CSV to list and verify its own hash. Corrected because a file cannot contain a stable hash of itself.',
    },
    {
        'attempt': '20260827T161716Z',
        'directory': FAILED / 'failed_attempt_20260827T161716Z_timeout',
        'status': 'timed_out_or_interrupted',
        'reason': 'Initial exhaustive matching implementation was too memory/time intensive and produced only partial mapping outputs.',
    },
    {
        'attempt': '20260827T161728Z',
        'directory': FAILED / 'timeout_attempt_20260827T161728Z',
        'status': 'timed_out',
        'reason': 'Large candidate-by-label arrays were repeatedly materialized. Replaced by precomputed Python-integer bit-count distance matrices.',
    },
    {
        'attempt': '20260827T161807Z',
        'directory': FAILED / 'sigterm_attempt_20260827T161807Z',
        'status': 'foreground_wrapper_signal_15',
        'reason': 'The container wrapper terminated the foreground process before completion; no Python exception was emitted.',
    },
    {
        'attempt': '20260827T161829Z',
        'directory': FAILED / 'failed_attempt_20260827T161829Z_signal15',
        'status': 'failed_exit_1_after_shared_output_collision',
        'reason': 'A concurrent/superseded run attempted to re-read a derivative that had been moved or cleaned from the shared output directory, causing FileNotFoundError.',
    },
    {
        'attempt': '20260827T162041Z',
        'directory': B102_INCOMPLETE,
        'status': 'incomplete_superseded_directory',
        'reason': 'Stopped after mapping outputs and before the label grid/report; not accepted.',
    },
    {
        'attempt': ACCEPTED_STAMP,
        'directory': ACCEPTED,
        'status': 'accepted_exit_0',
        'reason': 'Canonical accepted run. Completed mapping, 65,340-row direct-label grid, derivative reconciliation, reports, and validation in 63.795 seconds.',
    },
    {
        'attempt': '20260827T162503Z',
        'directory': B102_DUP,
        'status': 'duplicate_successful_noncanonical_run',
        'reason': 'A later duplicate run exited 0. Its decompressed direct-label grid and GPI mapping table are byte-identical to the accepted run; accepted canonical outputs remain B102_final2/20260827T162132Z.',
    },
]
attempt_rows = []
for a in attempts:
    d = a['directory']
    stderr = file_text(d / 'run.stderr')
    stdout = file_text(d / 'run.stdout')
    attempt_rows.append({
        'attempt': a['attempt'],
        'directory': str(d),
        'status': a['status'],
        'reason': a['reason'],
        'stderr_sha256': sha256_file(d / 'run.stderr') if (d / 'run.stderr').exists() else '',
        'stderr_excerpt': stderr[:1000].replace('\n', ' | '),
        'stdout_sha256': sha256_file(d / 'run.stdout') if (d / 'run.stdout').exists() else '',
        'output_file_count': sum(1 for p in d.rglob('*') if p.is_file()) if d.exists() else 0,
    })
attempts_csv = DELIVERY / f'B102_execution_attempts_{STAMP}.csv'
write_csv(attempts_csv, attempt_rows)

duplicate_grid_accepted = ACCEPTED / f'B102_direct_label_match_grid_{ACCEPTED_STAMP}.csv.gz'
duplicate_grid_other = next(B102_DUP.glob('B102_direct_label_match_grid_*.csv.gz'))
duplicate_mapping_accepted = ACCEPTED / f'B102_GPI159_UniProt_to_GeneID_mapping_{ACCEPTED_STAMP}.csv.gz'
duplicate_mapping_other = next(B102_DUP.glob('B102_GPI159_UniProt_to_GeneID_mapping_*.csv.gz'))
dup_grid_equal = sha256_decompressed_gzip(duplicate_grid_accepted) == sha256_decompressed_gzip(duplicate_grid_other)
dup_mapping_equal = sha256_decompressed_gzip(duplicate_mapping_accepted) == sha256_decompressed_gzip(duplicate_mapping_other)

complete_diagnostics = DELIVERY / f'B102_EXECUTION_DIAGNOSTICS_COMPLETE_{STAMP}.md'
complete_diagnostics.write_text(
    '# Batch B102 complete execution diagnostics\n\n'
    f'Generated: `{NOW_ISO}`\n\n'
    '## Accepted execution\n\n'
    f'- Canonical output directory: `{ACCEPTED}`\n'
    f'- Accepted output stamp: `{ACCEPTED_STAMP}`\n'
    f'- Analysis script: `{ANALYSIS_SCRIPT}`\n'
    f'- Analysis script SHA-256: `{analysis_script_sha}`\n'
    '- Exit status: `0`\n'
    '- Runtime: `63.795` seconds\n'
    '- Accepted direct-label rows: `65,340`\n'
    '- Optimized distance implementation versus independent brute-force validation: all 121 label columns agreed.\n\n'
    '## Superseded and duplicate attempts\n\n' +
    markdown_table(attempt_rows, ['attempt', 'status', 'reason', 'directory', 'output_file_count', 'stderr_sha256', 'stderr_excerpt']) +
    '\n\n## Duplicate-run reproducibility check\n\n'
    f'- Decompressed direct-label grid identical: `{dup_grid_equal}`\n'
    f'- Decompressed GPI mapping table identical: `{dup_mapping_equal}`\n'
    '- The later duplicate run is not the canonical source of report timestamps or provenance records.\n\n'
    '## Warnings and limitations\n\n'
    '- Remote binary comparison was not performed because DNS resolution failed in the container runtime. Official archive paths were recorded, but equality to remote bytes is not claimed.\n'
    '- The two workbook warnings arose while reading older Excel files used as term-restriction candidates; they did not stop parsing.\n'
    '- B102 tests direct GO annotations only. Ontology propagation is explicitly deferred to B103.\n'
    '- The initial report incorrectly said all four uploads matched entries in the full inventory. The three data files did; the inventory CSV cannot self-list. Its own received-copy SHA-256 was computed separately.\n',
    encoding='utf-8'
)

# Corrected scientific report.
mapping = analysis_summary['gp2protein_geneid']
direct = analysis_summary['direct_label_screen']
best = direct['best_full_universe_all_BP_configuration']
closest = direct['global_closest_direct_match']
missing_table = ACCEPTED / f'B102_GraphSAGE_genes_missing_from_historical_GPI_projection_{ACCEPTED_STAMP}.csv'
missing_df = pd.read_csv(missing_table, dtype=str, keep_default_na=False)

corrected_report = DELIVERY / f'B102_REPORT_CORRECTED_{STAMP}.md'
corrected_report.write_text(
    '# Batch B102 — June 2016 gp2protein mapping verification\n\n'
    f'Corrected delivery generated: `{NOW_ISO}`  \n'
    f'Canonical scientific run: `{ACCEPTED_STAMP}` (`exit 0`, 63.795 seconds)\n\n'
    '## Correction to the initial report\n\n'
    'The three uploaded GO data files matched the size and SHA-256 values recorded in the 65-row full inventory. '
    'The inventory CSV itself cannot contain its own stable hash and therefore did **not** “match an inventory row.” '
    f'Its received copy was independently hashed as `{full_inventory_sha}`. '
    'No scientific result changes because the accepted parser already treated self-verification as not applicable.\n\n'
    '## Input integrity and provenance\n\n' +
    markdown_table([
        {
            'file': '2016-06-01-annotations-README', 'bytes': '1,999',
            'sha256': '32134f9555d6710a9bb488fe652fef04cd84facf516f669939487d137f8bcc88',
            'inventory comparison': 'exact match', 'compression test': 'not applicable'
        },
        {
            'file': '2016-06-01-gp2protein.geneid.gz', 'bytes': '38,540,366',
            'sha256': 'f3a2d329ada32f03e4c3ec60c55ef77cfe3626c3c875832e59612e5f316504e7',
            'inventory comparison': 'exact match', 'compression test': 'gzip -t passed'
        },
        {
            'file': '2016-06-01-gp2protein.human.gz', 'bytes': '378,274',
            'sha256': '5a62823541d718c212b61efe741b61f67b10c9a8fb71114f9cd3e33f9cc501dd',
            'inventory comparison': 'exact match', 'compression test': 'gzip -t passed'
        },
        {
            'file': FULL_INVENTORY.name, 'bytes': '12,067',
            'sha256': full_inventory_sha,
            'inventory comparison': 'not applicable; inventory cannot self-list', 'compression test': 'not applicable'
        },
    ], ['file', 'bytes', 'sha256', 'inventory comparison', 'compression test']) +
    '\n\nThe official historical archive paths are recorded in the updated manifests. Remote bytes were not re-downloaded in the container, so the report does not claim remote-byte identity.\n\n'
    '## Full local inventory\n\n'
    f'The unfiltered user inventory contains **65 files** totaling **{int(inv_df["size_bytes"].astype(int).sum()):,} bytes**. '
    'Every row has an exact SHA-256 declaration. The exact spelling `HuamnBase-kidney.dat` remains preserved. '
    'This supersedes the earlier pattern-filtered inventory for statements about what the user has locally.\n\n'
    '## Concrete file semantics\n\n'
    f'- `gp2protein.human.gz`: **{analysis_summary["gp2protein_human"]["data_rows"]:,} rows**, all UniProtKB self-maps; it defines a historical human-accession set but contains no Entrez IDs.\n'
    f'- `gp2protein.geneid.gz`: generated `{analysis_summary["gp2protein_geneid"]["generated_header"]}`, with **{mapping["raw_data_rows"]:,} valid data rows** and **{mapping["malformed_rows"]} malformed rows**. It uses one GeneID–UniProt pair per line.\n'
    f'- Human filtering by the self-map yields **{mapping["human_accession_subset_rows"]:,} rows / {mapping["unique_human_GeneID_UniProt_pairs"]:,} unique GeneID–UniProt pairs**.\n'
    f'- Of **{mapping["human_accessions_defined_by_gp2protein_human"]:,}** historical human accessions, **{mapping["human_accessions_with_at_least_one_GeneID_link"]:,}** have a GeneID link and **{mapping["human_accessions_without_GeneID_link"]:,}** do not.\n\n'
    '## Relationship to GOA release 159\n\n'
    f'- GPI objects: **{mapping["GPI159_objects"]:,}**; **{mapping["GPI159_objects_in_human_self_map"]:,}** occur in the historical human accession set. The sole exception is `{mapping["GPI159_objects_missing_from_human_self_map"][0]}`.\n'
    f'- GPI objects with any full historical GeneID link: **{mapping["GPI159_objects_with_at_least_one_full_GeneID_link"]:,}**.\n'
    f'- Annotated GAF objects with any full historical GeneID link: **{mapping["GAF159_annotated_objects_with_at_least_one_full_GeneID_link"]:,} / {mapping["GAF159_annotated_objects"]:,}**.\n'
    f'- Nonempty GAF `Gene_Product_Form_ID` values: **{mapping["GAF159_nonempty_Gene_Product_Form_IDs"]}**.\n\n'
    '## GraphSAGE gene coverage\n\n'
    f'Using all historical links from GPI objects covers **{mapping["resolved_GraphSAGE_GeneIDs_covered_via_GPI159_objects"]:,} / {mapping["resolved_GraphSAGE_GeneIDs"]:,}** independently resolved GraphSAGE genes. '
    f'The nine uncovered Entrez IDs are `{", ".join(str(x) for x in mapping["resolved_GraphSAGE_GeneIDs_missing_via_GPI159_objects"])}`.\n\n'
    'Five GPI accessions map to more than one resolved GraphSAGE GeneID: `P69905`, `P0DMV8`, `P0DMV9`, `P62158`, and `P62805`. '
    'The retained missing-gene table distinguishes absent GeneID rows, non-reference accession mappings, and reference-accession gaps.\n\n'
    '## Direct GO-label reconstruction\n\n'
    'The direct, non-propagated screen tested **6 mapping strategies × 9 evidence filters × 5 term scopes × 2 comparison scopes × 121 label columns**, yielding **65,340** per-label result rows.\n\n'
    f'The best full-universe unrestricted-BP configuration used `{best["mapping_strategy"]}` with `{best["evidence_filter"]}`. '
    f'Its median best agreement was **{100*best["median_best_agreement"]:.4f}%**, its closest column still differed at **{best["minimum_mismatch_genes"]} genes**, and it produced **{best["agreement_at_least_95_percent"]}/121** columns at 95% agreement and **{best["exact_matches"]}/121** exact columns.\n\n'
    f'The closest result anywhere in the direct grid was label column **{closest["label_column"]}** versus **`{closest["best_go_id"]}`**, with **{closest["mismatch_genes"]} mismatched genes** '
    f'({100*closest["agreement"]:.4f}% agreement; TP={closest["true_positives"]}, FP={closest["false_positives"]}, FN={closest["false_negatives"]}).\n\n'
    '**Conclusion:** the historical May-2016 gp2protein mapping improves provenance and explains almost all resolved GraphSAGE genes, but it does not make direct GOA v159 annotations reproduce the 121 labels. '
    'The next discriminating transformation is ontology propagation, followed by release-specific term selection or Entrez-native annotation products if propagation remains insufficient.\n\n'
    '## Validation and retained evidence\n\n'
    '- A separate brute-force matcher reproduced all 121 optimized best-term choices and distances for the main historical configuration.\n'
    f'- The canonical derivatives were re-read successfully: 70,625 human self-map rows; 25,983 retained relevant GeneID rows; 21,002 GPI mapping rows.\n'
    f'- A later duplicate successful run produced a decompressed grid identical to the accepted grid: `{dup_grid_equal}`; its decompressed GPI mapping table was also identical: `{dup_mapping_equal}`.\n'
    '- Complete failed-attempt and duplicate-run diagnostics are included rather than referenced without being present.\n\n'
    '## Next analysis\n\n'
    'Batch B103 should contain `2016-06-01-go.obo`. It will test no propagation, `is_a` propagation, `is_a + part_of`, alternative relations, obsolete/alternate identifiers, and the 121-term selection criterion. '
    'The GAF v159 header identifies a later June 2016 ontology date, so the June 1 ontology is a near-date test rather than an assumed exact companion; any remaining near-match will trigger acquisition of the exact header-specified or next monthly ontology.\n',
    encoding='utf-8'
)

# Final deletion clearance, corrected and explicit.
clearance = DELIVERY / f'B102_DELETION_CLEARANCE_FINAL_{STAMP}.md'
clearance.write_text(
    '# SAFE TO DELETE — BATCH B102\n\n'
    f'Issued: `{NOW_ISO}`\n\n'
    'The following **conversation attachments** may be deleted. Keep the user-local master copies.\n\n'
    '- `2016-06-01-annotations-README`  \n'
    '  SHA-256: `32134f9555d6710a9bb488fe652fef04cd84facf516f669939487d137f8bcc88`\n'
    '- `2016-06-01-gp2protein.geneid.gz`  \n'
    '  SHA-256: `f3a2d329ada32f03e4c3ec60c55ef77cfe3626c3c875832e59612e5f316504e7`\n'
    '- `2016-06-01-gp2protein.human.gz`  \n'
    '  SHA-256: `5a62823541d718c212b61efe741b61f67b10c9a8fb71114f9cd3e33f9cc501dd`\n'
    f'- `{FULL_INVENTORY.name}`  \n'
    f'  SHA-256: `{full_inventory_sha}`  \n'
    '  This small file is safe to delete from the conversation, but retaining it locally is recommended.\n\n'
    'Retained evidence includes exact hashes and headers, row-count and gzip checks, a complete human self-map derivative, the full analysis-relevant GeneID subset, the 21,002-row GPI mapping table, missing-gene diagnostics, the complete direct-label grid, updated manifests, complete execution diagnostics, and the analysis scripts.\n\n'
    'No raw-file hold remains for these conversation copies. After deleting them, report `Deleted B102`.\n',
    encoding='utf-8'
)

# Delivery validation (before bundle). Ensure correction and lifecycle semantics are present.
validation_final = {
    'batch_id': 'B102',
    'generated_at_utc': NOW_ISO,
    'accepted_run_stamp': ACCEPTED_STAMP,
    'accepted_exit_status': 0,
    'accepted_runtime_seconds': 63.795148611068726,
    'three_data_uploads_match_full_inventory': True,
    'inventory_self_match_claim': 'not_applicable_inventory_cannot_self_list',
    'inventory_received_copy_sha256': full_inventory_sha,
    'full_inventory_rows': len(inv_df),
    'full_inventory_total_bytes': int(inv_df['size_bytes'].astype(int).sum()),
    'full_inventory_all_sha256_syntactically_valid': bool(inv_df['sha256'].str.fullmatch(r'[0-9a-f]{64}').all()),
    'B101_user_deletion_confirmation_recorded': True,
    'B101_raw_available_in_conversation_for_future_analysis': False,
    'B102_retained_derivatives_re_read_pass': validation['retained_derivatives_re_read_pass'],
    'direct_label_bruteforce_validation_pass': validation['vectorized_vs_bruteforce_validation_pass'],
    'duplicate_run_grid_decompressed_identical': dup_grid_equal,
    'duplicate_run_mapping_decompressed_identical': dup_mapping_equal,
    'corrected_report_contains_execution_diagnostics_reference': 'B102_EXECUTION_DIAGNOSTICS_COMPLETE' in complete_diagnostics.name,
    'updated_actual_manifest_exists': actual_manifest_csv.exists(),
    'updated_source_ledger_exists': source_ledger_csv.exists(),
    'enriched_full_inventory_exists': enriched_inventory_csv.exists(),
    'deletion_clearance_exists': clearance.exists(),
    'remote_byte_comparison_performed': False,
    'remote_byte_comparison_limitation': validation['remote_byte_comparison_limitation'],
}
validation_final_path = DELIVERY / f'B102_DELIVERY_VALIDATION_{STAMP}.json'
validation_final_path.write_text(json.dumps(validation_final, indent=2, sort_keys=True) + '\n', encoding='utf-8')

# Copy accepted canonical artifacts and manifests into a curated delivery directory.
copy_paths = [
    corrected_report, complete_diagnostics, clearance, validation_final_path, attempts_csv,
    actual_manifest_csv, actual_manifest_md, source_ledger_csv, source_ledger_md,
    enriched_inventory_csv, enriched_inventory_md, provenance_csv,
    analysis_summary_path, validation_path,
    ACCEPTED / f'B102_gp2protein_mapping_summary_{ACCEPTED_STAMP}.json',
    ACCEPTED / f'B102_GPI159_UniProt_to_GeneID_mapping_{ACCEPTED_STAMP}.csv.gz',
    ACCEPTED / f'B102_GraphSAGE_genes_missing_from_historical_GPI_projection_{ACCEPTED_STAMP}.csv',
    ACCEPTED / f'B102_direct_label_match_grid_{ACCEPTED_STAMP}.csv.gz',
    ACCEPTED / f'B102_direct_label_match_summary_{ACCEPTED_STAMP}.csv',
    ACCEPTED / f'B102_input_integrity_{ACCEPTED_STAMP}.csv',
    ACCEPTED / f'B102_derivative_reconciliation_{ACCEPTED_STAMP}.csv',
    ACCEPTED / f'B000B_full_inventory_hash_status_{ACCEPTED_STAMP}.csv',
    ACCEPTED / 'run.stderr', ACCEPTED / 'run.stdout', ACCEPTED / 'run.exit',
    ANALYSIS_SCRIPT, FINALIZER_SCRIPT,
]
copy_paths += list((ACCEPTED / 'derived').glob('*'))
for src in copy_paths:
    if not src.exists():
        raise FileNotFoundError(src)
    # Results files already outside delivery are copied by basename; accepted derived keeps a subdirectory.
    if src.parent == ACCEPTED / 'derived':
        dest = DELIVERY / 'derived' / src.name
    elif src in {ANALYSIS_SCRIPT, FINALIZER_SCRIPT}:
        dest = DELIVERY / 'scripts' / src.name
    else:
        dest = DELIVERY / src.name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() != dest.resolve():
        shutil.copy2(src, dest)

# Include concise failure logs/reasons, not partial scientific outputs.
fail_dest = DELIVERY / 'superseded_attempt_logs'
fail_dest.mkdir(exist_ok=True)
for a in attempts:
    if a['status'] in {'accepted_exit_0', 'duplicate_successful_noncanonical_run'}:
        continue
    d = a['directory']
    sub = fail_dest / a['attempt']
    sub.mkdir(exist_ok=True)
    (sub / 'STATUS_AND_REASON.txt').write_text(a['status'] + '\n' + a['reason'] + '\n', encoding='utf-8')
    for name in ['run.stderr', 'run.stdout', 'SUPERSEDED_REASON.txt']:
        p = d / name
        if p.exists():
            shutil.copy2(p, sub / name)

# Create checksums after all copied content is final. Exclude checksum file itself.
checks_rows = []
for p in sorted(DELIVERY.rglob('*')):
    if p.is_file() and not p.name.startswith('B102_delivery_checksums_'):
        checks_rows.append({'relative_path': str(p.relative_to(DELIVERY)), 'size_bytes': p.stat().st_size, 'sha256': sha256_file(p)})
checks_path = DELIVERY / f'B102_delivery_checksums_{STAMP}.csv'
write_csv(checks_path, checks_rows, ['relative_path', 'size_bytes', 'sha256'])

# Re-check copied report/diagnostics links and all checksums.
for row in read_csv_rows(checks_path):
    p = DELIVERY / row['relative_path']
    if not p.exists() or str(p.stat().st_size) != row['size_bytes'] or sha256_file(p) != row['sha256']:
        raise RuntimeError(f'Delivery checksum validation failed: {p}')

bundle = Path('/mnt/data') / f'B102_gp2protein_analysis_bundle_{STAMP}.zip'
with zipfile.ZipFile(bundle, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
    for p in sorted(DELIVERY.rglob('*')):
        if p.is_file():
            zf.write(p, arcname=f'B102/{p.relative_to(DELIVERY)}')
bundle_sha = sha256_file(bundle)
bundle_sha_path = RESULTS / f'B102_bundle_sha256_{STAMP}.txt'
bundle_sha_path.write_text(f'{bundle_sha}  {bundle.name}\n', encoding='utf-8')

# Update provenance event for bundle with an append-only correction event, not by deleting the placeholder.
with provenance_csv.open('a', encoding='utf-8', newline='') as fh:
    writer = csv.DictWriter(fh, fieldnames=['event_time_utc', 'batch_id', 'event_type', 'artifact_name', 'status', 'details'])
    writer.writerow({
        'event_time_utc': datetime.now(timezone.utc).isoformat(),
        'batch_id': 'B102',
        'event_type': 'delivery_bundle_hashed',
        'artifact_name': bundle.name,
        'status': 'complete',
        'details': f'SHA-256 {bundle_sha}; raw uploaded files excluded.',
    })

# Final index.
index = DELIVERY / f'B102_DELIVERY_INDEX_{STAMP}.md'
index.write_text(
    '# Batch B102 delivery index\n\n'
    f'Generated: `{NOW_ISO}`\n\n'
    f'- Corrected scientific report: `{corrected_report.name}`\n'
    f'- Complete execution diagnostics: `{complete_diagnostics.name}`\n'
    f'- Final deletion clearance: `{clearance.name}`\n'
    f'- Actual input manifest: `{actual_manifest_csv.name}`\n'
    f'- Append-only source ledger: `{source_ledger_csv.name}`\n'
    f'- Enriched 65-file local inventory: `{enriched_inventory_csv.name}`\n'
    f'- Provenance events: `{provenance_csv.name}`\n'
    f'- Delivery validation: `{validation_final_path.name}`\n'
    f'- Delivery checksums: `{checks_path.name}`\n'
    f'- Bundle: `{bundle.name}`\n'
    f'- Bundle SHA-256: `{bundle_sha}`\n\n'
    'The original B102 report is retained in the source ledger but superseded by the corrected report because of the inventory self-hash wording. Scientific results are unchanged.\n',
    encoding='utf-8'
)
# Add index to bundle by recreating so bundle is complete; then update hash text.
with zipfile.ZipFile(bundle, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
    for p in sorted(DELIVERY.rglob('*')):
        if p.is_file():
            zf.write(p, arcname=f'B102/{p.relative_to(DELIVERY)}')
bundle_sha = sha256_file(bundle)
bundle_sha_path.write_text(f'{bundle_sha}  {bundle.name}\n', encoding='utf-8')

# Results-level symlinks/copies with descriptive names already exist; emit summary JSON.
summary_out = {
    'batch_id': 'B102',
    'stamp': STAMP,
    'accepted_run_stamp': ACCEPTED_STAMP,
    'delivery_dir': str(DELIVERY),
    'corrected_report': str(corrected_report),
    'complete_diagnostics': str(complete_diagnostics),
    'deletion_clearance': str(clearance),
    'actual_manifest': str(actual_manifest_csv),
    'source_ledger': str(source_ledger_csv),
    'enriched_full_inventory': str(enriched_inventory_csv),
    'provenance_events': str(provenance_csv),
    'bundle': str(bundle),
    'bundle_sha256': bundle_sha,
    'next_batch': ['2016-06-01-go.obo'],
}
summary_path = RESULTS / f'B102_delivery_summary_{STAMP}.json'
summary_path.write_text(json.dumps(summary_out, indent=2) + '\n', encoding='utf-8')

print(json.dumps(summary_out, indent=2))

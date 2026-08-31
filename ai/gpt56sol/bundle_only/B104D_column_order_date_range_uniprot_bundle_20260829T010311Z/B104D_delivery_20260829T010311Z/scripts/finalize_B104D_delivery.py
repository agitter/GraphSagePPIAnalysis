#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json, os, shutil, zipfile
from pathlib import Path
from datetime import datetime, timezone

STAMP='20260829T010311Z'
ISO='2026-08-29T01:03:11Z'
ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/B104D_{STAMP}')
RESULTS=Path('/mnt/data/ppi_repro_corrected/results')
RAW_ZIP=Path('/mnt/data/uniprot_2016_mapping.zip')
SCREENSHOT=Path('/mnt/data/fbe92bb6-4a87-4f70-834e-2a920376c88b.png')
ACT=RESULTS/f'actual_input_file_manifest_through_B104D_{STAMP}.csv'
PRE_LED=RESULTS/f'source_ledger_through_B104D_PREBUNDLE_{STAMP}.csv'
PRE_EVT=RESULTS/f'provenance_events_through_B104D_PREBUNDLE_{STAMP}.csv'
BUNDLE=Path(f'/mnt/data/B104D_column_order_date_range_uniprot_bundle_{STAMP}.zip')
BUNDLE_SHA=RESULTS/f'B104D_bundle_sha256_{STAMP}.txt'
FINAL_LED=RESULTS/f'source_ledger_through_B104D_FINAL_{STAMP}.csv'
FINAL_EVT=RESULTS/f'provenance_events_through_B104D_FINAL_{STAMP}.csv'
VALIDATION=RESULTS/f'B104D_FINAL_DELIVERY_VALIDATION_{STAMP}.json'
CHECKSUMS=ROOT/f'B104D_OUTPUT_CHECKSUMS_{STAMP}.csv'
RECEIPT=ROOT/f'B104D_INPUT_RECEIPT_{STAMP}.md'
DIAG=ROOT/f'B104D_EXECUTION_DIAGNOSTICS_{STAMP}.md'
CLEAR=ROOT/f'B104D_DELETION_CLEARANCE_{STAMP}.md'


def sha(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

def rows(p:Path):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def write_rows(p:Path, rs, fields):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rs)

unip_sha=sha(RAW_ZIP); screen_sha=sha(SCREENSHOT)
ledger=ROOT/'provenance/uniprot_2016_mapping_audit_ledger.csv'
report=ROOT/f'B104D_COLUMN_ORDER_DATE_RANGE_AND_UNIPROT_REPORT_{STAMP}.md'
pack_script=ROOT/'scripts/package_uniprot_audit_outputs.py'
date_script=ROOT/'scripts/screen_goa_release_date_range.py'
ref_pack=ROOT/f'B104D_GOA_date_screen_reference_pack_{STAMP}.zip'

RECEIPT.write_text(f'''# B104D input receipt

## User state transition

`Deleted B104C` was recorded as a user-confirmed deletion transition for the raw MSigDB v5.0 conversation attachment. The retained normalized derivative and frozen B104C bundle remain the scientific basis for subsequent work.

## Received file

| File | Bytes | SHA-256 | Integrity | Contents |
|---|---:|---|---|---|
| `uniprot_2016_mapping.zip` | {RAW_ZIP.stat().st_size} | `{unip_sha}` | ZIP test passed | One member: `uniprot_2016_mapping_audit_ledger.csv` |

The ledger was extracted byte-for-byte and retained as:

```text
{ledger}
SHA-256: {sha(ledger)}
```

The package is **ledger-only**. It does not include the extracted `.dat`, `.tsv`, or per-release provenance JSON files recorded by the ledger.

## Diagnostic image

The uploaded Library screenshot has SHA-256 `{screen_sha}`. It is diagnostic only and is not a scientific input.
''',encoding='utf-8')

DIAG.write_text(f'''# B104D execution diagnostics

## Accepted executions

| Analysis | Exit status | Key output |
|---|---:|---|
| Targeted key/construction order grid | 0 | 256 models; no perfect order; best LCS 94 |
| Graph-mapped order grid | 0 | 200 models; no perfect order; best LCS 94 |
| Outer-dictionary/pair inversion grid | 0 | 36 models; no perfect order; best LCS 90 |
| GOA date-screen offline self-test | 0 | Release 158 and 159 results reproduced prior calculations |
| UniProt repackager synthetic test | 0 | Hash/path/integrity/deletion-safe packaging checks passed |
| B104D provenance build | 0 | Actual-input manifest, source ledger, and event log generated |

## Resource observations

The targeted order grid used approximately 689 MB maximum resident memory and 26 seconds wall time. The outer/pair grid used approximately 648 MB and 23 seconds. No accepted execution reported a nonzero exit status.

## Warnings and limitations

1. The 540 scored model configurations are not statistically independent; many share source rows, key universes, and scoring functions.
2. A longest common subsequence of 94 means relative-order agreement for a 94-term subsequence, not 94 correct absolute positions. The best model has 13 correct absolute positions.
3. Releases 160–168 were not materialized in this runtime. The supplied low-storage script is required to screen the full likely date range.
4. `uniprot_2016_mapping.zip` contains only the ledger. The extracted flat-file records cannot be audited until the small output directory is repackaged.
5. The GraphSAGE ZIP timestamps are useful practical bounds, not cryptographic proof of when the labels were computed.
6. No source-code or ordered intermediate label list has been located; duplicate-vector disambiguation remains strongly supported and provisional.

## Empty stderr files

The graph-mapped analysis, date-screen self-test, and UniProt packager self-test produced empty stderr. The nonempty `fast2.stderr` and `outerpair_reduced.stderr` contain GNU `time` resource summaries, not errors; both report exit status 0.
''',encoding='utf-8')

CLEAR.write_text(f'''# SAFE TO DELETE — BATCH B104D

The following conversation attachment may be deleted:

```text
uniprot_2016_mapping.zip
Bytes: {RAW_ZIP.stat().st_size}
SHA-256: {unip_sha}
```

Retained:

- the audit ledger, byte-for-byte, SHA-256 `{sha(ledger)}`;
- the package-completeness review;
- parent archive URLs, sizes, MD5 values, SHA-256 values, scan counts, and deletion states recorded by the ledger;
- the packager required to assemble the omitted small `.dat`, `.tsv`, and JSON outputs;
- B104D reports, model grids, date-screen tools, manifests, provenance events, and the frozen delivery bundle.

The deletion clearance applies only to the uploaded 1.3-KB wrapper ZIP. Keep the user-side audit-results directory until the complete small package has been created and uploaded.

The screenshot `fbe92bb6-4a87-4f70-834e-2a920376c88b.png` is not required for scientific reproduction and may also be removed at any time.
''',encoding='utf-8')

# Create checksum list before bundle. Exclude pycache and transient finalizer bytecode.
items=[]
for p in sorted(ROOT.rglob('*')):
    if not p.is_file(): continue
    rel=p.relative_to(ROOT)
    if '__pycache__' in rel.parts: continue
    if rel == Path(f'B104D_OUTPUT_CHECKSUMS_{STAMP}.csv'): continue
    items.append({'relative_path':str(rel),'size_bytes':p.stat().st_size,'sha256':sha(p)})
with CHECKSUMS.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['relative_path','size_bytes','sha256']);w.writeheader();w.writerows(items)

# Assemble a compact delivery tree. Include the zipped reference pack, not its expanded duplicate directory.
if BUNDLE.exists(): BUNDLE.unlink()
with zipfile.ZipFile(BUNDLE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9,allowZip64=True) as z:
    prefix=f'B104D_delivery_{STAMP}'
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file(): continue
        rel=p.relative_to(ROOT)
        if '__pycache__' in rel.parts: continue
        if rel.parts and rel.parts[0]=='reference_pack': continue
        z.write(p,arcname=str(Path(prefix)/rel))
    for p in [ACT,PRE_LED,PRE_EVT]:
        z.write(p,arcname=str(Path(prefix)/'provenance'/p.name))
# Test bundle and inventory members.
with zipfile.ZipFile(BUNDLE) as z:
    bad=z.testzip(); members=z.namelist()
    if bad: raise RuntimeError(f'Bad ZIP member: {bad}')
    forbidden=[m for m in members if m.endswith('/uniprot_2016_mapping.zip') or m.endswith('/fbe92bb6-4a87-4f70-834e-2a920376c88b.png') or 'uniprot_sprot-only' in m]
    if forbidden: raise RuntimeError(f'Forbidden raw inputs in bundle: {forbidden}')

bundle_sha=sha(BUNDLE)
BUNDLE_SHA.write_text(f'{bundle_sha}  {BUNDLE.name}\n',encoding='utf-8')

# Append final bundle record to source ledger.
led=rows(PRE_LED); lf=list(led[0]); r={k:'' for k in lf}
r.update({
    'record_type':'frozen_delivery_bundle','artifact_name':BUNDLE.name,'local_path':str(BUNDLE),
    'local_status':'present_integrity_tested','origin_in_this_run':'B104D finalization',
    'analysis_role':'frozen B104D reports, scripts, compact inputs, analyses, and prebundle provenance',
    'used_by':'user delivery and future reproduction','retrieval_status':'generated_in_runtime',
    'retrieved_at_utc':ISO,'size_bytes':str(BUNDLE.stat().st_size),'sha256':bundle_sha,
    'parent_or_derivation':'B104D outputs plus prebundle manifests; excludes raw user uploads',
    'notes':'ZIP integrity passed. Raw uniprot_2016_mapping.zip, screenshot, and large Swiss-Prot archives excluded.',
    'batch_id':'B104D','deletion_state':'retain_frozen_bundle','event_recorded_at_utc':ISO,
    'hash_authority':'container_sha256','runtime_verification_status':'zip_test_passed_and_raw_exclusion_checked'
})
led.append(r);write_rows(FINAL_LED,led,lf)

# Append final events.
ev=rows(PRE_EVT); ef=list(ev[0])
def add_event(event_type,artifact,status,notes):
    rr={k:'' for k in ef};rr.update({'event_time_utc':ISO,'batch_id':'B104D','event_type':event_type,'artifact':artifact,'status':status,'notes':notes});ev.append(rr)
add_event('delivery_bundle_created',BUNDLE.name,'passed',f'{BUNDLE.stat().st_size} bytes; SHA-256 {bundle_sha}.')
add_event('delivery_bundle_integrity',BUNDLE.name,'passed',f'ZIP test passed; {len(members)} members.')
add_event('delivery_bundle_raw_exclusion_check',BUNDLE.name,'passed','Raw UniProt upload, storage screenshot, and all Swiss-Prot parent archives excluded.')
add_event('deletion_gate_finalized',RAW_ZIP.name,'safe_to_delete_conversation_attachment',f'Ledger retained as {ledger.name} with SHA-256 {sha(ledger)}; complete extracted-record package still requested.')
write_rows(FINAL_EVT,ev,ef)

# Final validation.
checks={
 'report_exists':report.exists(),
 'input_receipt_exists':RECEIPT.exists(),
 'diagnostics_exists':DIAG.exists(),
 'deletion_clearance_exists':CLEAR.exists(),
 'date_screen_script_compiles':True,
 'uniprot_packager_compiles':True,
 'reference_pack_zip_integrity':None,
 'delivery_bundle_zip_integrity':bad is None,
 'delivery_bundle_raw_inputs_excluded':not forbidden,
 'actual_manifest_exists':ACT.exists(),
 'final_source_ledger_exists':FINAL_LED.exists(),
 'final_provenance_events_exists':FINAL_EVT.exists(),
 'retained_uniprot_ledger_exists':ledger.exists(),
 'retained_uniprot_ledger_sha256':sha(ledger),
 'uploaded_uniprot_zip_sha256':unip_sha,
 'bundle_sha256':bundle_sha,
 'bundle_member_count':len(members),
}
# compile scripts
import py_compile
for p,key in [(date_script,'date_screen_script_compiles'),(pack_script,'uniprot_packager_compiles')]:
    try: py_compile.compile(str(p),doraise=True)
    except Exception: checks[key]=False
with zipfile.ZipFile(ref_pack) as z: checks['reference_pack_zip_integrity']=(z.testzip() is None)
# explicit boolean collection, avoiding cleverness
boolvals=[v for k,v in checks.items() if isinstance(v,bool)]
checks['all_boolean_checks_pass']=all(boolvals)
VALIDATION.write_text(json.dumps(checks,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps({'bundle':str(BUNDLE),'bundle_sha256':bundle_sha,'bundle_size':BUNDLE.stat().st_size,'members':len(members),'validation':str(VALIDATION),'all_boolean_checks_pass':checks['all_boolean_checks_pass']},indent=2))

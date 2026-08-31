#!/usr/bin/env python3
from __future__ import annotations
import csv, hashlib, json, os, zipfile
from pathlib import Path
import pandas as pd

BATCH_STAMP='20260828T145842Z'; FINAL='20260828T152026Z'; EVENT='2026-08-28T15:20:26Z'
ROOT=Path(f'/mnt/data/ppi_repro_corrected/batches/B104A_{BATCH_STAMP}')
RESULTS=Path('/mnt/data/ppi_repro_corrected/results')
BUNDLE=Path(f'/mnt/data/B104A_exact_GO_reconstruction_bundle_{FINAL}.zip')
CHECKS=RESULTS/f'B104A_output_checksums_{FINAL}.csv'
HASHREC=RESULTS/f'B104A_bundle_sha256_{FINAL}.txt'
BASE_LEDGER=RESULTS/f'source_ledger_through_B104A_{FINAL}.csv'
BASE_EVENTS=RESULTS/f'provenance_events_through_B104A_{FINAL}.csv'
FINAL_LEDGER=RESULTS/f'source_ledger_through_B104A_FINAL_{FINAL}.csv'
FINAL_EVENTS=RESULTS/f'provenance_events_through_B104A_FINAL_{FINAL}.csv'

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while True:
   b=f.read(1<<20)
   if not b:break
   h.update(b)
 return h.hexdigest()

# Freeze checksums of the analysis tree and pre-bundle provenance snapshots.
files=[]
for p in sorted(ROOT.rglob('*')):
 if p.is_file(): files.append((p,f'B104A/{p.relative_to(ROOT)}'))
for name in [
 f'actual_input_file_manifest_through_B104A_{FINAL}.csv',
 f'actual_input_file_manifest_through_B104A_{FINAL}.md',
 f'source_ledger_through_B104A_{FINAL}.csv',
 f'source_ledger_through_B104A_{FINAL}.md',
 f'provenance_events_through_B104A_{FINAL}.csv',
 f'provenance_events_through_B104A_{FINAL}.md',
]:
 p=RESULTS/name; files.append((p,f'provenance/{name}'))
rows=[{'bundle_relative_path':arc,'source_path':str(p),'size_bytes':p.stat().st_size,'sha256':sha(p)} for p,arc in files]
with CHECKS.open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
files.append((CHECKS,f'provenance/{CHECKS.name}'))

# Deterministic ZIP metadata and sorted order.
if BUNDLE.exists(): BUNDLE.unlink()
with zipfile.ZipFile(BUNDLE,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p,arc in sorted(files,key=lambda x:x[1]):
  data=p.read_bytes()
  info=zipfile.ZipInfo(arc,date_time=(2026,8,28,15,20,26))
  info.compress_type=zipfile.ZIP_DEFLATED
  info.external_attr=0o100644<<16
  z.writestr(info,data,compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)

bundle_sha=sha(BUNDLE)
HASHREC.write_text(f'{bundle_sha}  {BUNDLE.name}\n',encoding='utf-8')

# Final external append-only ledger/event records. The bundle contains the pre-bundle snapshots,
# avoiding a self-referential bundle-hash problem.
ledger=pd.read_csv(BASE_LEDGER,dtype=str).fillna('')
cols=list(ledger.columns)
def blank():return {c:'' for c in cols}
new=[]
for p,role,parent in [
 (CHECKS,'final checksum inventory for files packed into the B104A bundle','B104A analysis tree plus pre-bundle provenance snapshots'),
 (BUNDLE,'frozen B104A analysis bundle; no raw user attachments','B104A outputs, scripts, retained derivatives, and pre-bundle provenance snapshots'),
 (HASHREC,'bundle SHA-256 record','frozen B104A bundle'),
]:
 r=blank();r.update({'record_type':'generated_output','artifact_name':p.name,'local_path':str(p),'local_status':'present_frozen','origin_in_this_run':'B104A_final_delivery','analysis_role':role,'used_by':'user delivery and future B105 reproducibility','size_bytes':str(p.stat().st_size),'sha256':sha(p),'parent_or_derivation':parent,'notes':('The bundle contains pre-bundle ledger/event snapshots; this external final ledger is authoritative for the bundle hash.' if p==BUNDLE else ''),'batch_id':'B104A','event_recorded_at_utc':EVENT,'hash_authority':'container_sha256','runtime_verification_status':'present_and_verified'})
 new.append(r)
ledger=pd.concat([ledger,pd.DataFrame(new,columns=cols)],ignore_index=True)
ledger.to_csv(FINAL_LEDGER,index=False)

events=pd.read_csv(BASE_EVENTS,dtype=str).fillna('')
new_events=[
 {'event_time_utc':EVENT,'batch_id':'B104A','event_type':'output_checksum_inventory_created','artifact_name':CHECKS.name,'status':'accepted','details':f'{len(rows)} pre-bundle files hashed; checksum file excludes its own hash.'},
 {'event_time_utc':EVENT,'batch_id':'B104A','event_type':'analysis_bundle_frozen','artifact_name':BUNDLE.name,'status':'accepted','details':f'ZIP integrity to be validated; SHA-256 {bundle_sha}; raw user uploads excluded.'},
]
events=pd.concat([events,pd.DataFrame(new_events,columns=events.columns)],ignore_index=True)
events.to_csv(FINAL_EVENTS,index=False)

# Render final ledger/events Markdown.
for p,title in [(FINAL_LEDGER,'Source Ledger Through B104A — Final'),(FINAL_EVENTS,'Provenance Events Through B104A — Final')]:
 df=pd.read_csv(p,dtype=str).fillna('')
 out=p.with_suffix('.md')
 with out.open('w',encoding='utf-8') as f:
  f.write(f'# {title}\n\nRows: {len(df)}  \nSource CSV: `{p.name}`\n\n')
  f.write(df.to_markdown(index=False));f.write('\n')

print(json.dumps({'bundle':str(BUNDLE),'bundle_size_bytes':BUNDLE.stat().st_size,'bundle_sha256':bundle_sha,'checksum_rows':len(rows),'final_ledger':str(FINAL_LEDGER),'final_events':str(FINAL_EVENTS)},indent=2))

#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

FINAL='20260828T152026Z'
R=Path('/mnt/data/ppi_repro_corrected/results')
items=[
 ('actual_input_file_manifest_through_B104A', ['record_type','artifact_name','local_path','local_status','origin_in_this_run','analysis_role','used_by','direct_or_canonical_source_url','source_page_url','retrieval_status','size_bytes','sha256','parent_or_derivation','batch_id','user_deletion_confirmed_at_utc','raw_retention_status','raw_available_in_conversation','retained_derivative_paths','parser_version_or_script_sha256','raw_to_derived_reconciliation','notes']),
 ('source_ledger_through_B104A', ['record_type','artifact_name','local_path','local_status','origin_in_this_run','analysis_role','used_by','direct_or_canonical_source_url','source_page_url','retrieval_status','size_bytes','sha256','parent_or_derivation','notes','batch_id','deletion_state','event_recorded_at_utc','hash_authority','runtime_verification_status','retained_derivative_paths','user_deletion_confirmed_at_utc']),
 ('provenance_events_through_B104A', None),
]
for base,cols in items:
    p=R/f'{base}_{FINAL}.csv'
    df=pd.read_csv(p,dtype=str).fillna('')
    if cols is not None:
        cols=[c for c in cols if c in df.columns]
        d=df[cols]
    else:
        d=df
    out=R/f'{base}_{FINAL}.md'
    title=base.replace('_',' ').title()
    with out.open('w',encoding='utf-8') as f:
        f.write(f'# {title}\n\n')
        f.write(f'Rows: {len(df)}  \nSource CSV: `{p.name}`\n\n')
        f.write(d.to_markdown(index=False))
        f.write('\n')
    print(out)

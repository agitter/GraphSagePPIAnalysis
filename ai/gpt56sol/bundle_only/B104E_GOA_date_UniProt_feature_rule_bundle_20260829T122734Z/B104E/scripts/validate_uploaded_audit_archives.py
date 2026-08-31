#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,hashlib,json,zipfile
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--uniprot-dir',type=Path,required=True);ap.add_argument('--goa-dir',type=Path,required=True);ap.add_argument('--uniprot-zip',type=Path,required=True);ap.add_argument('--goa-zip',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
 a=ap.parse_args(); checks=[]
 def ck(name,ok,details=''):
  checks.append({'check':name,'passed':bool(ok),'details':details})
 # ZIP integrity
 for label,p in [('uniprot',a.uniprot_zip),('goa',a.goa_zip)]:
  with zipfile.ZipFile(p) as z: bad=z.testzip()
  ck(label+'_zip_integrity',bad is None,str(bad))
 # UniProt ledger and retained file hashes
 ledger=list(csv.DictReader((a.uniprot_dir/'uniprot_2016_mapping_audit_ledger.csv').open(newline='',encoding='utf-8')))
 ck('uniprot_three_successful_releases',len(ledger)==3 and all(r['status']=='success' for r in ledger),f'rows={len(ledger)}')
 flat_summary=[]
 for r in ledger:
  rel=r['release']
  for kind,pathcol,hashcol in [('dat','record_dat_path','record_dat_sha256'),('tsv','summary_tsv_path','summary_tsv_sha256')]:
   p=a.uniprot_dir/Path(r[pathcol].replace('\\','/'))
   ck(f'uniprot_{rel}_{kind}_exists',p.exists(),str(p))
   if p.exists():ck(f'uniprot_{rel}_{kind}_sha256',sha(p)==r[hashcol],f'observed={sha(p)} expected={r[hashcol]}')
  prov=a.uniprot_dir/rel/f'O95073_Q9Y620_{rel}_provenance.json'
  d=json.load(prov.open())
  ck(f'uniprot_{rel}_archive_size',str(d['archive']['size_bytes'])==r['archive_observed_size'])
  ck(f'uniprot_{rel}_archive_md5',d['archive']['md5']==r['archive_observed_md5'])
  ck(f'uniprot_{rel}_archive_sha256',d['archive']['sha256']==r['archive_sha256'])
  ck(f'uniprot_{rel}_archive_deleted',bool(d['archive_deletion_policy']) and r['archive_deleted_after_success']=='1')
  tsv=a.uniprot_dir/rel/f'O95073_Q9Y620_{rel}.tsv'
  trs=list(csv.DictReader(tsv.open(newline='',encoding='utf-8'),delimiter='\t'))
  by={x['matched_target_accessions']:x for x in trs}
  ck(f'uniprot_{rel}_both_targets',set(by)=={'O95073','Q9Y620'},str(set(by)))
  ck(f'uniprot_{rel}_O95073_symbol',by.get('O95073',{}).get('GN')=='Name=FSBP;')
  ck(f'uniprot_{rel}_Q9Y620_symbol',by.get('Q9Y620',{}).get('GN')=='Name=RAD54B;')
  ck(f'uniprot_{rel}_O95073_geneids',by.get('O95073',{}).get('GeneID_cross_references')=='100861412|25788')
  ck(f'uniprot_{rel}_Q9Y620_geneids',by.get('Q9Y620',{}).get('GeneID_cross_references')=='25788')
  for x in trs: flat_summary.append({'release':rel,'release_date':r['release_date'],**x})
 # GOA summary vs JSONs
 summary_path=next(a.goa_dir.glob('goa_release_date_screen_summary_*.csv'))
 summary=list(csv.DictReader(summary_path.open(newline='',encoding='utf-8')))
 ck('goa_twelve_releases', [int(r['release']) for r in summary]==list(range(158,170)),str([r['release'] for r in summary]))
 exact_releases=[]; orientations=set(); lcs_by={}
 for r in summary:
  rel=int(r['release']); p=a.goa_dir/f'goa_release_{rel}_screen.json';d=json.load(p.open())
  checks_map={
   'exact_label_columns':d['label_comparison']['exact_columns'],
   'label_false_positives':d['label_comparison']['false_positives'],
   'label_false_negatives':d['label_comparison']['false_negatives'],
   'label_total_mismatches':d['label_comparison']['total_mismatches'],
   'terms_ge1000':d['term_selection']['terms_at_least_1000'],
   'candidate_overlap_ge1000':d['term_selection']['candidate_overlap_ge1000'],
   'candidate_overlap_top121':d['term_selection']['candidate_overlap_top121'],
   'order_lcs':d['column_order']['lcs'],
   'order_exact_positions':d['column_order']['exact_positions'],
   'order_exact_prefix':d['column_order']['exact_prefix'],
   'dictionary_terms':d['column_order']['unique_dictionary_terms'],
   'dictionary_table_size':d['column_order']['dictionary_table_size'],
  }
  ok=all(float(r[k])==float(v) for k,v in checks_map.items())
  ck(f'goa_release_{rel}_summary_matches_json',ok)
  if d['label_comparison']['total_mismatches']==0:exact_releases.append(rel)
  orientations.add(d['column_order']['duplicate_orientation']);lcs_by[rel]=d['column_order']['lcs']
 ck('goa_release159_uniquely_exact',exact_releases==[159],str(exact_releases))
 ck('goa_duplicate_orientation_stable',orientations=={'001'},str(orientations))
 # events include successful deletion of downloaded raw files for all releases
 events=list(csv.DictReader((a.goa_dir/'goa_date_screen_events.csv').open(newline='',encoding='utf-8')))
 for rel in range(158,170):
  rr=[e for e in events if str(rel) in (e.get('artifact_name','')+e.get('details','')+e.get('batch_id',''))]
  text=' '.join(str(v) for e in rr for v in e.values()).lower()
  ck(f'goa_release_{rel}_events_recorded',bool(rr),f'events={len(rr)}')
  ck(f'goa_release_{rel}_source_deleted_after_success','deleted' in text or 'delete' in text,text[:300])
 # Write normalized UniProt summary as sibling
 outdir=a.output.parent; outdir.mkdir(parents=True,exist_ok=True)
 upout=outdir/'B104E_UniProt_2016_O95073_Q9Y620_mapping_audit_summary.csv'
 with upout.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(flat_summary[0]));w.writeheader();w.writerows(flat_summary)
 result={'created_from':{'uniprot_zip':{'path':str(a.uniprot_zip),'bytes':a.uniprot_zip.stat().st_size,'sha256':sha(a.uniprot_zip)},'goa_zip':{'path':str(a.goa_zip),'bytes':a.goa_zip.stat().st_size,'sha256':sha(a.goa_zip)}},
         'checks':checks,'all_checks_passed':all(x['passed'] for x in checks),'uniprot_normalized_summary':str(upout),'goa_exact_releases':exact_releases,'goa_lcs_by_release':lcs_by}
 a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'all_checks_passed':result['all_checks_passed'],'checks':len(checks),'failed':[x for x in checks if not x['passed']]},indent=2))
 if not result['all_checks_passed']:raise SystemExit(1)
if __name__=='__main__':main()

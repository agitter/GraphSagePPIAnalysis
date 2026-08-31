#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,re
from pathlib import Path

def header_value(headers,prefix):
 for line in headers:
  if line.startswith(prefix):return line[len(prefix):].strip()
 return ''
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--summary-output',type=Path,required=True);ap.add_argument('--trajectory-output',type=Path,required=True);ap.add_argument('--analysis-json',type=Path,required=True)
 a=ap.parse_args(); data={}
 for p in sorted(a.input_dir.glob('goa_release_*_screen.json')):
  d=json.load(p.open());data[int(d['release'])]=d
 rows=[]
 for rel,d in sorted(data.items()):
  gafsrc=d['source_files']['gaf'];gpisrc=d['source_files']['gpi']
  row={
   'release':rel,'listed_release_date':d['listed_release_date'],
   'gaf_generated_header':header_value(d['gaf']['headers'],'!Generated:'),
   'gaf_GO_version_header':header_value(d['gaf']['headers'],'!GO-version:'),
   'gpi_generated_header':header_value(d['gpi']['headers'],'!Generated:'),
   'gaf_url':gafsrc.get('url_or_local',gafsrc.get('url','')) if isinstance(gafsrc,dict) else '',
   'gaf_size_bytes':gafsrc.get('size_bytes','') if isinstance(gafsrc,dict) else '',
   'gaf_sha256':gafsrc.get('sha256','') if isinstance(gafsrc,dict) else '',
   'gpi_url':gpisrc.get('url_or_local',gpisrc.get('url','')) if isinstance(gpisrc,dict) else '',
   'gpi_size_bytes':gpisrc.get('size_bytes','') if isinstance(gpisrc,dict) else '',
   'gpi_sha256':gpisrc.get('sha256','') if isinstance(gpisrc,dict) else '',
   'gaf_total_rows':d['gaf']['total_rows'],'gaf_accepted_rows':d['gaf']['accepted_rows'],
   'mapped_geneids':d['mapping']['mapped_geneids'],'symbol_fallback_edges':d['mapping']['unique_symbol_fallback_edges'],
   'exact_label_columns':d['label_comparison']['exact_columns'],'false_positives':d['label_comparison']['false_positives'],'false_negatives':d['label_comparison']['false_negatives'],'total_mismatches':d['label_comparison']['total_mismatches'],
   'terms_at_least_1000':d['term_selection']['terms_at_least_1000'],'candidate_overlap_ge1000':d['term_selection']['candidate_overlap_ge1000'],'candidate_overlap_top121':d['term_selection']['candidate_overlap_top121'],
   'rank121_count':d['term_selection']['rank121_count'],'rank122_count':d['term_selection']['rank122_count'],
   'order_lcs':d['column_order']['lcs'],'order_kendall_tau':d['column_order']['kendall_tau'],'order_exact_positions':d['column_order']['exact_positions'],'order_exact_prefix':d['column_order']['exact_prefix'],'duplicate_orientation':d['column_order']['duplicate_orientation'],
   'dictionary_unique_terms':d['column_order']['unique_dictionary_terms'],'dictionary_table_size':d['column_order']['dictionary_table_size'],
  };rows.append(row)
 a.summary_output.parent.mkdir(parents=True,exist_ok=True)
 with a.summary_output.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 # Per-column trajectory
 traj=[]
 for col in range(121):
  first=data[min(data)]['label_comparison']['per_column'][col]
  r={'label_column':col,'candidate_GO_ID':first['best_GO_ID'],'observed_positive_genes':first['observed_positive_genes']}
  for rel,d in sorted(data.items()):
   x=d['label_comparison']['per_column'][col]
   r[f'r{rel}_mismatches']=x['mismatches'];r[f'r{rel}_FP']=x['false_positives'];r[f'r{rel}_FN']=x['false_negatives'];r[f'r{rel}_predicted_positive']=x['predicted_positive_genes']
  exact=[rel for rel,d in sorted(data.items()) if d['label_comparison']['per_column'][col]['mismatches']==0]
  r['exact_releases']='|'.join(map(str,exact));r['release159_exact']=int(159 in exact);r['minimum_mismatches']=min(d['label_comparison']['per_column'][col]['mismatches'] for d in data.values())
  traj.append(r)
 with a.trajectory_output.open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(traj[0]));w.writeheader();w.writerows(traj)
 # conclusions strictly from result files
 exact_releases=[rel for rel,d in data.items() if d['label_comparison']['total_mismatches']==0]
 top121_exact=[rel for rel,d in data.items() if d['term_selection']['candidate_overlap_top121']==121]
 threshold_exact=[rel for rel,d in data.items() if d['term_selection']['terms_at_least_1000']==121 and d['term_selection']['candidate_overlap_ge1000']==121]
 best_other=sorted((d['label_comparison']['total_mismatches'],rel) for rel,d in data.items() if rel!=159)
 result={
  'release_range':[min(data),max(data)],'exact_membership_releases':exact_releases,'top121_candidate_set_exact_releases':top121_exact,'ge1000_candidate_set_exact_releases':threshold_exact,
  'closest_non159_release':{'release':best_other[0][1],'total_mismatches':best_other[0][0]},
  'column_order_best_lcs':max(d['column_order']['lcs'] for d in data.values()),
  'column_order_releases_with_best_lcs':[rel for rel,d in data.items() if d['column_order']['lcs']==max(x['column_order']['lcs'] for x in data.values())],
  'column_order_perfect_in_tested_range':any(d['column_order']['lcs']==121 and d['column_order']['exact_positions']==121 for d in data.values()),
  'duplicate_orientation_values':sorted(set(d['column_order']['duplicate_orientation'] for d in data.values())),
  'interpretive_limits':['Identifier mapping and ontology were fixed across releases; this isolates annotation-release effects but is not a joint search over all historical inputs.','The screen results were produced by the user-side audited script; source hashes and deletion events are retained, but raw releases 160-169 were not re-downloaded in this runtime.']
 }
 a.analysis_json.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()

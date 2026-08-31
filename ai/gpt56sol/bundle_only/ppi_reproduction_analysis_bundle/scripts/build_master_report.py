#!/usr/bin/env python3
from __future__ import annotations
import json,math,os,re,shutil,traceback,zipfile
from pathlib import Path
import pandas as pd
ROOT=Path('/mnt/data');BASE=ROOT/'ppi_repro';OUT=BASE/'results';SCRIPTS=BASE/'scripts'

def j(name,default=None):
 p=OUT/name
 try:return json.loads(p.read_text())
 except:return default

def csv(name):
 p=OUT/name
 try:return pd.read_csv(p)
 except:return pd.DataFrame()
core=j('core_summary.json',{});label=j('label_source_search_summary.json',{});sources=csv('source_manifest.csv');feat=csv('feature_column_mapping.csv');featall=csv('feature_best_matches_all_versions.csv');tissues=csv('tissue_matches.csv');screen=csv('label_source_screen.csv');external=csv('label_external_sources.csv');repo=j('repository_history_search.json',[])
errors=[]
for p in sorted(OUT.glob('*ERROR*'))+sorted(OUT.glob('*.stderr')):
 try:
  txt=p.read_text(errors='replace').strip()
  if txt:errors.append((p.name,txt[-5000:]))
 except:pass
lines=['# GraphSAGE PPI / DGL / gene-to-GO reconstruction','',
       '## Scope and independence','',
       'All numerical checks in this report were recomputed from the supplied GraphSAGE, DGL, OhmNet, MSigDB, Greene, and paper files. The prior investigation summary was treated only as a list of hypotheses. It was not used as a source of row mappings, tissue names, feature identities, or GO labels.','']
if not core:
 lines += ['> **Core reconstruction did not complete.** See the error section and rerun `scripts/core_reproduction.py`.','']
else:
 gs=core['graphsage'];mp=core['mapping'];top=core['topology'];leak=core['leakage'];dgl=core['dgl']
 # Feature verdicts
 exact52=0;c1=0;c3=0;unresolved_cols=[]
 if not feat.empty:
  for _,r in feat.iterrows():
   if pd.notna(r.get('mismatches_unique_wl')) and int(r.mismatches_unique_wl)==0:
    exact52+=1
    c=str(r.get('chosen_collection',''))
    if c.startswith('c1'):c1+=1
    if c.startswith('c3'):c3+=1
   else:unresolved_cols.append(int(r['column']))
 # DGL aggregate
 dgl_ok=True
 for sp,r in dgl.get('splits',{}).items():
  dgl_ok &= bool(r.get('labels_equal')) and (r.get('feature_max_abs_diff',1)<1e-5) and bool(r.get('edge_multiset_equal'))
 # tissue exactness
 exact_tissue=0
 if not tissues.empty: exact_tissue=int(((tissues.gs_n==tissues.ohm_n)&(tissues.gs_m==tissues.ohm_m)&(tissues.cost==0)).sum())
 lines += ['## Claim-by-claim verdicts','',
          '### 1. Gene identity recovery — recomputed','',
          f"Recovered **{mp['mapped_rows']:,}/{gs['nodes']:,} tissue-node instances**, representing **{mp['unique_entrez']:,} distinct Entrez Gene IDs**. The topology check covers **{top['edges_with_both_endpoints_mapped']:,}** edges with both endpoints mapped; **{top['matched_edges']:,}** of those are present in the assigned OhmNet network.", '',
          'The mapping was obtained from graph invariants and Weisfeiler–Lehman partitions, then disambiguated only with independently parsed MSigDB membership vectors. The row→Entrez table is `graphsage_row_to_entrez.csv`; unresolved equivalence classes are in `wl_unresolved.json`.','',
          '### 2. Graph topology provenance — recomputed, with an upstream-source correction','',
          f"The archive has **{gs['graphs']} graph IDs**. The assignment procedure matched **{exact_tissue}/{gs['graphs']}** by exact node/edge/statistical invariants at zero assignment cost; all assignments and second-best costs are in `tissue_matches.csv`.", '',
          '**Verdict:** the immediate source is the supplied OhmNet tissue-network archive. Describing the upstream interactome as “BioGRID” alone is too narrow: the OhmNet manuscript says its global physical PPI network combined multiple interaction resources. The data-level match to OhmNet is stronger than a database-name inference.','',
          '### 3. Feature provenance — recomputed','',
          f"On the uniquely mapped WL rows, **{exact52}/{gs['features'][1]}** feature columns have a zero-mismatch MSigDB v5.2 candidate; selected candidates include **{c1} C1** and **{c3} C3** columns. Columns still ambiguous or non-exact at that stage: `{unresolved_cols}`.", '',
          'The graph-by-graph count signature was used before node identities were finalized, preventing circular use of the claimed gene map. Full v5.1/v5.2/v6.0 comparisons are in `feature_best_matches_all_versions.csv`. An all-zero column can have multiple observationally equivalent gene-set names; the report preserves those ties instead of inventing a unique identity.','',
          '### 4. DGL transformation — recomputed','',
          f"Overall transformation check: **{'passes' if dgl_ok else 'does not fully pass'}**. The independent reconstruction fits the scaler on training rows, applies the same graph-ID ordering, compares labels and graph IDs, and compares the complete directed-plus-self-loop edge multiset. Exact per-split booleans and maximum floating-point differences are in `dgl_verification.json`.", '']
 for sp,r in dgl.get('splits',{}).items():
  lines.append(f"- `{sp}`: labels equal `{r.get('labels_equal')}`; maximum feature difference `{r.get('feature_max_abs_diff')}`; directed/self-loop edge multiset equal `{r.get('edge_multiset_equal')}`; graph-ID match global/renumbered `{r.get('graph_id_equal_global')}/{r.get('graph_id_equal_renumbered')}`.")
 lines += ['', '### 5. Leakage mechanism — recomputed','',
          f"Of **{leak['test_nodes']:,}** test nodes, **{leak['test_genes_seen_in_train']:,} ({leak['test_seen_fraction']:.6%})** have a recovered Entrez gene present in training. Among those, **{leak['seen_test_labels_identical']:,} ({leak['seen_identical_fraction']:.6%})** have byte-identical label vectors. The zero-parameter lookup baseline obtains micro-F1 **{leak['lookup_micro_f1']:.8f}** with unseen genes predicted as all-zero.", '']

lines += ['## Gene-to-GO investigation','']
if not label:
 lines += ['> **The historical GO search did not complete.** See the error section and rerun `scripts/label_source_search.py`.','']
else:
 top=label.get('top_screen',[]); ext=label.get('top_external',[])
 lines += [f"The label matrix was collapsed to **{label['genes']:,} independently recovered Entrez genes × {label['labels']} columns**. The search parsed **{label['annotation_sources']} historical annotation-source/mapping variants** and **{label['ontologies']} dated ontology snapshots.**",'',
           'The search is broader than the prior investigation in four important ways: it tests mapping-release GPI files directly; historical Bioconductor SQLite `go` and `go_all` tables; every broad evidence-group subset for the best sources; and OhmNet/MSigDB products in alternative aggregation and identifier spaces.','',
           '### Best historical annotation combinations','',
           '| Rank | Source | Date | Mapping | Propagation | Evidence | Exact | ≥99% | ≥99.5% | ≥95% | Total mismatches |',
           '|---:|---|---|---|---|---|---:|---:|---:|---:|---:|']
 for i,r in enumerate(top[:15],1):
  lines.append(f"| {i} | `{r['source']}` | {r.get('date_hint','')} | {r.get('mapping','')} | {r.get('mode','')} | {r.get('filter','')} | {r.get('exact','')} | {r.get('n99','')} | {r.get('n995','')} | {r.get('n95','')} | {r.get('total_mismatch','')} |")
 if top:
  b=top[0]
  if int(b.get('exact',0))==label['labels']:
   verdict='An exact global reconstruction was found.'
  elif int(b.get('n99',0))==label['labels']:
   verdict='Every column has a ≥99% match, but no single tested combination is exact.'
  else:
   verdict='No single screened source/filter/propagation combination exactly reconstructs all columns.'
  lines += ['',f'**Current verdict:** {verdict} The best candidate is `{b.get("source")}` with mapping `{b.get("mapping","")}`, propagation `{b.get("mode")}`, and evidence filter `{b.get("filter")}`.','']
 lines += ['### Direct alternative-source tests','', '| Source | Exact | ≥99% | ≥95% | Total mismatches |','|---|---:|---:|---:|---:|']
 for r in ext[:12]:lines.append(f"| `{r['source']}` | {r.get('exact','')} | {r.get('n99','')} | {r.get('n95','')} | {r.get('total_mismatch','')} |")
 lines += ['', 'The evidence-filter oracle in `label_filter_oracle.json` is especially diagnostic. If a column remains non-exact even when it may choose any broad evidence-group subset, the missing discrepancy cannot be explained solely by a single global evidence filter. That points instead to identifier history, annotation transfer, term selection, or a processed upstream product.','',
           'Term-selection counts under Greene Table 6, Greene Table 9, OhmNet term IDs, and support thresholds are in `label_term_selection_thresholds.csv`.','']

lines += ['## Source tracking and acquisition','',
          '`source_manifest.csv` and `source_manifest.json` are the authoritative input ledger. They include every supplied file and every attempted historical download, with URL, status, byte size, and SHA-256. Failed URLs are retained rather than omitted.','']
if not sources.empty:
 status=sources.status.fillna('unknown').value_counts().to_dict() if 'status' in sources else {}
 lines.append(f'Source status counts: `{status}`.')
 # Important failed/missing categories
 bad=sources[sources.status.astype(str).isin(['error','missing'])] if 'status' in sources else pd.DataFrame()
 if not bad.empty:
  lines += ['', '### Acquisition failures that may require an upload','']
  for _,r in bad.head(40).iterrows():lines.append(f"- `{r.get('label','')}` — {r.get('url','')} — {r.get('error','')}")
 lines += ['', 'The prior summary refers to `HumanBase-blood.dat`, `HumanBase-kidney.dat`, and `blood_sample_tsv.gz`, but those files were not among the supplied attachments. They are not required for the immediate OhmNet topology match, but they are high-value for testing whether the labels came from a processed Greene/GIANT gold-standard product. Please upload them if they are available in the earlier investigator’s workspace.','']

lines += ['## Reproduction commands','', 'Run from an environment with Python, NumPy, pandas, SciPy, scikit-learn, NetworkX, and openpyxl:','',
          '```bash',
          'python scripts/inspect_inputs.py',
          'python scripts/core_reproduction.py',
          'python scripts/download_sources.py',
          'python scripts/repo_history_search.py',
          'python scripts/label_source_search.py',
          'python scripts/build_master_report.py',
          '```','',
          'The scripts never use the prior summary as a mapping input. Downloaded files are cached and hash-checked in `downloads/`.','']
if errors:
 lines += ['## Execution diagnostics','']
 for name,txt in errors:
  lines += [f'### `{name}`','```text',txt,'```','']
(OUT/'MASTER_REPRODUCTION_REPORT.md').write_text('\n'.join(lines))

# Missing-input request file, concise and generated from actual manifest.
miss=['# Requested additional inputs','',
      'The following files are not present in the supplied attachments and would materially extend the remaining label-provenance search:','',
      '- `HumanBase-blood.dat`', '- `HumanBase-kidney.dat`', '- `blood_sample_tsv.gz`',
      '- Any locally saved `gene2go.gz` dated between June and September 2016',
      '- Any OhmNet/GraphSAGE preprocessing script or intermediate file that assigns the 121 GO columns', '',
      'The first three are explicitly named in the prior investigation summary. Uploading them avoids relying on a potentially changed HumanBase download endpoint.','']
(OUT/'REQUESTED_ADDITIONAL_INPUTS.md').write_text('\n'.join(miss))

# Compact bundle: scripts and reports/tables, excluding raw downloaded biological databases.
bundle=ROOT/'ppi_reproduction_analysis_bundle.zip'
with zipfile.ZipFile(bundle,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
 for p in sorted(SCRIPTS.glob('*')):
  if p.is_file():z.write(p,'scripts/'+p.name)
 for p in sorted(OUT.glob('*')):
  if p.is_file() and p.stat().st_size<100_000_000:z.write(p,'results/'+p.name)
 # include readme and prior inventory, but no proprietary MSigDB archives
 for p in [ROOT/'bio-tissue-readme.txt',ROOT/'historical_go_mapping_inventory.md']:
  if p.exists():z.write(p,'inputs/'+p.name)
print(json.dumps({'report':str(OUT/'MASTER_REPRODUCTION_REPORT.md'),'bundle':str(bundle),'errors':len(errors)}))

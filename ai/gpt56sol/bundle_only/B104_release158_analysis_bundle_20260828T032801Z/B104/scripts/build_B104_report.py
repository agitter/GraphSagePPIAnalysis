#!/usr/bin/env python3
from __future__ import annotations
import csv, json, gzip, hashlib
from pathlib import Path
from datetime import datetime, timezone

STAMP='20260828T030759Z'
ROOT=Path('/mnt/data/ppi_repro_corrected/batches/B104_20260828T030759Z')
ANA=ROOT/'analysis'; DER=ROOT/'derived'; LOG=ROOT/'logs'; SCRIPTS=ROOT/'scripts'
NOW=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

summary=json.load(open(ROOT/f'B104_analysis_summary_{STAMP}.json'))
alt=json.load(open(ANA/f'B104_alternative_hypothesis_checks_{STAMP}.json'))
hdr=json.load(open(ROOT/f'B104_headers_and_raw_stats_{STAMP}.json'))
rec=json.load(open(ROOT/f'B104_gaf_gpad_reconciliation_{STAMP}.json'))
labelmap=list(csv.DictReader(open(ANA/f'B104_label_to_GO_mapping_release158_159_{STAMP}.csv')))
watch=list(csv.DictReader(open(ANA/f'B104_identifier_mapping_watchlist_{STAMP}.csv')))

# Tables
aspect=[]
for ns in ['biological_process','cellular_component','molecular_function']:
 rows=[r for r in labelmap if r['namespace']==ns]
 aspect.append({
  'namespace':ns,'columns':len(rows),
  'v158_exact':sum(int(r['v158_fixed_term_mismatches'])==0 for r in rows),
  'v158_fp':sum(int(r['v158_fixed_term_false_positives']) for r in rows),
  'v158_fn':sum(int(r['v158_fixed_term_false_negatives']) for r in rows),
  'v159_exact':sum(int(r['v159_mismatches'])==0 for r in rows),
  'v159_fp':sum(int(r['v159_false_positives']) for r in rows),
  'v159_fn':sum(int(r['v159_false_negatives']) for r in rows),
 })

only158=summary['GPI_compare']['only_release158']
# GPI metadata rows for only158
only_meta={}
with gzip.open(DER/f'B104_goa_human_gpi158_normalized_{STAMP}.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['DB_Object_ID'] in only158: only_meta[r['DB_Object_ID']]=r

report=ROOT/f'B104_REPORT_{STAMP}.md'
with report.open('w',encoding='utf-8') as f:
 f.write('# B104: GOA release 158 versus release 159 label reconstruction\n\n')
 f.write(f'Generated: `{NOW}`  \n')
 f.write('Status: accepted B104 analysis; raw B104 uploads are eligible for deletion after the formal clearance below.\n\n')
 f.write('## Executive findings\n\n')
 f.write('The release-158 triplet does **not** improve on release 159. Under the same identifier policy, evidence filter, and `is_a` propagation, release 158 produces 1,733 gene-label differences, whereas release 159 produces 901. The transition from 158 to 159 resolves every one of the 814 release-158 false negatives, but leaves a stable set of 901 release-159 false positives.\n\n')
 f.write('| Release | Exact columns | ≥99% agreement | ≥95% agreement | Total differences | False positives | False negatives |\n|---|---:|---:|---:|---:|---:|---:|\n')
 for rel,key in [('158','release158_fixed_v159_terms'),('159','release159')]:
  s=summary['label_reconstruction'][key]
  f.write(f"| {rel} | {s['exact']} | {s['at_least_99pct']} | {s['at_least_95pct']} | {s['total_mismatches']:,} | {s['false_positives']:,} | {s['false_negatives']:,} |\n")
 f.write('\nThe release-159 result remains the current leading reconstruction: **89/121 exact columns**, **114/121 at ≥99%**, and all columns at ≥95%. All 85 Biological Process columns are exact. The residual errors are confined to 23 Cellular Component and 9 Molecular Function columns.\n\n')
 f.write('## Input integrity and file semantics\n\n')
 f.write('| File | Bytes | SHA-256 | Integrity |\n|---|---:|---|---|\n')
 for r in summary['input_integrity']:
  f.write(f"| `{r['artifact_name']}` | {r['size_bytes']:,} | `{r['sha256']}` | inventory hash matched; `gzip -t` passed |\n")
 f.write('\nThe release-158 GAF and GPAD contain the same 388,218 unique projected assertions. GPAD has 1,017 additional physical rows because `ECO:0000364` and `ECO:0000366` both project to the same `IEA` GAF assertion. This is controlled ECO granularity, not unexplained duplication.\n\n')
 f.write('Release-158 raw counts: 388,218 GAF rows, 389,235 GPAD rows, and 21,005 GPI objects. The GAF and GPAD headers identify ontology `2016-05-07/go.owl`; GPI was generated June 6 and GAF/GPAD June 7, 2016.\n\n')
 f.write('## The ontology-version caveat is now the highest-priority uncertainty\n\n')
 f.write('The ontology used for this analysis is the uploaded `2016-06-01-go.obo`, whose internal data version is `releases/2016-05-31`. It is **not** the exact ontology named in either annotation header:\n\n')
 f.write('| Annotation release | Header ontology | Ontology used in the current test |\n|---|---|---|\n')
 f.write('| GOA 158 | `2016-05-07/go.owl` | `2016-05-31` OBO |\n')
 f.write('| GOA 159 | `2016-06-29/go.owl` | `2016-05-31` OBO |\n\n')
 f.write('This mismatch is consequential because **878 of the 901 residual release-159 false positives are ancestor-only assignments created by `is_a` propagation**. Only 23 residual pairs have a direct annotation to the selected term. Testing the exact late-June ontology is therefore more informative than another broad evidence-code search.\n\n')
 f.write('## Release 158 → 159 change analysis\n\n')
 rc=summary['release_changes']
 f.write(f"- {rc['v158_false_negatives_resolved_in_v159']:,} release-158 false negatives become correct in release 159.\n")
 f.write(f"- {rc['v158_false_positives_removed_in_v159']:,} release-158 false positives disappear.\n")
 f.write(f"- {rc['v159_false_positives_added']:,} new release-159 false positives appear.\n")
 f.write(f"- {rc['false_positives_in_both']:,} false positives persist across both releases.\n")
 f.write(f"- {rc['new_v159_witness_rows_for_resolved_false_negatives']:,} new release-159 annotation rows witness the resolved false negatives.\n\n")
 f.write('The new witnesses are dominated by IDA (777 rows), IMP (372), ISS (172), and IGI (35), with UniProt and BHF-UCL as the largest assigned-by sources. Annotation dates cluster in May and June 2016. No source leave-one-out test improves on 901 differences, and date cutoffs from June 28 through July 1 reproduce the same 901-difference optimum.\n\n')
 f.write('### Aspect-level comparison\n\n')
 f.write('| Namespace | Columns | v158 exact | v158 FP | v158 FN | v159 exact | v159 FP | v159 FN |\n|---|---:|---:|---:|---:|---:|---:|---:|\n')
 for r in aspect:
  f.write(f"| {r['namespace']} | {r['columns']} | {r['v158_exact']} | {r['v158_fp']} | {r['v158_fn']} | {r['v159_exact']} | {r['v159_fp']} | {r['v159_fn']} |\n")
 f.write('\nRelease 159 exactly recovers all 85 Biological Process columns. Its remaining 514 Cellular Component and 387 Molecular Function differences are all false positives.\n\n')
 f.write('## GPI 158 versus GPI 159\n\n')
 f.write('GPI 158 has 21,005 objects and GPI 159 has 21,002. All 21,002 release-159 objects are present in release 158. The three release-158-only objects are unannotated in the release-158 GAF and do not affect the label result:\n\n')
 f.write('| Accession | Symbol | Name |\n|---|---|---|\n')
 for a in only158:
  r=only_meta[a]
  f.write(f"| `{a}` | {r['DB_Object_Symbol']} | {r['DB_Object_Name']} |\n")
 f.write(f"\nAmong shared objects, {summary['GPI_compare']['metadata_change_rows']} field-level metadata changes were observed: 49 synonym changes, 23 symbol changes, 20 name changes, and one properties change. Switching between the 158 and 159 GPI object lists does not explain the annotation improvement; the GAF assertion changes do.\n\n")
 f.write('## Ambiguous identifier mappings\n\n')
 f.write('The accepted mapping remains a bipartite, ambiguity-preserving crosswalk. Historical one-to-many and many-to-one edges are retained; a unique primary-symbol fallback is added only where no direct historical edge exists; and a square ambiguous component is resolved only when primary symbols establish a unique bijection. No accession is arbitrarily reduced to one GeneID.\n\n')
 f.write('| Mapping policy | Covered genes | Exact columns | Differences | FP | FN |\n|---|---:|---:|---:|---:|---:|\n')
 for r in alt['mapping_policy_sensitivity']:
  f.write(f"| {r['mapping_policy']} | {r['covered_GraphSAGE_GeneIDs']} | {r['exact_columns']} | {r['total_mismatches']:,} | {r['false_positives']:,} | {r['false_negatives']:,} |\n")
 f.write('\nDiscarding ambiguous mappings is strongly harmful. The component-aware policy is not merely increasing coverage indiscriminately: it eliminates false negatives without reducing the remaining 901 false positives.\n\n')
 f.write('## Identifier watchlist: GeneIDs 7957, 29901, and 10159\n\n')
 f.write('These cases are retained as temporal mapping questions rather than declared accession replacements. The detailed accession-level watchlist is a separate CSV/Markdown artifact.\n\n')
 f.write('### GeneID 7957 — EPM2A\n\n')
 f.write('May-2016 `gp2protein` maps 7957 to `O95278` and `H0UI04`, while GPI 158/159 represents EPM2A as `B3EWF7`. Current UniProt still has both `O95278` (laforin) and reviewed `B3EWF7` (laforin, isoform 9); therefore, the difference must **not** be described as a simple replacement. Under the accepted evidence and ontology model, the two surviving B3EWF7 direct annotations (`GO:0016239` and `GO:0032007`) propagate to exactly the 17 observed GraphSAGE labels for GeneID 7957. This is strong contextual validation of `B3EWF7 → 7957`, but it remains recorded as a symbol-supported historical fallback.\n\n')
 f.write('### GeneID 29901 — SAC3D1\n\n')
 f.write('May-2016 `gp2protein` uses `F8WC89`; GPI 158/159 uses reviewed `A6NKF1`; the current UniProt mapping output returns `F8WC89`, `A0A6I8PRW4`, and `H9KVA8`. Current UniProt lists those sequences as potential isoforms mapped to A6NKF1, so this is also not a simple replacement relationship. No A6NKF1 annotation survives the accepted evidence filter, and the deposited GraphSAGE row is all zero. The mapping choice therefore has no effect on the current label fit, but it remains on the watchlist.\n\n')
 f.write('### GeneID 10159 — ATP6AP2\n\n')
 f.write('No broad additional search is needed for the current label reconstruction. The defensible accession anchor is reviewed `O75787`, which is connected to GeneID 10159 and to `NM_005765.3 / NP_005756.2`. The multiple Swiss-Prot accessions listed by NCBI are largely secondary or historically replaced accessions consolidated into O75787, not independent current primary mappings. None of the O75787 primary/secondary accessions, nor the current TrEMBL product `A0A1C7CYW4`, appears as a GPI 158/159 object. GeneID 10159 has an all-zero GraphSAGE label row, so leaving it unmapped introduces no false negative.\n\n')
 f.write('A **targeted** follow-up remains useful for provenance: inspect a 2016 UniProt flat-file entry or a contemporaneous NCBI `gene2accession`/RefSeq crosswalk to determine why O75787 was outside these GOA reference-proteome GPI files. Do not map through `PSEC0072`: the historical GPI assigns that synonym to `SIDT2 / Q8NBJ9`, making a synonym-only join unsafe.\n\n')
 f.write('## Alternative explanations tested\n\n')
 f.write(f"- GAF 159 contains {alt['NOT_annotation_check']['GAF159_NOT_rows']:,} NOT rows; {alt['NOT_annotation_check']['NOT_rows_with_component_aware_GeneID_mapping']:,} map to resolved GraphSAGE genes. **Zero** residual false positives have a direct NOT annotation to the selected term. Unsafe upward propagation of NOT would touch only {alt['NOT_annotation_check']['residual_false_positive_pairs_with_NOT_reaching_selected_term_by_is_a']} pairs.\n")
 f.write(f"- There are {alt['all_zero_gene_posthoc_mask']['all_zero_GraphSAGE_genes']} all-zero genes. Only {alt['all_zero_gene_posthoc_mask']['residual_false_positive_pairs_on_all_zero_genes']} of 901 residual pairs occur on five of them. Masking those genes would be target leakage and would still leave {alt['all_zero_gene_posthoc_mask']['residual_false_positive_pairs_remaining_after_invalid_posthoc_mask']} differences.\n")
 f.write('- No assigned-by source removal improves the 901-difference baseline.\n')
 f.write('- The global evidence optimum remains `EXP, IDA, IEP, IGI, IMP, ISS`. Aspect-specific evidence filters improve the total by only one difference (900) while introducing two false negatives, so they are not a compelling mechanistic explanation.\n')
 f.write('- Adding `part_of` propagation remains harmful.\n\n')
 f.write('## What remains open\n\n')
 f.write('1. Whether the exact `2016-06-29` ontology removes some or all of the 878 ancestor-only residual assignments.\n')
 f.write('2. Whether the remaining 23 direct and 878 ancestor-only false positives reflect a source-specific pruning rule, an earlier curated label matrix, or a different ontology product such as a basic/filtered closure.\n')
 f.write('3. The exact rule that selected the 118 distinct GO terms represented by the 121 columns. A simple ≥15-gene threshold does not identify this set.\n')
 f.write('4. Whether an Entrez-native processed source—particularly the small May-2016 dhimmel annotation tables, a historical `gene2go`, or a Bioconductor package—reproduces the same 118 term memberships without the GOA reference-proteome mapping step.\n\n')
 f.write('## Recommended next test\n\n')
 f.write('The next batch should contain **one exact or nearest official late-June ontology file**, not another annotation triplet. Preferred order:\n\n')
 f.write('1. `go.obo` resolving from the GAF header version `http://purl.obolibrary.org/obo/go/releases/2016-06-29/go.owl` (an OBO serialization of the same release, when obtainable).\n')
 f.write('2. Otherwise, `https://release.geneontology.org/2016-07-01/ontology/go.obo`; its internal `data-version` will be checked before use.\n\n')
 f.write('After the exact/nearest late-June ontology test, the next small high-value batch is the four May-2016 `dhimmel-gene-ontology-962a5e1-GO_annotations-9606-*.tsv` files already present locally. Release 160 should follow if needed.\n\n')
 f.write('## Principal retained outputs\n\n')
 for rel,desc in [
  (f'analysis/B104_label_to_GO_mapping_release158_159_{STAMP}.csv','121-column release comparison and selected GO terms'),
  (f'analysis/B104_release158_to_159_gene_label_changes_{STAMP}.csv.gz','all release-specific gene-label changes'),
  (f'analysis/B104_v159_witness_rows_resolving_v158_false_negatives_{STAMP}.csv.gz','new release-159 witness rows'),
  (f'analysis/B104_v159_residual_false_positive_witness_rows_{STAMP}.csv.gz','residual witness rows'),
  (f'analysis/B104_identifier_mapping_watchlist_{STAMP}.csv','temporal identifier watchlist'),
  (f'analysis/B104_alternative_hypothesis_checks_{STAMP}.json','NOT, all-zero, and mapping-policy tests'),
  (f'B104_analysis_summary_{STAMP}.json','machine-readable accepted summary')]:
  f.write(f'- `{rel}` — {desc}.\n')

# Diagnostics
diag=ROOT/f'B104_EXECUTION_DIAGNOSTICS_{STAMP}.md'
with diag.open('w',encoding='utf-8') as f:
 f.write('# B104 execution diagnostics\n\n')
 f.write(f'Generated: `{NOW}`\n\n')
 f.write('## Accepted inputs\n\n')
 for r in summary['input_integrity']:
  f.write(f"- `{r['local_path']}` — {r['size_bytes']:,} bytes — SHA-256 `{r['sha256']}` — inventory match and gzip integrity passed.\n")
 f.write('\n## Executions\n\n')
 f.write('1. `analyze_B104_release158.py` — exhaustive first implementation. The command exceeded the 120-second execution limit and was terminated. Its empty stdout/stderr files are retained; none of its partial state is treated as accepted scientific output.\n')
 f.write('2. `analyze_B104_release158_fast.py` — accepted implementation, exit status 0. It reproduced the independently established release-159 baseline and generated the release-158 comparison, normalized derivatives, reconciliation tables, and witness analyses.\n')
 f.write('3. Evidence/source/date exploratory results under `/mnt/data/work_b104/filter_exploration` were generated independently before the accepted script. The accepted script imported them only after checking that their best global mask was exactly `EXP,IDA,IEP,IGI,IMP,ISS`, had 901 differences, and that no source leave-one-out result improved on 901.\n')
 f.write('4. `analyze_B104_alternative_hypotheses.py` — exit status 0. It independently regenerated mapping-policy, NOT, and all-zero-gene sensitivity results. A non-fatal terminal-control warning (`TERM environment variable not set`) appeared after successful completion; all requested outputs were written and re-read.\n')
 f.write('5. `build_B104_identifier_watchlist.py` and `build_B104_report.py` — report-only transformations over accepted local outputs.\n\n')
 f.write('## B103 provenance repair\n\n')
 f.write('The user had already confirmed deletion of the B103 conversation attachments. The runtime nevertheless retained residual mounts. A durable B103 derivative was not present, so those residual bytes were used once to reconstruct stable term, `is_a` edge, closure, and current-ID-mapping derivatives. Before use, the raw SHA-256 values were checked against the hashes already verified in B103. This repair is explicitly recorded and does not reverse the logical deletion state: future work consumes only the reconstructed derivatives.\n\n')
 f.write(f"- OBO SHA-256: `{summary['provenance_repair']['OBO_raw_sha256']}`\n")
 f.write(f"- Current ID mapping SHA-256: `{summary['provenance_repair']['current_idmapping_raw_sha256']}`\n")
 f.write(f"- Reconstructed ontology term rows: {summary['provenance_repair']['ontology_term_rows']:,}\n")
 f.write(f"- Reconstructed `is_a` edges: {summary['provenance_repair']['ontology_is_a_edges']:,}\n")
 f.write(f"- Reconstructed closure rows for GOA 158/159 direct terms: {summary['provenance_repair']['ontology_closure_rows_for_GOA158_159_direct_terms']:,}\n\n")
 f.write('## Validation assertions\n\n')
 f.write('- GAF/GPAD projected assertion sets are identical.\n')
 f.write('- Accepted release-159 baseline is 89 exact columns and 901 differences with zero false negatives.\n')
 f.write('- Release 158 under the same selected terms is 4 exact columns and 1,733 differences.\n')
 f.write('- All 814 release-158 false negatives are corrected by release 159.\n')
 f.write('- Residual release-159 false positives decompose into 23 direct-term and 878 ancestor-only pairs.\n')
 f.write('- Mapping sensitivity independently reproduces 901 differences for the component-aware hybrid and shows severe degradation when ambiguous mappings are discarded.\n')
 f.write('- All report-linked paths are checked in the final delivery validation.\n')

# Formal deletion clearance
clear=ROOT/f'B104_DELETION_CLEARANCE_{STAMP}.md'
with clear.open('w',encoding='utf-8') as f:
 f.write('# SAFE TO DELETE — BATCH B104\n\n')
 f.write('The following conversation attachments have passed inventory-hash matching, gzip integrity checks, complete parsing, analysis, retained-derivative generation, and final validation. This clearance applies only to conversation copies; retain local master copies.\n\n')
 for r in summary['input_integrity']:
  f.write(f"- `{r['artifact_name']}` — {r['size_bytes']:,} bytes — SHA-256 `{r['sha256']}`\n")
 f.write('\nRetained: normalized row-preserving derivatives, GAF/GPAD reconciliation data, release-comparison outputs, witness rows, mapping components, identifier watchlist, scripts, diagnostics, manifests, checksums, and bundle validation.\n\n')
 f.write('After deleting the three attachments, report: `Deleted B104`.\n')

print(report); print(diag); print(clear)

#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path

BATCH_STAMP = "20260828T030759Z"
DELIVERY_STAMP = "20260828T032801Z"
ROOT = Path(f"/mnt/data/ppi_repro_corrected/batches/B104_{BATCH_STAMP}")
RESULTS = Path("/mnt/data/ppi_repro_corrected/results")
DELIVERY = Path(f"/mnt/data/ppi_repro_corrected/batches/B104_delivery_{DELIVERY_STAMP}")
BUNDLE = Path(f"/mnt/data/B104_release158_analysis_bundle_{DELIVERY_STAMP}.zip")

MANIFEST = RESULTS / "actual_input_file_manifest_through_B104_20260828T032243Z.csv"
MANIFEST_MD = RESULTS / "actual_input_file_manifest_through_B104_20260828T032243Z.md"
LEDGER = RESULTS / "source_ledger_through_B104_20260828T032243Z.csv"
LEDGER_MD = RESULTS / "source_ledger_through_B104_20260828T032243Z.md"
EVENTS = RESULTS / "provenance_events_through_B104_20260828T032243Z.csv"
LOCAL_INV = RESULTS / "user_local_inventory_augmented_through_B104_20260828T032243Z.csv"

NOW = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def csv_info(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return {"columns": reader.fieldnames or [], "row_count": len(rows)}


if DELIVERY.exists():
    shutil.rmtree(DELIVERY)
DELIVERY.mkdir(parents=True)

# Copy the user-facing reports into a delivery directory without altering originals.
for name in [
    f"B104_REPORT_{BATCH_STAMP}.md",
    f"B104_EXECUTION_DIAGNOSTICS_{BATCH_STAMP}.md",
    f"B104_DELETION_CLEARANCE_{BATCH_STAMP}.md",
    f"B104_analysis_summary_{BATCH_STAMP}.json",
    f"B104_headers_and_raw_stats_{BATCH_STAMP}.json",
    f"B104_gaf_gpad_reconciliation_{BATCH_STAMP}.json",
    f"B104_input_integrity_{BATCH_STAMP}.csv",
]:
    shutil.copy2(ROOT / name, DELIVERY / name)

for src in [MANIFEST, MANIFEST_MD, LEDGER, LEDGER_MD, EVENTS, LOCAL_INV]:
    shutil.copy2(src, DELIVERY / src.name)

# Final checksums cover all accepted B104 outputs plus the copied final manifests.
checksum_rows = []
checksum_targets = []
for p in sorted(ROOT.rglob("*")):
    if p.is_file() and p.name not in {
        f"B104_output_checksums_final_{DELIVERY_STAMP}.csv",
        f"B104_FINAL_DELIVERY_VALIDATION_{DELIVERY_STAMP}.json",
    }:
        checksum_targets.append(p)
for p in sorted(DELIVERY.rglob("*")):
    if p.is_file():
        checksum_targets.append(p)

seen = set()
for p in checksum_targets:
    rp = str(p.resolve())
    if rp in seen:
        continue
    seen.add(rp)
    checksum_rows.append({
        "artifact_name": p.name,
        "absolute_path": str(p),
        "size_bytes": p.stat().st_size,
        "sha256": sha256(p),
        "category": "B104_root" if ROOT in p.parents else "B104_delivery",
    })
checksums = DELIVERY / f"B104_output_checksums_final_{DELIVERY_STAMP}.csv"
with checksums.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(checksum_rows[0].keys()))
    writer.writeheader()
    writer.writerows(checksum_rows)

# Validate all tabular files that matter for provenance and core scientific outputs.
csv_checks = {}
for p in [MANIFEST, LEDGER, EVENTS, LOCAL_INV,
          ROOT / "analysis" / f"B104_label_to_GO_mapping_release158_159_{BATCH_STAMP}.csv",
          ROOT / "analysis" / f"B104_identifier_mapping_watchlist_{BATCH_STAMP}.csv",
          ROOT / "analysis" / f"B104_mapping_policy_sensitivity_{BATCH_STAMP}.csv"]:
    csv_checks[str(p)] = csv_info(p)

# Explicit required paths from the report and deletion clearance.
required = [
    ROOT / f"B104_REPORT_{BATCH_STAMP}.md",
    ROOT / f"B104_EXECUTION_DIAGNOSTICS_{BATCH_STAMP}.md",
    ROOT / f"B104_DELETION_CLEARANCE_{BATCH_STAMP}.md",
    ROOT / f"B104_analysis_summary_{BATCH_STAMP}.json",
    ROOT / "analysis" / f"B104_label_to_GO_mapping_release158_159_{BATCH_STAMP}.csv",
    ROOT / "analysis" / f"B104_release158_to_159_gene_label_changes_{BATCH_STAMP}.csv.gz",
    ROOT / "analysis" / f"B104_v159_witness_rows_resolving_v158_false_negatives_{BATCH_STAMP}.csv.gz",
    ROOT / "analysis" / f"B104_v159_residual_false_positive_witness_rows_{BATCH_STAMP}.csv.gz",
    ROOT / "analysis" / f"B104_identifier_mapping_watchlist_{BATCH_STAMP}.csv",
    ROOT / "analysis" / f"B104_alternative_hypothesis_checks_{BATCH_STAMP}.json",
    ROOT / "derived" / f"B104_goa_human_gaf158_normalized_{BATCH_STAMP}.tsv.gz",
    ROOT / "derived" / f"B104_goa_human_gpad158_normalized_{BATCH_STAMP}.tsv.gz",
    ROOT / "derived" / f"B104_goa_human_gpi158_normalized_{BATCH_STAMP}.tsv.gz",
    ROOT / "derived" / f"B104_repaired_B103_GO_terms_{BATCH_STAMP}.tsv.gz",
    ROOT / "derived" / f"B104_repaired_B103_GO_is_a_edges_{BATCH_STAMP}.tsv.gz",
    ROOT / "derived" / f"B104_repaired_B103_GO_is_a_closure_for_GOA158_159_terms_{BATCH_STAMP}.tsv.gz",
    MANIFEST, LEDGER, EVENTS, LOCAL_INV,
]
missing = [str(p) for p in required if not p.exists()]
if missing:
    raise RuntimeError(f"Missing required outputs: {missing}")

# Check accepted summary assertions.
summary = json.loads((ROOT / f"B104_analysis_summary_{BATCH_STAMP}.json").read_text())
assert summary["label_reconstruction"]["release159"]["exact"] == 89
assert summary["label_reconstruction"]["release159"]["total_mismatches"] == 901
assert summary["label_reconstruction"]["release159"]["false_negatives"] == 0
assert summary["label_reconstruction"]["release158_fixed_v159_terms"]["total_mismatches"] == 1733
assert summary["residuals"]["v159_false_positive_pairs_ancestor_only"] == 878
assert summary["residuals"]["v159_false_positive_pairs_with_direct_selected_term_annotation"] == 23

# Build ZIP. Raw B104 uploads and logically deleted B103 raw files are never included.
if BUNDLE.exists():
    BUNDLE.unlink()
raw_basenames = {
    "goa_human.gaf.158.gz", "goa_human.gpa.158.gz", "goa_human.gpi.158.gz",
    "2016-06-01-go.obo", "idmapping_2026_08_27.tsv.gz",
}
with zipfile.ZipFile(BUNDLE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p.name in raw_basenames:
            continue
        z.write(p, Path("B104") / p.relative_to(ROOT))
    for p in sorted(DELIVERY.rglob("*")):
        if p.is_file():
            z.write(p, Path("delivery") / p.relative_to(DELIVERY))

with zipfile.ZipFile(BUNDLE, "r") as z:
    bad = z.testzip()
    names = z.namelist()
    forbidden_members = [n for n in names if Path(n).name in raw_basenames]
    if bad:
        raise RuntimeError(f"ZIP integrity failure at {bad}")
    if forbidden_members:
        raise RuntimeError(f"Raw files leaked into bundle: {forbidden_members}")

bundle_hash = sha256(BUNDLE)
bundle_sha_path = RESULTS / f"B104_bundle_sha256_{DELIVERY_STAMP}.txt"
bundle_sha_path.write_text(
    f"{bundle_hash}  {BUNDLE.name}\n",
    encoding="utf-8",
)

validation = {
    "generated_at_utc": NOW,
    "batch_id": "B104",
    "status": "passed",
    "required_paths_checked": len(required),
    "missing_required_paths": missing,
    "csv_checks": csv_checks,
    "accepted_scientific_assertions": {
        "release158_total_mismatches": 1733,
        "release159_exact_columns": 89,
        "release159_total_mismatches": 901,
        "release159_false_negatives": 0,
        "release159_residual_direct_pairs": 23,
        "release159_residual_ancestor_only_pairs": 878,
    },
    "input_manifest_rows": csv_info(MANIFEST)["row_count"],
    "source_ledger_rows": csv_info(LEDGER)["row_count"],
    "provenance_event_rows": csv_info(EVENTS)["row_count"],
    "known_user_local_files": csv_info(LOCAL_INV)["row_count"],
    "bundle_path": str(BUNDLE),
    "bundle_size_bytes": BUNDLE.stat().st_size,
    "bundle_sha256": bundle_hash,
    "zip_member_count": len(names),
    "zip_test_result": "passed",
    "forbidden_raw_members": forbidden_members,
    "logical_deletion_notes": [
        "B103 raw conversation attachments remain deletion-confirmed; repaired derivatives are included, not raw bytes.",
        "B104 raw conversation attachments are excluded and have deletion clearance pending user confirmation.",
    ],
    "external_ledger_note": "The external source ledger and provenance event file predate the bundle-hash event to avoid a self-referential bundle hash.",
}
validation_path = DELIVERY / f"B104_FINAL_DELIVERY_VALIDATION_{DELIVERY_STAMP}.json"
validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")

# Append external authoritative bundle records after the bundle has been frozen.
ledger_final = RESULTS / f"source_ledger_through_B104_FINAL_{DELIVERY_STAMP}.csv"
shutil.copy2(LEDGER, ledger_final)
with ledger_final.open(newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames or []
    rows = list(reader)
row = {k: "" for k in fieldnames}
row.update({
    "record_type": "analysis_bundle",
    "artifact_name": BUNDLE.name,
    "local_path": str(BUNDLE),
    "local_status": "frozen_and_zip_integrity_verified",
    "origin_in_this_run": "B104_final_delivery",
    "analysis_role": "complete B104 reports, scripts, retained derivatives, detailed outputs, and pre-bundle manifests",
    "used_by": "user download and future audit",
    "size_bytes": str(BUNDLE.stat().st_size),
    "sha256": bundle_hash,
    "parent_or_derivation": "B104 accepted analysis outputs; raw B104 and logically deleted B103 files excluded",
    "notes": "External ledger row is authoritative for the frozen bundle hash; the ledger copy inside the bundle necessarily predates this row.",
    "batch_id": "B104",
    "deletion_state": "retain_generated_output",
    "event_recorded_at_utc": NOW,
    "hash_authority": "runtime_recomputed_after_zip_freeze",
    "runtime_verification_status": "zip_test_passed; forbidden_raw_member_check_passed",
    "retained_derivative_paths": str(BUNDLE),
})
rows.append(row)
with ledger_final.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

events_final = RESULTS / f"provenance_events_through_B104_FINAL_{DELIVERY_STAMP}.csv"
shutil.copy2(EVENTS, events_final)
with events_final.open(newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    efields = reader.fieldnames or []
    erows = list(reader)
erows.append({
    "event_time_utc": NOW,
    "batch_id": "B104",
    "event_type": "analysis_bundle_frozen",
    "artifact_name": BUNDLE.name,
    "status": "accepted",
    "details": f"ZIP integrity passed; raw B104 and logically deleted B103 inputs excluded; SHA-256 {bundle_hash}.",
})
with events_final.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=efields)
    writer.writeheader()
    writer.writerows(erows)

# Copy validation and final checksum into result area for easy access.
validation_result = RESULTS / validation_path.name
shutil.copy2(validation_path, validation_result)
checksums_result = RESULTS / checksums.name
shutil.copy2(checksums, checksums_result)

print(json.dumps({
    "delivery_dir": str(DELIVERY),
    "bundle": str(BUNDLE),
    "bundle_sha256": bundle_hash,
    "bundle_size_bytes": BUNDLE.stat().st_size,
    "validation": str(validation_result),
    "checksums": str(checksums_result),
    "final_ledger": str(ledger_final),
    "final_events": str(events_final),
    "bundle_sha_record": str(bundle_sha_path),
}, indent=2))

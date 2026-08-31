#!/usr/bin/env python3
"""Package the small retained outputs from the UniProt 2016 mapping audit.

The downloader/audit script deletes each ~1.5 GB Swiss-Prot archive after
successful extraction.  This packager verifies the hashes recorded in the
ledger and creates a compact ZIP containing only the ledger, extracted records,
summary tables, provenance JSON files, and optional metalinks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PureWindowsPath
from typing import Iterable


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def normalize_ledger_path(value: str) -> Path:
    # The audit may have run on Windows, so its relative paths can contain '\\'.
    value = value.strip()
    if not value:
        return Path()
    return Path(*PureWindowsPath(value).parts)


def resolve_output(root: Path, ledger_value: str) -> Path:
    rel = normalize_ledger_path(ledger_value)
    candidates = []
    if rel:
        candidates.extend([root / rel, root / rel.name])
    # Last resort: find a unique basename recursively.  This is safe only when
    # unique; ambiguity is treated as an error rather than silently selecting.
    if rel and rel.name:
        candidates.extend(root.rglob(rel.name))
    existing = []
    seen = set()
    for candidate in candidates:
        try:
            key = candidate.resolve()
        except OSError:
            key = candidate
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            existing.append(candidate)
    if not existing:
        raise FileNotFoundError(
            f"Could not locate ledger output {ledger_value!r} below {root}"
        )
    if len(existing) > 1:
        exact = root / rel
        if exact.is_file():
            return exact
        raise RuntimeError(
            f"Ambiguous ledger output {ledger_value!r}: "
            + ", ".join(str(p) for p in existing)
        )
    return existing[0]


def expected_provenance_sha(notes: str) -> str:
    match = re.search(r"(?:^|[; ])provenance_json_sha256=([0-9a-fA-F]{64})(?:$|[; ])", notes)
    return match.group(1).lower() if match else ""


def arcname_for(root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return f"external/{path.name}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True,
                        help="Root directory containing release subdirectories")
    parser.add_argument("--ledger", type=Path, required=True,
                        help="uniprot_2016_mapping_audit_ledger.csv")
    parser.add_argument("--output", type=Path, default=Path("uniprot_2016_mapping_complete.zip"))
    parser.add_argument("--include-metalinks", action="store_true",
                        help="Include downloaded RELEASE.metalink files when present")
    args = parser.parse_args()

    root = args.root.expanduser().resolve()
    ledger = args.ledger.expanduser().resolve()
    if not root.is_dir():
        parser.error(f"--root is not a directory: {root}")
    if not ledger.is_file():
        parser.error(f"--ledger is not a file: {ledger}")

    rows = list(csv.DictReader(ledger.open("r", encoding="utf-8-sig", newline="")))
    if not rows:
        raise RuntimeError("Ledger has no data rows")

    required_fields = {
        "release", "status", "record_dat_path", "record_dat_sha256",
        "summary_tsv_path", "summary_tsv_sha256", "provenance_json_path", "notes"
    }
    missing_fields = sorted(required_fields - set(rows[0]))
    if missing_fields:
        raise RuntimeError(f"Ledger is missing required columns: {missing_fields}")

    files_to_add: list[tuple[Path, str, str, str]] = []
    manifest_rows: list[dict[str, str | int]] = []

    ledger_sha = sha256(ledger)
    files_to_add.append((ledger, ledger.name, ledger_sha, "audit_ledger"))

    for row in rows:
        if row.get("status", "").strip().lower() != "success":
            raise RuntimeError(
                f"Release {row.get('release', '?')} has non-success status: {row.get('status')!r}"
            )
        release = row["release"].strip()
        specs = [
            ("record_dat_path", "record_dat_sha256", "extracted_swissprot_records"),
            ("summary_tsv_path", "summary_tsv_sha256", "extracted_record_summary"),
        ]
        for path_field, hash_field, role in specs:
            path = resolve_output(root, row[path_field])
            observed = sha256(path)
            expected = row[hash_field].strip().lower()
            if observed != expected:
                raise RuntimeError(
                    f"SHA-256 mismatch for {release} {path_field}: expected {expected}, observed {observed}"
                )
            files_to_add.append((path, arcname_for(root, path), observed, role))

        prov = resolve_output(root, row["provenance_json_path"])
        prov_observed = sha256(prov)
        prov_expected = expected_provenance_sha(row.get("notes", ""))
        if not prov_expected:
            raise RuntimeError(f"No provenance_json_sha256 found in ledger notes for {release}")
        if prov_observed != prov_expected:
            raise RuntimeError(
                f"SHA-256 mismatch for {release} provenance JSON: expected {prov_expected}, observed {prov_observed}"
            )
        files_to_add.append((prov, arcname_for(root, prov), prov_observed, "release_provenance"))

        if args.include_metalinks and row.get("metalink_path", "").strip():
            try:
                metalink = resolve_output(root, row["metalink_path"])
            except FileNotFoundError:
                metalink = None
            if metalink is not None:
                files_to_add.append((metalink, arcname_for(root, metalink), sha256(metalink), "official_metalink"))

    # Deduplicate exact paths without erasing role information from the manifest.
    unique: dict[Path, tuple[Path, str, str, str]] = {}
    for item in files_to_add:
        unique[item[0].resolve()] = item
    files_to_add = list(unique.values())

    for path, arcname, digest, role in sorted(files_to_add, key=lambda x: x[1]):
        manifest_rows.append({
            "archive_path": arcname,
            "source_path": str(path),
            "size_bytes": path.stat().st_size,
            "sha256": digest,
            "role": role,
        })

    output = args.output.expanduser()
    stem = output.stem
    suffix = output.suffix or ".zip"
    if "{datestamp}" in output.name:
        output = output.with_name(output.name.replace("{datestamp}", utc_stamp()))
    elif not re.search(r"\d{8}T\d{6}Z", output.name):
        output = output.with_name(f"{stem}_{utc_stamp()}{suffix}")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".partial")
    if tmp.exists():
        tmp.unlink()

    manifest_csv = "package_manifest.csv"
    manifest_json = "package_manifest.json"
    readme = "README.txt"
    package_meta = {
        "created_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "ledger_sha256": ledger_sha,
        "release_count": len(rows),
        "file_count_excluding_generated_manifest": len(manifest_rows),
        "raw_swissprot_archives_included": False,
        "purpose": "Small verifiable outputs for auditing O95073 and Q9Y620 across UniProt releases 2016_04 through 2016_06",
    }

    with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for path, arcname, _, _ in files_to_add:
            zf.write(path, arcname)
        from io import StringIO
        sio = StringIO()
        writer = csv.DictWriter(sio, fieldnames=["archive_path", "source_path", "size_bytes", "sha256", "role"])
        writer.writeheader()
        writer.writerows(manifest_rows)
        zf.writestr(manifest_csv, sio.getvalue())
        zf.writestr(manifest_json, json.dumps({"package": package_meta, "files": manifest_rows}, indent=2, sort_keys=True) + "\n")
        zf.writestr(readme,
            "This archive intentionally excludes the ~1.5 GB UniProt Swiss-Prot source archives.\n"
            "Every retained extraction is verified against the SHA-256 values recorded by the audit ledger.\n"
            "The ledger also records the official archive URL, byte count, MD5, and source-archive SHA-256.\n")

    # Verify ZIP integrity before making it visible under its final name.
    with zipfile.ZipFile(tmp, "r") as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"ZIP integrity test failed at member {bad}")
    os.replace(tmp, output)
    print(json.dumps({
        "output": str(output),
        "size_bytes": output.stat().st_size,
        "sha256": sha256(output),
        "included_files": len(manifest_rows),
        "releases": [r["release"] for r in rows],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

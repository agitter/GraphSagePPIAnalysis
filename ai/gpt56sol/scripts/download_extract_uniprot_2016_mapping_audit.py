#!/usr/bin/env python3
"""Sequentially download, verify, extract, and delete large UniProt 2016 archives.

The script is intentionally conservative:
  * downloads one archive at a time;
  * verifies the official size and MD5 from UniProt RELEASE.metalink;
  * also records SHA-256 locally;
  * scans the complete Swiss-Prot flat file for O95073 and Q9Y620;
  * writes complete records, a compact TSV, JSON provenance, and an append-only CSV ledger;
  * deletes the large archive only after all checks and retained outputs succeed;
  * leaves the archive in place after any failure.

Requirements: Python 3.9+ and curl on PATH.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

TARGET_ACCESSIONS = ("O95073", "Q9Y620")
BASE_URL = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases"
RELEASES = {
    "2016_04": {
        "release_date": "2016-04-13",
        "archive_name": "uniprot_sprot-only2016_04.tar.gz",
        "expected_size": 1516525310,
        "expected_md5": "e607b83de1ac87e6f63b13715c049a3f",
    },
    "2016_05": {
        "release_date": "2016-05-11",
        "archive_name": "uniprot_sprot-only2016_05.tar.gz",
        "expected_size": 1504161063,
        "expected_md5": "fe9525832026b03ab34f0971b43c0c81",
    },
    "2016_06": {
        "release_date": "2016-06-08",
        "archive_name": "uniprot_sprot-only2016_06.tar.gz",
        "expected_size": 1504963399,
        "expected_md5": "e3a5ac5a166efc95e9ad06465d5bd2c4",
    },
}


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def hash_file(path: Path, algorithms=("md5", "sha256")) -> dict[str, str]:
    hs = {name: hashlib.new(name) for name in algorithms}
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            for h in hs.values():
                h.update(block)
    return {name: h.hexdigest() for name, h in hs.items()}


def sha256_file(path: Path) -> str:
    return hash_file(path, ("sha256",))["sha256"]


def run_curl(url: str, destination: Path, resume: bool = False) -> None:
    cmd = [
        "curl", "-fL", "--retry", "5", "--retry-delay", "5",
        "--connect-timeout", "30", "--speed-time", "120", "--speed-limit", "1024",
    ]
    if resume:
        cmd += ["--continue-at", "-"]
    cmd += ["--output", str(destination), url]
    subprocess.run(cmd, check=True)


def parse_metalink(path: Path, archive_name: str) -> dict[str, object]:
    root = ET.parse(path).getroot()
    ns = {"m": root.tag.split("}")[0].strip("{")} if "}" in root.tag else {}
    def q(tag: str) -> str:
        return f"m:{tag}" if ns else tag
    version = root.findtext(q("version"), namespaces=ns) or ""
    for file_element in root.findall(f".//{q('file')}", ns):
        if file_element.attrib.get("name") != archive_name:
            continue
        size_text = file_element.findtext(q("size"), namespaces=ns)
        md5 = ""
        for h in file_element.findall(f".//{q('hash')}", ns):
            if h.attrib.get("type", "").lower() == "md5":
                md5 = (h.text or "").strip().lower()
        urls = [(u.text or "").strip() for u in file_element.findall(f".//{q('url')}", ns)]
        return {"version": version, "size": int(size_text or 0), "md5": md5, "resource_urls": urls}
    raise RuntimeError(f"{archive_name} not found in {path}")


def accession_lines_from_record(lines: list[str]) -> set[str]:
    accessions: set[str] = set()
    for line in lines:
        if line.startswith("AC   "):
            payload = line[5:].strip()
            for value in payload.split(";"):
                value = value.strip()
                if value:
                    accessions.add(value)
    return accessions


def summarize_record(lines: list[str], release: str, member_name: str) -> dict[str, str]:
    accessions = accession_lines_from_record(lines)
    row: dict[str, str] = {
        "release": release,
        "archive_member": member_name,
        "matched_target_accessions": "|".join(sorted(accessions & set(TARGET_ACCESSIONS))),
        "all_accessions": "|".join(sorted(accessions)),
        "ID": "",
        "DT": "",
        "DE": "",
        "GN": "",
        "GeneID_cross_references": "",
    }
    dt_lines: list[str] = []
    de_lines: list[str] = []
    gn_lines: list[str] = []
    geneids: list[str] = []
    for line in lines:
        if line.startswith("ID   "):
            row["ID"] = line[5:].strip()
        elif line.startswith("DT   "):
            dt_lines.append(line[5:].strip())
        elif line.startswith("DE   "):
            de_lines.append(line[5:].strip())
        elif line.startswith("GN   "):
            gn_lines.append(line[5:].strip())
        elif line.startswith("DR   GeneID;"):
            parts = [x.strip() for x in line.split(";")]
            if len(parts) >= 2:
                geneids.append(parts[1])
    row["DT"] = " | ".join(dt_lines)
    row["DE"] = " ".join(de_lines)
    row["GN"] = " ".join(gn_lines)
    row["GeneID_cross_references"] = "|".join(sorted(set(geneids)))
    return row


def scan_swissprot_archive(archive: Path, release: str, output_dat: Path, output_tsv: Path) -> dict[str, object]:
    candidates = []
    found_records: list[tuple[list[str], set[str]]] = []
    member_name = ""
    member_declared_size: int | None = None
    scanned_records = 0
    # Streaming mode avoids extracting the full multi-gigabyte uncompressed flat file.
    with tarfile.open(archive, mode="r|gz") as tf:
        for member in tf:
            basename = Path(member.name).name.lower()
            if not member.isfile() or not (basename.endswith("uniprot_sprot.dat") or basename.endswith("uniprot_sprot.dat.gz")):
                continue
            candidates.append(member.name)
            if member_name:
                raise RuntimeError(f"multiple Swiss-Prot data members found: {candidates}")
            member_name = member.name
            member_declared_size = member.size
            raw = tf.extractfile(member)
            if raw is None:
                raise RuntimeError(f"unable to extract {member.name}")
            stream = gzip.GzipFile(fileobj=raw, mode="rb") if basename.endswith(".gz") else raw
            record_lines: list[str] = []
            for raw_line in stream:
                line = raw_line.decode("utf-8", "replace").rstrip("\n")
                record_lines.append(line)
                if line == "//":
                    scanned_records += 1
                    accessions = accession_lines_from_record(record_lines)
                    if accessions & set(TARGET_ACCESSIONS):
                        found_records.append((record_lines.copy(), accessions))
                    record_lines.clear()
            if record_lines:
                raise RuntimeError(f"unterminated final record in {member.name}")
    if not member_name:
        raise RuntimeError(f"no uniprot_sprot.dat[.gz] member found in {archive}")
    if scanned_records < 100_000:
        raise RuntimeError(f"implausibly few Swiss-Prot records scanned: {scanned_records}")

    output_dat.parent.mkdir(parents=True, exist_ok=True)
    with output_dat.open("w", encoding="utf-8", newline="\n") as fh:
        for lines, _ in found_records:
            fh.write("\n".join(lines) + "\n")

    summary_rows = [summarize_record(lines, release, member_name) for lines, _ in found_records]
    fields = ["release", "archive_member", "matched_target_accessions", "all_accessions", "ID", "DT", "DE", "GN", "GeneID_cross_references"]
    with output_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    found_targets = sorted(set().union(*(accessions & set(TARGET_ACCESSIONS) for _, accessions in found_records)) if found_records else set())
    return {
        "archive_member": member_name,
        "archive_member_declared_size": member_declared_size,
        "swissprot_records_scanned": scanned_records,
        "matching_records": len(found_records),
        "target_accessions_requested": list(TARGET_ACCESSIONS),
        "target_accessions_found": found_targets,
        "target_accessions_absent": sorted(set(TARGET_ACCESSIONS) - set(found_targets)),
    }


def append_ledger(path: Path, row: dict[str, object]) -> None:
    fields = [
        "event_time_utc", "release", "release_date", "archive_url", "metalink_url",
        "archive_expected_size", "archive_observed_size", "archive_expected_md5", "archive_observed_md5",
        "archive_sha256", "archive_deleted_after_success", "archive_path", "metalink_path",
        "data_member", "swissprot_records_scanned", "target_accessions_found", "target_accessions_absent",
        "record_dat_path", "record_dat_sha256", "summary_tsv_path", "summary_tsv_sha256",
        "provenance_json_path", "script_sha256", "status", "notes",
    ]
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def audit_release(release: str, work_dir: Path, output_dir: Path, keep_archives: bool, dry_run: bool) -> None:
    meta = RELEASES[release]
    archive_name = str(meta["archive_name"])
    archive_url = f"{BASE_URL}/release-{release}/knowledgebase/{archive_name}"
    metalink_url = f"{BASE_URL}/release-{release}/knowledgebase/RELEASE.metalink"
    release_work = work_dir / release
    release_out = output_dir / release
    release_work.mkdir(parents=True, exist_ok=True)
    release_out.mkdir(parents=True, exist_ok=True)
    archive_path = release_work / archive_name
    partial_path = release_work / f"{archive_name}.part"
    metalink_path = release_out / f"RELEASE_{release}.metalink"
    output_dat = release_out / f"O95073_Q9Y620_{release}.dat"
    output_tsv = release_out / f"O95073_Q9Y620_{release}.tsv"
    provenance_path = release_out / f"O95073_Q9Y620_{release}_provenance.json"
    ledger_path = output_dir / "uniprot_2016_mapping_audit_ledger.csv"
    script_path = Path(__file__).resolve()
    script_hash = sha256_file(script_path)

    print(f"\n[{release}] {archive_url}")
    print(f"  expected size: {meta['expected_size']} bytes")
    print(f"  expected MD5:  {meta['expected_md5']}")
    print(f"  scratch:       {archive_path}")
    if dry_run:
        return

    started = utc_now()
    deleted = False
    try:
        run_curl(metalink_url, metalink_path, resume=False)
        remote = parse_metalink(metalink_path, archive_name)
        if remote["version"] != release:
            raise RuntimeError(f"metalink version mismatch: {remote['version']} != {release}")
        if remote["size"] != meta["expected_size"] or str(remote["md5"]).lower() != meta["expected_md5"]:
            raise RuntimeError(f"metalink metadata changed or does not match hardcoded audit values: {remote}")

        free = shutil.disk_usage(release_work).free
        required = int(meta["expected_size"]) + 100_000_000
        if free < required:
            raise RuntimeError(f"insufficient free space: {free} bytes available; at least {required} required")

        # Resume into .part. If an already verified final archive exists, reuse it.
        if not archive_path.exists():
            run_curl(archive_url, partial_path, resume=True)
            if partial_path.stat().st_size != meta["expected_size"]:
                raise RuntimeError(f"download size mismatch: {partial_path.stat().st_size} != {meta['expected_size']}")
            partial_path.replace(archive_path)

        observed_size = archive_path.stat().st_size
        archive_hashes = hash_file(archive_path)
        if observed_size != meta["expected_size"]:
            raise RuntimeError(f"archive size mismatch: {observed_size} != {meta['expected_size']}")
        if archive_hashes["md5"].lower() != meta["expected_md5"]:
            raise RuntimeError(f"archive MD5 mismatch: {archive_hashes['md5']} != {meta['expected_md5']}")

        scan = scan_swissprot_archive(archive_path, release, output_dat, output_tsv)
        retained_hashes = {
            "record_dat_sha256": sha256_file(output_dat),
            "summary_tsv_sha256": sha256_file(output_tsv),
            "metalink_sha256": sha256_file(metalink_path),
        }
        provenance = {
            "audit_started_at_utc": started,
            "audit_completed_at_utc": utc_now(),
            "release": release,
            "release_date": meta["release_date"],
            "source_archive_url": archive_url,
            "official_metalink_url": metalink_url,
            "official_metalink_local_path": str(metalink_path),
            "official_metalink_parsed": remote,
            "archive": {
                "filename": archive_name,
                "local_path_before_deletion": str(archive_path),
                "size_bytes": observed_size,
                "md5": archive_hashes["md5"],
                "sha256": archive_hashes["sha256"],
            },
            "extraction": scan,
            "retained_outputs": {
                "record_dat": {"path": str(output_dat), "size_bytes": output_dat.stat().st_size, "sha256": retained_hashes["record_dat_sha256"]},
                "summary_tsv": {"path": str(output_tsv), "size_bytes": output_tsv.stat().st_size, "sha256": retained_hashes["summary_tsv_sha256"]},
                "metalink": {"path": str(metalink_path), "size_bytes": metalink_path.stat().st_size, "sha256": retained_hashes["metalink_sha256"]},
            },
            "script": {"path": str(script_path), "sha256": script_hash},
            "archive_deletion_policy": "delete only after size, MD5, SHA-256, complete record scan, and retained-output hashing succeed",
            "archive_kept_by_request": keep_archives,
        }
        # Write provenance before deletion, then hash it after writing.
        with provenance_path.open("w", encoding="utf-8") as fh:
            json.dump(provenance, fh, indent=2, sort_keys=True)
        provenance_hash = sha256_file(provenance_path)

        if not keep_archives:
            archive_path.unlink()
            deleted = True

        row = {
            "event_time_utc": utc_now(),
            "release": release,
            "release_date": meta["release_date"],
            "archive_url": archive_url,
            "metalink_url": metalink_url,
            "archive_expected_size": meta["expected_size"],
            "archive_observed_size": observed_size,
            "archive_expected_md5": meta["expected_md5"],
            "archive_observed_md5": archive_hashes["md5"],
            "archive_sha256": archive_hashes["sha256"],
            "archive_deleted_after_success": int(deleted),
            "archive_path": str(archive_path),
            "metalink_path": str(metalink_path),
            "data_member": scan["archive_member"],
            "swissprot_records_scanned": scan["swissprot_records_scanned"],
            "target_accessions_found": "|".join(scan["target_accessions_found"]),
            "target_accessions_absent": "|".join(scan["target_accessions_absent"]),
            "record_dat_path": str(output_dat),
            "record_dat_sha256": retained_hashes["record_dat_sha256"],
            "summary_tsv_path": str(output_tsv),
            "summary_tsv_sha256": retained_hashes["summary_tsv_sha256"],
            "provenance_json_path": str(provenance_path),
            "script_sha256": script_hash,
            "status": "success",
            "notes": f"provenance_json_sha256={provenance_hash}; archive deletion={deleted}",
        }
        append_ledger(ledger_path, row)
        print(f"  success: scanned {scan['swissprot_records_scanned']} records")
        print(f"  found:   {scan['target_accessions_found']}")
        print(f"  deleted large archive: {deleted}")
    except Exception as exc:
        # Never delete the archive on failure.
        row = {
            "event_time_utc": utc_now(), "release": release, "release_date": meta["release_date"],
            "archive_url": archive_url, "metalink_url": metalink_url,
            "archive_expected_size": meta["expected_size"], "archive_observed_size": archive_path.stat().st_size if archive_path.exists() else partial_path.stat().st_size if partial_path.exists() else "",
            "archive_expected_md5": meta["expected_md5"], "archive_observed_md5": "", "archive_sha256": "",
            "archive_deleted_after_success": 0, "archive_path": str(archive_path), "metalink_path": str(metalink_path),
            "data_member": "", "swissprot_records_scanned": "", "target_accessions_found": "", "target_accessions_absent": "",
            "record_dat_path": str(output_dat), "record_dat_sha256": "", "summary_tsv_path": str(output_tsv), "summary_tsv_sha256": "",
            "provenance_json_path": str(provenance_path), "script_sha256": script_hash,
            "status": "failed_archive_retained", "notes": repr(exc),
        }
        append_ledger(ledger_path, row)
        print(f"  FAILED; archive/partial download retained for diagnosis or resume: {exc}", file=sys.stderr)
        raise


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="uniprot_audit_test_") as tmp:
        root = Path(tmp)
        source_dat = root / "uniprot_sprot.dat"
        source_dat.write_text(
            "ID   TEST1_HUMAN             Reviewed;          10 AA.\n"
            "AC   O95073; X11111;\n"
            "DT   01-JAN-2000, integrated into UniProtKB/Swiss-Prot.\n"
            "GN   Name=FSBP;\n"
            "DR   GeneID; 100861412; - .\n"
            "//\n"
            "ID   TEST2_HUMAN             Reviewed;          10 AA.\n"
            "AC   Q9Y620;\n"
            "GN   Name=RAD54B;\n"
            "DR   GeneID; 25788; - .\n"
            "//\n"
            "ID   OTHER_HUMAN             Reviewed;          10 AA.\n"
            "AC   P00001;\n"
            "//\n",
            encoding="utf-8",
        )
        archive = root / "synthetic.tar.gz"
        with tarfile.open(archive, "w:gz") as tf:
            tf.add(source_dat, arcname="synthetic/uniprot_sprot.dat")
        out_dat = root / "out.dat"
        out_tsv = root / "out.tsv"
        # The production guard requires >100k records. Patch with a direct tiny-test implementation check by duplicating records.
        big_dat = root / "big_uniprot_sprot.dat"
        with big_dat.open("w", encoding="utf-8") as fh:
            fh.write(source_dat.read_text(encoding="utf-8"))
            for i in range(100_001):
                fh.write(f"ID   X{i}_HUMAN Reviewed; 1 AA.\nAC   Z{i:06d};\n//\n")
        with tarfile.open(archive, "w:gz") as tf:
            tf.add(big_dat, arcname="synthetic/uniprot_sprot.dat")
        result = scan_swissprot_archive(archive, "SELF_TEST", out_dat, out_tsv)
        if result["target_accessions_found"] != ["O95073", "Q9Y620"]:
            raise AssertionError(result)
        text = out_dat.read_text(encoding="utf-8")
        if "O95073" not in text or "Q9Y620" not in text or "P00001" in text:
            raise AssertionError("self-test extraction content is incorrect")
        print(json.dumps({"self_test": "passed", "scan": result, "output_dat_sha256": sha256_file(out_dat), "output_tsv_sha256": sha256_file(out_tsv)}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--releases", nargs="+", choices=sorted(RELEASES), default=sorted(RELEASES), help="releases to process sequentially")
    parser.add_argument("--work-dir", type=Path, default=Path("uniprot_audit_work"), help="scratch directory holding at most one large archive at a time")
    parser.add_argument("--output-dir", type=Path, default=Path("uniprot_audit_results"), help="small retained records, summaries, provenance, and ledger")
    parser.add_argument("--keep-archives", action="store_true", help="do not delete verified archives after successful extraction")
    parser.add_argument("--dry-run", action="store_true", help="print the sequential plan without downloading")
    parser.add_argument("--self-test", action="store_true", help="run a local synthetic archive test and exit")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if shutil.which("curl") is None:
        raise SystemExit("curl is required but was not found on PATH")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for release in args.releases:
        audit_release(release, args.work_dir, args.output_dir, args.keep_archives, args.dry_run)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Verify supplied inputs or download public sources recorded in the corrected manifest.

Examples:
  python download_or_verify_sources.py --manifest results/actual_input_file_manifest_v2.csv --dest inputs --verify-only
  python download_or_verify_sources.py --manifest results/actual_input_file_manifest_v2.csv --dest inputs --download-missing
  python download_or_verify_sources.py --manifest results/actual_input_file_manifest_v2.csv --dest historical --download-missing --include-historical

MSigDB archives are intentionally not downloaded automatically because the official site
requires an authenticated account. Place those files in --dest and rerun verification.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def is_direct_download_url(url: str) -> bool:
    if not url.startswith(('http://', 'https://')):
        return False
    lower = url.lower().split('?', 1)[0]
    return lower.endswith(('.zip', '.tar.gz', '.gz', '.obo', '.owl', '.txt', '.json', '.csv', '.tsv'))


def download(url: str, dest: Path, retries: int = 3) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(1, retries + 1):
        tmp = dest.with_suffix(dest.suffix + '.part')
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'GraphSAGE-PPI-reproduction/2.0'})
            with urllib.request.urlopen(req, timeout=120) as response, tmp.open('wb') as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            os.replace(tmp, dest)
            return
        except Exception as exc:
            last = exc
            if tmp.exists():
                tmp.unlink()
            if attempt < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f'Failed to download {url}: {last!r}')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--dest', type=Path, required=True)
    ap.add_argument('--verify-only', action='store_true')
    ap.add_argument('--download-missing', action='store_true')
    ap.add_argument('--include-historical', action='store_true')
    ap.add_argument('--artifact', action='append', default=[], help='Process only the named artifact; repeat for multiple files')
    ap.add_argument('--log', type=Path, default=None)
    args = ap.parse_args()
    if not args.verify_only and not args.download_missing:
        ap.error('Choose --verify-only or --download-missing')
    args.dest.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.manifest.open()))
    results: List[Dict[str, object]] = []
    failures = 0
    for row in rows:
        record_type = row['record_type']
        if record_type == 'actual_input':
            pass
        elif record_type == 'historical_candidate_not_materialized' and args.include_historical:
            pass
        else:
            continue
        name = row['artifact_name']
        if args.artifact and name not in set(args.artifact):
            continue
        expected = row.get('sha256', '')
        original_path = Path(row['local_path']) if row.get('local_path') else None
        target = args.dest / name
        source_url = row.get('direct_or_canonical_source_url', '')
        status = ''
        detail = ''
        used_path = None
        if original_path and original_path.exists():
            used_path = original_path
        elif target.exists():
            used_path = target
        elif args.download_missing and is_direct_download_url(source_url):
            try:
                download(source_url, target)
                used_path = target
                status = 'downloaded'
            except Exception as exc:
                status = 'download_failed'
                detail = repr(exc)
                failures += 1
        elif args.download_missing and source_url:
            status = 'manual_acquisition_required'
            detail = 'URL is a source page/DOI/authenticated endpoint rather than a direct public file URL.'
        else:
            status = 'missing'
            detail = 'File not present at manifest local_path or destination.'
            if record_type == 'actual_input':
                failures += 1
        if used_path is not None:
            actual = sha256(used_path)
            if expected and actual != expected:
                status = 'checksum_mismatch'
                detail = f'expected {expected}, got {actual}'
                failures += 1
            else:
                status = status or 'verified_present'
                detail = f'sha256={actual}; bytes={used_path.stat().st_size}'
        results.append({
            'record_type': record_type,
            'artifact_name': name,
            'status': status,
            'path': str(used_path or target),
            'source_url': source_url,
            'detail': detail,
        })
    log_path = args.log or (args.dest / 'download_or_verify_log.csv')
    with log_path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()) if results else ['record_type'])
        w.writeheader()
        w.writerows(results)
    counts = {}
    for row in results:
        counts[row['status']] = counts.get(row['status'], 0) + 1
    print(f'Checked {len(results)} records; status counts: {counts}; log: {log_path}')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())

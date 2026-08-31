#!/usr/bin/env python3
import csv, hashlib, pathlib, sys
root = pathlib.Path(__file__).resolve().parents[1]
manifest = root / 'metadata' / 'SNAPSHOT_FILE_MANIFEST.csv'
failures = []
with manifest.open(newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        path = root / row['relative_path']
        if not path.is_file():
            failures.append(f"missing: {row['relative_path']}")
            continue
        h = hashlib.sha256()
        with path.open('rb') as fh:
            for block in iter(lambda: fh.read(1024 * 1024), b''):
                h.update(block)
        if h.hexdigest() != row['sha256']:
            failures.append(f"hash mismatch: {row['relative_path']}")
if failures:
    print('\n'.join(failures), file=sys.stderr)
    raise SystemExit(1)
print('All snapshot files verified.')

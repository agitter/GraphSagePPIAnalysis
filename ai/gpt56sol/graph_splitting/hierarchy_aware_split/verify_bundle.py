#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import hashlib
import sys

root = Path(__file__).resolve().parent
manifest = root / 'SHA256SUMS.txt'
failures = []
checked = 0
for line in manifest.read_text(encoding='utf-8').splitlines():
    if not line.strip():
        continue
    expected, rel = line.split('  ', 1)
    path = root / rel
    if not path.is_file():
        failures.append(f'MISSING {rel}')
        continue
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    actual = h.hexdigest()
    checked += 1
    if actual != expected:
        failures.append(f'HASH {rel} expected={expected} actual={actual}')
if failures:
    print('\n'.join(failures))
    print(f'FAIL: checked={checked} failures={len(failures)}')
    sys.exit(1)
print(f'PASS: checked {checked} files')

#!/usr/bin/env python3
"""Create a provenance inventory for files that may be uploaded in batches.

The script never opens archives beyond reading their bytes for SHA-256. It does
not modify source files. The resulting CSV is intended to be uploaded before
large files so batches can be planned without consuming attachment storage.

Every output filename receives a UTC datestamp in YYYYMMDDTHHMMSSZ format.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

OUTPUT_FIELDS = [
    "artifact_name",
    "relative_path",
    "size_bytes",
    "sha256",
    "mtime_utc",
    "file_suffixes",
    "direct_or_canonical_source_url",
    "source_page_url",
    "obtained_from",
    "acquired_at_utc",
    "provenance_notes",
]

SOURCE_MAP_FIELDS = {
    "artifact_name",
    "direct_or_canonical_source_url",
    "source_page_url",
    "obtained_from",
    "acquired_at_utc",
    "provenance_notes",
}


def utc_iso(timestamp: float) -> str:
    return datetime.fromtimestamp(timestamp, tz=timezone.utc).isoformat()


def utc_datestamp(now: datetime | None = None) -> str:
    """Return a filesystem-safe UTC datestamp such as 20260827T154233Z."""
    value = now or datetime.now(tz=timezone.utc)
    return value.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def datestamped_output_path(requested: Path, datestamp: str) -> Path:
    """Return an output path containing *datestamp*.

    If ``{datestamp}`` occurs anywhere in the requested path, it is replaced.
    Otherwise, the datestamp is inserted before the final suffix. A path with no
    suffix receives ``.csv``.

    Examples:
        local_upload_inventory.csv
          -> local_upload_inventory_20260827T154233Z.csv
        inventories/{datestamp}_inventory.csv
          -> inventories/20260827T154233Z_inventory.csv
        local_upload_inventory
          -> local_upload_inventory_20260827T154233Z.csv
    """
    requested_text = str(requested)
    if "{datestamp}" in requested_text:
        return Path(requested_text.replace("{datestamp}", datestamp))

    if requested.suffix:
        stamped_name = f"{requested.stem}_{datestamp}{requested.suffix}"
    else:
        stamped_name = f"{requested.name}_{datestamp}.csv"
    return requested.with_name(stamped_name)


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def load_patterns(path: Path | None) -> List[str]:
    if path is None:
        return []
    patterns: List[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            value = raw.strip()
            if value and not value.startswith("#"):
                patterns.append(value)
    return patterns


def load_source_map(path: Path | None) -> Mapping[str, Dict[str, str]]:
    if path is None:
        return {}
    result: Dict[str, Dict[str, str]] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "artifact_name" not in reader.fieldnames:
            raise ValueError("source-map CSV must contain an artifact_name column")
        for row in reader:
            name = (row.get("artifact_name") or "").strip()
            if not name:
                continue
            result[name] = {
                key: (row.get(key) or "").strip()
                for key in SOURCE_MAP_FIELDS
                if key != "artifact_name"
            }
    return result


def matches_any(value: str, patterns: Sequence[str]) -> bool:
    return any(fnmatch.fnmatch(value, pattern) for pattern in patterns)


def iter_files(
    roots: Sequence[Path],
    include_patterns: Sequence[str],
    exclude_patterns: Sequence[str],
) -> Iterable[tuple[Path, Path]]:
    seen: set[Path] = set()
    for root in roots:
        resolved_root = root.expanduser().resolve()
        if not resolved_root.exists():
            raise FileNotFoundError(f"Input path does not exist: {root}")

        candidates = (
            [resolved_root]
            if resolved_root.is_file()
            else resolved_root.rglob("*")
        )
        for candidate in candidates:
            if not candidate.is_file() or candidate.is_symlink():
                continue
            resolved = candidate.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)

            rel = (
                candidate.name
                if resolved_root.is_file()
                else str(candidate.relative_to(resolved_root))
            )
            target_values = (candidate.name, rel)
            if include_patterns and not any(
                matches_any(value, include_patterns) for value in target_values
            ):
                continue
            if exclude_patterns and any(
                matches_any(value, exclude_patterns) for value in target_values
            ):
                continue
            yield resolved_root, candidate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Hash local files and write a datestamped provenance inventory CSV "
            "for staged uploads."
        )
    )
    parser.add_argument(
        "paths",
        nargs="+",
        type=Path,
        help="One or more files or directories to inventory recursively.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("local_upload_inventory.csv"),
        help=(
            "Base output CSV path. A UTC datestamp is always inserted before "
            "the suffix. Use {datestamp} to choose its position "
            "(default base: local_upload_inventory.csv)."
        ),
    )
    parser.add_argument(
        "--datestamp",
        help=(
            "Optional explicit datestamp for reproducible naming. Default: "
            "current UTC time in YYYYMMDDTHHMMSSZ format."
        ),
    )
    parser.add_argument(
        "--patterns-file",
        type=Path,
        help=(
            "Optional text file containing filename/glob patterns to include, "
            "one per line."
        ),
    )
    parser.add_argument(
        "--include",
        action="append",
        default=[],
        help="Additional filename/glob pattern to include; may be repeated.",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=[".DS_Store", "Thumbs.db", "*.tmp", "*.part"],
        help="Filename/glob pattern to exclude; may be repeated.",
    )
    parser.add_argument(
        "--source-map",
        type=Path,
        help=(
            "Optional CSV keyed by artifact_name with URL/provenance columns. "
            "Blank values are acceptable."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inventory_datestamp = args.datestamp or utc_datestamp()
    output_path = datestamped_output_path(args.output, inventory_datestamp)

    include_patterns = load_patterns(args.patterns_file) + list(args.include)
    source_map = load_source_map(args.source_map)
    rows: List[Dict[str, str]] = []

    for root, path in iter_files(args.paths, include_patterns, args.exclude):
        stat = path.stat()
        rel = path.name if root.is_file() else str(path.relative_to(root))
        provenance = source_map.get(path.name, {})
        rows.append(
            {
                "artifact_name": path.name,
                "relative_path": rel,
                "size_bytes": str(stat.st_size),
                "sha256": sha256_file(path),
                "mtime_utc": utc_iso(stat.st_mtime),
                "file_suffixes": "".join(path.suffixes),
                "direct_or_canonical_source_url": provenance.get(
                    "direct_or_canonical_source_url", ""
                ),
                "source_page_url": provenance.get("source_page_url", ""),
                "obtained_from": provenance.get("obtained_from", ""),
                "acquired_at_utc": provenance.get("acquired_at_utc", ""),
                "provenance_notes": provenance.get("provenance_notes", ""),
            }
        )

    rows.sort(
        key=lambda row: (
            row["artifact_name"].lower(),
            row["relative_path"].lower(),
        )
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    total_bytes = sum(int(row["size_bytes"]) for row in rows)
    print(f"Wrote {len(rows)} records to {output_path}")
    print(f"Inventory datestamp (UTC): {inventory_datestamp}")
    print(f"Total inventoried bytes: {total_bytes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

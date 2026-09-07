"""Source acquisition, integrity checks, and run-level provenance.

The source inventory deliberately separates *upstream* files, which may affect
reconstruction, from released *reference* files, which validation alone may
read.  Every file is accepted only when both its byte size and SHA-256 match the
committed specification.  Expected identities are never learned from whatever
happens to be available at a URL.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import platform
import shutil
import subprocess
import tarfile
import tempfile
import urllib.error
import urllib.request
import zipfile
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import yaml

SOURCE_FIELDS = (
    "source_id",
    "role",
    "acquisition",
    "filename",
    "primary_url",
    "mirror_url",
    "sha256",
    "size_bytes",
    "license",
    "description",
)
VALID_ROLES = frozenset({"upstream", "reference"})
VALID_ACQUISITION = frozenset({"download", "manual"})


class SourceError(RuntimeError):
    """Raised when a source record or file violates the frozen specification."""


@dataclass(frozen=True)
class SourceRecord:
    """One row from :file:`spec/sources.tsv`."""

    source_id: str
    role: str
    acquisition: str
    filename: str
    primary_url: str
    mirror_url: str
    sha256: str
    size_bytes: int
    license: str
    description: str

    @property
    def urls(self) -> tuple[str, ...]:
        """Return configured acquisition locations in priority order."""

        return tuple(url for url in (self.primary_url, self.mirror_url) if url)


@dataclass(frozen=True)
class FileVerification:
    """Observed identity and structural check for one source file."""

    source_id: str
    role: str
    path: str
    size_bytes: int
    sha256: str
    archive_or_format_check: str
    acquisition_status: str
    resolved_url: str | None = None


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp with second precision."""

    return datetime.now(UTC).isoformat(timespec="seconds")


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of *path* without loading it into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_atomic(path: Path, value: object) -> None:
    """Write stable UTF-8 JSON through a temporary file beside *path*."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_text_atomic(path: Path, text: str) -> None:
    """Write UTF-8 text with LF newlines through an adjacent temporary file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(text)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_source_manifest(path: Path) -> list[SourceRecord]:
    """Parse and strictly validate the flat source inventory."""

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SOURCE_FIELDS:
            raise SourceError(
                f"Unexpected columns in {path}: {reader.fieldnames!r}; "
                f"expected {list(SOURCE_FIELDS)!r}"
            )

        records: list[SourceRecord] = []
        seen_ids: set[str] = set()
        seen_filenames: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            try:
                record = SourceRecord(
                    source_id=row["source_id"].strip(),
                    role=row["role"].strip(),
                    acquisition=row["acquisition"].strip(),
                    filename=row["filename"].strip(),
                    primary_url=row["primary_url"].strip(),
                    mirror_url=row["mirror_url"].strip(),
                    sha256=row["sha256"].strip().lower(),
                    size_bytes=int(row["size_bytes"]),
                    license=row["license"].strip(),
                    description=row["description"].strip(),
                )
            except (KeyError, ValueError) as exc:
                raise SourceError(
                    f"Invalid source row {line_number} in {path}: {exc}"
                ) from exc

            _validate_source_record(record, path=path, line_number=line_number)
            if record.source_id in seen_ids:
                raise SourceError(f"Duplicate source_id {record.source_id!r} in {path}")
            if record.filename in seen_filenames:
                raise SourceError(f"Duplicate filename {record.filename!r} in {path}")
            seen_ids.add(record.source_id)
            seen_filenames.add(record.filename)
            records.append(record)
    return records


def _validate_source_record(record: SourceRecord, *, path: Path, line_number: int) -> None:
    location = f"{path}:{line_number}"
    if not record.source_id:
        raise SourceError(f"Blank source_id at {location}")
    if record.role not in VALID_ROLES:
        raise SourceError(f"Invalid role {record.role!r} at {location}")
    if record.acquisition not in VALID_ACQUISITION:
        raise SourceError(f"Invalid acquisition {record.acquisition!r} at {location}")
    if not record.filename or Path(record.filename).name != record.filename:
        raise SourceError(
            f"filename must be a plain basename, got {record.filename!r} at {location}"
        )
    if len(record.sha256) != 64 or any(
        character not in "0123456789abcdef" for character in record.sha256
    ):
        raise SourceError(f"Invalid SHA-256 at {location}")
    if record.size_bytes < 1:
        raise SourceError(f"Invalid size_bytes at {location}")
    if record.acquisition == "download" and not record.urls:
        raise SourceError(f"Download source has no URL at {location}")


def records_by_id(records: Iterable[SourceRecord]) -> dict[str, SourceRecord]:
    """Index source records by their unique IDs."""

    return {record.source_id: record for record in records}


def _looks_like_html(path: Path) -> bool:
    with path.open("rb") as handle:
        start = handle.read(4096).lstrip().lower()
    return (
        start.startswith(b"<!doctype html")
        or start.startswith(b"<html")
        or b"<title>login" in start
    )


def _verify_zip(path: Path) -> str:
    try:
        with zipfile.ZipFile(path) as archive:
            members = archive.infolist()
            if not members:
                raise SourceError(f"ZIP archive is empty: {path}")
            bad_member = archive.testzip()
            if bad_member is not None:
                raise SourceError(f"ZIP CRC failure in {path}: {bad_member}")
    except zipfile.BadZipFile as exc:
        raise SourceError(f"Invalid ZIP archive: {path}") from exc
    return f"zip_ok:{len(members)}_members"


def _verify_tar(path: Path) -> str:
    try:
        with tarfile.open(path, mode="r:*") as archive:
            members = archive.getmembers()
            if not members:
                raise SourceError(f"TAR archive is empty: {path}")
    except tarfile.TarError as exc:
        raise SourceError(f"Invalid TAR archive: {path}") from exc
    return f"tar_ok:{len(members)}_members"


def _verify_gzip(path: Path) -> str:
    uncompressed = 0
    try:
        with gzip.open(path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                uncompressed += len(block)
    except (OSError, EOFError) as exc:
        raise SourceError(f"Invalid gzip stream: {path}") from exc
    if uncompressed == 0:
        raise SourceError(f"Gzip stream is empty: {path}")
    return f"gzip_ok:{uncompressed}_uncompressed_bytes"


def _verify_text_format(path: Path, expected_filename: str) -> str:
    lower = expected_filename.lower()
    if lower.endswith(".gmt"):
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    if len(line.rstrip("\r\n").split("\t")) < 3:
                        raise SourceError(
                            f"First GMT row has fewer than three fields: {path}"
                        )
                    return "gmt_ok"
        raise SourceError(f"GMT file is empty: {path}")

    if lower.endswith(".obo"):
        with path.open(encoding="utf-8") as handle:
            head = "".join(handle.readline() for _ in range(20))
        if "format-version:" not in head:
            raise SourceError(f"OBO header not found near start of {path}")
        return "obo_ok"

    return "plain_file_ok"


def _verify_named_format(path: Path, expected_filename: str) -> str:
    """Check structure according to the committed filename, not temp suffixes."""

    lower = expected_filename.lower()
    if lower.endswith(".zip"):
        return _verify_zip(path)
    if lower.endswith((".tar.gz", ".tgz", ".tar")):
        return _verify_tar(path)
    if lower.endswith(".gz"):
        return _verify_gzip(path)
    return _verify_text_format(path, expected_filename)


def validate_source_file(path: Path, record: SourceRecord) -> FileVerification:
    """Verify byte identity and a lightweight structure or format check."""

    if not path.is_file():
        raise SourceError(f"Required source is missing: {path}")
    observed_size = path.stat().st_size
    if observed_size != record.size_bytes:
        raise SourceError(
            f"Size mismatch for {record.source_id}: expected {record.size_bytes}, "
            f"observed {observed_size} at {path}"
        )
    observed_hash = sha256_file(path)
    if observed_hash != record.sha256:
        raise SourceError(
            f"SHA-256 mismatch for {record.source_id}: expected {record.sha256}, "
            f"observed {observed_hash} at {path}"
        )
    if _looks_like_html(path):
        raise SourceError(f"Source looks like an HTML error or login page: {path}")

    structure = _verify_named_format(path, record.filename)
    return FileVerification(
        source_id=record.source_id,
        role=record.role,
        path=str(path.resolve()),
        size_bytes=observed_size,
        sha256=observed_hash,
        archive_or_format_check=structure,
        acquisition_status="verified_existing",
    )


def _download_to(url: str, temporary: Path, timeout_seconds: int) -> str:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "GraphSAGE-PPI-reproduction/0.1 (+public research workflow)"
        },
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        resolved_url = response.geturl()
        content_type = (response.headers.get("Content-Type") or "").lower()
        if "text/html" in content_type:
            raise SourceError(f"Server returned HTML for data URL {url}")
        with temporary.open("wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
            output.flush()
            os.fsync(output.fileno())
    return resolved_url


def acquire_source(
    record: SourceRecord,
    data_dir: Path,
    *,
    timeout_seconds: int = 120,
) -> FileVerification:
    """Verify an existing cache entry or atomically install a missing one.

    Existing paths are never replaced: a valid file is returned unchanged and
    an invalid file raises :class:`SourceError`. Only an absent destination can
    be populated by a verified download.
    """

    data_dir.mkdir(parents=True, exist_ok=True)
    destination = data_dir / record.filename
    if destination.exists():
        return validate_source_file(destination, record)
    if record.acquisition == "manual":
        raise SourceError(
            f"{record.source_id} must be supplied manually as {destination}; "
            "the expected size and SHA-256 are recorded in spec/sources.tsv"
        )

    errors: list[str] = []
    for attempt, url in enumerate(record.urls, start=1):
        temporary = data_dir / f".{record.filename}.{os.getpid()}.{attempt}.part"
        temporary.unlink(missing_ok=True)
        try:
            resolved_url = _download_to(url, temporary, timeout_seconds)
            verification = validate_source_file(temporary, record)
            os.replace(temporary, destination)
            return FileVerification(
                source_id=verification.source_id,
                role=verification.role,
                path=str(destination.resolve()),
                size_bytes=verification.size_bytes,
                sha256=verification.sha256,
                archive_or_format_check=verification.archive_or_format_check,
                acquisition_status="downloaded_and_verified",
                resolved_url=resolved_url,
            )
        except (OSError, urllib.error.URLError, SourceError) as exc:
            errors.append(f"{url}: {exc}")
        finally:
            temporary.unlink(missing_ok=True)

    raise SourceError(
        f"Could not acquire {record.source_id} ({record.filename}). Attempts:\n  "
        + "\n  ".join(errors)
    )


def ensure_sources(
    manifest_path: Path,
    data_dir: Path,
    *,
    source_ids: Sequence[str] | None = None,
    roles: Sequence[str] | None = None,
    report_path: Path | None = None,
    timeout_seconds: int = 120,
) -> list[FileVerification]:
    """Acquire and verify an explicit subset of the source inventory."""

    records = read_source_manifest(manifest_path)
    by_id = records_by_id(records)
    if source_ids:
        missing = sorted(set(source_ids) - set(by_id))
        if missing:
            raise SourceError(f"Unknown source IDs: {', '.join(missing)}")
        selected = [by_id[source_id] for source_id in source_ids]
    else:
        selected = records

    if roles:
        role_set = set(roles)
        invalid_roles = sorted(role_set - VALID_ROLES)
        if invalid_roles:
            raise SourceError(f"Unknown roles: {', '.join(invalid_roles)}")
        selected = [record for record in selected if record.role in role_set]

    verified = [
        acquire_source(record, data_dir, timeout_seconds=timeout_seconds)
        for record in selected
    ]
    if report_path is not None:
        write_json_atomic(
            report_path,
            {
                "schema_version": 1,
                "generated_at_utc": utc_now(),
                "source_manifest": str(manifest_path.resolve()),
                "source_manifest_sha256": sha256_file(manifest_path),
                "data_directory": str(data_dir.resolve()),
                "records": [asdict(item) for item in verified],
            },
        )
    return verified


def _read_json_object(path: Path) -> dict[str, Any]:
    """Read a JSON object and reject other top-level JSON values."""

    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceError(f"Could not read JSON object from {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SourceError(f"Expected a JSON object in {path}")
    return value


def _read_yaml_object(path: Path) -> dict[str, Any]:
    """Read a YAML mapping used as the frozen reconstruction specification."""

    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise SourceError(f"Could not read YAML object from {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SourceError(f"Expected a YAML mapping in {path}")
    return value


def _require_mapping(value: object, description: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise SourceError(f"Expected {description} to be a mapping")
    return value


def check_reconstruction_milestone(
    *,
    specification_path: Path,
    topology_summary_path: Path,
    feature_summary_path: Path,
    output_path: Path,
) -> dict[str, object]:
    """Check target-independent topology and feature invariants.

    This check reads only committed expectations and summaries generated from
    upstream sources.  It does not open the released GraphSAGE reference ZIP,
    so it can be part of the reconstruction dependency graph without crossing
    the project's non-circularity boundary.
    """

    specification = _read_yaml_object(specification_path)
    validation = _require_mapping(specification.get("validation"), "validation")
    expected = _require_mapping(validation.get("expected"), "validation.expected")
    topology = _read_json_object(topology_summary_path)
    features = _read_json_object(feature_summary_path)

    topology_counts = _require_mapping(topology.get("counts"), "topology counts")
    topology_hashes = _require_mapping(
        topology.get("content_hashes"), "topology content hashes"
    )
    feature_counts = _require_mapping(features.get("counts"), "feature counts")
    feature_hashes = _require_mapping(features.get("hashes"), "feature hashes")

    expected_topology_hashes = _require_mapping(
        expected.get("topology_content_hashes"),
        "validation.expected.topology_content_hashes",
    )
    expected_shape = expected.get("feature_shape")
    observed_shape = [feature_counts.get("rows"), feature_counts.get("columns")]

    checks: dict[str, bool] = {
        "graph_count": topology_counts.get("graphs") == expected.get("graphs"),
        "row_count": topology_counts.get("rows") == expected.get("rows"),
        "distinct_entrez_gene_ids": (
            topology_counts.get("distinct_entrez_gene_ids")
            == expected.get("distinct_entrez_gene_ids")
        ),
        "edge_record_count": (
            topology_counts.get("edge_records") == expected.get("graphsage_edge_records")
        ),
        "split_graph_counts": (
            topology_counts.get("split_graph_counts") == expected.get("split_graph_counts")
        ),
        "split_row_counts": (
            topology_counts.get("split_row_counts") == expected.get("split_row_counts")
        ),
        "topology_content_hashes": topology_hashes == expected_topology_hashes,
        "feature_shape": observed_shape == expected_shape,
        "feature_collection_counts": (
            feature_counts.get("selected_by_collection")
            == expected.get("feature_collection_counts")
        ),
        "feature_all_zero_columns": (
            feature_counts.get("all_zero_columns_0based")
            == expected.get("feature_all_zero_columns_0based")
        ),
        "feature_float64_data_hash": (
            feature_hashes.get("float64_c_order_data_sha256")
            == expected.get("feature_float64_c_order_data_sha256")
        ),
        "feature_uint8_data_hash": (
            feature_hashes.get("uint8_c_order_data_sha256")
            == expected.get("feature_uint8_c_order_data_sha256")
        ),
        "graphsage_feature_npy_hash": (
            feature_hashes.get("npy_file_sha256")
            == expected.get("graphsage_feature_npy_sha256")
        ),
    }

    result: dict[str, object] = {
        "schema_version": 1,
        "scope": "target-independent topology and feature invariants",
        "generated_at_utc": utc_now(),
        "inputs": {
            "specification": str(specification_path.resolve()),
            "specification_sha256": sha256_file(specification_path),
            "topology_summary": str(topology_summary_path.resolve()),
            "topology_summary_sha256": sha256_file(topology_summary_path),
            "feature_summary": str(feature_summary_path.resolve()),
            "feature_summary_sha256": sha256_file(feature_summary_path),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    write_json_atomic(output_path, result)
    if not result["all_checks_pass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise SourceError(
            "Target-independent reconstruction checks failed: " + ", ".join(failed)
        )
    return result


def _git_value(project_root: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(project_root), *arguments],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def _manifest_artifact(path: Path, project_root: Path) -> dict[str, object]:
    if not path.is_file():
        raise SourceError(f"Cannot manifest missing artifact: {path}")
    try:
        relative = path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError:
        relative = None
    return {
        "path": str(path.resolve()),
        "project_relative_path": relative,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def write_run_manifest(
    output_path: Path,
    *,
    project_root: Path,
    reproduction_root: Path,
    source_report: Path,
    artifact_paths: Sequence[Path],
    specification_paths: Sequence[Path] = (),
    stage: str,
) -> None:
    """Write a compact manifest for one completed reconstruction stage."""

    if not source_report.is_file():
        raise SourceError(f"Source-verification report is missing: {source_report}")
    source_value = json.loads(source_report.read_text(encoding="utf-8"))
    lock_path = reproduction_root / "pixi.lock"
    status = _git_value(project_root, "status", "--short")
    write_json_atomic(
        output_path,
        {
            "schema_version": 1,
            "generated_at_utc": utc_now(),
            "stage": stage,
            "git": {
                "commit": _git_value(project_root, "rev-parse", "HEAD"),
                "branch": _git_value(project_root, "branch", "--show-current"),
                "dirty": bool(status),
                "status_short": status,
            },
            "environment": {
                "python": platform.python_version(),
                "implementation": platform.python_implementation(),
                "platform": platform.platform(),
                "pixi_lock_sha256": (
                    sha256_file(lock_path) if lock_path.is_file() else None
                ),
            },
            "source_verification": source_value,
            "specifications": [
                _manifest_artifact(path, project_root)
                for path in sorted(specification_paths, key=lambda item: str(item))
            ],
            "artifacts": [
                _manifest_artifact(path, project_root)
                for path in sorted(artifact_paths, key=lambda item: str(item))
            ],
        },
    )


def write_tsv_atomic(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[dict[str, object]],
) -> None:
    """Write an optionally gzip-compressed TSV with deterministic newlines.

    For ``.gz`` outputs the gzip member timestamp and embedded filename are
    fixed, so identical rows yield identical bytes across runs.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(file_descriptor)
    temporary = Path(temporary_name)
    try:
        if path.suffix == ".gz":
            with temporary.open("wb") as raw:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw, mtime=0
                ) as compressed:
                    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
                    writer = csv.DictWriter(
                        text,
                        fieldnames=list(fieldnames),
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(rows)
                    text.flush()
                    text.detach()
        else:
            with temporary.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=list(fieldnames),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)

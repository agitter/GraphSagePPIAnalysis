#!/usr/bin/env python3
"""Validate MSigDB v6.1 C1/C3 as a canonical GraphSAGE PPI feature source."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def membership_sha256(members: set[int]) -> str:
    payload = "\n".join(str(value) for value in sorted(members)) + "\n"
    return sha256_bytes(payload.encode("ascii"))


def membership_sequence_sha256(records: list[dict[str, Any]]) -> str:
    payload = "\n".join(record["membership_sha256"] for record in records) + "\n"
    return sha256_bytes(payload.encode("ascii"))


def detect_newline_facts(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    return {
        "ends_with_newline": data.endswith(b"\n"),
        "contains_crlf": b"\r\n" in data,
        "contains_bare_cr": b"\r" in data.replace(b"\r\n", b""),
    }


def read_gmt(path: Path, collection: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for source_row_0based, raw_line in enumerate(handle):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{source_row_0based + 1}: fewer than three GMT fields")
            name, description = fields[:2]
            raw_members = fields[2:]
            members: list[int] = []
            invalid_tokens: list[str] = []
            for token in raw_members:
                try:
                    members.append(int(token))
                except ValueError:
                    invalid_tokens.append(token)
            if invalid_tokens:
                raise ValueError(
                    f"{path}:{source_row_0based + 1}: non-integer member tokens: {invalid_tokens[:5]}"
                )
            member_set = set(members)
            rows.append(
                {
                    "collection": collection,
                    "source_filename": path.name,
                    "source_row_index_0based": source_row_0based,
                    "source_row_index_1based": source_row_0based + 1,
                    "name": name,
                    "description": description,
                    "raw_member_count": len(members),
                    "unique_entrez_count": len(member_set),
                    "duplicate_member_fields": len(members) - len(member_set),
                    "members": member_set,
                    "membership_sha256": membership_sha256(member_set),
                }
            )
    return rows


def select_features(
    collection_rows: dict[str, list[dict[str, Any]]],
    *,
    minimum: int,
    inclusive: bool,
    cap: int,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    selected: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    predicate = (lambda n: n >= minimum) if inclusive else (lambda n: n > minimum)
    for collection in ("C1", "C3"):
        qualifying = [row for row in collection_rows[collection] if predicate(row["unique_entrez_count"])]
        counts[collection] = len(qualifying)
        for row in qualifying:
            if len(selected) >= cap:
                break
            selected.append(row)
        if len(selected) >= cap:
            break
    return selected, counts


def read_row_map(path: Path) -> np.ndarray:
    opener = gzip.open if path.suffix == ".gz" else open
    row_to_gene: dict[int, int] = {}
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        for record in csv.DictReader(handle):
            row = int(record["graphsage_row"])
            gene = int(record["entrez_gene_id"])
            if row in row_to_gene:
                raise ValueError(f"duplicate GraphSAGE row {row}")
            row_to_gene[row] = gene
    expected_rows = list(range(len(row_to_gene)))
    if sorted(row_to_gene) != expected_rows:
        raise ValueError("row map is not complete and contiguous from zero")
    return np.asarray([row_to_gene[row] for row in expected_rows], dtype=np.int64)


def read_graphsage_features(path: Path) -> tuple[np.ndarray, bytes]:
    with zipfile.ZipFile(path) as archive:
        member_bytes = archive.read("ppi/ppi-feats.npy")
    matrix = np.load(io.BytesIO(member_bytes), allow_pickle=False)
    return matrix, member_bytes


def read_prior_details(path: Path) -> list[dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    details = payload["details"]["6.0"]
    details.sort(key=lambda record: int(record["column"]))
    return details


def write_csv(path: Path, records: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c1", type=Path, required=True)
    parser.add_argument("--c3", type=Path, required=True)
    parser.add_argument("--graphsage-zip", type=Path, required=True)
    parser.add_argument("--row-map", type=Path, required=True)
    parser.add_argument("--prior-v60-details", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    source_paths = {"C1": args.c1, "C3": args.c3}
    rows_by_collection = {
        collection: read_gmt(path, collection) for collection, path in source_paths.items()
    }

    selected_ge, qualifying_ge = select_features(
        rows_by_collection, minimum=200, inclusive=True, cap=50
    )
    selected_gt, qualifying_gt = select_features(
        rows_by_collection, minimum=200, inclusive=False, cap=50
    )

    if len(selected_ge) != 50 or len(selected_gt) != 50:
        raise SystemExit("v6.1 did not yield 50 selected features under both threshold operators")

    sequence_ge = [(row["collection"], row["membership_sha256"]) for row in selected_ge]
    sequence_gt = [(row["collection"], row["membership_sha256"]) for row in selected_gt]
    threshold_selected_sequences_identical = sequence_ge == sequence_gt

    prior_v60 = read_prior_details(args.prior_v60_details)
    if len(prior_v60) != 50:
        raise ValueError("prior v6.0 details do not contain 50 features")

    genes = read_row_map(args.row_map)
    observed, observed_member_bytes = read_graphsage_features(args.graphsage_zip)
    if observed.shape != (len(genes), 50):
        raise ValueError(f"unexpected GraphSAGE feature shape: {observed.shape}")
    if not np.all((observed == 0) | (observed == 1)):
        raise ValueError("GraphSAGE feature matrix is not binary")

    expected = np.zeros(observed.shape, dtype=observed.dtype)
    selected_records: list[dict[str, Any]] = []
    changed_name_records: list[dict[str, Any]] = []
    column_mismatch_counts: list[int] = []

    for column, (row, old) in enumerate(zip(selected_ge, prior_v60, strict=True)):
        expected[:, column] = np.fromiter(
            (int(gene) in row["members"] for gene in genes),
            dtype=observed.dtype,
            count=len(genes),
        )
        mismatch_count = int(np.count_nonzero(expected[:, column] != observed[:, column]))
        column_mismatch_counts.append(mismatch_count)
        membership_matches_v60 = row["membership_sha256"] == old["membership_sha256"]
        collection_matches_v60 = row["collection"] == old["collection"]
        old_name = old["name"]
        name_changed = old_name != row["name"]
        record = {
            "graphsage_feature_column_0based": column,
            "graphsage_feature_column_1based": column + 1,
            "collection": row["collection"],
            "source_filename": row["source_filename"],
            "source_row_index_0based": row["source_row_index_0based"],
            "source_row_index_1based": row["source_row_index_1based"],
            "v6_1_name": row["name"],
            "v6_0_name": old_name,
            "name_changed_from_v6_0": int(name_changed),
            "description": row["description"],
            "unique_entrez_count": row["unique_entrez_count"],
            "raw_member_count": row["raw_member_count"],
            "duplicate_member_fields": row["duplicate_member_fields"],
            "membership_sha256": row["membership_sha256"],
            "membership_matches_v6_0": int(membership_matches_v60),
            "collection_matches_v6_0": int(collection_matches_v60),
            "observed_positive_rows": int(np.count_nonzero(observed[:, column])),
            "expected_positive_rows": int(np.count_nonzero(expected[:, column])),
            "mismatch_rows": mismatch_count,
            "exact_graphsage_column": int(mismatch_count == 0),
            "all_zero_graphsage_column": int(np.count_nonzero(observed[:, column]) == 0),
        }
        selected_records.append(record)
        if name_changed:
            changed_name_records.append(record)

    total_mismatches = int(np.count_nonzero(expected != observed))
    observed_raw_matrix_sha = sha256_bytes(np.ascontiguousarray(observed).tobytes())
    expected_raw_matrix_sha = sha256_bytes(np.ascontiguousarray(expected).tobytes())

    gene_to_first_row: dict[int, int] = {}
    repeated_gene_conflicts: list[dict[str, int]] = []
    for row_index, gene in enumerate(genes.tolist()):
        if gene in gene_to_first_row:
            first = gene_to_first_row[gene]
            if not np.array_equal(observed[first], observed[row_index]):
                repeated_gene_conflicts.append(
                    {"entrez_gene_id": int(gene), "first_row": first, "conflicting_row": row_index}
                )
        else:
            gene_to_first_row[gene] = row_index

    selected_sequence_hash = membership_sequence_sha256(selected_ge)
    prior_sequence_hash = membership_sequence_sha256(
        [
            {"membership_sha256": record["membership_sha256"]}
            for record in prior_v60
        ]
    )

    exact_boundary_rows = {
        collection: [
            {
                "source_row_index_1based": row["source_row_index_1based"],
                "name": row["name"],
            }
            for row in rows
            if row["unique_entrez_count"] == 200
        ]
        for collection, rows in rows_by_collection.items()
    }

    source_manifest = []
    for source_id, path in (
        ("msigdb_c1_v6_1", args.c1),
        ("msigdb_c3_v6_1", args.c3),
        ("graphsage_ppi_reference", args.graphsage_zip),
        ("full_row_map", args.row_map),
        ("prior_v6_0_feature_details", args.prior_v60_details),
    ):
        source_manifest.append(
            {
                "source_id": source_id,
                "filename": path.name,
                "size_bytes": path.stat().st_size,
                "sha256": sha256_path(path),
            }
        )

    summary = {
        "result": "PASS" if total_mismatches == 0 else "FAIL",
        "canonical_source_recommendation": (
            "MSigDB v6.1 C1/C3 is acceptable as the canonical reproduction source"
            if total_mismatches == 0
            and all(r["membership_matches_v6_0"] for r in selected_records)
            and threshold_selected_sequences_identical
            else "MSigDB v6.1 should not be adopted without resolving failed checks"
        ),
        "inputs": {record["source_id"]: record for record in source_manifest},
        "gmt": {
            collection: {
                "row_count": len(rows_by_collection[collection]),
                "unique_name_count": len({row["name"] for row in rows_by_collection[collection]}),
                "duplicate_name_count": len(rows_by_collection[collection])
                - len({row["name"] for row in rows_by_collection[collection]}),
                "newline": detect_newline_facts(source_paths[collection]),
            }
            for collection in ("C1", "C3")
        },
        "selection_rule": {
            "collection_order": ["C1", "C3"],
            "minimum_source_unique_entrez_members": 200,
            "global_cap": 50,
            "preserve_source_row_order": True,
            "greater_than_or_equal_qualifying_counts": qualifying_ge,
            "strictly_greater_qualifying_counts": qualifying_gt,
            "sets_with_exactly_200_unique_members": exact_boundary_rows,
            "selected_sequences_identical_for_ge_200_and_gt_200": threshold_selected_sequences_identical,
            "selected_counts": dict(Counter(row["collection"] for row in selected_ge)),
        },
        "cross_version": {
            "prior_v6_0_selected_membership_sequence_sha256": prior_sequence_hash,
            "v6_1_selected_membership_sequence_sha256": selected_sequence_hash,
            "expected_previously_observed_sequence_sha256": "41c5d821e1b706ec4c8dceb47ab25c5dbad689998483e34e9e41d08094448101",
            "all_50_memberships_match_v6_0_by_column": all(
                r["membership_matches_v6_0"] for r in selected_records
            ),
            "all_50_collections_match_v6_0_by_column": all(
                r["collection_matches_v6_0"] for r in selected_records
            ),
            "changed_selected_names_count": len(changed_name_records),
            "changed_selected_names": [
                {
                    "column_0based": r["graphsage_feature_column_0based"],
                    "column_1based": r["graphsage_feature_column_1based"],
                    "v6_0_name": r["v6_0_name"],
                    "v6_1_name": r["v6_1_name"],
                    "membership_sha256": r["membership_sha256"],
                }
                for r in changed_name_records
            ],
        },
        "graphsage_comparison": {
            "matrix_shape": list(observed.shape),
            "matrix_dtype": str(observed.dtype),
            "cells_compared": int(observed.size),
            "total_mismatches": total_mismatches,
            "exact_columns": int(sum(count == 0 for count in column_mismatch_counts)),
            "all_zero_columns_0based": [
                int(column)
                for column in range(observed.shape[1])
                if np.count_nonzero(observed[:, column]) == 0
            ],
            "observed_matrix_raw_c_order_sha256": observed_raw_matrix_sha,
            "reconstructed_matrix_raw_c_order_sha256": expected_raw_matrix_sha,
            "ppi_feats_npy_member_sha256": sha256_bytes(observed_member_bytes),
            "row_count": len(genes),
            "distinct_entrez_gene_count": len(set(int(gene) for gene in genes)),
            "repeated_gene_feature_conflicts": repeated_gene_conflicts,
        },
        "acceptance_checks": {
            "exact_30_C1_and_20_C3": Counter(row["collection"] for row in selected_ge)
            == Counter({"C1": 30, "C3": 20}),
            "same_selected_sequence_under_ge_200_and_gt_200": threshold_selected_sequences_identical,
            "all_50_memberships_match_v6_0": all(
                r["membership_matches_v6_0"] for r in selected_records
            ),
            "membership_sequence_hash_matches_prior": selected_sequence_hash
            == "41c5d821e1b706ec4c8dceb47ab25c5dbad689998483e34e9e41d08094448101",
            "all_2847200_graphsage_cells_exact": total_mismatches == 0,
            "all_50_graphsage_columns_exact": all(count == 0 for count in column_mismatch_counts),
            "no_repeated_gene_feature_conflicts": not repeated_gene_conflicts,
        },
    }

    summary_path = args.output_dir / "validation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    write_csv(
        args.output_dir / "selected_features_v6.1.csv",
        selected_records,
        [
            "graphsage_feature_column_0based",
            "graphsage_feature_column_1based",
            "collection",
            "source_filename",
            "source_row_index_0based",
            "source_row_index_1based",
            "v6_1_name",
            "v6_0_name",
            "name_changed_from_v6_0",
            "description",
            "unique_entrez_count",
            "raw_member_count",
            "duplicate_member_fields",
            "membership_sha256",
            "membership_matches_v6_0",
            "collection_matches_v6_0",
            "observed_positive_rows",
            "expected_positive_rows",
            "mismatch_rows",
            "exact_graphsage_column",
            "all_zero_graphsage_column",
        ],
    )
    write_csv(
        args.output_dir / "v6.0_to_v6.1_selected_name_changes.csv",
        changed_name_records,
        [
            "graphsage_feature_column_0based",
            "graphsage_feature_column_1based",
            "collection",
            "v6_0_name",
            "v6_1_name",
            "membership_sha256",
            "unique_entrez_count",
            "source_row_index_1based",
        ],
    )
    write_csv(
        args.output_dir / "source_manifest.csv",
        source_manifest,
        ["source_id", "filename", "size_bytes", "sha256"],
    )

    failed = [name for name, passed in summary["acceptance_checks"].items() if not passed]
    print(json.dumps({"result": summary["result"], "failed_checks": failed, **summary["cross_version"], **summary["graphsage_comparison"]}, indent=2, sort_keys=True))
    if failed:
        raise SystemExit("Validation failed: " + ", ".join(failed))


if __name__ == "__main__":
    main()

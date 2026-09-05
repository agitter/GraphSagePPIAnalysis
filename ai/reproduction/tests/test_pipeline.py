"""Small end-to-end and provenance tests without external data downloads."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import tarfile
import zipfile
from pathlib import Path

import pytest
import yaml

from graphsage_ppi_repro.features import reconstruct_features
from graphsage_ppi_repro.legacy_order import ordered_string_keys
from graphsage_ppi_repro.provenance import (
    SourceError,
    check_reconstruction_milestone,
    ensure_sources,
)
from graphsage_ppi_repro.topology import reconstruct_topology
from graphsage_ppi_repro.validate import validate_topology_and_features


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_existing_source_is_verified_without_downloading(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    source = data_dir / "example.txt"
    source.write_text("verified\n", encoding="utf-8", newline="\n")
    manifest = tmp_path / "sources.tsv"
    manifest.write_text(
        "source_id\trole\tacquisition\tfilename\tprimary_url\tmirror_url\tsha256\t"
        "size_bytes\tlicense\tdescription\n"
        f"example\tupstream\tmanual\texample.txt\t\t\t{_sha256(source)}\t"
        f"{source.stat().st_size}\ttest\tSynthetic source\n",
        encoding="utf-8",
        newline="\n",
    )
    report = tmp_path / "report.json"
    verified = ensure_sources(
        manifest,
        data_dir,
        source_ids=["example"],
        report_path=report,
    )
    assert verified[0].acquisition_status == "verified_existing"
    assert json.loads(report.read_text())["records"][0]["sha256"] == _sha256(source)


def test_source_hash_mismatch_is_rejected(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    source = data_dir / "example.txt"
    source.write_text("changed\n", encoding="utf-8")
    manifest = tmp_path / "sources.tsv"
    manifest.write_text(
        "source_id\trole\tacquisition\tfilename\tprimary_url\tmirror_url\tsha256\t"
        "size_bytes\tlicense\tdescription\n"
        f"example\tupstream\tmanual\texample.txt\t\t\t{'0' * 64}\t"
        f"{source.stat().st_size}\ttest\tSynthetic source\n",
        encoding="utf-8",
    )
    with pytest.raises(SourceError, match="SHA-256 mismatch"):
        ensure_sources(manifest, data_dir, source_ids=["example"])


def _write_graph_inputs(tmp_path: Path) -> tuple[Path, Path, bytes]:
    member = "bio-tissue-networks/tiny.edgelist"
    payload = b"1 2\n2 3\n3 3\n"
    archive_path = tmp_path / "networks.tar.gz"
    with tarfile.open(archive_path, "w:gz") as archive:
        info = tarfile.TarInfo(member)
        info.size = len(payload)
        info.mtime = 0
        archive.addfile(info, io.BytesIO(payload))
    spec_path = tmp_path / "selected_graphs.tsv"
    spec_path.write_text(
        "graph_index_1based\ttissue\tsplit\tohmnet_member\tohmnet_member_sha256\t"
        "node_count\tedge_count\tevidence_status\tevidence_reference\n"
        f"1\ttiny\ttrain\t{member}\t{hashlib.sha256(payload).hexdigest()}\t"
        "3\t3\tsynthetic\ttest\n",
        encoding="utf-8",
    )
    return archive_path, spec_path, payload


def _write_feature_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    c1 = tmp_path / "c1.gmt"
    c3 = tmp_path / "c3.gmt"
    c1.write_text("set_a\tdesc\t1\t2\n", encoding="utf-8")
    c3.write_text("set_b\tdesc\t2\t3\n", encoding="utf-8")
    fields = [
        "column_index_0based",
        "collection",
        "source_file",
        "source_row_1based",
        "canonical_name",
        "historical_alias",
        "membership_sha256",
        "source_unique_entrez_count",
        "expected_positive_rows",
        "all_zero_after_projection",
        "evidence_status",
        "evidence_reference",
    ]
    spec = tmp_path / "feature_columns.tsv"
    with spec.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for index, (collection, filename, name, members) in enumerate(
            [
                ("C1", "c1.gmt", "set_a", {1, 2}),
                ("C3", "c3.gmt", "set_b", {2, 3}),
            ]
        ):
            membership = hashlib.sha256(
                "".join(f"{member}\n" for member in sorted(members)).encode("ascii")
            ).hexdigest()
            writer.writerow(
                {
                    "column_index_0based": index,
                    "collection": collection,
                    "source_file": filename,
                    "source_row_1based": 1,
                    "canonical_name": name,
                    "historical_alias": "",
                    "membership_sha256": membership,
                    "source_unique_entrez_count": 2,
                    "expected_positive_rows": 2,
                    "all_zero_after_projection": 0,
                    "evidence_status": "synthetic",
                    "evidence_reference": "test",
                }
            )
    return c1, c3, spec


def _read_mapping(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_edges(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_reconstruct_without_target_then_validate_independently(tmp_path: Path) -> None:
    archive, graph_spec, _ = _write_graph_inputs(tmp_path)
    mapping = tmp_path / "mapping.tsv.gz"
    edges = tmp_path / "edges.tsv.gz"
    topology_summary_path = tmp_path / "topology.json"

    # At this point no released target exists. Reconstruction nevertheless completes.
    missing_reference = tmp_path / "reference.zip"
    assert not missing_reference.exists()
    topology_summary = reconstruct_topology(
        archive_path=archive,
        graph_spec_path=graph_spec,
        mapping_output=mapping,
        edge_output=edges,
        summary_output=topology_summary_path,
    )

    c1, c3, feature_spec = _write_feature_inputs(tmp_path)
    features = tmp_path / "features.npy"
    selected = tmp_path / "selected.tsv"
    feature_summary_path = tmp_path / "features.json"
    feature_summary = reconstruct_features(
        mapping_path=mapping,
        c1_path=c1,
        c3_path=c3,
        feature_spec_path=feature_spec,
        matrix_output=features,
        selected_output=selected,
        summary_output=feature_summary_path,
        minimum_source_members=2,
        maximum_columns=2,
    )

    specification = tmp_path / "specification.yaml"
    specification.write_text(
        yaml.safe_dump(
            {
                "validation": {
                    "expected": {
                        "graphs": 1,
                        "rows": 3,
                        "distinct_entrez_gene_ids": 3,
                        "graphsage_edge_records": 3,
                        "split_graph_counts": {"train": 1},
                        "split_row_counts": {"train": 3},
                        "topology_content_hashes": topology_summary["content_hashes"],
                        "feature_shape": [3, 2],
                        "feature_collection_counts": {"C1": 1, "C3": 1},
                        "feature_all_zero_columns_0based": [],
                        "feature_float64_c_order_data_sha256": feature_summary["hashes"][
                            "float64_c_order_data_sha256"
                        ],
                        "feature_uint8_c_order_data_sha256": feature_summary["hashes"][
                            "uint8_c_order_data_sha256"
                        ],
                        "graphsage_feature_npy_sha256": feature_summary["hashes"][
                            "npy_file_sha256"
                        ],
                    }
                }
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    checks = check_reconstruction_milestone(
        specification_path=specification,
        topology_summary_path=topology_summary_path,
        feature_summary_path=feature_summary_path,
        output_path=tmp_path / "checks.json",
    )
    assert checks["all_checks_pass"]

    mapping_rows = _read_mapping(mapping)
    edge_rows = _read_edges(edges)
    graph_nodes = [
        {"id": int(row["graphsage_node_id"]), "val": False, "test": False}
        for row in mapping_rows
    ]
    graph_links = [
        {"source": int(row["source_node_id"]), "target": int(row["target_node_id"])}
        for row in edge_rows
    ]
    id_map = {str(index): index for index in range(len(mapping_rows))}
    with zipfile.ZipFile(missing_reference, "w") as archive_file:
        archive_file.writestr(
            "ppi/ppi-G.json",
            json.dumps({"nodes": graph_nodes, "links": graph_links}),
        )
        archive_file.writestr("ppi/ppi-id_map.json", json.dumps(id_map))
        archive_file.writestr("ppi/ppi-feats.npy", features.read_bytes())

    validation = validate_topology_and_features(
        graphsage_reference_zip=missing_reference,
        mapping_path=mapping,
        edge_path=edges,
        feature_matrix_path=features,
        output_json=tmp_path / "validation.json",
        output_markdown=tmp_path / "validation.md",
    )
    assert validation["all_checks_pass"]
    assert [row["entrez_gene_id"] for row in mapping_rows] == ordered_string_keys(
        ["1", "2", "2", "3", "3", "3"]
    )

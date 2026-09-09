"""Small end-to-end and provenance tests without external data downloads."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
import tomllib
import zipfile
from pathlib import Path

import numpy as np
import pytest
import yaml

from graphsage_ppi_repro import provenance
from graphsage_ppi_repro.cli import main as cli_main
from graphsage_ppi_repro.features import reconstruct_features
from graphsage_ppi_repro.graphsage import assemble_graphsage
from graphsage_ppi_repro.legacy_order import ordered_string_keys
from graphsage_ppi_repro.provenance import (
    SourceError,
    check_reconstruction,
    ensure_sources,
)
from graphsage_ppi_repro.topology import reconstruct_topology
from graphsage_ppi_repro.validate import validate_graphsage_dataset


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_existing_source_is_verified_without_modification(tmp_path: Path) -> None:
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
    original_bytes = source.read_bytes()
    original_mtime_ns = source.stat().st_mtime_ns

    report = tmp_path / "report.json"
    verified = ensure_sources(
        manifest,
        data_dir,
        source_ids=["example"],
        report_path=report,
    )

    assert verified[0].acquisition_status == "verified_existing"
    assert json.loads(report.read_text())["records"][0]["sha256"] == _sha256(source)
    assert source.read_bytes() == original_bytes
    assert source.stat().st_mtime_ns == original_mtime_ns


def test_source_hash_mismatch_is_rejected_without_replacement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    source = data_dir / "example.txt"
    source.write_text("changed\n", encoding="utf-8")
    manifest = tmp_path / "sources.tsv"
    manifest.write_text(
        "source_id\trole\tacquisition\tfilename\tprimary_url\tmirror_url\tsha256\t"
        "size_bytes\tlicense\tdescription\n"
        f"example\tupstream\tdownload\texample.txt\thttps://example.test/file\t\t"
        f"{'0' * 64}\t{source.stat().st_size}\ttest\tSynthetic source\n",
        encoding="utf-8",
    )
    original_bytes = source.read_bytes()
    original_mtime_ns = source.stat().st_mtime_ns

    def unexpected_download(url: str, temporary: Path, timeout_seconds: int) -> str:
        raise AssertionError(
            f"download attempted for existing invalid source: {url}, "
            f"{temporary}, {timeout_seconds}"
        )

    monkeypatch.setattr(provenance, "_download_to", unexpected_download)
    with pytest.raises(SourceError, match="SHA-256 mismatch"):
        ensure_sources(manifest, data_dir, source_ids=["example"])

    assert source.read_bytes() == original_bytes
    assert source.stat().st_mtime_ns == original_mtime_ns


def test_missing_source_is_downloaded_through_a_temporary_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_dir = tmp_path / "data"
    payload = b"downloaded\n"
    expected_hash = hashlib.sha256(payload).hexdigest()
    manifest = tmp_path / "sources.tsv"
    manifest.write_text(
        "source_id\trole\tacquisition\tfilename\tprimary_url\tmirror_url\tsha256\t"
        "size_bytes\tlicense\tdescription\n"
        "example\tupstream\tdownload\texample.txt\thttps://example.test/file\t\t"
        f"{expected_hash}\t{len(payload)}\ttest\tSynthetic source\n",
        encoding="utf-8",
        newline="\n",
    )

    def fake_download(url: str, temporary: Path, timeout_seconds: int) -> str:
        assert url == "https://example.test/file"
        assert timeout_seconds == 120
        assert temporary.name.endswith(".part")
        assert not (data_dir / "example.txt").exists()
        temporary.parent.mkdir(parents=True, exist_ok=True)
        temporary.write_bytes(payload)
        return url

    monkeypatch.setattr(provenance, "_download_to", fake_download)
    verified = ensure_sources(manifest, data_dir, source_ids=["example"])

    assert verified[0].acquisition_status == "downloaded_and_verified"
    assert (data_dir / "example.txt").read_bytes() == payload
    assert not list(data_dir.glob("*.part"))


def test_clean_refuses_to_remove_the_source_cache(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    marker = data_dir / "source.txt"
    marker.write_text("keep\n", encoding="utf-8")

    assert cli_main(["clean", "--directory", str(data_dir)]) == 2
    assert marker.read_text(encoding="utf-8") == "keep\n"


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

    # No released target exists while the complete synthetic dataset is built.
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

    labels = tmp_path / "labels.npy"
    label_array = np.asarray([[1, 0], [0, 1], [1, 1]], dtype=np.uint8)
    np.save(labels, label_array, allow_pickle=False)
    label_summary_path = tmp_path / "labels.json"
    label_summary = {
        "counts": {
            "rows": 3,
            "columns": 2,
            "positive_cells": 4,
            "selected_by_namespace": {"biological_process": 2},
            "distinct_column_membership_vectors": 2,
            "repeated_gene_vector_conflicts": 0,
            "unmapped_graph_gene_ids": [],
        },
        "hashes": {
            "uint8_c_order_data_sha256": hashlib.sha256(
                label_array.tobytes(order="C")
            ).hexdigest()
        },
        "term_selection": {
            "derived_set_matches_column_specification": True,
        },
    }
    label_summary_path.write_text(json.dumps(label_summary), encoding="utf-8", newline="\n")

    graphsage_root = tmp_path / "graphsage"
    graphsage_summary_path = tmp_path / "graphsage-summary.json"
    graphsage_summary = assemble_graphsage(
        mapping_path=mapping,
        edge_path=edges,
        feature_matrix_path=features,
        label_matrix_path=labels,
        output_directory=graphsage_root / "ppi",
        summary_output=graphsage_summary_path,
    )
    assert not missing_reference.exists()

    artifact_hashes = {
        name: metadata["sha256"]
        for name, metadata in graphsage_summary["outputs"]["files"].items()
    }
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
                        "label_shape": [3, 2],
                        "label_positive_cells": 4,
                        "label_namespace_counts": {"biological_process": 2},
                        "label_distinct_membership_vectors": 2,
                        "label_unmapped_graph_gene_ids": [],
                        "label_uint8_c_order_data_sha256": label_summary["hashes"][
                            "uint8_c_order_data_sha256"
                        ],
                        "graphsage_content_hashes": graphsage_summary["content_hashes"],
                        "graphsage_artifact_hashes": artifact_hashes,
                        "graphsage_checksums_sha256": graphsage_summary["outputs"][
                            "checksums_sha256"
                        ],
                    }
                }
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    checks = check_reconstruction(
        specification_path=specification,
        topology_summary_path=topology_summary_path,
        feature_summary_path=feature_summary_path,
        label_summary_path=label_summary_path,
        graphsage_summary_path=graphsage_summary_path,
        output_path=tmp_path / "checks.json",
    )
    assert checks["all_checks_pass"]

    graph_path = graphsage_root / "ppi" / "ppi-G.json"
    id_map_path = graphsage_root / "ppi" / "ppi-id_map.json"
    class_map_path = graphsage_root / "ppi" / "ppi-class_map.json"
    feature_path = graphsage_root / "ppi" / "ppi-feats.npy"
    reference_graph = json.loads(graph_path.read_text())
    reference_graph["links"] = list(reversed(reference_graph["links"]))
    with zipfile.ZipFile(missing_reference, "w") as archive_file:
        archive_file.writestr("ppi/ppi-G.json", json.dumps(reference_graph))
        archive_file.writestr("ppi/ppi-id_map.json", id_map_path.read_bytes())
        archive_file.writestr("ppi/ppi-class_map.json", class_map_path.read_bytes())
        archive_file.writestr("ppi/ppi-feats.npy", feature_path.read_bytes())

    validation = validate_graphsage_dataset(
        reference_archive=missing_reference,
        graph_path=graph_path,
        id_map_path=id_map_path,
        class_map_path=class_map_path,
        feature_path=feature_path,
        json_output=tmp_path / "validation.json",
        markdown_output=tmp_path / "validation.md",
    )
    assert validation["all_checks_pass"]
    assert not validation["byte_equality"]["ppi-G.json"]
    assert [row["entrez_gene_id"] for row in _read_mapping(mapping)] == (
        ordered_string_keys(["1", "2", "2", "3", "3", "3"])
    )


def test_workflow_treats_cached_sources_as_inputs_only() -> None:
    reproduction_root = Path(__file__).resolve().parents[1]
    snakefile_text = (reproduction_root / "Snakefile").read_text(encoding="utf-8")
    output_blocks = re.findall(r"(?m)^    output:\n((?:        .*\n)+)", snakefile_text)

    assert "rule acquire_upstream_sources:" not in snakefile_text
    assert "rule acquire_graphsage_reference:" not in snakefile_text
    assert "rule acquire_dgl_reference:" not in snakefile_text
    assert "rule verify_upstream_sources:" in snakefile_text
    assert "rule verify_graphsage_reference:" in snakefile_text
    assert "rule verify_dgl_reference:" in snakefile_text
    assert output_blocks
    for block in output_blocks:
        assert "OHMNET" not in block
        assert "OHMNET_README" not in block
        assert "MSIGDB_C1" not in block
        assert "MSIGDB_C3" not in block
        assert "GOA_GAF" not in block
        assert "GOA_GPI" not in block
        assert "GP2PROTEIN" not in block
        assert "GO_ONTOLOGY" not in block
        assert "GRAPHSAGE_REFERENCE" not in block
        assert re.search(r"\bDGL_REFERENCE\b", block) is None

    with (reproduction_root / "pixi.toml").open("rb") as handle:
        tasks = tomllib.load(handle)["tasks"]
    assert tasks["reproduce"]["depends-on"] == ["acquire-upstream"]
    assert tasks["validate"]["depends-on"] == [
        "acquire-upstream",
        "acquire-graphsage-reference",
        "acquire-dgl-reference",
    ]
    upstream_command = tasks["acquire-upstream"]
    for source_id in (
        "ohmnet_networks",
        "ohmnet_readme",
        "msigdb_c1_v61",
        "msigdb_c3_v61",
        "goa_human_gaf_159",
        "goa_human_gpi_159",
        "gp2protein_geneid_20160601",
        "go_ontology_20160601",
    ):
        assert f"--source-id {source_id}" in upstream_command
    assert "rule reconstruct_labels:" in snakefile_text
    assert "rule assemble_graphsage:" in snakefile_text
    assert "rule validate_graphsage:" in snakefile_text

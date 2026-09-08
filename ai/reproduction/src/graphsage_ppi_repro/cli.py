"""Thin command-line adapters used by Pixi and the Snakemake workflow.

The CLI contains no scientific transformations.  Each subcommand validates
arguments and delegates to an importable, unit-tested function in the module
named for that stage.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from .features import FeatureError, reconstruct_features
from .graphsage import GraphSAGEError, assemble_graphsage
from .identifiers import IdentifierError
from .labels import LabelError, reconstruct_labels
from .provenance import (
    SourceError,
    check_reconstruction_milestone,
    ensure_sources,
    write_run_manifest,
)
from .topology import TopologyError, reconstruct_topology
from .validate import (
    ValidationError,
    validate_graphsage_dataset,
    validate_labels_against_graphsage,
    validate_topology_and_features,
)


def _path(value: str) -> Path:
    return Path(value).expanduser()


def _add_path_argument(parser: argparse.ArgumentParser, name: str, help_text: str) -> None:
    parser.add_argument(name, type=_path, help=help_text)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="graphsage-ppi-repro",
        description="Reconstruct and validate the GraphSAGE PPI data workflow.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sources = subparsers.add_parser(
        "ensure-sources", help="Acquire or verify checksum-locked source files."
    )
    sources.add_argument("--manifest", type=_path, required=True)
    sources.add_argument("--data-dir", type=_path, required=True)
    sources.add_argument(
        "--source-id",
        action="append",
        dest="source_ids",
        help="Source ID to acquire; repeat for multiple records.",
    )
    sources.add_argument(
        "--role",
        action="append",
        dest="roles",
        choices=("upstream", "reference"),
        help="Optional role filter; repeat to allow both roles.",
    )
    sources.add_argument("--report", type=_path)
    sources.add_argument("--timeout-seconds", type=int, default=120)

    topology = subparsers.add_parser(
        "reconstruct-topology",
        help="Reconstruct graph blocks, row order, and Entrez identities.",
    )
    for name, text in (
        ("--archive", "OhmNet tissue-network tar archive."),
        ("--graph-spec", "Accepted selected_graphs.tsv specification."),
        ("--mapping-output", "Output node-to-Entrez TSV or TSV.GZ."),
        ("--edge-output", "Output logical edge TSV or TSV.GZ."),
        ("--summary-output", "Output JSON summary."),
    ):
        topology.add_argument(name, type=_path, required=True, help=text)
    topology.add_argument("--word-size-bits", type=int, choices=(32, 64), default=64)

    features = subparsers.add_parser(
        "reconstruct-features",
        help="Select and project the 50 MSigDB feature columns.",
    )
    for name, text in (
        ("--mapping", "Reconstructed node-to-Entrez mapping."),
        ("--c1", "MSigDB C1 Entrez GMT file."),
        ("--c3", "MSigDB C3 Entrez GMT file."),
        ("--feature-spec", "Accepted feature_columns.tsv specification."),
        ("--matrix-output", "Output NumPy matrix."),
        ("--selected-output", "Output selected-feature TSV."),
        ("--summary-output", "Output JSON summary."),
    ):
        features.add_argument(name, type=_path, required=True, help=text)
    features.add_argument("--minimum-source-members", type=int, default=200)
    features.add_argument("--maximum-columns", type=int, default=50)

    labels = subparsers.add_parser(
        "reconstruct-labels",
        help="Reconstruct the 121 GO-label columns from dated upstream sources.",
    )
    for name, text in (
        ("--mapping", "Reconstructed node-to-Entrez mapping."),
        ("--gaf", "GOA human release-159 GAF file."),
        ("--gpi", "GOA human release-159 GPI file."),
        ("--gp2protein", "Historical GeneID-to-UniProt mapping."),
        ("--ontology", "June 2016 GO OBO ontology."),
        ("--label-columns", "Accepted label_columns.tsv specification."),
        ("--identifier-decisions", "Accepted identifier decisions."),
        ("--matrix-output", "Output NumPy label matrix."),
        ("--class-map-output", "Output GraphSAGE-style class map."),
        ("--selected-labels-output", "Output selected-label metadata TSV."),
        ("--summary-output", "Output JSON reconstruction summary."),
    ):
        labels.add_argument(name, type=_path, required=True, help=text)
    labels.add_argument(
        "--evidence-code",
        action="append",
        dest="evidence_codes",
        required=True,
        help="Accepted GAF evidence code; repeat for each code.",
    )
    labels.add_argument(
        "--allowed-relation",
        action="append",
        dest="allowed_relations",
        required=True,
        help="Accepted normalized relation; repeat for each relation.",
    )
    labels.add_argument("--selected-term-count", type=int, default=121)

    graphsage = subparsers.add_parser(
        "assemble-graphsage",
        help="Assemble the four deterministic files used by supervised GraphSAGE.",
    )
    for name, text in (
        ("--mapping", "Reconstructed node-to-Entrez mapping."),
        ("--edges", "Reconstructed undirected edge table."),
        ("--features", "Reconstructed GraphSAGE feature matrix."),
        ("--labels", "Reconstructed GraphSAGE label matrix."),
        ("--output-directory", "Directory for the four ppi-* files."),
        ("--summary-output", "Output JSON assembly summary."),
    ):
        graphsage.add_argument(name, type=_path, required=True, help=text)
    graphsage.add_argument("--word-size-bits", type=int, choices=(32, 64), default=64)

    manifest = subparsers.add_parser(
        "write-manifest", help="Write a reconstruction-stage provenance manifest."
    )
    for name, text in (
        ("--output", "Output JSON manifest."),
        ("--project-root", "Logical project root (the ai directory)."),
        ("--reproduction-root", "Reproduction package root."),
        ("--source-report", "Source-verification JSON report."),
    ):
        manifest.add_argument(name, type=_path, required=True, help=text)
    manifest.add_argument("--stage", required=True)
    manifest.add_argument(
        "--artifact",
        action="append",
        type=_path,
        required=True,
        help="Generated artifact to hash; repeat as needed.",
    )
    manifest.add_argument(
        "--specification",
        action="append",
        type=_path,
        default=[],
        help="Committed specification to hash; repeat as needed.",
    )

    checks = subparsers.add_parser(
        "check-milestone",
        help="Check target-independent topology and feature invariants.",
    )
    for name, text in (
        ("--specification", "Frozen specification.yaml file."),
        ("--topology-summary", "Generated topology summary JSON."),
        ("--feature-summary", "Generated feature summary JSON."),
        ("--label-summary", "Generated GO-label summary JSON."),
        ("--graphsage-summary", "Generated GraphSAGE assembly summary JSON."),
        ("--output", "Output target-independent check report."),
    ):
        checks.add_argument(name, type=_path, required=True, help=text)

    validation = subparsers.add_parser(
        "validate-topology-features",
        help="Compare topology and features with the released GraphSAGE target.",
    )
    for name, text in (
        ("--reference", "Released GraphSAGE PPI ZIP."),
        ("--mapping", "Reconstructed node-to-Entrez mapping."),
        ("--edges", "Reconstructed logical edge table."),
        ("--features", "Reconstructed ppi-feats.npy."),
        ("--json-output", "Machine-readable validation report."),
        ("--markdown-output", "Human-readable validation report."),
    ):
        validation.add_argument(name, type=_path, required=True, help=text)

    graphsage_validation = subparsers.add_parser(
        "validate-graphsage",
        help="Validate the complete supervised GraphSAGE artifact set.",
    )
    for name, text in (
        ("--reference", "Released GraphSAGE PPI ZIP."),
        ("--graph", "Reconstructed ppi-G.json."),
        ("--id-map", "Reconstructed ppi-id_map.json."),
        ("--class-map", "Reconstructed ppi-class_map.json."),
        ("--features", "Reconstructed ppi-feats.npy."),
        ("--json-output", "Machine-readable validation report."),
        ("--markdown-output", "Human-readable validation report."),
    ):
        graphsage_validation.add_argument(name, type=_path, required=True, help=text)

    label_validation = subparsers.add_parser(
        "validate-labels",
        help="Compare reconstructed labels with the released GraphSAGE target.",
    )
    for name, text in (
        ("--reference", "Released GraphSAGE PPI ZIP."),
        ("--matrix", "Reconstructed ppi-labels.npy."),
        ("--class-map", "Reconstructed ppi-class_map.json."),
        ("--json-output", "Machine-readable validation report."),
        ("--markdown-output", "Human-readable validation report."),
    ):
        label_validation.add_argument(name, type=_path, required=True, help=text)

    clean = subparsers.add_parser(
        "clean", help="Remove generated build and result directories only."
    )
    clean.add_argument(
        "--directory",
        action="append",
        type=_path,
        required=True,
        help="Generated directory to remove; repeat as needed.",
    )
    return parser


def _remove_generated_directory(path: Path) -> None:
    resolved = path.resolve()
    if path.name not in {"build", "results"}:
        raise SourceError(
            f"Refusing to remove {path}: only directories named build or results "
            "are allowed"
        )
    if resolved == Path(resolved.anchor):
        raise SourceError(f"Refusing to remove filesystem root: {resolved}")
    if resolved.exists() and not resolved.is_dir():
        raise SourceError(f"Refusing to remove non-directory path: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved, ignore_errors=False)


def _run(arguments: argparse.Namespace) -> None:
    if arguments.command == "ensure-sources":
        ensure_sources(
            arguments.manifest,
            arguments.data_dir,
            source_ids=arguments.source_ids,
            roles=arguments.roles,
            report_path=arguments.report,
            timeout_seconds=arguments.timeout_seconds,
        )
        return

    if arguments.command == "reconstruct-topology":
        reconstruct_topology(
            archive_path=arguments.archive,
            graph_spec_path=arguments.graph_spec,
            mapping_output=arguments.mapping_output,
            edge_output=arguments.edge_output,
            summary_output=arguments.summary_output,
            word_size_bits=arguments.word_size_bits,
        )
        return

    if arguments.command == "reconstruct-features":
        reconstruct_features(
            mapping_path=arguments.mapping,
            c1_path=arguments.c1,
            c3_path=arguments.c3,
            feature_spec_path=arguments.feature_spec,
            matrix_output=arguments.matrix_output,
            selected_output=arguments.selected_output,
            summary_output=arguments.summary_output,
            minimum_source_members=arguments.minimum_source_members,
            maximum_columns=arguments.maximum_columns,
        )
        return

    if arguments.command == "reconstruct-labels":
        reconstruct_labels(
            mapping_path=arguments.mapping,
            gaf_path=arguments.gaf,
            gpi_path=arguments.gpi,
            gp2protein_path=arguments.gp2protein,
            ontology_path=arguments.ontology,
            label_columns_path=arguments.label_columns,
            identifier_decisions_path=arguments.identifier_decisions,
            matrix_output=arguments.matrix_output,
            class_map_output=arguments.class_map_output,
            selected_labels_output=arguments.selected_labels_output,
            summary_output=arguments.summary_output,
            evidence_codes=arguments.evidence_codes,
            allowed_relations=arguments.allowed_relations,
            selected_term_count=arguments.selected_term_count,
        )
        return

    if arguments.command == "assemble-graphsage":
        assemble_graphsage(
            mapping_path=arguments.mapping,
            edge_path=arguments.edges,
            feature_matrix_path=arguments.features,
            label_matrix_path=arguments.labels,
            output_directory=arguments.output_directory,
            summary_output=arguments.summary_output,
            word_size_bits=arguments.word_size_bits,
        )
        return

    if arguments.command == "check-milestone":
        check_reconstruction_milestone(
            specification_path=arguments.specification,
            topology_summary_path=arguments.topology_summary,
            feature_summary_path=arguments.feature_summary,
            label_summary_path=arguments.label_summary,
            graphsage_summary_path=arguments.graphsage_summary,
            output_path=arguments.output,
        )
        return

    if arguments.command == "write-manifest":
        write_run_manifest(
            arguments.output,
            project_root=arguments.project_root,
            reproduction_root=arguments.reproduction_root,
            source_report=arguments.source_report,
            artifact_paths=arguments.artifact,
            specification_paths=arguments.specification,
            stage=arguments.stage,
        )
        return

    if arguments.command == "validate-topology-features":
        validate_topology_and_features(
            graphsage_reference_zip=arguments.reference,
            mapping_path=arguments.mapping,
            edge_path=arguments.edges,
            feature_matrix_path=arguments.features,
            output_json=arguments.json_output,
            output_markdown=arguments.markdown_output,
        )
        return

    if arguments.command == "validate-graphsage":
        validate_graphsage_dataset(
            reference_archive=arguments.reference,
            graph_path=arguments.graph,
            id_map_path=arguments.id_map,
            class_map_path=arguments.class_map,
            feature_path=arguments.features,
            json_output=arguments.json_output,
            markdown_output=arguments.markdown_output,
        )
        return

    if arguments.command == "validate-labels":
        validate_labels_against_graphsage(
            reference_archive=arguments.reference,
            reconstructed_matrix=arguments.matrix,
            reconstructed_class_map=arguments.class_map,
            json_output=arguments.json_output,
            markdown_output=arguments.markdown_output,
        )
        return

    if arguments.command == "clean":
        for directory in arguments.directory:
            _remove_generated_directory(directory)
        return

    raise AssertionError(f"Unhandled command: {arguments.command}")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line interface and return a process status code."""

    parser = _build_parser()
    arguments = parser.parse_args(argv)
    try:
        _run(arguments)
    except (
        FeatureError,
        GraphSAGEError,
        IdentifierError,
        LabelError,
        SourceError,
        TopologyError,
        ValidationError,
        OSError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

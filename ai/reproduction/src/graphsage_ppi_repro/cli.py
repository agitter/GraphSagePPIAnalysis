"""Thin command-line adapters used by the Snakemake workflow.

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
from .provenance import (
    SourceError,
    check_reconstruction_milestone,
    ensure_sources,
    write_run_manifest,
)
from .topology import TopologyError, reconstruct_topology
from .validate import ValidationError, validate_topology_and_features


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

    if arguments.command == "check-milestone":
        check_reconstruction_milestone(
            specification_path=arguments.specification,
            topology_summary_path=arguments.topology_summary,
            feature_summary_path=arguments.feature_summary,
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
    except (FeatureError, SourceError, TopologyError, ValidationError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

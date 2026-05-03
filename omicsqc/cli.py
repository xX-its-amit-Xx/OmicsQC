"""Command-line entry point for OmicsQC.

Examples:
    python -m omicsqc --input examples/example.fastq --outdir outputs/fastq_report
    python -m omicsqc --input examples/example.fasta --outdir outputs/fasta_report
    zcat sample.fq.gz | python -m omicsqc -i - --format fastq -o outputs/sample
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .metrics import (
    DEFAULT_Q_THRESHOLDS,
    PHRED33_OFFSET,
    compute_fasta_metrics,
    compute_fastq_metrics,
    summarize_fasta,
    summarize_fastq,
)
from .parsers import parse_fasta, parse_fastq
from .plots import plot_gc_histogram, plot_length_histogram, plot_per_base_quality
from .report import render_markdown_report
from .utils import (
    EmptyFileError,
    MalformedRecordError,
    OmicsQCError,
    UnsupportedExtensionError,
    configure_logging,
    detect_file_kind,
    ensure_outdir,
    is_stdin,
)


def _parse_thresholds(text: str) -> tuple[int, ...]:
    try:
        values = tuple(int(x) for x in text.split(",") if x.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"--q-thresholds must be a comma-separated list of integers, got {text!r}"
        ) from exc
    if not values or any(v < 0 or v > 93 for v in values):
        raise argparse.ArgumentTypeError(
            "--q-thresholds values must be in [0, 93]"
        )
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="omicsqc",
        description="Run quality control on a FASTQ or FASTA file and emit a report.",
        epilog=(
            "Examples:\n"
            "  omicsqc -i reads.fastq.gz -o out/\n"
            "  zcat reads.fq.gz | omicsqc -i - --format fastq -o out/\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to a .fastq, .fq, .fasta, .fa file (gz/bgz allowed). Use '-' for stdin.",
    )
    parser.add_argument(
        "--outdir", "-o", required=True,
        help="Directory where the report and plots will be written.",
    )
    parser.add_argument(
        "--format", "-f", choices=("fastq", "fasta"), default=None,
        help="Override file-type detection (required when reading from stdin).",
    )
    parser.add_argument(
        "--phred-offset", type=int, default=PHRED33_OFFSET, choices=(33, 64),
        help="Phred quality offset (default 33; use 64 for legacy Illumina <1.8).",
    )
    parser.add_argument(
        "--q-thresholds", type=_parse_thresholds, default=DEFAULT_Q_THRESHOLDS,
        metavar="N[,N,...]",
        help="Comma-separated Phred thresholds for 'reads-below' counts (default 20,30).",
    )
    parser.add_argument(
        "--max-reads", type=int, default=None, metavar="N",
        help="Stop after the first N records (useful for sampling huge files).",
    )
    parser.add_argument(
        "--no-plots", action="store_true",
        help="Skip plot generation (writes summary metrics and report only).",
    )

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("--verbose", "-v", action="store_true", help="Enable debug logging.")
    verbosity.add_argument("--quiet", "-q", action="store_true", help="Suppress info logging.")

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def run(
    input_path: str,
    outdir: str,
    *,
    format: str | None = None,
    phred_offset: int = PHRED33_OFFSET,
    q_thresholds: tuple[int, ...] = DEFAULT_Q_THRESHOLDS,
    max_reads: int | None = None,
    no_plots: bool = False,
    verbose: bool = False,
    quiet: bool = False,
) -> Path:
    """Top-level pipeline. Returns the path to the generated Markdown report."""
    log = configure_logging(verbose=verbose, quiet=quiet)
    in_path = Path(input_path)
    out_path = ensure_outdir(outdir)

    log.info("OmicsQC v%s", __version__)
    log.info("Input: %s%s", in_path, " (stdin)" if is_stdin(input_path) else "")
    log.info("Output directory: %s", out_path)

    kind = detect_file_kind(input_path, hint=format)  # type: ignore[arg-type]
    log.info("File type: %s", kind)
    if max_reads is not None:
        log.info("Sampling at most %d records", max_reads)

    plots: list[Path] = []

    if kind == "fastq":
        log.info("Parsing FASTQ and computing metrics ...")
        metrics = compute_fastq_metrics(
            parse_fastq(input_path),
            phred_offset=phred_offset,
            q_thresholds=q_thresholds,
            max_reads=max_reads,
        )
        if metrics.total_reads == 0:
            raise OmicsQCError(f"No FASTQ records found in {in_path}.")
        summary = summarize_fastq(metrics)

        if not no_plots:
            log.info("Generating plots ...")
            plots.append(plot_length_histogram(metrics.read_lengths, out_path,
                                               title="FASTQ read length distribution"))
            plots.append(plot_gc_histogram(metrics.gc_per_read, out_path,
                                           title="FASTQ GC content distribution"))
            plots.append(plot_per_base_quality(metrics.per_base_mean_quality(), out_path))
    else:
        log.info("Parsing FASTA and computing metrics ...")
        metrics = compute_fasta_metrics(
            parse_fasta(input_path),
            max_records=max_reads,
        )
        if metrics.total_sequences == 0:
            raise OmicsQCError(f"No FASTA records found in {in_path}.")
        summary = summarize_fasta(metrics)

        if not no_plots:
            log.info("Generating plots ...")
            plots.append(plot_length_histogram(metrics.sequence_lengths, out_path,
                                               title="FASTA sequence length distribution"))
            plots.append(plot_gc_histogram(metrics.gc_per_sequence, out_path,
                                           title="FASTA GC content distribution"))

    log.info("Writing summary_metrics.json and summary_metrics.csv ...")
    (out_path / "summary_metrics.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.DataFrame([summary]).to_csv(out_path / "summary_metrics.csv", index=False)

    log.info("Rendering Markdown report ...")
    report_path = render_markdown_report(summary, plots, in_path, out_path)

    log.info("Done. Report: %s", report_path)
    return report_path


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    log = configure_logging(verbose=args.verbose, quiet=args.quiet)
    try:
        run(
            args.input,
            args.outdir,
            format=args.format,
            phred_offset=args.phred_offset,
            q_thresholds=args.q_thresholds,
            max_reads=args.max_reads,
            no_plots=args.no_plots,
            verbose=args.verbose,
            quiet=args.quiet,
        )
        return 0
    except FileNotFoundError as exc:
        log.error("File not found: %s", exc)
        return 2
    except UnsupportedExtensionError as exc:
        log.error("%s", exc)
        return 3
    except EmptyFileError as exc:
        log.error("%s", exc)
        return 4
    except MalformedRecordError as exc:
        log.error("Malformed input: %s", exc)
        return 5
    except OmicsQCError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.error("Interrupted by user.")
        return 130
    except Exception as exc:
        log.exception("Unexpected error: %s", exc)
        return 99


if __name__ == "__main__":
    sys.exit(main())

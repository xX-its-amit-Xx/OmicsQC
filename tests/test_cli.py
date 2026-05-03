"""End-to-end tests for the CLI pipeline (cli.run and cli.main)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from omicsqc.cli import build_parser, main, run

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def test_run_fastq_end_to_end(tmp_path: Path) -> None:
    report = run(str(EXAMPLES / "example.fastq"), str(tmp_path), quiet=True)
    assert report.exists()
    summary = json.loads((tmp_path / "summary_metrics.json").read_text())
    assert summary["file_kind"] == "fastq"
    assert summary["total_reads"] == 10
    assert (tmp_path / "summary_metrics.csv").exists()
    assert (tmp_path / "length_histogram.png").exists()
    assert (tmp_path / "gc_histogram.png").exists()
    assert (tmp_path / "per_base_quality.png").exists()


def test_run_fasta_end_to_end(tmp_path: Path) -> None:
    report = run(str(EXAMPLES / "example.fasta"), str(tmp_path), quiet=True)
    assert report.exists()
    summary = json.loads((tmp_path / "summary_metrics.json").read_text())
    assert summary["file_kind"] == "fasta"
    assert summary["total_sequences"] == 6
    assert "n50" in summary
    # per_base_quality.png should NOT be produced for FASTA input
    assert not (tmp_path / "per_base_quality.png").exists()


def test_run_no_plots_skips_pngs(tmp_path: Path) -> None:
    run(str(EXAMPLES / "example.fastq"), str(tmp_path), no_plots=True, quiet=True)
    assert (tmp_path / "summary_metrics.json").exists()
    assert not (tmp_path / "length_histogram.png").exists()


def test_run_max_reads_caps_record_count(tmp_path: Path) -> None:
    run(str(EXAMPLES / "example.fastq"), str(tmp_path), max_reads=3, quiet=True)
    summary = json.loads((tmp_path / "summary_metrics.json").read_text())
    assert summary["total_reads"] == 3


def test_run_custom_q_thresholds(tmp_path: Path) -> None:
    run(
        str(EXAMPLES / "example.fastq"),
        str(tmp_path),
        q_thresholds=(15, 35),
        quiet=True,
    )
    summary = json.loads((tmp_path / "summary_metrics.json").read_text())
    assert "reads_below_q15" in summary
    assert "reads_below_q35" in summary
    assert "reads_below_q20" not in summary


def test_main_returns_zero_on_success(tmp_path: Path) -> None:
    rc = main(["-i", str(EXAMPLES / "example.fastq"), "-o", str(tmp_path), "--quiet"])
    assert rc == 0


def test_main_returns_2_for_missing_file(tmp_path: Path) -> None:
    rc = main(["-i", str(tmp_path / "nope.fastq"), "-o", str(tmp_path / "out"), "--quiet"])
    assert rc == 2


def test_main_returns_3_for_unsupported_extension(tmp_path: Path) -> None:
    bad = tmp_path / "data.txt"
    bad.write_text("hello\n")
    rc = main(["-i", str(bad), "-o", str(tmp_path / "out"), "--quiet"])
    assert rc == 3


def test_main_returns_4_for_empty_file(tmp_path: Path) -> None:
    empty = tmp_path / "empty.fastq"
    empty.write_text("")
    rc = main(["-i", str(empty), "-o", str(tmp_path / "out"), "--quiet"])
    assert rc == 4


def test_main_returns_5_for_malformed_fastq(tmp_path: Path) -> None:
    bad = tmp_path / "bad.fastq"
    bad.write_text("@r\nACGT\n+\nII\n")  # quality length mismatch
    rc = main(["-i", str(bad), "-o", str(tmp_path / "out"), "--quiet"])
    assert rc == 5


def test_q_thresholds_arg_validation_rejects_garbage() -> None:
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["-i", "x.fq", "-o", "y", "--q-thresholds", "abc"])
    with pytest.raises(SystemExit):
        parser.parse_args(["-i", "x.fq", "-o", "y", "--q-thresholds", "200"])

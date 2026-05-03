"""Tests for FASTQ and FASTA parsers, including malformed-input handling."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from omicsqc.parsers import FastaRecord, FastqRecord, parse_fasta, parse_fastq
from omicsqc.utils import (
    EmptyFileError,
    MalformedRecordError,
    UnsupportedExtensionError,
    detect_file_kind,
)


VALID_FASTQ = (
    "@read1\n"
    "ACGTACGT\n"
    "+\n"
    "IIIIIIII\n"
    "@read2\n"
    "GGCCAATT\n"
    "+\n"
    "########\n"
)

VALID_FASTA = (
    ">seq1 description\n"
    "ACGTACGT\n"
    "GGCCAATT\n"
    ">seq2\n"
    "NNNNAAAA\n"
)


def test_parse_fastq_yields_two_records(tmp_path: Path) -> None:
    p = tmp_path / "ok.fastq"
    p.write_text(VALID_FASTQ, encoding="utf-8")

    records = list(parse_fastq(p))
    assert len(records) == 2
    assert records[0] == FastqRecord(header="read1", sequence="ACGTACGT", quality="IIIIIIII")
    assert records[1].sequence == "GGCCAATT"
    assert all(len(r.sequence) == len(r.quality) for r in records)


def test_parse_fastq_handles_gzip(tmp_path: Path) -> None:
    p = tmp_path / "ok.fastq.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(VALID_FASTQ)

    records = list(parse_fastq(p))
    assert len(records) == 2
    assert records[0].header == "read1"


def test_parse_fastq_truncated_record_raises(tmp_path: Path) -> None:
    truncated = "@read1\nACGT\n+\n"  # missing quality line
    p = tmp_path / "bad.fastq"
    p.write_text(truncated, encoding="utf-8")
    with pytest.raises(MalformedRecordError):
        list(parse_fastq(p))


def test_parse_fastq_missing_at_raises(tmp_path: Path) -> None:
    bad = "read1\nACGT\n+\nIIII\n"
    p = tmp_path / "bad.fastq"
    p.write_text(bad, encoding="utf-8")
    with pytest.raises(MalformedRecordError):
        list(parse_fastq(p))


def test_parse_fastq_seq_quality_length_mismatch_raises(tmp_path: Path) -> None:
    bad = "@r\nACGTAA\n+\nIII\n"
    p = tmp_path / "bad.fastq"
    p.write_text(bad, encoding="utf-8")
    with pytest.raises(MalformedRecordError):
        list(parse_fastq(p))


def test_parse_fasta_multiline_sequence(tmp_path: Path) -> None:
    p = tmp_path / "ok.fasta"
    p.write_text(VALID_FASTA, encoding="utf-8")

    records = list(parse_fasta(p))
    assert len(records) == 2
    assert records[0] == FastaRecord(header="seq1 description", sequence="ACGTACGTGGCCAATT")
    assert records[1].sequence == "NNNNAAAA"


def test_parse_fasta_handles_gzip(tmp_path: Path) -> None:
    p = tmp_path / "ok.fa.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(VALID_FASTA)
    records = list(parse_fasta(p))
    assert [r.header for r in records] == ["seq1 description", "seq2"]


def test_parse_fasta_sequence_before_header_raises(tmp_path: Path) -> None:
    bad = "ACGT\n>seq1\nACGT\n"
    p = tmp_path / "bad.fasta"
    p.write_text(bad, encoding="utf-8")
    with pytest.raises(MalformedRecordError):
        list(parse_fasta(p))


def test_parse_fastq_tolerates_trailing_blank_lines(tmp_path: Path) -> None:
    p = tmp_path / "trailing.fastq"
    p.write_text(VALID_FASTQ + "\n\n  \n", encoding="utf-8")
    assert len(list(parse_fastq(p))) == 2


def test_parse_fasta_tolerates_blank_lines_between_records(tmp_path: Path) -> None:
    p = tmp_path / "spaced.fasta"
    p.write_text(">s1\nACGT\n\n\n>s2\nGGCC\n\n", encoding="utf-8")
    records = list(parse_fasta(p))
    assert [r.sequence for r in records] == ["ACGT", "GGCC"]


def test_parse_fasta_strips_internal_whitespace(tmp_path: Path) -> None:
    p = tmp_path / "wrapped.fasta"
    # Some pretty-printed FASTAs use spaces every N bases.
    p.write_text(">s1\nACGT ACGT\nGGCC AATT\n", encoding="utf-8")
    records = list(parse_fasta(p))
    assert records[0].sequence == "ACGTACGTGGCCAATT"


def test_empty_file_raises(tmp_path: Path) -> None:
    p = tmp_path / "empty.fastq"
    p.write_text("", encoding="utf-8")
    with pytest.raises(EmptyFileError):
        list(parse_fastq(p))


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        list(parse_fastq(tmp_path / "does_not_exist.fastq"))


def test_detect_file_kind() -> None:
    assert detect_file_kind("a.fastq") == "fastq"
    assert detect_file_kind("a.fq") == "fastq"
    assert detect_file_kind("a.fastq.gz") == "fastq"
    assert detect_file_kind("a.fq.gz") == "fastq"
    assert detect_file_kind("a.fasta") == "fasta"
    assert detect_file_kind("a.fa") == "fasta"
    assert detect_file_kind("a.fasta.gz") == "fasta"
    assert detect_file_kind("a.fa.gz") == "fasta"


def test_detect_file_kind_unsupported_raises() -> None:
    with pytest.raises(UnsupportedExtensionError):
        detect_file_kind("a.txt")
    with pytest.raises(UnsupportedExtensionError):
        detect_file_kind("a.bam")


def test_detect_file_kind_with_hint() -> None:
    # Hint overrides extension, including for stdin.
    assert detect_file_kind("-", hint="fastq") == "fastq"
    assert detect_file_kind("/dev/stdin", hint="fasta") == "fasta"
    assert detect_file_kind("foo.txt", hint="fastq") == "fastq"


def test_detect_file_kind_stdin_without_hint_raises() -> None:
    with pytest.raises(UnsupportedExtensionError):
        detect_file_kind("-")


def test_open_text_detects_gzip_by_magic_bytes(tmp_path: Path) -> None:
    """A .gz extension is not required if the file is actually gzipped."""
    from omicsqc.utils import open_text
    p = tmp_path / "mislabeled.fastq"  # no .gz extension, gzipped content
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write(VALID_FASTQ)
    with open_text(p) as fh:
        assert fh.read().startswith("@read1")

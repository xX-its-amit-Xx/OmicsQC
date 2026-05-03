"""Tests for GC content, Phred scoring, N50, and FASTQ/FASTA metric aggregation."""

from __future__ import annotations

import math

from omicsqc.metrics import (
    PHRED64_OFFSET,
    compute_fasta_metrics,
    compute_fastq_metrics,
    gc_content,
    mean_quality,
    n50,
    n_content,
    phred_scores,
    summarize_fasta,
    summarize_fastq,
)
from omicsqc.parsers import FastaRecord, FastqRecord


def test_gc_content_basic() -> None:
    assert gc_content("GCGC") == 1.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ACGT") == 0.5


def test_gc_content_case_insensitive() -> None:
    assert gc_content("gcgc") == 1.0
    assert gc_content("AcGt") == 0.5


def test_gc_content_excludes_ambiguous_bases() -> None:
    # N excluded from denominator
    assert gc_content("NNAA") == 0.0
    assert gc_content("NNGC") == 1.0
    # Other IUPAC ambiguity codes are also excluded
    assert gc_content("RYACGT") == 0.5


def test_gc_content_handles_rna() -> None:
    # U should be treated as T
    assert gc_content("AUGC") == 0.5
    assert gc_content("UUUU") == 0.0


def test_gc_content_empty_returns_zero() -> None:
    assert gc_content("") == 0.0
    assert gc_content("NNNN") == 0.0


def test_phred_scores_phred33() -> None:
    assert phred_scores("!I") == [0, 40]
    assert phred_scores("") == []


def test_phred_scores_phred64() -> None:
    # '@' is 64 -> Q0; 'h' is 104 -> Q40 in Phred+64
    assert phred_scores("@h", offset=PHRED64_OFFSET) == [0, 40]


def test_mean_quality() -> None:
    assert mean_quality("IIII") == 40.0
    assert mean_quality("!I") == 20.0
    assert mean_quality("") == 0.0


def test_n_content_counts_case_insensitive() -> None:
    assert n_content("NNNN") == 4
    assert n_content("nNAa") == 2
    assert n_content("ACGT") == 0


def test_n50_basic() -> None:
    # Cumulative: [9,8,7,6,5] -> sum 35, half=17.5; running 9, 17, 24 — first to cross 17.5 is 7
    assert n50([9, 8, 7, 6, 5]) == 7
    # Single contig
    assert n50([100]) == 100


def test_n50_empty() -> None:
    assert n50([]) == 0


def test_compute_fastq_metrics_aggregates_correctly() -> None:
    records = [
        FastqRecord(header="r1", sequence="ACGT", quality="IIII"),    # GC=0.5, mean Q=40
        FastqRecord(header="r2", sequence="GGCC", quality="!!!!"),    # GC=1.0, mean Q=0
        FastqRecord(header="r3", sequence="ATATAT", quality="555555"),  # GC=0.0, mean Q=20
    ]
    m = compute_fastq_metrics(records)
    assert m.total_reads == 3
    assert m.read_lengths == [4, 4, 6]
    assert math.isclose(m.gc_per_read[0], 0.5)
    assert math.isclose(m.gc_per_read[1], 1.0)
    assert math.isclose(m.gc_per_read[2], 0.0)
    # Strict <: r1 (40) no, r2 (0) yes, r3 (20) no
    assert m.reads_below_threshold[20] == 1
    # r1 (40) no, r2 (0) yes, r3 (20) yes
    assert m.reads_below_threshold[30] == 2
    assert len(m.per_base_quality_sums) == 6
    assert len(m.per_base_counts) == 6
    assert m.per_base_counts[4] == 1
    assert m.per_base_counts[5] == 1


def test_compute_fastq_custom_thresholds() -> None:
    records = [
        FastqRecord(header="r1", sequence="ACGT", quality="IIII"),  # mean Q=40
        FastqRecord(header="r2", sequence="ACGT", quality="5555"),  # mean Q=20
    ]
    m = compute_fastq_metrics(records, q_thresholds=(15, 25, 35))
    assert m.reads_below_threshold == {15: 0, 25: 1, 35: 1}


def test_compute_fastq_max_reads_caps_iteration() -> None:
    records = [FastqRecord(header=f"r{i}", sequence="A", quality="!") for i in range(100)]
    m = compute_fastq_metrics(records, max_reads=10)
    assert m.total_reads == 10


def test_compute_fastq_phred64_offset() -> None:
    # All-'h' (Q40 in Phred+64) should give mean Q=40, no reads below Q20/Q30.
    rec = [FastqRecord(header="r", sequence="ACGT", quality="hhhh")]
    m = compute_fastq_metrics(rec, phred_offset=PHRED64_OFFSET)
    assert m.mean_quality_per_read == [40.0]
    assert m.reads_below_threshold[20] == 0


def test_summarize_fastq_shape() -> None:
    records = [FastqRecord(header="r1", sequence="ACGT", quality="IIII")]
    summary = summarize_fastq(compute_fastq_metrics(records))
    expected_keys = {
        "file_kind", "total_reads", "total_bases", "min_read_length", "max_read_length",
        "mean_read_length", "mean_gc_content", "mean_phred_quality",
        "reads_below_q20", "reads_below_q20_pct",
        "reads_below_q30", "reads_below_q30_pct",
    }
    assert expected_keys.issubset(summary.keys())
    assert summary["total_reads"] == 1
    assert summary["total_bases"] == 4
    assert summary["mean_phred_quality"] == 40.0


def test_summarize_fastq_empty_metrics_does_not_divide_by_zero() -> None:
    summary = summarize_fastq(compute_fastq_metrics(iter([])))
    assert summary["total_reads"] == 0
    assert summary["mean_read_length"] == 0.0
    assert summary["mean_gc_content"] == 0.0
    assert summary["mean_phred_quality"] == 0.0
    assert summary["reads_below_q20_pct"] == 0.0


def test_compute_fasta_metrics_aggregates_correctly() -> None:
    records = [
        FastaRecord(header="s1", sequence="ACGT"),
        FastaRecord(header="s2", sequence="NNGGCC"),
        FastaRecord(header="s3", sequence="AAAA"),
    ]
    m = compute_fasta_metrics(records)
    assert m.total_sequences == 3
    assert m.sequence_lengths == [4, 6, 4]
    assert math.isclose(m.gc_per_sequence[0], 0.5)
    assert math.isclose(m.gc_per_sequence[1], 1.0)
    assert math.isclose(m.gc_per_sequence[2], 0.0)
    assert m.n_count == 2
    assert m.total_bases == 14


def test_summarize_fasta_includes_n50() -> None:
    records = [FastaRecord(header=f"s{i}", sequence="A" * length)
               for i, length in enumerate([9, 8, 7, 6, 5])]
    summary = summarize_fasta(compute_fasta_metrics(records))
    assert summary["n50"] == 7


def test_summarize_fasta_shape() -> None:
    records = [FastaRecord(header="s1", sequence="NNGC")]
    summary = summarize_fasta(compute_fasta_metrics(records))
    expected_keys = {
        "file_kind", "total_sequences", "min_sequence_length",
        "max_sequence_length", "mean_sequence_length", "n50", "mean_gc_content",
        "total_bases", "n_base_count", "n_base_pct",
    }
    assert expected_keys.issubset(summary.keys())
    assert summary["n_base_count"] == 2
    assert summary["total_bases"] == 4
    assert math.isclose(summary["n_base_pct"], 50.0)

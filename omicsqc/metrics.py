"""QC metric calculations for FASTQ and FASTA records.

All functions are pure and operate on plain strings or in-memory record streams,
which keeps them easy to unit test.

Note on Phred averaging: ``mean_quality`` returns the arithmetic mean of Phred
scores, matching what samtools/seqkit report. The probabilistically-correct
mean (convert Phred -> error probability, mean, convert back) gives slightly
different numbers; we choose the arithmetic version for cross-tool consistency.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from .parsers import FastaRecord, FastqRecord

PHRED33_OFFSET = 33
PHRED64_OFFSET = 64
DEFAULT_Q_THRESHOLDS = (20, 30)


@dataclass
class FastqMetrics:
    total_reads: int = 0
    read_lengths: list[int] = field(default_factory=list)
    gc_per_read: list[float] = field(default_factory=list)
    mean_quality_per_read: list[float] = field(default_factory=list)
    per_base_quality_sums: list[int] = field(default_factory=list)
    per_base_counts: list[int] = field(default_factory=list)
    reads_below_threshold: dict[int, int] = field(default_factory=dict)

    def per_base_mean_quality(self) -> list[float]:
        """Average Phred score at each base position across all reads."""
        return [
            (s / c) if c else 0.0
            for s, c in zip(self.per_base_quality_sums, self.per_base_counts)
        ]


@dataclass
class FastaMetrics:
    total_sequences: int = 0
    sequence_lengths: list[int] = field(default_factory=list)
    gc_per_sequence: list[float] = field(default_factory=list)
    n_count: int = 0
    total_bases: int = 0


def gc_content(sequence: str) -> float:
    """Return GC fraction in [0.0, 1.0]. Empty / all-ambiguous sequence -> 0.0.

    Counts only A/C/G/T/U toward the denominator (so N, R, Y, etc. are
    excluded). U is treated as T for GC computation, allowing the same
    function to handle RNA.
    """
    if not sequence:
        return 0.0
    upper = sequence.upper()
    a = upper.count("A")
    c = upper.count("C")
    g = upper.count("G")
    t = upper.count("T") + upper.count("U")
    valid = a + c + g + t
    if valid == 0:
        return 0.0
    return (c + g) / valid


def phred_scores(quality_string: str, offset: int = PHRED33_OFFSET) -> list[int]:
    """Convert an ASCII quality string into integer Q-scores using ``offset``."""
    return [ord(ch) - offset for ch in quality_string]


def mean_quality(quality_string: str, offset: int = PHRED33_OFFSET) -> float:
    """Mean Phred score for a single read (arithmetic mean). Empty string -> 0.0."""
    if not quality_string:
        return 0.0
    total = 0
    for ch in quality_string:
        total += ord(ch) - offset
    return total / len(quality_string)


def n_content(sequence: str) -> int:
    """Count of 'N' bases (case-insensitive) in a sequence."""
    return sequence.upper().count("N")


def n50(lengths: Iterable[int]) -> int:
    """Return N50: the length L such that contigs >= L cover >= 50% of total bases.

    Returns 0 for empty input.
    """
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    if total == 0:
        return 0
    half = total / 2.0
    running = 0
    for length in sorted_lengths:
        running += length
        if running >= half:
            return length
    return sorted_lengths[-1]


def compute_fastq_metrics(
    records: Iterable[FastqRecord],
    *,
    phred_offset: int = PHRED33_OFFSET,
    q_thresholds: tuple[int, ...] = DEFAULT_Q_THRESHOLDS,
    max_reads: int | None = None,
) -> FastqMetrics:
    """Walk a stream of FASTQ records and accumulate QC metrics in one pass.

    ``max_reads`` stops parsing after that many records (useful for very large
    files); ``None`` means no cap.
    """
    m = FastqMetrics(reads_below_threshold={t: 0 for t in q_thresholds})
    for rec in records:
        if max_reads is not None and m.total_reads >= max_reads:
            break

        m.total_reads += 1
        m.read_lengths.append(len(rec.sequence))
        m.gc_per_read.append(gc_content(rec.sequence))

        scores = phred_scores(rec.quality, offset=phred_offset)
        avg = (sum(scores) / len(scores)) if scores else 0.0
        m.mean_quality_per_read.append(avg)

        for t in q_thresholds:
            if avg < t:
                m.reads_below_threshold[t] += 1

        if len(scores) > len(m.per_base_quality_sums):
            extra = len(scores) - len(m.per_base_quality_sums)
            m.per_base_quality_sums.extend([0] * extra)
            m.per_base_counts.extend([0] * extra)
        for i, q in enumerate(scores):
            m.per_base_quality_sums[i] += q
            m.per_base_counts[i] += 1

    return m


def compute_fasta_metrics(
    records: Iterable[FastaRecord],
    *,
    max_records: int | None = None,
) -> FastaMetrics:
    """Walk a stream of FASTA records and accumulate QC metrics in one pass."""
    m = FastaMetrics()
    for rec in records:
        if max_records is not None and m.total_sequences >= max_records:
            break
        m.total_sequences += 1
        m.sequence_lengths.append(len(rec.sequence))
        m.gc_per_sequence.append(gc_content(rec.sequence))
        m.n_count += n_content(rec.sequence)
        m.total_bases += len(rec.sequence)
    return m


def summarize_fastq(m: FastqMetrics) -> dict:
    """Build a JSON/CSV-friendly summary dict from accumulated FASTQ metrics."""
    n = max(1, m.total_reads)
    avg_len = sum(m.read_lengths) / n
    avg_gc = sum(m.gc_per_read) / n
    avg_q = sum(m.mean_quality_per_read) / n

    summary: dict = {
        "file_kind": "fastq",
        "total_reads": m.total_reads,
        "total_bases": sum(m.read_lengths),
        "min_read_length": min(m.read_lengths) if m.read_lengths else 0,
        "max_read_length": max(m.read_lengths) if m.read_lengths else 0,
        "mean_read_length": round(avg_len, 3) if m.total_reads else 0.0,
        "mean_gc_content": round(avg_gc, 4) if m.total_reads else 0.0,
        "mean_phred_quality": round(avg_q, 3) if m.total_reads else 0.0,
    }
    for threshold, count in sorted(m.reads_below_threshold.items()):
        pct = (100.0 * count / n) if m.total_reads else 0.0
        summary[f"reads_below_q{threshold}"] = count
        summary[f"reads_below_q{threshold}_pct"] = round(pct, 3)
    return summary


def summarize_fasta(m: FastaMetrics) -> dict:
    """Build a JSON/CSV-friendly summary dict from accumulated FASTA metrics."""
    n = max(1, m.total_sequences)
    avg_len = sum(m.sequence_lengths) / n
    avg_gc = sum(m.gc_per_sequence) / n
    n_pct = (100.0 * m.n_count / m.total_bases) if m.total_bases else 0.0
    return {
        "file_kind": "fasta",
        "total_sequences": m.total_sequences,
        "min_sequence_length": min(m.sequence_lengths) if m.sequence_lengths else 0,
        "max_sequence_length": max(m.sequence_lengths) if m.sequence_lengths else 0,
        "mean_sequence_length": round(avg_len, 3) if m.total_sequences else 0.0,
        "n50": n50(m.sequence_lengths),
        "mean_gc_content": round(avg_gc, 4) if m.total_sequences else 0.0,
        "total_bases": m.total_bases,
        "n_base_count": m.n_count,
        "n_base_pct": round(n_pct, 4),
    }

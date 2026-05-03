"""Minimal FASTQ and FASTA parsers.

These parsers are intentionally implemented from scratch (no Biopython) so the
record-level logic is transparent and easy to audit. They yield small dataclass
records and raise :class:`MalformedRecordError` on structural problems.

Assumptions:
  * FASTQ is the strict 4-line-per-record form (header / sequence / '+' /
    quality). Multi-line sequence/quality variants are rare in practice and
    are not supported. Tools like FastQC and bcl2fastq emit the 4-line form.
  * Phred encoding is Phred+33. Phred+64 (Illumina <1.8) is not auto-detected;
    pass ``--phred-offset 64`` if you have legacy data.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Union

from .utils import MalformedRecordError, open_text

PathLike = Union[str, Path]


@dataclass(frozen=True)
class FastqRecord:
    """A single FASTQ entry."""
    header: str
    sequence: str
    quality: str

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.quality):
            raise MalformedRecordError(
                f"Sequence and quality lengths differ for record '{self.header}': "
                f"{len(self.sequence)} vs {len(self.quality)}"
            )


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA entry."""
    header: str
    sequence: str


def parse_fastq(path: PathLike) -> Iterator[FastqRecord]:
    """Yield :class:`FastqRecord` entries from a FASTQ file (gz-aware).

    Tolerates trailing whitespace at the end of the file. Raises
    :class:`MalformedRecordError` on truncated records, missing markers, or
    sequence/quality length mismatches.
    """
    with open_text(path) as fh:
        line_no = 0
        while True:
            header = fh.readline()
            line_no += 1
            if not header:
                return
            if not header.strip():
                # Tolerate blank lines between records or at EOF.
                continue

            sequence = fh.readline()
            plus = fh.readline()
            quality = fh.readline()
            line_no += 3

            if not (sequence and plus and quality):
                raise MalformedRecordError(
                    f"Truncated FASTQ record near line {line_no} in {path}."
                )

            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")

            if not header.startswith("@"):
                raise MalformedRecordError(
                    f"Expected FASTQ header starting with '@' at line {line_no - 3} "
                    f"in {path}, got: {header[:40]!r}"
                )
            if not plus.startswith("+"):
                raise MalformedRecordError(
                    f"Expected '+' separator at line {line_no - 1} in {path}, "
                    f"got: {plus[:40]!r}"
                )

            yield FastqRecord(header=header[1:], sequence=sequence, quality=quality)


def parse_fasta(path: PathLike) -> Iterator[FastaRecord]:
    """Yield :class:`FastaRecord` entries from a FASTA file (gz-aware).

    Multi-line sequences are concatenated. Whitespace inside a sequence line
    (occasionally present in pretty-printed FASTA) is stripped. Raises
    :class:`MalformedRecordError` if a sequence appears before any '>' header.
    """
    with open_text(path) as fh:
        header: str | None = None
        chunks: list[str] = []
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if header is not None:
                    yield FastaRecord(header=header, sequence="".join(chunks))
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise MalformedRecordError(
                        f"FASTA sequence appears before any '>' header in {path}."
                    )
                # Drop internal whitespace defensively (some tools wrap with spaces).
                chunks.append("".join(line.split()))

        if header is not None:
            yield FastaRecord(header=header, sequence="".join(chunks))

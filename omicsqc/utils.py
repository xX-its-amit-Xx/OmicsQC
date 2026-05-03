"""Utility helpers: file-type detection, opening (gzipped or plain), and logging setup."""

from __future__ import annotations

import gzip
import logging
import os
import sys
from pathlib import Path
from typing import IO, Literal

FileKind = Literal["fastq", "fasta"]

FASTQ_SUFFIXES = (".fastq", ".fq")
FASTA_SUFFIXES = (".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn")
COMPRESSED_SUFFIXES = (".gz", ".bgz")

STDIN_SENTINELS = ("-", "/dev/stdin")


class OmicsQCError(Exception):
    """Base exception for OmicsQC errors."""


class UnsupportedExtensionError(OmicsQCError):
    """Raised when an input file has an unsupported extension."""


class EmptyFileError(OmicsQCError):
    """Raised when an input file is empty."""


class MalformedRecordError(OmicsQCError):
    """Raised when a sequencing record cannot be parsed."""


def configure_logging(verbose: bool = False, quiet: bool = False) -> logging.Logger:
    """Configure and return the package logger.

    Verbose enables DEBUG; quiet drops below WARNING. Idempotent: re-calling
    will adjust level without duplicating handlers.
    """
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logger = logging.getLogger("omicsqc")
    logger.setLevel(level)

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                "[%(asctime)s] %(levelname)s %(name)s: %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        logger.addHandler(handler)

    logger.propagate = False
    return logger


def is_stdin(path: str | os.PathLike[str]) -> bool:
    """Return True if the path refers to stdin ('-' or /dev/stdin)."""
    return str(path) in STDIN_SENTINELS


def detect_file_kind(path: str | os.PathLike[str], hint: FileKind | None = None) -> FileKind:
    """Detect 'fastq' or 'fasta' from extension; honor explicit ``hint`` if given.

    Stdin (``-``) always requires an explicit hint, since there's no extension
    to look at.
    """
    if hint is not None:
        if hint not in ("fastq", "fasta"):
            raise UnsupportedExtensionError(
                f"hint must be 'fastq' or 'fasta', got {hint!r}"
            )
        return hint

    if is_stdin(path):
        raise UnsupportedExtensionError(
            "Cannot infer file type from stdin; pass --format fastq|fasta."
        )

    name = str(path).lower()
    for comp in COMPRESSED_SUFFIXES:
        if name.endswith(comp):
            name = name[: -len(comp)]
            break
    suffix = Path(name).suffix
    if suffix in FASTQ_SUFFIXES:
        return "fastq"
    if suffix in FASTA_SUFFIXES:
        return "fasta"
    raise UnsupportedExtensionError(
        f"Unsupported file extension for {path!s}. "
        f"Expected one of: {FASTQ_SUFFIXES + FASTA_SUFFIXES} "
        f"(optionally with {COMPRESSED_SUFFIXES})."
    )


def _is_gzipped(path: Path) -> bool:
    """Return True if the file starts with the gzip magic number 0x1f8b."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def open_text(path: str | os.PathLike[str]) -> IO[str]:
    """Open a sequencing file in text mode.

    Supports stdin via ``-``, gzipped files (detected by magic number, not just
    extension), and plain text. Decoding is lenient (``errors='replace'``) so a
    stray non-ASCII byte in a header doesn't blow up an entire run.
    """
    if is_stdin(path):
        # sys.stdin is already text-mode; strip any encoding strictness.
        return sys.stdin

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file does not exist: {p}")
    if p.stat().st_size == 0:
        raise EmptyFileError(f"Input file is empty: {p}")

    if _is_gzipped(p):
        return gzip.open(p, mode="rt", encoding="utf-8", errors="replace")
    return open(p, mode="rt", encoding="utf-8", errors="replace")


def ensure_outdir(outdir: str | os.PathLike[str]) -> Path:
    """Create the output directory (and parents) if needed; return as Path."""
    p = Path(outdir)
    p.mkdir(parents=True, exist_ok=True)
    return p

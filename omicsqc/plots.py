"""Matplotlib plotting helpers. Each function writes a PNG and returns its path.

The matplotlib backend is selected lazily inside each function rather than at
module import, so importing :mod:`omicsqc.plots` from a Jupyter notebook (or
any environment with a configured backend) won't override the user's choice.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path


def _ensure_headless_backend() -> None:
    """Switch matplotlib to the Agg backend if no GUI backend has been chosen."""
    import matplotlib

    current = matplotlib.get_backend().lower()
    if current in ("agg", "module://matplotlib_inline.backend_inline"):
        return
    # Only force Agg when there's no display configured.
    import os
    if os.environ.get("DISPLAY") or os.environ.get("MPLBACKEND"):
        return
    matplotlib.use("Agg", force=True)


def _save(fig, path: Path) -> Path:
    import matplotlib.pyplot as plt

    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def plot_length_histogram(
    lengths: Sequence[int],
    outdir: Path,
    *,
    title: str = "Read length distribution",
    filename: str = "length_histogram.png",
) -> Path:
    """Histogram of sequence/read lengths."""
    _ensure_headless_backend()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4))
    if lengths:
        bins = min(50, max(10, len(set(lengths))))
        ax.hist(lengths, bins=bins, color="#4C72B0", edgecolor="black", linewidth=0.4)
    ax.set_xlabel("Length (bp)")
    ax.set_ylabel("Count")
    ax.set_title(title)
    return _save(fig, outdir / filename)


def plot_gc_histogram(
    gc_values: Sequence[float],
    outdir: Path,
    *,
    title: str = "GC content distribution",
    filename: str = "gc_histogram.png",
) -> Path:
    """Histogram of per-read/per-sequence GC fractions (0-1)."""
    _ensure_headless_backend()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4))
    if gc_values:
        ax.hist(
            [v * 100 for v in gc_values],
            bins=40, range=(0, 100),
            color="#55A868", edgecolor="black", linewidth=0.4,
        )
    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("Count")
    ax.set_title(title)
    return _save(fig, outdir / filename)


def plot_per_base_quality(
    per_base_mean: Sequence[float],
    outdir: Path,
    *,
    filename: str = "per_base_quality.png",
    y_max: float | None = None,
) -> Path:
    """Line plot of mean Phred quality vs. base position.

    ``y_max`` defaults to ``max(observed, 42)`` so the plot frames Q20/Q30
    references comfortably for Illumina-scale data, but auto-extends for
    long-read platforms (PacBio HiFi, Q40+).
    """
    _ensure_headless_backend()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4))
    if per_base_mean:
        positions = list(range(1, len(per_base_mean) + 1))
        ax.plot(positions, per_base_mean, color="#C44E52", linewidth=1.6)
        ax.axhline(20, color="#888", linestyle="--", linewidth=0.8, label="Q20")
        ax.axhline(30, color="#444", linestyle="--", linewidth=0.8, label="Q30")
        ax.legend(loc="lower left")
        if y_max is None:
            y_max = max(42.0, max(per_base_mean) * 1.1)
    ax.set_xlabel("Base position")
    ax.set_ylabel("Mean Phred quality")
    ax.set_ylim(0, y_max if y_max is not None else 42)
    ax.set_title("Per-base mean quality")
    return _save(fig, outdir / filename)

# OmicsQC

A lightweight, dependency-light command-line toolkit for **basic quality control of FASTQ and FASTA files**. OmicsQC parses sequencing files, computes standard QC metrics (read length, GC content, Phred quality), generates plots, and produces a clean Markdown report — all with a single command.

The project is intentionally implemented in plain Python (no Biopython) so the parsing and metric logic is transparent and easy to read. It is meant as a small, well-tested starting point for a bioinformatics QC pipeline.

---

## Why sequencing QC matters

Before any downstream analysis — alignment, variant calling, expression quantification, assembly — sequencing data should be checked for problems that, if missed, will quietly bias every result that follows. Common issues include:

- Reads with very low Phred quality, especially toward the 3' end.
- Unexpected GC distributions, which can indicate contamination or library-prep bias.
- An excess of `N` (ambiguous) bases.
- Truncated, malformed, or empty files from a failed transfer or run.

OmicsQC surfaces these signals early so they can be addressed before they propagate.

---

## Features

- Accepts `.fastq`, `.fq`, `.fasta`, `.fa`, `.fna` and their gzipped variants (`.gz`).
- **FASTQ metrics:** total reads, length distribution, GC% per read and overall mean, per-base Phred quality, mean Phred quality, and the count and percentage of reads below Q20 and Q30.
- **FASTA metrics:** total sequences, length distribution, GC% per sequence and overall mean, total bases, and N-base count and percentage.
- Plots written as PNGs:
  - read/sequence length histogram
  - GC content histogram
  - per-base mean quality (FASTQ only)
- Outputs `summary_metrics.json`, `summary_metrics.csv`, plot PNGs, and a `final_report.md` that embeds everything.
- Friendly error messages for missing files, unsupported extensions, malformed FASTQ records, and empty files.
- Logging at every step.
- Unit tests with `pytest`.

---

## Installation

Clone the repository and install in a virtual environment:

```bash
git clone https://github.com/your-org/OmicsQC.git
cd OmicsQC

python -m venv .venv
source .venv/bin/activate          # macOS / Linux
.venv\Scripts\activate             # Windows

pip install -r requirements.txt
```

Or install as a package (gives you the `omicsqc` console script):

```bash
pip install -e .
```

Python 3.9 or newer is required.

---

## Usage

Run the toolkit on a FASTQ file:

```bash
python -m omicsqc --input examples/example.fastq --outdir outputs/fastq_report
```

Run it on a FASTA file:

```bash
python -m omicsqc --input examples/example.fasta --outdir outputs/fasta_report
```

Gzipped inputs work transparently:

```bash
python -m omicsqc -i sample.fastq.gz -o outputs/sample_report
```

Pipe from stdin (handy for streaming large or remote files):

```bash
zcat sample.fq.gz | python -m omicsqc -i - --format fastq -o outputs/sample_report
```

If you installed the package via `pip install -e .`, the same thing works as:

```bash
omicsqc --input examples/example.fastq --outdir outputs/fastq_report
```

### CLI flags

| Flag | Description |
| --- | --- |
| `--input`, `-i` | Path to a FASTQ or FASTA file (gz/bgz allowed). Use `-` for stdin. |
| `--outdir`, `-o` | Directory where outputs will be written (created if needed). |
| `--format`, `-f` | Override file-type detection (`fastq` or `fasta`). Required with stdin. |
| `--phred-offset` | 33 (default) or 64 (legacy Illumina <1.8). |
| `--q-thresholds` | Comma-separated Phred thresholds for "reads below" counts (default `20,30`). |
| `--max-reads` | Stop after the first N records — useful for sampling huge files. |
| `--no-plots` | Skip plot generation (writes summary metrics + report only). |
| `--verbose`, `-v` | Enable debug-level logging. |
| `--quiet`, `-q` | Suppress info-level logging. |
| `--version` | Print the OmicsQC version and exit. |

### Exit codes

`0` success · `2` input file not found · `3` unsupported extension · `4` empty file · `5` malformed record · `1` other OmicsQC error · `99` unexpected error · `130` interrupted.

---

## Output

After a run, the output directory contains:

```
outputs/fastq_report/
├── summary_metrics.json     # All summary metrics, machine-readable
├── summary_metrics.csv      # Same metrics as a single-row CSV
├── length_histogram.png     # Read/sequence length distribution
├── gc_histogram.png         # GC content distribution (%)
├── per_base_quality.png     # Per-base mean Phred quality (FASTQ only)
└── final_report.md          # Human-readable Markdown summary
```

Open `final_report.md` in any Markdown viewer (GitHub, VS Code, `grip`, etc.) to see the metrics table alongside the plots.

### What the metrics mean

| Metric | Notes |
| --- | --- |
| Total reads / sequences | Record count in the file. |
| Read/sequence length | Min, max, and mean length in bp. |
| Mean GC content | Fraction of `G`+`C` among unambiguous A/C/G/T bases. Reported as a percentage in the report. |
| Mean Phred quality | Arithmetic mean of Phred scores across all bases of all reads. (Note: this matches samtools/seqkit. The probabilistic mean — average error probability, then convert back — is a slightly different number; we use the arithmetic mean for cross-tool consistency.) |
| Reads below Q*N* | Reads whose **mean** quality is strictly below the threshold *N*. Configurable via `--q-thresholds`. |
| N50 | FASTA only. The contig length L such that contigs ≥ L cover ≥ 50% of total bases — a standard assembly-contiguity statistic. |
| N-base count / % | FASTA only. High values can indicate masked or low-coverage assembly regions. |

---

## Running the tests

```bash
pytest
```

The suite covers the FASTQ and FASTA parsers (including malformed-record handling), GC computation, Phred conversion, and the end-to-end metric aggregation.

---

## Project layout

```
OmicsQC/
├── omicsqc/
│   ├── __init__.py
│   ├── __main__.py     # enables `python -m omicsqc`
│   ├── cli.py          # argparse entry point
│   ├── parsers.py      # FASTQ / FASTA parsing
│   ├── metrics.py      # GC, Phred, N50, aggregation
│   ├── plots.py        # matplotlib helpers
│   ├── report.py       # Markdown report
│   └── utils.py        # logging, file detection, errors
├── tests/
│   ├── test_parsers.py
│   ├── test_metrics.py
│   └── test_cli.py     # end-to-end pipeline tests
├── .github/
│   └── workflows/ci.yml
├── examples/
│   ├── example.fastq
│   └── example.fasta
├── outputs/            # populated by CLI runs (git-ignored)
├── requirements.txt
├── pyproject.toml
├── .gitignore
└── README.md
```

---

## Future improvements

- Adapter / overrepresented-sequence detection (k-mer scan).
- Per-base composition plot (A/C/G/T/N proportions across positions).
- Paired-end FASTQ awareness (R1/R2 cross-checks, insert size hints).
- HTML report with interactive Plotly figures alongside the Markdown one.
- Multi-file batch mode that produces a combined cohort summary.
- Optional Biopython backend for spec-edge-case FASTQ variants.
- Sampling/streaming flags for very large files (e.g., first N reads, reservoir sample).

---

## License

GPL-3.0 — see [LICENSE](LICENSE).

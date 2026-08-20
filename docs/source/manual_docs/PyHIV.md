# 🧬 PyHIV: A Python Package for Local HIV‑1 Sequence Alignment, Subtyping and Gene Splitting

![CI](https://github.com/anaapspereira/PyHIV/actions/workflows/ci.yml/badge.svg)
[![codecov](https://codecov.io/gh/anaapspereira/PyHIV/branch/main/graph/badge.svg)](https://codecov.io/gh/anaapspereira/PyHIV)
![Python Version](https://img.shields.io/pypi/pyversions/pyhiv-tools)
![OS Supported](https://img.shields.io/badge/OS-Windows%20%7C%20Linux%20%7C%20macOS-blue)

![PyPI version](https://img.shields.io/pypi/v/pyhiv-tools)
![Documentation Status](https://readthedocs.org/projects/pyhiv/badge/?version=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub issues](https://img.shields.io/github/issues/anaapspereira/PyHIV)

---

## 📖 Overview

**PyHIV** is a Python tool that aligns HIV nucleotide sequences against reference genomes to determine the **most similar subtype** and optionally **split the aligned sequences into gene regions**.

It produces:
- Best reference alignment per sequence  
- Subtype and reference metadata  
- Ranked top 3 closest HIV-1 subtypes
- Gene-region–specific FASTA files (optional)  
- A final summary table (`final_table.tsv`)  

---

## ⚙️ How It Works
```
┌─────────────────────────────────────────────┐
│  User FASTA sequences                       │
└─────────────────────────────────────────────┘
                │
                ▼
       Read and preprocess input
                │
                ▼
 Align sequences against reference genomes
                │
                ▼
    Identify best matching reference
                │
                ▼
     (Optional) Split by gene region
                │
                ▼
  Save results and summary table (.tsv)

```

---

## 📦 Installation

You can install PyHIV using pip:

```bash
pip install pyhiv-tools
```

Alternatively, you can clone the repository and install it manually:

```bash
git clone https://github.com/anaapspereira/PyHIV.git
cd PyHIV
python setup.py install
```
## 🚀 Getting Started

Basic usage:

```python
from pyhiv import PyHIV

PyHIV(
    fastas_dir="path/to/fasta/files",
    subtyping=True,
    splitting=True,  # True/"subtype", "hxb2"/"reference", or False/"none"
    output_dir="results_folder",
    n_jobs=4,
    reporting=True,
    alignment_tool="edlib-HW",
    kmer_size=15,
    reference_top_k=30,
    reference_groups="M",
    show_progress=False,
)
```

### Parameters:

| Parameter    | Type   | Default           | Description                                                                |
| ------------ | ------ | ----------------- | -------------------------------------------------------------------------- |
| `fastas_dir` | `str`  | *Required*        | Directory containing user FASTA files.                                     |
| `subtyping`  | `bool` | `True`            | Aligns against subtype reference genomes. If `False`, aligns only to HXB2. |
| `splitting`  | `bool` or `str` | `True` | Splits aligned sequences into gene regions. Use `True`/`"subtype"` to split by the best subtype reference when it has annotated features, with automatic HXB2 fallback when features are missing; use `"hxb2"`/`"reference"` to subtype but always split against HXB2, or `False`/`"none"` to skip splitting. |
| `output_dir` | `str`  | `"PyHIV_results"` | Output directory for results.                                              |
| `n_jobs`     | `int`  | `None`            | Number of parallel worker processes, created once per run and reused for every sequence. `None` uses all available CPU cores. |
| `reporting`  | `bool` | `True`            | Generates a PDF report (`PyHIV_report_all_sequences.pdf`) with per-sequence visualizations. |
| `alignment_tool` | `str` | `"edlib-HW"` | Alignment backend: `edlib-HW`, `parasail-NW`/`parasail`, `MAFFT`, or `PyFamsa`. |
| `kmer_size` | `int` | `15` | K-mer size used to prefilter candidate references. |
| `reference_top_k` | `int` | `30` | Number of top k-mer ranked references to align. Use `0` to align all references. |
| `reference_groups` | `str` or iterable | `"M"` | HIV-1 reference groups used for subtyping. Use `"M,N,O,P"` or `"all"` to include groups N, O, and P. |
| `show_progress` | `bool` | `False` | Displays a terminal progress bar for processed input sequences.            |

`edlib-HW` is the default and projects alignments onto full-reference genome coordinates. `parasail-NW`/`parasail`, `PyFamsa`, and `MAFFT` remain available as alternatives. Before final alignment, PyHIV ranks references using query/reference k-mer containment and aligns only the top candidates by default. Use `reference_top_k=0` to keep the original all-reference strategy. By default, subtyping uses group M references from `reference_fastas`, selected through the `group` column in `sequences_with_locations.tsv`; set `reference_groups="M,N,O,P"` (or `"all"`) to include groups N, O, and P. `edlib` is installed with PyHIV. `parasail` and `PyFamsa` are optional extras: `parasail` — install with `pip install pyhiv-tools[parasail]` — since it has no prebuilt wheel on some platforms (e.g. macOS on Apple Silicon); `PyFamsa` — install with `pip install pyhiv-tools[famsa]` — since `pyfamsa` is GPL-3.0 licensed and is kept out of the default (MIT) install. `MAFFT` requires an external `mafft` executable. PyHIV resolves MAFFT from `PYHIV_MAFFT_BIN`, then `mafft` on `PATH`. PyHIV validates that the selected `alignment_tool` is actually available before processing starts and raises immediately with an install hint if it isn't, instead of failing per-reference partway through a run.

When `subtyping=False`, any active splitting mode is treated as HXB2-based splitting, even if `splitting="subtype"` is provided.

When `subtyping=True` and `splitting=True`, PyHIV keeps the selected subtype reference in the `Reference` column. If that reference has no annotated `features` in `sequences_with_locations.tsv`, PyHIV splits against HXB2 instead and records `K03455` in `Splitting Reference`; the PDF report also shows the splitting reference.

Input sequences longer than 12000 nucleotides are skipped with this warning: `The submitted sequence is longer than the HIV-1 genome.`

`fastas_dir` is searched recursively, but `output_dir` is always excluded from that search — so setting `output_dir` to a subfolder of `fastas_dir` is safe and won't cause a rerun to re-ingest files from a previous run.

The reference dataset root (`reference_fastas/`, `HXB2_fasta/`, `sequences_with_locations.tsv`) is resolved from the `REFERENCE_GENOMES_DIR` environment variable, falling back to the packaged `reference_genomes` directory. Set it to point PyHIV at a custom or offline reference dataset, e.g. one prepared with `pyhiv update reference-dataset --reference-dir <path>`.

### 📂 Output Structure

After running PyHIV, your output directory (default: PyHIV_results/) will contain:

```
PyHIV_results/
│
├── best_alignment_<file_stem>_<sequence>.fasta      # Alignment to best reference
├── splitting_alignment_<file_stem>_<sequence>.fasta # HXB2 alignment when splitting uses HXB2 with subtyping
├── final_table.tsv                     # Summary of results
├── PyHIV_report_all_sequences.pdf      # Per-sequence PDF report (when reporting=True)
│
├── gag/
│   ├── <file_stem>_<sequence>_gag.fasta
│   └── ...
├── pol/
│   ├── <sequence>_pol.fasta
│   └── ...
└── env/
    ├── <sequence>_env.fasta
    └── ...
```

### Final Table Columns:

| Column                        | Description                                     |
| ----------------------------- | ----------------------------------------------- |
| **File Name**                 | Source FASTA file name                          |
| **Sequence**                  | Input sequence name                             |
| **Reference**                 | Best matching reference accession               |
| **Group**                     | Predicted HIV-1 reference group                 |
| **Subtype**                   | Predicted HIV-1 subtype                         |
| **Closest Subtypes**          | Top 3 closest group/subtype calls by alignment score |
| **Subtype Score Warning**     | `Low score margin: review top 3 subtype matches` when the top two subtype scores are within 1% of each other, otherwise empty |
| **Splitting Reference**       | Reference accession used for gene splitting     |
| **Most Matching Gene Region** | Region with highest similarity                  |
| **Present Gene Regions**      | All detected gene regions with valid alignments |


---

## 📟 Command Line Interface

PyHIV provides a user-friendly CLI for HIV-1 sequence analysis.

### 🚀 Getting Started

```bash
# Basic usage
pyhiv run sequences/

# With custom options
pyhiv run sequences/ -o results/ -j 4 -v

# Validate inputs first
pyhiv validate sequences/

# Refresh the packaged reference dataset
pyhiv update reference-dataset --yes
```

### ⚙️ Main Options

| Option | Description |
|--------|-------------|
| `--subtyping BOOL` | Enable/disable HIV-1 subtyping (default: `true`) |
| `--splitting TEXT` | Splitting mode: `true`/`subtype`, `hxb2`/`reference`, or `false`/`none` (default: `true`). If `--subtyping false`, active splitting uses HXB2. |
| `-o`, `--output-dir PATH` | Output directory (default: `PyHIV_results`) |
| `-j`, `--n-jobs INTEGER` | Number of parallel jobs (default: all CPUs) |
| `--progress` / `--no-progress` | Show a terminal progress bar (default: `--progress`) |
| `--reporting` / `--no-reporting` | Generate a PDF report (default: `--reporting`) |
| `-v`, `--verbose` | Detailed output |
| `-q`, `--quiet` | Suppress non-error output |

See [CLI.md](CLI.md) for the full option list, including `--alignment-tool`, `--kmer-size`, `--reference-top-k`, `--reference-groups`, and the `pyhiv update reference-dataset` command.

### 💼 Common Use Cases

**Full analysis with subtyping and splitting:**
```bash
pyhiv run data/sequences/
```

**Alignment only:**
```bash
pyhiv run data/sequences/ --subtyping false --splitting false
```

**Subtyping without gene splitting:**
```bash
pyhiv run data/sequences/ --splitting false
```

**Subtyping with HXB2-based splitting:**
```bash
pyhiv run data/sequences/ --splitting hxb2
```

**Parallel processing:**
```bash
pyhiv run data/sequences/ -j 8 -o results/batch1/
```

**Validation:**
```bash
pyhiv validate data/sequences/
```

### 📤 Output

PyHIV generates:
- `final_table.tsv` - Summary with sequence IDs, references, group/subtype calls, closest subtypes, and gene regions
- `best_alignment_*.fasta` - Best alignment for each sequence
- Gene-specific folders (when `--splitting` is enabled) with extracted regions
- `PyHIV_report_all_sequences.pdf` - Per-sequence PDF report (when `--reporting` is enabled, the default)

### 🆘 Getting Help

```bash
pyhiv --help                          # Show all commands
pyhiv run --help                      # Show options for run command
pyhiv update reference-dataset --help # Show options for the reference dataset updater
pyhiv --version                       # Show version
```

For comprehensive documentation, see [CLI.md](CLI.md).

---

## 🗂️ Citation

Manuscript in preparation. Please cite this repository if you use PyHIV in your research.

---

## 🧾 License

This project is licensed under the MIT License — see the LICENSE
 file for details.

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
    splitting=True,
    output_dir="results_folder",
    n_jobs=4,
    alignment_tool="parasail-NW"
)
```

### Parameters:

| Parameter    | Type   | Default           | Description                                                                |
| ------------ | ------ | ----------------- | -------------------------------------------------------------------------- |
| `fastas_dir` | `str`  | *Required*        | Directory containing user FASTA files.                                     |
| `subtyping`  | `bool` | `True`            | Aligns against subtype reference genomes. If `False`, aligns only to HXB2. |
| `splitting`  | `bool` | `True`            | Splits aligned sequences into gene regions.                                |
| `output_dir` | `str`  | `"PyHIV_results"` | Output directory for results.                                              |
| `n_jobs`     | `int`  | `None`            | Number of parallel jobs for alignment.                                     |
| `alignment_tool` | `str` | `"parasail-NW"` | Alignment backend: `parasail-NW`, `edlib-HW`, `edlib-NW`, `MAFFT`, or `PyFamsa`. |

`edlib`, `parasail`, and `PyFamsa` are installed with PyHIV. `MAFFT` requires an external `mafft` executable. PyHIV resolves MAFFT from `PYHIV_MAFFT_BIN`, then `mafft` on `PATH`.

Input sequences longer than 12000 nucleotides are skipped with this warning: `The submitted sequence is longer than the HIV-1 genome.`

### 📂 Output Structure

After running PyHIV, your output directory (default: PyHIV_results/) will contain:

```
PyHIV_results/
│
├── best_alignment_<sequence>.fasta     # Alignment to best reference
├── final_table.tsv                     # Summary of results
│
├── gag/
│   ├── <sequence>_gag.fasta
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
| **Sequence**                  | Input sequence name                             |
| **Reference**                 | Best matching reference accession               |
| **Subtype**                   | Predicted HIV-1 subtype                         |
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
```

### ⚙️ Main Options

| Option | Description |
|--------|-------------|
| `--subtyping` / `--no-subtyping` | Enable/disable HIV-1 subtyping (default: enabled) |
| `--splitting` / `--no-splitting` | Enable/disable gene region splitting (default: enabled) |
| `-o`, `--output-dir PATH` | Output directory (default: `PyHIV_results`) |
| `-j`, `--n-jobs INTEGER` | Number of parallel jobs (default: all CPUs) |
| `-v`, `--verbose` | Detailed output |
| `-q`, `--quiet` | Suppress non-error output |

### 💼 Common Use Cases

**Full analysis with subtyping and splitting:**
```bash
pyhiv run data/sequences/
```

**Alignment only:**
```bash
pyhiv run data/sequences/ --no-subtyping --no-splitting
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
- `final_table.tsv` - Summary with sequence IDs, references, subtypes, and gene regions
- `best_alignment_*.fasta` - Best alignment for each sequence
- Gene-specific folders (when `--splitting` is enabled) with extracted regions

### 🆘 Getting Help

```bash
pyhiv --help           # Show all commands
pyhiv run --help       # Show options for run command
pyhiv --version        # Show version
```

For comprehensive documentation, see [CLI_README.md](CLI.md).

---

## 🗂️ Citation

Manuscript in preparation. Please cite this repository if you use PyHIV in your research.

---

## 🧾 License

This project is licensed under the MIT License — see the LICENSE
 file for details.

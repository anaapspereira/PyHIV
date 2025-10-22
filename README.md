# 🧬 PyHIV: A Python Package for Local HIV‑1 Sequence Alignment, Subtyping and Gene Splitting

<div align="center">

[![CI](https://github.com/anaapspereira/PyHIV/actions/workflows/ci.yml/badge.svg)](https://github.com/anaapspereira/PyHIV/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/anaapspereira/PyHIV/branch/main/graph/badge.svg)](https://codecov.io/gh/anaapspereira/PyHIV)
[![Python Version](https://img.shields.io/pypi/pyversions/pyhiv-tools)](https://pypi.org/project/pyhiv-tools/)
[![OS Supported](https://img.shields.io/badge/OS-Windows%20%7C%20Linux%20%7C%20macOS-blue)](https://pypi.org/project/pyhiv-tools/)

[![PyPI version](https://img.shields.io/pypi/v/pyhiv-tools)](https://pypi.org/project/pyhiv-tools/)
[![Documentation Status](https://readthedocs.org/projects/pyhiv/badge/?version=latest)](https://pyhiv.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub issues](https://img.shields.io/github/issues/anaapspereira/PyHIV)](https://github.com/anaapspereira/PyHIV/issues)

</div>

---

## 📋 Table of Contents

- [Overview](#-overview)
- [How It Works](#-how-it-works)
- [Installation](#-installation)
- [Getting Started](#-getting-started)
- [Command Line Interface](#-command-line-interface)
- [Output Structure](#-output-structure)
- [Citation](#-citation)
- [License](#-license)

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
```pgsql
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

### From PyPI (Recommended)

```bash
pip install pyhiv-tools
```

### From Source

```bash
git clone https://github.com/anaapspereira/PyHIV.git
cd PyHIV
pip install -e .
```

### Development Installation

```bash
git clone https://github.com/anaapspereira/PyHIV.git
cd PyHIV
pip install -e ".[dev]"
```

### Requirements

- Python 3.10+
- pandas
- biopython
- pyfamsa
- click
- matplotlib
## 🚀 Getting Started

### Quick Start (CLI)

The easiest way to use PyHIV is through the command line:

```bash
# Install PyHIV
pip install pyhiv-tools

# Run analysis on your sequences
pyhiv run /path/to/your/fasta/files

# Check results
ls PyHIV_results/
```

### Python API Usage

```python
from pyhiv import PyHIV

PyHIV(
    fastas_dir="path/to/fasta/files",
    subtyping=True,
    splitting=True,
    output_dir="results_folder",
    n_jobs=4,
    reporting=True
)
```

### Parameters:

| Parameter    | Type   | Default           | Description                                                                |
| ------------ | ------ |-------------------| -------------------------------------------------------------------------- |
| `fastas_dir` | `str`  | *Required*        | Directory containing user FASTA files.                                     |
| `subtyping`  | `bool` | `True`            | Aligns against subtype reference genomes. If `False`, aligns only to HXB2. |
| `splitting`  | `bool` | `True`            | Splits aligned sequences into gene regions.                                |
| `output_dir` | `str`  | `"PyHIV_results"` | Output directory for results.                                              |
| `n_jobs`     | `int`  | `None`            | Number of parallel jobs for alignment.                                     |
| `reporting`  | `bool` | `True`            | Generates PDF report with sequence visualizations.                         |

### 📂 Output Structure

After running PyHIV, your output directory (default: PyHIV_results/) will contain:

```
PyHIV_results/
│
├── best_alignment_<sequence>.fasta     # Alignment to best reference
├── final_table.tsv                     # Summary of results
├── PyHIV_report_all_sequences.pdf     # PDF report (if reporting=True)
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

PyHIV provides a comprehensive command-line interface for HIV-1 sequence analysis.

### 🚀 Basic Commands

```bash
# Basic usage - process all FASTA files in a directory
pyhiv run sequences/

# With custom output directory
pyhiv run sequences/ -o my_results/

# Parallel processing with 8 jobs
pyhiv run sequences/ -j 8

# Validate input files before processing
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

**Alignment only (no subtyping or splitting):**
```bash
pyhiv run data/sequences/ --no-subtyping --no-splitting
```

**Subtyping without gene splitting:**
```bash
pyhiv run data/sequences/ --no-splitting
```

**Parallel processing for large datasets:**
```bash
pyhiv run data/sequences/ -j 8 -o results/batch1/
```

**Validation before processing:**
```bash
pyhiv validate data/sequences/
```

### 🆘 Getting Help

```bash
pyhiv --help           # Show all commands
pyhiv run --help       # Show options for run command
pyhiv validate --help # Show validation options
pyhiv --version        # Show version
```

For comprehensive CLI documentation, see [CLI_README.md](CLI_README.md).

---

## 🗂️ Citation

If you use PyHIV in your research, please cite:

```bibtex
@software{pyhiv2024,
  title={PyHIV: A Python Package for Local HIV-1 Sequence Alignment, Subtyping and Gene Splitting},
  author={Santos-Pereira, Ana},
  year={2024},
  url={https://github.com/anaapspereira/PyHIV},
  license={MIT}
}
```

**Note:** Manuscript in preparation. Please cite this repository if you use PyHIV in your research.

---

## 🤝 Contributing

### Reporting Issues

Please report bugs and request features through [GitHub Issues](https://github.com/anaapspereira/PyHIV/issues).

---

## 📚 Documentation

- **Full Documentation**: [https://pyhiv.readthedocs.io/](https://pyhiv.readthedocs.io/)
- **CLI Reference**: [CLI_README.md](CLI_README.md)
- **API Reference**: Available in the documentation

---

## 🧾 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

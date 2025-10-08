# 🧬 PyHIV: A Python Package for Local HIV‑1 Sequence Alignment, Subtyping and Gene Splitting

![CI](https://github.com/anaapspereira/PyHIV/actions/workflows/ci.yml/badge.svg)
[![Coverage](https://img.shields.io/codecov/c/github/anaapspereira/PyHIV/main.svg)](https://codecov.io/gh/anaapspereira/PyHIV)
![Python Version](https://img.shields.io/pypi/pyversions/PyHIV)


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

You can install PyHIV using pip:

```bash
pip install pihiv-tools
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
    n_jobs=4
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

## 🗂️ Citation

Manuscript in preparation. Please cite this repository if you use PyHIV in your research.

---

## 🧾 License

This project is licensed under the MIT License — see the LICENSE
 file for details.

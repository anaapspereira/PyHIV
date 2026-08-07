# PyHIV: A Python Package for Local HIV-1 Sequence Alignment, Subtyping, and Gene Splitting

<div align="center">

[![PyPI version](https://img.shields.io/pypi/v/pyhiv-tools)](https://pypi.org/project/pyhiv-tools/)
[![Python Version](https://img.shields.io/pypi/pyversions/pyhiv-tools)](https://pypi.org/project/pyhiv-tools/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation Status](https://readthedocs.org/projects/pyhiv/badge/?version=latest)](https://pyhiv.readthedocs.io/)

</div>

---

## Overview

PyHIV is a Python package that aligns HIV nucleotide sequences against reference genomes to determine the most similar subtype and optionally split the aligned sequences into gene regions.

**Key Features:**
- 🧬 **Local HIV-1 sequence alignment** against reference genomes
- 🏷️ **Automated subtyping** with comprehensive reference database
- ✂️ **Gene region splitting** (gag, pol, env, etc.)
- 📊 **Detailed reporting** with summary tables and visualizations
- ⚡ **Parallel processing** for efficient analysis
- 🖥️ **Command-line interface** for easy integration

## Installation

```bash
pip install pyhiv-tools
```

## Quick Start

### Command Line Interface

```bash
# Basic usage
pyhiv run /path/to/fasta/files

# With custom options
pyhiv run /path/to/fasta/files -o results/ -j 4 -v

# Validate inputs first
pyhiv validate /path/to/fasta/files
```

### Python API

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
    reference_groups="M"
)
```

When `subtyping=False`, any active splitting mode uses HXB2 coordinates.

## What PyHIV Produces

- **Best reference alignment** per sequence
- **Subtype and reference metadata**
- **Ranked top 3 closest HIV-1 subtypes**
- **Gene-region–specific FASTA files** (optional)
- **Final summary table** (`final_table.tsv`)
- **PDF reports** with sequence visualizations (optional)

## Output Structure

```
PyHIV_results/
├── final_table.tsv                     # Summary of results
├── best_alignment_<sequence>.fasta     # Alignment to best reference
├── PyHIV_report_all_sequences.pdf     # PDF report (if enabled)
├── gag/                               # Gene regions (if splitting enabled)
│   ├── <sequence>_gag.fasta
│   └── ...
├── pol/
│   ├── <sequence>_pol.fasta
│   └── ...
└── env/
    ├── <sequence>_env.fasta
    └── ...
```

## Requirements

- Python 3.10+
- pandas
- biopython
- edlib
- pyfamsa
- click
- matplotlib
- parasail (optional, only for `--alignment-tool parasail-NW`; install with `pip install pyhiv-tools[parasail]`)

PyHIV supports `edlib-HW` (default), `parasail-NW`/`parasail`, `PyFamsa`, and `MAFFT` as alignment tool names. `parasail-NW` requires the optional `parasail` extra. `MAFFT` requires an external `mafft` executable. Before final alignment, PyHIV ranks references using query/reference k-mer containment and aligns only the top candidates by default. Use `reference_top_k=0` to keep the original all-reference strategy. By default, subtyping uses group M references from `reference_fastas`, selected through the `group` column in `sequences_with_locations.tsv`; set `reference_groups="M,N,O,P"` to include groups N, O, and P. The final table and PDF report include the best `Group`/`Subtype` call plus a ranked `Closest Subtypes` field with the top 3 closest unique group/subtype calls.

PyHIV resolves MAFFT from `PYHIV_MAFFT_BIN`, then `mafft` on `PATH`.

Input sequences longer than 12000 nucleotides are skipped with this warning: `The submitted sequence is longer than the HIV-1 genome.`

## Documentation

- **Full Documentation**: [https://pyhiv.readthedocs.io/](https://pyhiv.readthedocs.io/)
- **CLI Reference**: Available in the package documentation
- **GitHub Repository**: [https://github.com/anaapspereira/PyHIV](https://github.com/anaapspereira/PyHIV)

## Citation

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

## License

This project is licensed under the MIT License. See [LICENSE](https://github.com/anaapspereira/PyHIV/blob/main/LICENSE) file for details.

## Project Links

- **Source Code**: [https://github.com/anaapspereira/PyHIV](https://github.com/anaapspereira/PyHIV)
- **Issues**: [https://github.com/anaapspereira/PyHIV/issues](https://github.com/anaapspereira/PyHIV/issues)
- **PyPI Package**: [https://pypi.org/project/pyhiv-tools/](https://pypi.org/project/pyhiv-tools/)

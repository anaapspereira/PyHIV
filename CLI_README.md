# 🧬 PyHIV Command Line Interface

<div align="center">

A comprehensive command-line interface for HIV-1 sequence alignment, subtyping, and gene region splitting.

[![PyPI version](https://img.shields.io/pypi/v/pyhiv-tools)](https://pypi.org/project/pyhiv-tools/)
[![Python Version](https://img.shields.io/pypi/pyversions/pyhiv-tools)](https://pypi.org/project/pyhiv-tools/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

</div>

---

## 📋 Table of Contents

- [Installation](#-installation)
- [Getting Started](#-getting-started)
- [Commands](#-commands)
- [Options](#-options)
- [Usage Examples](#-usage-examples)
- [Input Requirements](#-input-requirements)
- [Output Structure](#-output-structure)
- [Advanced Usage](#-advanced-usage)
- [Troubleshooting](#-troubleshooting)
- [Performance Tips](#-performance-tips)
- [Contributing](#-contributing)
- [License](#-license)
- [Citation](#-citation)

---

## 📦 Installation

### From PyPI (Recommended)

```bash
pip install pyhiv-tools
```

This installs the default `edlib-HW` alignment backend. Two other backends are
available as optional extras:

```bash
pip install pyhiv-tools[parasail]  # parasail-NW backend
pip install pyhiv-tools[famsa]     # PyFamsa backend (GPL-3.0 licensed, kept out of the default MIT install)
```

`MAFFT` is also available as a backend, but requires a separate `mafft` executable on `PATH` (not installed via pip).

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

### Verify Installation

```bash
pyhiv --version
```

**Expected output:**
```
pyhiv, version 0.1.0
```

## 🚀 Getting Started

### Quick Start

Run PyHIV with default settings:

```bash
pyhiv run /path/to/fastas/
```

This will:
- ✅ Align sequences with reference genomes
- ✅ Perform HIV-1 subtyping
- ✅ Split sequences into gene regions
- ✅ Save results to `PyHIV_results/`

### First Time Setup

1. **Prepare your sequences**: Place FASTA files in a directory
2. **Validate input**: Check your files before processing
   ```bash
   pyhiv validate /path/to/fastas/
   ```
3. **Run analysis**: Process your sequences
   ```bash
   pyhiv run /path/to/fastas/ -v
   ```
4. **Check results**: Review the output
   ```bash
   ls PyHIV_results/
   cat PyHIV_results/final_table.tsv
   ```

## 🧭 Commands

### `pyhiv run` - Main Processing Command

Processes HIV-1 sequences for alignment, subtyping, and gene splitting.

```bash
pyhiv run [OPTIONS] FASTAS_DIR
```

**Arguments:**
- `FASTAS_DIR`: Directory containing input FASTA files (required)

**What it does:**
- 🔍 Reads all FASTA files from the specified directory
- 🧬 Aligns sequences against reference genomes
- 🏷️ Determines HIV-1 subtypes (if enabled)
- ✂️ Splits sequences into gene regions (if enabled)
- 📊 Generates summary table and reports

### `pyhiv validate` - Input Validation

Validates input directory and files without processing.

```bash
pyhiv validate FASTAS_DIR
```

**What it checks:**
- ✅ Directory exists and is readable
- ✅ FASTA files are present
- ✅ File format validation
- 📋 Lists found files

## ⚙️ Options

### 🔧 Processing Options

Control the analysis pipeline behavior.

| Option | Default | Description |
|--------|---------|-------------|
| `--subtyping BOOL` | `true` | Enable/disable HIV-1 subtyping against reference genomes |
| `--splitting TEXT` | `true` | Splitting mode: `true`/`subtype`, `hxb2`/`reference`, or `false`/`none`. With `--subtyping true`, `true` uses the subtype reference when it has features and falls back to HXB2 when features are missing. If `--subtyping false`, active splitting uses HXB2 |

**Examples:**
```bash
# Full analysis (default)
pyhiv run sequences/

# Alignment only
pyhiv run sequences/ --subtyping false --splitting false

# Subtyping without splitting
pyhiv run sequences/ --splitting false

# Subtyping with HXB2-based splitting
pyhiv run sequences/ --splitting hxb2
```

### 📁 Output Options

Control where and how results are saved.

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output-dir PATH` | `PyHIV_results` | Output directory for results |

**Examples:**
```bash
# Custom output directory
pyhiv run sequences/ -o my_analysis_results/

# Timestamped results
pyhiv run sequences/ -o results_$(date +%Y%m%d)/
```

### ⚡ Performance Options

Optimize processing speed and resource usage.

| Option | Default | Description |
|--------|---------|-------------|
| `-j`, `--n-jobs INTEGER` | All CPUs | Number of parallel jobs for alignment |
| `--alignment-tool [edlib-HW|MAFFT|parasail-NW|PyFamsa|parasail]` | `edlib-HW` | Alignment backend |
| `--kmer-size INTEGER` | `15` | K-mer size used to prefilter candidate references |
| `--reference-top-k INTEGER` | `30` | Number of top k-mer ranked references to align; use `0` to align all references |
| `--reference-groups TEXT` | `M` | Comma-separated HIV-1 reference groups used for subtyping; use `M,N,O,P` to include groups N, O, and P |

**Examples:**
```bash
# Use all available CPUs (default)
pyhiv run sequences/

# Limit to 4 cores
pyhiv run sequences/ -j 4

# Use PyFamsa instead of the default edlib-HW
pyhiv run sequences/ --alignment-tool PyFamsa

# Align only the 10 references with the highest k-mer containment
pyhiv run sequences/ --reference-top-k 10

# Include reference groups N, O, and P in addition to the default group M
pyhiv run sequences/ --reference-groups M,N,O,P

# Single-threaded processing
pyhiv run sequences/ -j 1
```

`edlib-HW` is the default and projects alignments onto full-reference genome coordinates. `parasail-NW`/`parasail`, `PyFamsa`, and `MAFFT` remain available through `--alignment-tool`. `edlib` is installed with PyHIV. `parasail` and `PyFamsa` are optional extras. `parasail` — install with `pip install pyhiv-tools[parasail]` — since it has no prebuilt wheel on some platforms (e.g. macOS on Apple Silicon). `PyFamsa` — install with `pip install pyhiv-tools[famsa]` — since `pyfamsa` is GPL-3.0 licensed and is kept out of the default (MIT) install. `MAFFT` requires an external `mafft` executable. PyHIV resolves MAFFT from `PYHIV_MAFFT_BIN`, then `mafft` on `PATH`.

Before final alignment, PyHIV ranks references using query/reference k-mer containment and aligns only the top candidates by default. Use `--reference-top-k 0` to keep the original all-reference strategy. By default, subtyping uses group M references from `reference_fastas`, selected through the `group` column in `sequences_with_locations.tsv`; use `--reference-groups M,N,O,P` to include groups N, O, and P.

Sequences longer than 12000 nucleotides are skipped with this warning: `The submitted sequence is longer than the HIV-1 genome.`

### 📺 Display Options

Control output verbosity and information display.

| Option | Description |
|--------|-------------|
| `-v`, `--verbose` | Enable detailed output with timing information |
| `-q`, `--quiet` | Suppress all non-error output |
| `--version` | Show version and exit |
| `--help` | Show help message and exit |

**Examples:**
```bash
# Verbose output for debugging
pyhiv run sequences/ -v

# Quiet mode for scripting
pyhiv run sequences/ -q

# Get help
pyhiv run --help
```

## 💼 Usage Examples

### Basic Usage

**Default processing:**
```bash
pyhiv run sequences/
```

**Custom output directory:**
```bash
pyhiv run sequences/ -o my_results/
```

**Parallel processing with 8 jobs:**
```bash
pyhiv run sequences/ -j 8
```

### Advanced Options

**Alignment only (no subtyping or splitting):**
```bash
pyhiv run sequences/ --subtyping false --splitting false
```

**Subtyping without gene splitting:**
```bash
pyhiv run sequences/ --splitting false
```

**Subtyping with HXB2-based splitting:**
```bash
pyhiv run sequences/ --splitting hxb2
```

**Verbose output with timing:**
```bash
pyhiv run sequences/ -v
```

**Quiet mode for scripting:**
```bash
pyhiv run sequences/ -q
```

### Validation

**Check input files before processing:**
```bash
pyhiv validate sequences/
```

### Pipeline Examples

**Complete workflow:**
```bash
# 1. Validate inputs
pyhiv validate data/raw_sequences/

# 2. Process with verbose output
pyhiv run data/raw_sequences/ -o results/run1/ -v

# 3. Process subset without splitting
pyhiv run data/subset/ -o results/run2/ --splitting false -j 4
```

## 📥 Input Requirements

### Supported Formats

PyHIV accepts FASTA files with the following extensions:
- `.fasta`
- `.fa`
- `.fna` (nucleic acid)
- `.ffn` (nucleotide coding regions)

### Directory Structure

```
sequences/
├── sample1.fasta
├── sample2.fa
├── sample3.fasta
└── subfolder/
    └── sample4.fasta
```

PyHIV recursively searches for FASTA files in all subdirectories.

## 📂 Output Structure

### Default Output (`PyHIV_results/`)

```
PyHIV_results/
├── final_table.tsv                    # Summary table
├── best_alignment_inputA_sample1.fasta # Best alignments
├── best_alignment_inputB_sample2.fasta
├── gag/                               # Gene regions (if --splitting)
│   ├── inputA_sample1_gag.fasta
│   └── inputB_sample2_gag.fasta
├── pol/
│   ├── inputA_sample1_pol.fasta
│   └── inputB_sample2_pol.fasta
├── env/
└── ...
```

### Output Files

#### `final_table.tsv`

Summary table with columns:

| Column | Description |
|--------|-------------|
| File Name | Source FASTA file name |
| Sequence | Input sequence ID |
| Reference | Best matching reference accession |
| Group | HIV-1 reference group |
| Subtype | HIV-1 subtype (if `--subtyping` enabled) |
| Closest Subtypes | Top 3 closest unique group/subtype calls by alignment score |
| Splitting Reference | Reference accession used for gene splitting |
| Most Matching Gene Region | Gene with most matches |
| Present Gene Regions | All detected gene regions |

**Example:**
```
File Name        Sequence    Reference    Group    Subtype    Closest Subtypes                                  Splitting Reference    Most Matching Gene Region    Present Gene Regions
sample_1.fasta   seq001      K03455       M        B          M:B (score=9120); M:C (score=8894); M:A1 (...)   K03455                 pol                          gag, pol, env
sample_2.fasta   seq002      AF004885     M        C          M:C (score=9015); M:B (score=8801); M:A1 (...)   AF004885               env                          pol, env
```

#### Alignment Files

- `best_alignment_<file_stem>_<sequence_id>.fasta`: Contains reference and query alignment
- `splitting_alignment_<file_stem>_<sequence_id>.fasta`: Contains the HXB2 alignment when `--splitting hxb2` is used with subtyping, or when `--splitting true` falls back to HXB2 because the selected subtype reference has no annotated features
- Format: Multi-FASTA with reference sequence and aligned query

#### Gene Region Files

When `--splitting` is enabled:
- Organized by gene (gag, pol, env, etc.)
- One file per sequence per gene
- Contains extracted gene region from alignment

## ⚡ Advanced Usage

### Performance Tuning

**Optimize for large datasets:**
```bash
# Use all CPUs
pyhiv run sequences/ -j -1

# Limit to 4 cores to avoid memory issues
pyhiv run sequences/ -j 4
```

**Memory considerations:**
- Each job loads reference sequences
- Reduce `-j` value if encountering memory errors
- Process in batches for very large datasets

### Batch Processing

```bash
# Process multiple directories
for dir in batch1/ batch2/ batch3/; do
    pyhiv run "$dir" -o "results_$(basename $dir)" -q
done
```

### Integration with Other Tools

**Export to CSV:**
```bash
pyhiv run sequences/ -o results/
python -c "import pandas as pd; df = pd.read_csv('results/final_table.tsv', sep='\t'); df.to_csv('results.csv', index=False)"
```

**Filter by subtype:**
```bash
pyhiv run sequences/ -o results/ -v
awk -F'\t' '$5 == "B"' results/final_table.tsv > subtype_B.tsv
```

## 🛠️ Troubleshooting

### Common Issues

**No FASTA files found:**
```
Error: No FASTA files found in the input directory.
```
- Check file extensions (must be .fasta, .fa, etc.)
- Verify directory path is correct
- Use `pyhiv validate` to diagnose

**Output directory exists:**
```
Warning: Output directory 'PyHIV_results' already exists. Files may be overwritten.
```
- Files will be overwritten
- Use `-o` to specify a different directory
- Or remove existing directory: `rm -rf PyHIV_results/`

**Memory errors:**
```
MemoryError: Unable to allocate array
```
- Reduce parallel jobs: `-j 2` or `-j 1`
- Process sequences in smaller batches
- Close other applications

**Import errors:**
```
Error: Could not import PyHIV module
```
- Verify installation: `pip list | grep pyhiv-tools`
- Reinstall: `pip install --force-reinstall pyhiv-tools`
- Check Python version compatibility

### Debug Mode

Enable verbose output for debugging:
```bash
pyhiv run sequences/ -v
```

This shows:
- Version information
- Number of input files
- Processing parameters
- Elapsed time
- Generated output files
- Full stack traces on errors

### Getting Help

```bash
# Show all commands
pyhiv --help

# Show help for specific command
pyhiv run --help
pyhiv validate --help
```

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error during processing |
| 130 | Interrupted by user (Ctrl+C) |

## 📈 Performance Tips

1. **Use validation first** - `pyhiv validate` is fast and catches input errors
2. **Adjust parallelism** - Start with default (all CPUs), reduce if memory is limited
3. **Disable unused features** - Use `--splitting false` if you only need alignments
4. **Batch processing** - For thousands of sequences, split into smaller batches
5. **SSD storage** - Use SSD for output directory to improve I/O performance

## 🤝 Contributing

Found a bug or have a feature request? Please open an issue on GitHub.

## 🧾 License

PyHIV is released under the MIT License. See LICENSE file for details.

## 🗂️ Citation

If you use PyHIV in your research, please cite:

```
Manuscript in preparation. Please cite this repository if you use PyHIV in your research.
```

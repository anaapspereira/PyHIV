# 🧬 PyHIV Command Line Interface

A comprehensive command-line interface for HIV-1 sequence alignment, subtyping, and gene region splitting.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Options](#options)
- [Usage Examples](#usage-examples)
- [Input Requirements](#input-requirements)
- [Output Structure](#output-structure)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)

## 📦 Installation

Install PyHIV using pip:

```bash
pip install pyhiv
```

Or install from source:

```bash
git clone https://github.com/yourusername/pyhiv.git
cd pyhiv
pip install -e .
```

Verify installation:

```bash
pyhiv --version
```

## 🚀 Getting Started

Run PyHIV with default settings:

```bash
pyhiv run /path/to/fastas/
```

This will:
- Align sequences with reference genomes
- Perform HIV-1 subtyping
- Split sequences into gene regions
- Save results to `PyHIV_results/`

## 🧭 Commands

### `pyhiv run`

Main command to process HIV-1 sequences.

```bash
pyhiv run [OPTIONS] FASTAS_DIR
```

**Arguments:**
- `FASTAS_DIR`: Directory containing input FASTA files (required)

### `pyhiv validate`

Validate input directory without processing.

```bash
pyhiv validate FASTAS_DIR
```

Checks:
- Directory exists and is readable
- FASTA files are present
- Lists found files (up to 10)

## ⚙️ Options

### Processing Options

| Option | Default | Description |
|--------|---------|-------------|
| `--subtyping` / `--no-subtyping` | `--subtyping` | Enable/disable HIV-1 subtyping |
| `--splitting` / `--no-splitting` | `--splitting` | Enable/disable gene region splitting |

### Output Options

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output-dir PATH` | `PyHIV_results` | Output directory for results |

### Performance Options

| Option | Default | Description |
|--------|---------|-------------|
| `-j`, `--n-jobs INTEGER` | All CPUs | Number of parallel jobs |

### Display Options

| Option | Description |
|--------|-------------|
| `-v`, `--verbose` | Enable detailed output |
| `-q`, `--quiet` | Suppress all non-error output |
| `--version` | Show version and exit |
| `--help` | Show help message and exit |

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
pyhiv run sequences/ --no-subtyping --no-splitting
```

**Subtyping without gene splitting:**
```bash
pyhiv run sequences/ --no-splitting
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

**Example output:**
```
✓ Found 15 FASTA file(s)

Files:
  • sequence1.fasta
  • sequence2.fasta
  • sequence3.fa
  ...
```

### Pipeline Examples

**Complete workflow:**
```bash
# 1. Validate inputs
pyhiv validate data/raw_sequences/

# 2. Process with verbose output
pyhiv run data/raw_sequences/ -o results/run1/ -v

# 3. Process subset without splitting
pyhiv run data/subset/ -o results/run2/ --no-splitting -j 4
```

**Integration with shell scripts:**
```bash
#!/bin/bash
INPUT_DIR="sequences/"
OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"

# Validate first
if pyhiv validate "$INPUT_DIR"; then
    echo "Validation passed, starting processing..."
    pyhiv run "$INPUT_DIR" -o "$OUTPUT_DIR" -j 8
else
    echo "Validation failed!"
    exit 1
fi
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

### File Requirements

- Valid FASTA format
- HIV-1 sequences (DNA or RNA)
- Sequence IDs should be unique

## 📂 Output Structure

### Default Output (`PyHIV_results/`)

```
PyHIV_results/
├── final_table.tsv                    # Summary table
├── best_alignment_sample1.fasta       # Best alignments
├── best_alignment_sample2.fasta
├── gag/                               # Gene regions (if --splitting)
│   ├── sample1_gag.fasta
│   └── sample2_gag.fasta
├── pol/
│   ├── sample1_pol.fasta
│   └── sample2_pol.fasta
├── env/
└── ...
```

### Output Files

#### `final_table.tsv`

Summary table with columns:

| Column | Description |
|--------|-------------|
| Sequence | Input sequence ID |
| Reference | Best matching reference accession |
| Subtype | HIV-1 subtype (if `--subtyping` enabled) |
| Most Matching Gene Region | Gene with most matches |
| Present Gene Regions | All detected gene regions |

**Example:**
```tsv
Sequence    Reference    Subtype    Most Matching Gene Region    Present Gene Regions
seq001      K03455       B          pol                          gag, pol, env
seq002      AF004885     C          env                          pol, env
```

#### Alignment Files

- `best_alignment_<sequence_id>.fasta`: Contains reference and query alignment
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
awk -F'\t' '$3 == "B"' results/final_table.tsv > subtype_B.tsv
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
- Verify installation: `pip list | grep pyhiv`
- Reinstall: `pip install --force-reinstall pyhiv`
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
3. **Disable unused features** - Use `--no-splitting` if you only need alignments
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

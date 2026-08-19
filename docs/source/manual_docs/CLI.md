# 🧬 PyHIV Command Line Interface

A comprehensive command-line interface for HIV-1 sequence alignment, subtyping, and gene region splitting.

## 📦 Installation

Install PyHIV using pip:

```bash
pip install pyhiv-tools
```

Or install from source:

```bash
git clone https://github.com/anaapspereira/pyhiv.git
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

### `pyhiv update reference-dataset`

Refresh the packaged HIV-1 reference dataset (`reference_fastas/` and `sequences_with_locations.tsv`) from LANL's latest reference alignment and CRF database, optionally refreshing GenBank feature annotations. Also available as the top-level alias `pyhiv update-reference-dataset`.

```bash
pyhiv update reference-dataset [OPTIONS]
```

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--reference-dir PATH` | `REFERENCE_GENOMES_DIR` or packaged directory | Reference dataset root to update |
| `--email TEXT` | `NCBI_EMAIL` env var | Email sent to NCBI E-utilities (required by `--refresh-features`) |
| `--ncbi-api-key TEXT` | `NCBI_API_KEY` env var | Optional NCBI API key |
| `--refresh-features` | off | Refetch GenBank feature locations for references whose features are `None` |
| `--dry-run` | off | Check for available updates without writing FASTA or TSV files |
| `-y`, `--yes` | off | Skip the confirmation prompt |

Without `--yes` or `--dry-run`, the command asks for confirmation before overwriting the reference dataset. `--refresh-features` requires an NCBI email, passed via `--email` or the `NCBI_EMAIL` environment variable.

```bash
# Preview available updates without changing anything
pyhiv update reference-dataset --dry-run

# Apply the update, skipping the confirmation prompt
pyhiv update reference-dataset --yes

# Also refresh missing GenBank feature annotations
pyhiv update reference-dataset --refresh-features --email you@example.com --yes
```

## ⚙️ Options

### Processing Options

| Option | Default | Description |
|--------|---------|-------------|
| `--subtyping BOOL` | `true` | Enable/disable HIV-1 subtyping |
| `--splitting TEXT` | `true` | Splitting mode: `true`/`subtype`, `hxb2`/`reference`, or `false`/`none`. With `--subtyping true`, `true` uses the subtype reference when it has features and falls back to HXB2 when features are missing. If `--subtyping false`, active splitting uses HXB2 |

### Output Options

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output-dir PATH` | `PyHIV_results` | Output directory for results |

### Performance Options

| Option | Default | Description |
|--------|---------|-------------|
| `-j`, `--n-jobs INTEGER` | All CPUs | Number of parallel jobs |
| `--alignment-tool [edlib-HW|MAFFT|parasail-NW|PyFamsa|parasail]` | `edlib-HW` | Alignment backend |
| `--kmer-size INTEGER` | `15` | K-mer size used to prefilter candidate references |
| `--reference-top-k INTEGER` | `30` | Number of top k-mer ranked references to align; use `0` to align all references |
| `--reference-groups TEXT` | `M` | Comma-separated HIV-1 reference groups used for subtyping; use `M,N,O,P` or `all` to include groups N, O, and P |

`edlib-HW` is the default alignment backend. `parasail-NW`/`parasail`, `PyFamsa`, and `MAFFT` are available through `--alignment-tool`. `parasail-NW` requires the optional `parasail` extra — install with `pip install pyhiv-tools[parasail]` — since it has no prebuilt wheel on some platforms (e.g. macOS on Apple Silicon). PyHIV checks the selected tool is actually usable before processing starts, and exits immediately with an install hint if it isn't (e.g. missing `parasail`, or no `mafft` executable on `PATH`) rather than failing per-reference partway through a run.

Before final alignment, PyHIV ranks references using query/reference k-mer containment and aligns only the top candidates by default. Use `--reference-top-k 0` to keep the original all-reference strategy. By default, subtyping uses group M references from `reference_fastas`, selected through the `group` column in `sequences_with_locations.tsv`; use `--reference-groups M,N,O,P` (or `--reference-groups all`) to include groups N, O, and P.

MAFFT can be configured with `PYHIV_MAFFT_BIN` or discovered as `mafft` on `PATH`.

Sequences longer than 12000 nucleotides are skipped with this warning: `The submitted sequence is longer than the HIV-1 genome.`

### Progress & Reporting Options

| Option | Default | Description |
|--------|---------|-------------|
| `--progress` / `--no-progress` | `--progress` | Show a terminal progress bar while processing input sequences |
| `--reporting` / `--no-reporting` | `--reporting` | Generate a PDF report (`PyHIV_report_all_sequences.pdf`) with per-sequence visualizations |

### Display Options

| Option | Description |
|--------|-------------|
| `-v`, `--verbose` | Enable detailed output |
| `-q`, `--quiet` | Suppress all non-error output |
| `--version` | Show version and exit |
| `--help` | Show help message and exit |

## 🔧 Environment Variables

| Variable | Used by | Description |
|----------|---------|-------------|
| `REFERENCE_GENOMES_DIR` | Every command (`run`, `validate`, `update reference-dataset`) | Overrides the reference dataset root for the whole tool, not just the update command. Defaults to the packaged `reference_genomes` directory. Point this at a directory populated by `pyhiv update reference-dataset --reference-dir <path>` to use a custom or offline reference dataset. |
| `NCBI_EMAIL` | `pyhiv update reference-dataset` | Email sent to NCBI E-utilities; equivalent to `--email`. Required when using `--refresh-features`. |
| `NCBI_API_KEY` | `pyhiv update reference-dataset` | Optional NCBI API key; equivalent to `--ncbi-api-key`. |
| `PYHIV_MAFFT_BIN` | `pyhiv run --alignment-tool MAFFT` | Path to the `mafft` executable. Falls back to `mafft` on `PATH` if unset. |

An explicit CLI flag always takes priority over its corresponding environment variable.

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
pyhiv run data/subset/ -o results/run2/ --splitting false -j 4
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

PyHIV recursively searches for FASTA files in all subdirectories. The run's own `--output-dir` is always excluded from this search, so pointing `--output-dir` at a subfolder of `FASTAS_DIR` is safe — a rerun won't re-ingest FASTA files it generated on a previous run (e.g. `best_alignment_*.fasta`, gene-region files).

### File Requirements

- Valid FASTA format
- HIV-1 sequences (DNA or RNA)
- Sequence IDs should be unique

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
| Subtype Score Warning | `Low score margin: review top 3 subtype matches` when the top two subtype scores are within 1% of each other, otherwise empty |
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

#### PDF Report

When `--reporting` is enabled (default), PyHIV writes `PyHIV_report_all_sequences.pdf` to the output directory with one page per processed sequence, summarizing its alignment, subtyping, and gene-splitting results. Disable it with `--no-reporting`.

## ⚡ Advanced Usage

### Performance Tuning

**Optimize for large datasets:**
```bash
# Use all CPUs (the default when -j is omitted)
pyhiv run sequences/

# Limit to 4 cores to avoid memory issues
pyhiv run sequences/ -j 4
```

**Memory considerations:**
- `-j` sizes one worker pool that's created once per run and reused for every sequence, not spawned per sequence, but each worker process still loads its own copy of the reference sequences
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
pyhiv update reference-dataset --help
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

"""
PyHIV Command Line Interface
"""
import click
from pathlib import Path
import sys
import time
import warnings
from pyhiv import __version__, normalize_reference_groups, normalize_splitting_mode
from pyhiv.align import (
    ALIGNMENT_TOOL_CHOICES,
    DEFAULT_ALIGNMENT_TOOL,
    DEFAULT_KMER_SIZE,
    DEFAULT_REFERENCE_TOP_K,
)
import logging
from pyhiv.config import get_reference_base_dir

SUPPORTED_FASTA_EXTENSIONS = {'.fasta', '.fa', '.fna', '.ffn'}


def validate_n_jobs(ctx, param, value):
    """Validate that n_jobs is positive if provided."""
    if value is not None and value < 1:
        raise click.BadParameter('must be at least 1')
    return value


def validate_positive(ctx, param, value):
    """Validate that an integer option is positive if provided."""
    if value is not None and value < 1:
        raise click.BadParameter('must be at least 1')
    return value


def validate_reference_top_k(ctx, param, value):
    """Validate candidate reference limit. Zero disables k-mer shortlisting."""
    if value is not None and value < 0:
        raise click.BadParameter('must be at least 0')
    return value


def validate_reference_groups(ctx, param, value):
    """Validate comma-separated HIV-1 reference groups."""
    if value is None:
        return None
    try:
        return ",".join(normalize_reference_groups(value))
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc


def validate_splitting(ctx, param, value):
    """Validate splitting mode while preserving true/false compatibility."""
    if isinstance(value, bool):
        return value

    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "on"}:
        return True
    if text in {"false", "0", "no", "off", "none"}:
        return False
    try:
        return normalize_splitting_mode(text)
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc


def count_fasta_files(directory):
    """Count FASTA files in the input directory."""
    return sum(1 for f in Path(directory).rglob('*') if f.is_file() and f.suffix.lower() in SUPPORTED_FASTA_EXTENSIONS)


@click.command()
@click.version_option(version=__version__, prog_name="PyHIV")
@click.argument(
    'fastas_dir',
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path, readable=True),
    required=True
)
@click.option(
    '--subtyping',
    type=click.BOOL,
    default=True,
    show_default=True,
    help='Enable or disable HIV-1 subtyping. Use true/false. When enabled, aligns with reference genomes for subtype identification.'
)
@click.option(
    '--splitting',
    type=str,
    default="true",
    show_default=True,
    callback=validate_splitting,
    help='Splitting mode: true/subtype, hxb2/reference, or false/none. With subtyping true, true falls back to HXB2 when the subtype reference has no features. If subtyping is false, active splitting uses HXB2.'
)
@click.option(
    '-o', '--output-dir',
    type=click.Path(path_type=Path),
    default=None,
    help='Output directory for results. Defaults to "PyHIV_results" in the current directory.'
)
@click.option(
    '-j', '--n-jobs',
    type=int,
    default=None,
    callback=validate_n_jobs,
    help='Number of parallel jobs to run. If not specified, uses all available CPU cores.'
)
@click.option(
    '-v', '--verbose',
    is_flag=True,
    help='Enable verbose output.'
)
@click.option(
    '-q', '--quiet',
    is_flag=True,
    help='Suppress all non-error output.'
)
@click.option(
    '--reporting/--no-reporting',
    default=True,
    show_default=True,
    help='Enable or disable PDF report generation. When enabled, generates a PDF report with sequence visualizations.'
)
@click.option(
    '--alignment-tool',
    type=click.Choice(ALIGNMENT_TOOL_CHOICES, case_sensitive=False),
    default=DEFAULT_ALIGNMENT_TOOL,
    show_default=True,
    help='Alignment tool to use.'
)
@click.option(
    '--kmer-size',
    type=int,
    default=DEFAULT_KMER_SIZE,
    show_default=True,
    callback=validate_positive,
    help='K-mer size used to prefilter candidate references before alignment.'
)
@click.option(
    '--reference-top-k',
    type=int,
    default=DEFAULT_REFERENCE_TOP_K,
    show_default=True,
    callback=validate_reference_top_k,
    help='Number of top k-mer ranked references to align. Use 0 to align all references.'
)
@click.option(
    '--reference-groups',
    default=None,
    show_default="M when splitting is true; all references when splitting is false",
    callback=validate_reference_groups,
    help='Comma-separated HIV-1 reference groups used for subtyping. Use M,N,O,P to include all groups.'
)
def main(
    fastas_dir,
    subtyping,
    splitting,
    output_dir,
    n_jobs,
    verbose,
    quiet,
    reporting,
    alignment_tool,
    kmer_size,
    reference_top_k,
    reference_groups,
):
    """
    PyHIV: HIV-1 sequence alignment, subtyping, and gene region splitting tool.

    FASTAS_DIR: Directory containing input FASTA files to process.

    \b
    Examples:
        # Basic usage with default settings
        pyhiv /path/to/fastas/

        # Disable subtyping
        pyhiv /path/to/fastas/ --subtyping false

        # Custom output directory with 4 parallel jobs
        pyhiv /path/to/fastas/ -o results/ -j 4

        # Only alignment, no splitting
        pyhiv /path/to/fastas/ --splitting false

        # Subtyping with HXB2-based splitting
        pyhiv /path/to/fastas/ --splitting hxb2

        # Quiet mode (only show errors)
        pyhiv /path/to/fastas/ -q

        # Select an alignment tool
        pyhiv /path/to/fastas/ --alignment-tool parasail-NW
    """

    # Handle conflicting flags
    if verbose and quiet:
        raise click.UsageError("Cannot use --verbose and --quiet together")

    # Configure logging based on flags
    if quiet:
        logging_level = logging.ERROR
    elif verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.WARNING

    if verbose:
        warnings.filterwarnings("default")
    else:
        warnings.filterwarnings("ignore")

    logging.basicConfig(
        level=logging_level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.getLogger().setLevel(logging_level)

    # Set output directory
    output_path = output_dir or Path('PyHIV_results')

    # Check if output directory exists and warn user
    if output_path.exists() and not quiet:
        click.secho(f"Warning: Output directory '{output_path}' already exists. Files may be overwritten.",
                    fg='yellow', err=True)

    # Count input files
    num_files = count_fasta_files(fastas_dir)
    if num_files == 0:
        click.secho("Error: No FASTA files found in the input directory.", fg='red', err=True)
        sys.exit(1)

    if verbose:
        splitting_mode = normalize_splitting_mode(splitting, subtyping=subtyping)
        click.echo(f"PyHIV v{__version__}")
        click.echo(f"Input directory: {fastas_dir}")
        click.echo(f"Found {num_files} FASTA file(s)")
        click.echo(f"Subtyping: {'enabled' if subtyping else 'disabled'}")
        click.echo(f"Splitting: {splitting_mode}")
        click.echo(f"Alignment tool: {alignment_tool}")
        click.echo(f"K-mer size: {kmer_size}")
        click.echo(f"Reference top-k: {reference_top_k or 'all'}")
        click.echo(f"Reference groups: {reference_groups}")
        click.echo(f"Output directory: {output_path}")
        click.echo(f"Parallel jobs: {n_jobs or 'auto (all CPUs)'}")
        click.echo()
    elif not quiet:
        click.echo(f"Processing {num_files} FASTA file(s)...")

    start_time = time.time()

    try:
        from pyhiv import PyHIV

        PyHIV(
            fastas_dir=str(fastas_dir),
            subtyping=subtyping,
            splitting=splitting,
            output_dir=str(output_dir) if output_dir else None,
            n_jobs=n_jobs,
            reporting=reporting,
            alignment_tool=alignment_tool,
            kmer_size=kmer_size,
            reference_top_k=reference_top_k,
            reference_groups=reference_groups,
        )

        elapsed_time = time.time() - start_time

        if not quiet:
            click.secho(f"\n✓ Processing complete!", fg='green', bold=True)
            click.echo(f"Results saved to: {output_path}")
            click.echo(f"Time elapsed: {elapsed_time:.2f}s")

        # Show key output files
        if verbose:
            click.echo("\nGenerated files:")
            final_table = output_path / 'final_table.tsv'
            if final_table.exists():
                click.echo(f"  • {final_table}")

            # List some alignment files
            alignment_files = list(output_path.glob('best_alignment_*.fasta'))
            for af in alignment_files[:3]:
                click.echo(f"  • {af}")
            if len(alignment_files) > 3:
                click.echo(f"  • ... and {len(alignment_files) - 3} more alignment file(s)")
            
            # Show PDF report if generated
            if reporting:
                pdf_report = output_path / 'PyHIV_report_all_sequences.pdf'
                if pdf_report.exists(): # pragma: no cover
                    click.echo(f"  • {pdf_report}")

    except ImportError as e:
        click.secho(f"Error: Could not import PyHIV module: {e}", fg='red', err=True)
        sys.exit(1)
    except KeyboardInterrupt:
        click.secho("\n\nProcessing interrupted by user.", fg='yellow', err=True)
        sys.exit(130)
    except Exception as e:
        click.secho(f"Error during processing: {e}", fg='red', err=True)
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


@click.command('validate')
@click.argument(
    'fastas_dir',
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
)
def validate(fastas_dir):
    """Validate FASTA files in the input directory without processing."""
    num_files = count_fasta_files(fastas_dir)

    if num_files == 0:
        click.secho("✗ No FASTA files found.", fg='red')
        sys.exit(1)

    click.secho(f"✓ Found {num_files} FASTA file(s)", fg='green')

    # List files if not too many
    if num_files <= 10:
        files = []
        for ext in SUPPORTED_FASTA_EXTENSIONS:
            files.extend(Path(fastas_dir).rglob(f'*{ext}'))
        files = list({f.resolve(): f for f in files}.values())  # Remove duplicates, preserve Path objects
        click.echo("\nFiles:")
        for f in files:
            click.echo(f"  • {f.name}")


@click.group('update')
def update():
    """Update packaged PyHIV datasets."""
    pass


@click.command('reference-dataset')
@click.option(
    '--reference-dir',
    type=click.Path(file_okay=False, dir_okay=True, path_type=Path),
    default=None,
    help='Reference dataset root. Defaults to REFERENCE_GENOMES_DIR or the packaged reference_genomes directory.',
)
@click.option(
    '--email',
    envvar='NCBI_EMAIL',
    default=None,
    help='Email sent to NCBI E-utilities. Can also be set with NCBI_EMAIL.',
)
@click.option(
    '--ncbi-api-key',
    envvar='NCBI_API_KEY',
    default=None,
    help='Optional NCBI API key. Can also be set with NCBI_API_KEY.',
)
@click.option(
    '--refresh-features',
    is_flag=True,
    help='Refetch GenBank feature locations for references whose features are None.',
)
@click.option(
    '--dry-run',
    is_flag=True,
    help='Check the update without writing FASTA or TSV files.',
)
@click.option(
    '-y',
    '--yes',
    is_flag=True,
    help='Run without asking for confirmation.',
)
def update_reference_dataset(reference_dir, email, ncbi_api_key, refresh_features, dry_run, yes):
    """Update HIV-1 reference FASTAs and sequences_with_locations.tsv."""
    from pyhiv.loading.reference_update import update_reference_dataset as run_update

    base_dir = reference_dir or get_reference_base_dir()

    if refresh_features and not email:
        raise click.UsageError(
            "--refresh-features needs an NCBI email. Pass --email or set NCBI_EMAIL."
        )

    if not yes and not dry_run:
        click.confirm(
            (
                f"This will update the reference dataset in '{base_dir}' and "
                "rewrite sequences_with_locations.tsv. Continue?"
            ),
            abort=True,
        )

    try:
        result = run_update(
            base_dir=base_dir,
            email=email,
            ncbi_api_key=ncbi_api_key,
            dry_run=dry_run,
            refresh_features=refresh_features,
        )
    except Exception as exc:
        click.secho(f"Error updating reference dataset: {exc}", fg='red', err=True)
        sys.exit(1)

    updated_rows = getattr(result, "updated_rows", 0)
    lanl_fasta_files_removed = getattr(result, "lanl_fasta_files_removed", 0)
    pending_ncbi_references = getattr(result, "pending_ncbi_references", 0)
    pending_ncbi_accessions = getattr(result, "pending_ncbi_accessions", ())

    if result.up_to_date:
        click.secho("Reference dataset is up to date.", fg='green')
    elif dry_run:
        click.secho(
            (
                "Reference dataset update available: "
                f"LANL alignment year {result.lanl_alignment_year or 'unknown'}; "
                f"LANL alignment {'would be refreshed' if result.lanl_alignment_updated else 'is current'}; "
                f"{result.lanl_fasta_files_added} LANL FASTA file(s) would be added; "
                f"{result.lanl_fasta_files_updated} LANL FASTA file(s) would be updated; "
                f"{lanl_fasta_files_removed} stale LANL FASTA file(s) would be removed; "
                f"{result.added_sequences} sequence(s) would be added; "
                f"{result.added_rows} TSV row(s) would be added; "
                f"{updated_rows} TSV row(s) would be updated; "
                f"{pending_ncbi_references} authoritative reference(s) pending NCBI; "
                f"{result.failed_sequences} sequence(s) would be skipped."
            ),
            fg='yellow',
        )
    else:
        click.secho(
            (
                "Reference dataset updated: "
                f"LANL alignment year {result.lanl_alignment_year or 'unknown'}; "
                f"LANL alignment {'refreshed' if result.lanl_alignment_updated else 'already current'}; "
                f"{result.lanl_fasta_files_added} LANL FASTA file(s) added; "
                f"{result.lanl_fasta_files_updated} LANL FASTA file(s) updated; "
                f"{lanl_fasta_files_removed} stale LANL FASTA file(s) removed; "
                f"{result.added_sequences} sequence(s) added; "
                f"{result.added_rows} TSV row(s) added; "
                f"{updated_rows} TSV row(s) updated; "
                f"{result.updated_feature_rows} feature row(s) updated; "
                f"{pending_ncbi_references} authoritative reference(s) pending NCBI; "
                f"{result.failed_sequences} sequence(s) skipped."
            ),
            fg='green',
        )
        if pending_ncbi_accessions:
            click.secho(
                "Authoritative reference pending NCBI: "
                + ", ".join(pending_ncbi_accessions),
                fg='yellow',
            )
        if result.failed_accessions:
            click.secho(
                "Skipped GenBank accession(s): "
                + ", ".join(result.failed_accessions),
                fg='yellow',
            )

    click.echo(f"References: {result.reference_fastas_dir}")
    click.echo(f"Sequences with locations: {result.sequences_with_locations}")


# Create a group to allow multiple commands
@click.group()
@click.version_option(version=__version__, prog_name="PyHIV")
def cli():
    """PyHIV: HIV-1 sequence analysis toolkit"""
    pass


cli.add_command(main, name='run')
cli.add_command(validate)
update.add_command(update_reference_dataset)
cli.add_command(update)
cli.add_command(update_reference_dataset, name='update-reference-dataset')

if __name__ == '__main__': # pragma: no cover
    cli()

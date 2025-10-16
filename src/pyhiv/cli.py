"""
PyHIV Command Line Interface
"""
import click
from pathlib import Path
import sys
import time
from pyhiv import __version__

SUPPORTED_FASTA_EXTENSIONS = {'.fasta', '.fa', '.fna', '.ffn'}


def validate_n_jobs(ctx, param, value):
    """Validate that n_jobs is positive if provided."""
    if value is not None and value < 1:
        raise click.BadParameter('must be at least 1')
    return value


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
    '--subtyping/--no-subtyping',
    default=True,
    show_default=True,
    help='Enable or disable HIV-1 subtyping. When enabled, aligns with reference genomes for subtype identification.'
)
@click.option(
    '--splitting/--no-splitting',
    default=True,
    show_default=True,
    help='Enable or disable gene region splitting. When enabled, splits sequences into gene regions.'
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
def main(fastas_dir, subtyping, splitting, output_dir, n_jobs, verbose, quiet):
    """
    PyHIV: HIV-1 sequence alignment, subtyping, and gene region splitting tool.

    FASTAS_DIR: Directory containing input FASTA files to process.

    \b
    Examples:
        # Basic usage with default settings
        pyhiv /path/to/fastas/

        # Disable subtyping
        pyhiv /path/to/fastas/ --no-subtyping

        # Custom output directory with 4 parallel jobs
        pyhiv /path/to/fastas/ -o results/ -j 4

        # Only alignment, no splitting
        pyhiv /path/to/fastas/ --no-splitting

        # Quiet mode (only show errors)
        pyhiv /path/to/fastas/ -q
    """

    # Handle conflicting flags
    if verbose and quiet:
        raise click.UsageError("Cannot use --verbose and --quiet together")

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
        click.echo(f"PyHIV v{__version__}")
        click.echo(f"Input directory: {fastas_dir}")
        click.echo(f"Found {num_files} FASTA file(s)")
        click.echo(f"Subtyping: {'enabled' if subtyping else 'disabled'}")
        click.echo(f"Splitting: {'enabled' if splitting else 'disabled'}")
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
            n_jobs=n_jobs
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


# Create a group to allow multiple commands
@click.group()
@click.version_option(version=__version__, prog_name="PyHIV")
def cli():
    """PyHIV: HIV-1 sequence analysis toolkit"""
    pass


cli.add_command(main, name='run')
cli.add_command(validate)

if __name__ == '__main__': # pragma: no cover
    cli()

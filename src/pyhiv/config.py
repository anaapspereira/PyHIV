import os
from pathlib import Path


def get_reference_base_dir() -> Path:
    """Return the root directory for reference genomes."""
    return Path(
        os.getenv(
            "REFERENCE_GENOMES_DIR",
            Path(__file__).parent / "reference_genomes"
        )
    )


def get_reference_paths(base_dir: Path | None = None) -> dict[str, Path]:
    """Return all key reference paths derived from the base directory."""
    base = Path(base_dir or get_reference_base_dir())
    return {
        "REFERENCE_GENOMES_DIR": base,
        "REFERENCE_GENOMES_FASTAS_DIR": base / "reference_fastas",
        "HXB2_GENOME_FASTA_DIR": base / "HXB2_fasta",
        "SEQUENCES_WITH_LOCATION": base / "sequences_with_locations.tsv",
    }


def validate_reference_paths(paths: dict[str, Path] | None = None):
    """Ensure all required reference files and directories exist."""
    paths = paths or get_reference_paths()
    for key, path in paths.items():
        if not path.exists():
            raise FileNotFoundError(f"{key} not found: {path}")

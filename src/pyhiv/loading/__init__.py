import os
from pathlib import Path

from .read_fastas import read_input_fastas


REFERENCE_GENOMES_DIR = Path(os.getenv('REFERENCE_GENOMES_DIR', str(Path(__file__).parent / 'reference_genomes')))
REFERENCE_GENOMES_FASTAS_DIR = REFERENCE_GENOMES_DIR / 'reference_fastas'
HXB2_GENOME_FASTA_DIR = REFERENCE_GENOMES_DIR / 'HXB2_fasta'
SEQUENCES_WITH_LOCATION = REFERENCE_GENOMES_DIR / 'sequences_with_locations.tsv'

if not REFERENCE_GENOMES_DIR.exists():
    raise FileNotFoundError(f"Directory not found: {REFERENCE_GENOMES_DIR}")
if not REFERENCE_GENOMES_FASTAS_DIR.exists():
    raise FileNotFoundError(f"Directory not found: {REFERENCE_GENOMES_FASTAS_DIR}")
if not HXB2_GENOME_FASTA_DIR.exists():
    raise FileNotFoundError(f"Directory not found: {HXB2_GENOME_FASTA_DIR}")
if not SEQUENCES_WITH_LOCATION.exists():
    raise FileNotFoundError(f"File not found: {SEQUENCES_WITH_LOCATION}")

# Explicitly declare the public API of the module
__all__ = [
    'read_input_fastas',
    'REFERENCE_GENOMES_DIR',
    'REFERENCE_GENOMES_FASTAS_DIR',
    'HXB2_GENOME_FASTA_DIR',
    'SEQUENCES_WITH_LOCATION',
]

from pyhiv.config import get_reference_paths, validate_reference_paths
from .read_fastas import read_input_fastas

paths = get_reference_paths()

REFERENCE_GENOMES_DIR = paths["REFERENCE_GENOMES_DIR"]
REFERENCE_GENOMES_FASTAS_DIR = paths["REFERENCE_GENOMES_FASTAS_DIR"]
HXB2_GENOME_FASTA_DIR = paths["HXB2_GENOME_FASTA_DIR"]
SEQUENCES_WITH_LOCATION = paths["SEQUENCES_WITH_LOCATION"]

__all__ = [
    "read_input_fastas",
    "REFERENCE_GENOMES_DIR",
    "REFERENCE_GENOMES_FASTAS_DIR",
    "HXB2_GENOME_FASTA_DIR",
    "SEQUENCES_WITH_LOCATION",
    "validate_reference_paths",
]

from pyhiv.config import get_reference_paths, validate_reference_paths
from .read_fastas import (
    INPUT_QC_TSV_COLUMNS,
    SUPPORTED_FASTA_EXTENSIONS,
    InputSequenceQC,
    discover_fasta_files,
    input_qc_rows_for_records,
    read_input_fastas,
    validate_input_sequence,
    write_input_qc_tsv,
)

paths = get_reference_paths()

REFERENCE_GENOMES_DIR = paths["REFERENCE_GENOMES_DIR"]
REFERENCE_GENOMES_FASTAS_DIR = paths["REFERENCE_GENOMES_FASTAS_DIR"]
HXB2_GENOME_FASTA_DIR = paths["HXB2_GENOME_FASTA_DIR"]
SEQUENCES_WITH_LOCATION = paths["SEQUENCES_WITH_LOCATION"]

__all__ = [
    "read_input_fastas",
    "discover_fasta_files",
    "SUPPORTED_FASTA_EXTENSIONS",
    "INPUT_QC_TSV_COLUMNS",
    "InputSequenceQC",
    "input_qc_rows_for_records",
    "validate_input_sequence",
    "write_input_qc_tsv",
    "REFERENCE_GENOMES_DIR",
    "REFERENCE_GENOMES_FASTAS_DIR",
    "HXB2_GENOME_FASTA_DIR",
    "SEQUENCES_WITH_LOCATION",
    "validate_reference_paths",
]

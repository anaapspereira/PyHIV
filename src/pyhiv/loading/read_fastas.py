import logging
from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_input_fastas(input_folder: Path) -> List[SeqRecord]:
    """
    Reads nucleotide FASTA files (.fasta, .fa, .fna, .ffn) from a specified input folder.

    Parameters
    ----------
    input_folder : Path
        Path to the folder containing the FASTA files.

    Returns
    -------
    List[SeqRecord]
        A list of BioPython SeqRecord objects containing sequence IDs and sequences.

    Raises
    ------
    NotADirectoryError
        If the input folder does not exist or is not a directory.
    """
    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input folder {input_folder} is not a directory.")

    supported_extensions = (".fasta", ".fa", ".fna", ".ffn")
    sequences = []

    fasta_files = [f for f in input_folder.iterdir() if f.suffix.lower() in supported_extensions]

    if not fasta_files:
        logging.warning(f"No FASTA files with supported extensions found in {input_folder}.")

    for fasta_file in fasta_files:
        try:
            with open(fasta_file, "r") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                logging.warning(f"File {fasta_file} contains no valid sequences.")
            else:
                sequences.extend(records)
                logging.info(f"Successfully read {len(records)} sequences from {fasta_file}")
        except Exception as e:
            logging.error(f"Error reading {fasta_file}: {e}")

    return sequences

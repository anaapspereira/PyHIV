import logging
from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_input_fastas(input_folder: Path) -> List[SeqRecord]:
    """
    Reads FASTA files from a specified input folder.

    Parameters
    ----------
    input_folder: Path
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

    sequences = []
    for fasta_file in input_folder.glob("*.fasta"):
        try:
            with open(fasta_file, "r") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                logging.warning(f"File {fasta_file} contains no valid sequences.")
            sequences.extend(records)
            logging.info(f"Successfully read {len(records)} sequences from {fasta_file}")
        except Exception as e:
            logging.error(f"Error reading {fasta_file}: {e}")
    return sequences

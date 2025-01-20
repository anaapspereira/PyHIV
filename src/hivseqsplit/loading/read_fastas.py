import logging
import os
from typing import List, Tuple

from Bio import SeqIO


def read_input_fastas(input_folder: str) -> List[Tuple[str, str]]:
    """
    Reads FASTA files from a specified input folder.

    Parameters
    ----------
    input_folder: str
        Path to the folder containing the FASTA files.

    Returns
    -------
    list
        List of tuples containing the sequence ID and the sequence itself.

    Raises
    ------
    FileNotFoundError
        If the input folder does not exist
    ValueError
        If a FASTA file cannot be parsed.
    """
    # TODO: handle logging properly
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder {input_folder} not found.")

    sequences = []
    fasta_files = [file for file in os.listdir(input_folder) if file.endswith('.fasta')]
    for file in fasta_files:
        file_path = os.path.join(input_folder, file)
        try:
            for record in SeqIO.parse(file_path, 'fasta'):
                sequences.append((record.id, str(record.seq)))
                logging.info(f"Read sequence {record.id}")
        except Exception as e:
            logging.error(f"Error reading {file_path}: {e}")
            raise ValueError(f"Failed to parse {file_path}") from e
    return sequences
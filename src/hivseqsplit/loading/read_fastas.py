import logging
import os
from pathlib import Path
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
    input_folder = Path(input_folder)
    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input folder {input_folder} is not a directory.")

    sequences = []
    fasta_files = [file for file in input_folder.iterdir() if file.suffix == '.fasta']
    for file in fasta_files:
        file_path = os.path.join(input_folder, file)
        try:
            for record in SeqIO.parse(file_path, 'fasta'):
                # sequences.append((record.id, str(record.seq)))
                sequences.append(record)
                logging.info(f"Read sequence {record.id}")
        except Exception as e:
            logging.error(f"Error reading {file_path}: {e}")
            raise ValueError(f"Failed to parse {file_path}") from e
    return sequences

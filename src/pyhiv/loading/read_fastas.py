import logging
from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


SUPPORTED_FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".ffn")


def discover_fasta_files(input_folder: Path, exclude_dirs: Iterable[Path] = ()) -> List[Path]:
    """
    Recursively find FASTA files with supported extensions under input_folder.

    Parameters
    ----------
    input_folder : Path
        Directory to search, including its subdirectories.
    exclude_dirs : Iterable[Path], optional
        Directories to skip entirely, e.g. a run's own output directory, so
        that previously generated alignment or gene-region FASTA files
        aren't picked up again as new input on a later run.

    Returns
    -------
    List[Path]
        Sorted list of matching FASTA file paths.
    """
    resolved_excludes = [Path(directory).resolve() for directory in exclude_dirs]

    found = []
    for path in input_folder.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in SUPPORTED_FASTA_EXTENSIONS:
            continue
        resolved = path.resolve()
        if any(excluded in resolved.parents for excluded in resolved_excludes):
            continue
        found.append(path)

    return sorted(found)


def read_input_fastas(input_folder: Path, exclude_dirs: Iterable[Path] = ()) -> List[SeqRecord]:
    """
    Reads nucleotide FASTA files (.fasta, .fa, .fna, .ffn) from a specified input
    folder and its subdirectories.

    Parameters
    ----------
    input_folder : Path
        Path to the folder containing the FASTA files.
    exclude_dirs : Iterable[Path], optional
        Directories to skip, e.g. this run's own output directory, so that
        previously generated result FASTAs aren't re-read as new input.

    Returns
    -------
    List[SeqRecord]
        A list of BioPython SeqRecord objects containing sequence IDs and sequences.
        Each record's annotations include "source_file" (the file's path relative
        to input_folder, using "/" separators) and "source_path" (its full path).

    Raises
    ------
    NotADirectoryError
        If the input folder does not exist or is not a directory.
    """
    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input folder {input_folder} is not a directory.")

    fasta_files = discover_fasta_files(input_folder, exclude_dirs=exclude_dirs)

    if not fasta_files:
        logging.warning(f"No FASTA files with supported extensions found in {input_folder}.")

    sequences = []
    for fasta_file in fasta_files:
        try:
            with open(fasta_file, "r") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                logging.warning(f"File {fasta_file} contains no valid sequences.")
            else:
                relative_source = fasta_file.relative_to(input_folder).as_posix()
                for record in records:
                    record.annotations["source_file"] = relative_source
                    record.annotations["source_path"] = str(fasta_file)
                sequences.extend(records)
                logging.info(f"Successfully read {len(records)} sequences from {fasta_file}")
        except Exception as e:
            logging.error(f"Error reading {fasta_file}: {e}")

    return sequences

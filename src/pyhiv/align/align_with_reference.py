import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pyhiv.align import pyfamsa_align
from pyhiv.loading import REFERENCE_GENOMES_FASTAS_DIR


def process_alignment(test_seq: SeqRecord, ref_seq: SeqRecord) -> Optional[Tuple[int, str, str, str]]:
    """
    Aligns test sequence with a reference sequence and calculates the score.

    Parameters
    ----------
    test_seq: SeqRecord
        A BioPython SeqRecord object representing the test sequence.
    ref_seq: SeqRecord
        A BioPython SeqRecord object representing the reference sequence.

    Returns
    -------
    Tuple[int, str, str, str]
        A tuple containing the alignment score, the aligned test sequence, the aligned reference sequence, and the reference sequence name.
    """
    try:
        test_aligned, ref_aligned = pyfamsa_align(test_seq, ref_seq)
        score = calculate_alignment_score(test_aligned, ref_aligned)
        return score, test_aligned, ref_aligned, ref_seq.name
    except Exception as e:
        logging.error(f"Failed to process {ref_seq.name}: {e}")
        return None

def align_with_references(test_sequence: SeqRecord,
                          references_dir: Optional[Path] = None,
                          n_jobs: Optional[int] = None) -> Optional[Tuple[str, str, str]]:
    """
    Aligns a test sequence with reference sequences in parallel and returns the best match.

    Parameters
    ----------
    test_sequence: SeqRecord
        A BioPython SeqRecord object representing the test sequence.
    references_dir: Path, optional
        Path to the directory containing reference sequences. Defaults to REFERENCE_GENOMES_FASTAS_DIR, containing
        reference genomes for HIV-1 subtyping.
    n_jobs: int, optional
        Number of worker processes to use for parallel processing. Defaults to using all available CPU cores.

    Returns
    -------
    Tuple[str, str, str]
        A tuple containing the test sequence, reference sequence, and the reference file name with the best alignment.
    """
    num_workers = n_jobs or 1
    references_dir = references_dir or REFERENCE_GENOMES_FASTAS_DIR

    if not isinstance(references_dir, Path) or not references_dir.exists():
        logging.error("Invalid reference directory provided.")
        return None

    # Load reference sequences efficiently
    ref_sequences: List[SeqRecord] = []
    for ref_file in references_dir.glob("*.fasta"):  # Only process FASTA files
        try:
            with open(ref_file, "r") as handle:
                ref_sequences.extend(list(SeqIO.parse(handle, "fasta")))
        except Exception as e:
            logging.error(f"Error reading {ref_file}: {e}")

    if not ref_sequences:
        logging.error("No valid reference sequences found.")
        return None

    best_alignment = None
    best_score = float('-inf')

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_alignment, test_sequence, ref): ref for ref in ref_sequences}

        for future in as_completed(futures):
            result = future.result()
            if result and result[0] > best_score:
                best_score = result[0]
                best_alignment = result[1:]

    return best_alignment


def calculate_alignment_score(seq1: str, seq2: str) -> int:
    """
    Calculate the alignment score between two sequences.

    Parameters
    ----------
    seq1: str
        First sequence to compare.
    seq2: str
        Second sequence to compare.

    Returns
    -------
    int
        Number of positions where the sequences are equal.
    """
    try:
        return sum(1 for seq1_nt, seq2_nt in zip(seq1, seq2)
               if seq1_nt.upper() == seq2_nt.upper() and seq1_nt != "-")
    except ValueError:
        logging.error("Sequences have different lengths, alignment might be incorrect.")
        return 0


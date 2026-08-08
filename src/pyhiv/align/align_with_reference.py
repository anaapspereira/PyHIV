from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional, Tuple, List

from pyhiv.align.tools import DEFAULT_ALIGNMENT_TOOL, align_sequences
from pyhiv.loading import REFERENCE_GENOMES_FASTAS_DIR

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError: # pragma: no cover
    raise ImportError("BioPython is required for this module. Please install it via 'pip install biopython'.")


DNA_BASES = set("ACGT")
DEFAULT_KMER_SIZE = 15
DEFAULT_REFERENCE_TOP_K = 30


def process_alignment(
    test_seq: SeqRecord,
    ref_seq: SeqRecord,
    alignment_tool: str = DEFAULT_ALIGNMENT_TOOL,
) -> Optional[Tuple[int, str, str, str]]:
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
        test_aligned, ref_aligned = align_sequences(test_seq, ref_seq, alignment_tool)
        score = calculate_alignment_score(test_aligned, ref_aligned)
        return score, test_aligned, ref_aligned, ref_seq.name
    except Exception as e:
        logging.error(f"Failed to process {ref_seq.name}: {e}")
        return None

def align_with_references(test_sequence: SeqRecord,
                          references_dir: Optional[Path] = None,
                          n_jobs: Optional[int] = None,
                          alignment_tool: str = DEFAULT_ALIGNMENT_TOOL,
                          kmer_size: int = DEFAULT_KMER_SIZE,
                          reference_top_k: int = DEFAULT_REFERENCE_TOP_K,
                          allowed_reference_accessions: Optional[set[str]] = None,
                          include_alignment_scores: bool = False) -> Optional[
                              Tuple[str, str, str] | Tuple[str, str, str, list[tuple[int, str]]]
                          ]:
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
    alignment_tool: str, optional
        Alignment tool to use. Defaults to edlib-HW.
    kmer_size: int, optional
        K-mer size used to prefilter candidate references before alignment.
    reference_top_k: int, optional
        Number of top k-mer ranked references to align. Use 0 or a negative value
        to align all references.
    allowed_reference_accessions: set[str], optional
        Reference accessions allowed for alignment. If omitted, all FASTA
        references in references_dir are eligible.
    include_alignment_scores: bool, optional
        If True, append a ranked list of (score, reference_name) tuples to the
        returned best alignment.

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

    # Load reference sequences efficiently.
    ref_sequences: List[SeqRecord] = []
    for ref_file in sorted(references_dir.glob("*.fasta")):  # Only process FASTA files
        if allowed_reference_accessions is not None:
            accession = reference_accession_from_name(ref_file.stem)
            if accession not in allowed_reference_accessions:
                continue

        try:
            with open(ref_file, "r") as handle:
                ref_sequences.extend(list(SeqIO.parse(handle, "fasta")))
        except Exception as e:
            logging.error(f"Error reading {ref_file}: {e}")

    if not ref_sequences:
        logging.error("No valid reference sequences found.")
        return None

    candidate_refs = kmer_shortlist(
        test_sequence,
        ref_sequences,
        size=kmer_size,
        top_k=reference_top_k,
    )

    if not candidate_refs:
        logging.error("No reference sequences selected by k-mer prefilter.")
        return None

    results: list[tuple[int, str, str, str]] = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(process_alignment, test_sequence, ref, alignment_tool): ref
            for ref in candidate_refs
        }

        for future in as_completed(futures):
            result = future.result()
            if result:
                results.append(result)

    if not results:
        return None

    results.sort(key=lambda item: (-item[0], item[3]))
    best_score, best_test_aligned, best_ref_aligned, best_reference_name = results[0]
    best_alignment = (best_test_aligned, best_ref_aligned, best_reference_name)

    if include_alignment_scores:
        alignment_scores = [(score, reference_name) for score, _, _, reference_name in results]
        return (*best_alignment, alignment_scores)

    return best_alignment


def reference_accession_from_name(reference_name: str) -> str:
    return Path(reference_name).stem.split("-", maxsplit=1)[0]


def kmers(sequence: str, size: int) -> set[str]:
    if size <= 0:
        raise ValueError("k-mer size must be positive.")

    sequence = sequence.upper()
    if len(sequence) < size:
        return set()

    found: set[str] = set()
    for index in range(0, len(sequence) - size + 1):
        kmer = sequence[index:index + size]
        if set(kmer) <= DNA_BASES:
            found.add(kmer)

    return found


def kmer_shortlist(
    test_sequence: SeqRecord,
    references: List[SeqRecord],
    size: int = DEFAULT_KMER_SIZE,
    top_k: int = DEFAULT_REFERENCE_TOP_K,
    reference_kmers: Optional[List[set[str]]] = None,
) -> List[SeqRecord]:
    if top_k <= 0 or top_k >= len(references):
        return references

    if reference_kmers is None:
        reference_kmers = build_reference_kmer_index(references, size)

    if len(reference_kmers) != len(references):
        raise ValueError("reference_kmers must have the same length as references.")

    query_kmers = kmers(str(test_sequence.seq), size)
    if not query_kmers:
        return references[:top_k]

    ranked: list[tuple[float, int, int, SeqRecord]] = []
    query_count = len(query_kmers)
    for index, (reference, ref_kmers) in enumerate(zip(references, reference_kmers)):
        if not ref_kmers:
            containment = 0.0
            shared = 0
        else:
            shared = len(query_kmers & ref_kmers)
            containment = shared / query_count
        ranked.append((containment, shared, index, reference))

    ranked.sort(key=lambda item: (-item[0], -item[1], item[2]))
    return [item[-1] for item in ranked[:top_k]]


def build_reference_kmer_index(
    references: List[SeqRecord],
    size: int = DEFAULT_KMER_SIZE,
) -> List[set[str]]:
    return [kmers(str(reference.seq), size) for reference in references]


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

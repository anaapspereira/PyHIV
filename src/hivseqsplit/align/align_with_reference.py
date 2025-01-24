import os

from Bio import SeqIO

from hivseqsplit.align import mafft_align
from hivseqsplit.loading import REFERENCE_GENOMES_FASTAS_DIR


# Align testing sequence with sequences from references folder to find best alignment with MAFFT
def align_with_references(test_sequence, references_dir=None):
    if references_dir is None:
        references_dir = REFERENCE_GENOMES_FASTAS_DIR
    # Read the test sequence
    #try:
    test_seq = test_sequence
    # Initialize the best alignment score to a high value
    best_score = -1
    best_alignment = None

    # Iterate over each reference sequence
    for ref_file in os.listdir(references_dir):
        ref_seq = SeqIO.read(f'{references_dir}/{ref_file}', "fasta")
        # Align the test sequence with the reference sequence
        test_aligned, ref_aligned = mafft_align(test_seq.seq, ref_seq.seq)
        # Calculate the alignment score
        score = calculate_alignment_score(test_aligned, ref_aligned)
        # Update the best alignment if the current score is better
        if score > best_score:
            best_score = score
            best_alignment = (test_aligned, ref_aligned, ref_file)
    return best_alignment
    # except Exception as e:
    #     print(f"Error reading the test sequence file: {e}")
    #     return None


# Calculate the alignment score between two sequences
def calculate_alignment_score(seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1
    return score

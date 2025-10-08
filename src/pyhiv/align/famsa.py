from pyfamsa import Aligner, Sequence


def pyfamsa_align(test_seq, ref_seq):
    """
    Aligns a test sequence to a reference sequence using the pyfamsa library.

    Parameters
    ----------
    test_seq : SeqRecord
        The test sequence to align. Should be a Biopython SeqRecord object containing the test sequence ID and sequence.
    ref_seq : SeqRecord
        The reference sequence to align against. Should be a Biopython SeqRecord object containing the reference
        sequence ID and sequence.

    Returns
    -------
    tuple
        test_aligned : str
            The aligned test sequence.
        ref_aligned : str
            The aligned reference sequence.
    """

    sequences = [
        Sequence(ref_seq.id.encode(), str(ref_seq.seq).encode()),
        Sequence(test_seq.id.encode(), str(test_seq.seq).encode())
    ]
    aligner = Aligner()
    alignment = aligner.align(sequences)
    test_aligned = alignment[1].sequence.decode()
    ref_aligned = alignment[0].sequence.decode()

    return test_aligned, ref_aligned

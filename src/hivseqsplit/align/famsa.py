from pyfamsa import Aligner, Sequence

def pyfamsa_align(test_seq, ref_seq):

    sequences = [
        Sequence(ref_seq.id.encode(), str(ref_seq.seq).encode()),
        Sequence(test_seq.id.encode(), str(test_seq.seq).encode())  # Fixed the typo here
    ]
    aligner = Aligner()
    alignment = aligner.align(sequences)
    test_aligned = alignment[1].sequence.decode()
    ref_aligned = alignment[0].sequence.decode()

    return test_aligned, ref_aligned
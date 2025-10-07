from unittest import TestCase, skip

from Bio import SeqIO

from pyhiv import align_with_references
from tests import TEST_DIR


@skip('mafft not installed in actions!')
class TestAlignments(TestCase):

    def test_align_with_references(self):
        test_sequence = TEST_DIR / 'data' / 'fastas' / '2636813482.fasta'
        references_dir = TEST_DIR / 'data' / 'references'
        test_sequence = SeqIO.read(test_sequence, "fasta")
        best_alignment = align_with_references(test_sequence, references_dir)
        self.assertEqual(len(best_alignment), 3)


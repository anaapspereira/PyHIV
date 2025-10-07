from unittest import TestCase, mock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pyhiv.align.famsa as pyfamsa_module  # replace with your actual module


class TestPyFamsaAlign(TestCase):

    @mock.patch("pyhiv.align.famsa.Aligner")
    @mock.patch("pyhiv.align.famsa.Sequence")
    def test_pyfamsa_align_returns_aligned_sequences(self, mock_sequence, mock_aligner):
        """Should return decoded aligned sequences in the correct order."""
        # Setup test sequences
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("TTT"), id="ref")

        # Mock Sequence constructor to just return the id (we don’t need actual objects)
        def sequence_side_effect(seq_id, seq_bytes):
            mock_obj = mock.MagicMock()
            mock_obj.sequence = seq_bytes  # store bytes so we can decode later
            return mock_obj
        mock_sequence.side_effect = sequence_side_effect

        # Mock alignment object returned by Aligner.align
        mock_alignment = [
            mock.MagicMock(sequence=b"TTT-aligned"),  # ref_aligned
            mock.MagicMock(sequence=b"AAA-aligned")   # test_aligned
        ]
        mock_aligner.return_value.align.return_value = mock_alignment

        # Call the function
        test_aligned, ref_aligned = pyfamsa_module.pyfamsa_align(test_seq, ref_seq)

        # Assertions
        self.assertEqual(test_aligned, "AAA-aligned")
        self.assertEqual(ref_aligned, "TTT-aligned")

        # Check that Sequence constructor was called correctly
        mock_sequence.assert_any_call(ref_seq.id.encode(), b"TTT")
        mock_sequence.assert_any_call(test_seq.id.encode(), b"AAA")

        # Check that Aligner.align was called with the sequence list
        mock_aligner.return_value.align.assert_called_once()

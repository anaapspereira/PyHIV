from pathlib import Path
from unittest import TestCase, mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pyhiv.align.align_with_reference as alignment


class TestCalculateAlignmentScore(TestCase):
    def test_valid_score(self):
        """Should count matching bases ignoring case and gaps."""
        seq1 = "ATGC"
        seq2 = "Atg-"
        result = alignment.calculate_alignment_score(seq1, seq2)
        self.assertEqual(result, 3)

    def test_handles_value_error(self):
        """Should log error and return 0 when zip() breaks (e.g. invalid input type)."""
        with mock.patch("logging.error") as mock_log:
            with mock.patch("builtins.zip", side_effect=ValueError("bad zip")):
                result = alignment.calculate_alignment_score("AAA", "AAA")
            mock_log.assert_called_once()
            self.assertEqual(result, 0)


class TestProcessAlignment(TestCase):
    @mock.patch("pyhiv.align.align_with_reference.calculate_alignment_score", return_value=42)
    @mock.patch("pyhiv.align.align_with_reference.pyfamsa_align", return_value=("AAA", "AAA"))
    def test_process_alignment_success(self, mock_align, mock_score):
        """Should return tuple with score, test, ref, and name."""
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("AAA"), id="ref")
        ref_seq.name = "REF1"

        result = alignment.process_alignment(test_seq, ref_seq)
        self.assertEqual(result, (42, "AAA", "AAA", "REF1"))
        mock_align.assert_called_once()
        mock_score.assert_called_once()

    @mock.patch("logging.error")
    @mock.patch("pyhiv.align.align_with_reference.pyfamsa_align", side_effect=Exception("boom"))
    def test_process_alignment_failure(self, mock_align, mock_log):
        """Should log error and return None when alignment fails."""
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("AAA"), id="ref")
        ref_seq.name = "REF1"

        result = alignment.process_alignment(test_seq, ref_seq)
        self.assertIsNone(result)
        mock_log.assert_called_once()
        self.assertIn("Failed to process REF1", mock_log.call_args[0][0])


class TestAlignWithReferences(TestCase):

    def setUp(self):
        self.ref_dir = Path("tmp_refs")
        self.ref_dir.mkdir(exist_ok=True)
        self.test_seq = SeqRecord(Seq("AAA"), id="test")

    def tearDown(self):
        for f in self.ref_dir.glob("*"):
            f.unlink()
        self.ref_dir.rmdir()

    @mock.patch("logging.error")
    def test_invalid_reference_dir(self, mock_log):
        """Should return None and log error if reference dir is invalid."""
        result = alignment.align_with_references(self.test_seq, references_dir="not_a_path")
        self.assertIsNone(result)
        mock_log.assert_called_once_with("Invalid reference directory provided.")

    @mock.patch("logging.error")
    def test_reference_dir_does_not_exist(self, mock_log):
        """Should return None if directory path does not exist."""
        fake_path = self.ref_dir / "missing"
        result = alignment.align_with_references(self.test_seq, references_dir=fake_path)
        self.assertIsNone(result)
        mock_log.assert_called_once_with("Invalid reference directory provided.")

    @mock.patch("logging.error")
    def test_no_valid_references(self, mock_log):
        """Should log error if no FASTA files found."""
        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir)
        self.assertIsNone(result)
        mock_log.assert_called_with("No valid reference sequences found.")

    @mock.patch("pyhiv.align.align_with_reference.as_completed")
    @mock.patch("pyhiv.align.align_with_reference.ProcessPoolExecutor")
    @mock.patch("pyhiv.align.align_with_reference.SeqIO.parse")
    def test_successful_alignment_with_best_score(self, mock_parse, mock_executor, mock_as_completed):
        """Test that the best alignment is selected correctly."""
        # Create fake FASTA files
        fasta1 = self.ref_dir / "ref1.fasta"
        fasta1.write_text(">ref1\nAAA\n")
        fasta2 = self.ref_dir / "ref2.fasta"
        fasta2.write_text(">ref2\nAAA\n")

        ref1 = SeqRecord(Seq("AAA"), id="ref1")
        ref2 = SeqRecord(Seq("AAA"), id="ref2")
        mock_parse.side_effect = [[ref1], [ref2]]

        # Dummy future objects
        class DummyFuture:
            def __init__(self, result_value):
                self._result_value = result_value

            def result(self):
                return self._result_value

        dummy_futures = [DummyFuture((5, "TEST", "REF", "refA")),
                         DummyFuture((10, "TEST2", "REF2", "refB"))]

        # Mock as_completed to just yield our dummy futures
        mock_as_completed.return_value = dummy_futures

        # Mock ProcessPoolExecutor to just return dummy futures mapping
        dummy_executor = mock.MagicMock()
        dummy_executor.__enter__.return_value = dummy_executor
        dummy_executor.__exit__.return_value = False
        dummy_executor.submit.side_effect = lambda *a, **kw: DummyFuture((10, "TEST2", "REF2", "refB"))
        mock_executor.return_value = dummy_executor

        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir, n_jobs=1)

        self.assertEqual(result, ("TEST2", "REF2", "refB"))

    @mock.patch("pyhiv.align.align_with_reference.SeqIO.parse", side_effect=Exception("bad parse"))
    @mock.patch("logging.error")
    def test_seqio_parse_error_logged(self, mock_log, mock_parse):
        """Should log error and skip bad FASTA file."""
        fasta_file = self.ref_dir / "bad.fasta"
        fasta_file.write_text(">bad\nAAA\n")

        result = alignment.align_with_references(self.test_seq, references_dir=self.ref_dir, n_jobs=1)
        self.assertIsNone(result)
        self.assertIn("No valid reference sequences found.", mock_log.call_args[0][0])

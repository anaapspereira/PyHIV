import tempfile
from pathlib import Path
from unittest import TestCase, mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from pyhiv.loading import read_input_fastas
from tests import TEST_DIR


class TestReadFastas(TestCase):

    def setUp(self):
        self.fastas_path = TEST_DIR / "data/fastas"
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def _write_fasta(self, filename: str, records: list[SeqRecord]):
        """Helper function to write FASTA files into temp directory"""
        fasta_path = self.input_path / filename
        SeqIO.write(records, fasta_path, "fasta")
        return fasta_path

    def test_read_fastas(self):
        """Should correctly read all FASTA files in the provided data folder."""
        fastas = read_input_fastas(self.fastas_path)
        self.assertEqual(len(fastas), 5, "Expected 5 total sequences from test FASTAs")

    def test_read_fastas_invalid_path(self):
        """Should raise NotADirectoryError when input path does not exist."""
        fastas_path = TEST_DIR / "data/invalid_path"
        with self.assertRaises(NotADirectoryError):
            read_input_fastas(fastas_path)

    def test_reads_valid_fasta_files(self):
        """Test that valid FASTA files are read correctly."""
        records = [
            SeqRecord(Seq("ATGC"), id="seq1"),
            SeqRecord(Seq("GATTACA"), id="seq2")
        ]
        self._write_fasta("sample.fasta", records)

        result = read_input_fastas(self.input_path)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].id, "seq1")
        self.assertEqual(result[1].seq, Seq("GATTACA"))

    def test_raises_error_if_not_directory(self):
        """Test that NotADirectoryError is raised if input path is not a directory."""
        fake_dir = self.input_path / "not_a_dir.fasta"
        fake_dir.touch()

        with self.assertRaises(NotADirectoryError):
            read_input_fastas(fake_dir)

    @mock.patch("logging.warning")
    def test_warns_if_empty_fasta(self, mock_warning):
        """Test that a warning is logged for empty FASTA files."""
        empty_file = self.input_path / "empty.fasta"
        empty_file.write_text("")  # Create empty file

        result = read_input_fastas(self.input_path)
        self.assertEqual(result, [])
        mock_warning.assert_called_once()
        self.assertIn("contains no valid sequences", mock_warning.call_args[0][0])

    @mock.patch("logging.error")
    def test_logs_error_on_parse_failure(self, mock_error):
        """Test that an error is logged if FASTA parsing fails."""
        bad_file = self.input_path / "corrupt.fasta"
        bad_file.write_text(">seq1\nATGC\n>broken_sequence")  # malformed FASTA

        # Force SeqIO.parse to raise an exception
        with mock.patch("Bio.SeqIO.parse", side_effect=ValueError("Parsing failed")):
            result = read_input_fastas(self.input_path)

        self.assertEqual(result, [])
        mock_error.assert_called_once()
        self.assertIn("Error reading", mock_error.call_args[0][0])

    @mock.patch("logging.info")
    def test_logs_info_on_successful_read(self, mock_info):
        """Test that info logging occurs on successful reads."""
        records = [SeqRecord(Seq("TTTT"), id="seqA")]
        self._write_fasta("one.fasta", records)

        result = read_input_fastas(self.input_path)
        self.assertEqual(len(result), 1)
        mock_info.assert_called_once()
        self.assertIn("Successfully read 1 sequences", mock_info.call_args[0][0])


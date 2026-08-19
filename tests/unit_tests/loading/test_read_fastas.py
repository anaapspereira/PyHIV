import tempfile
from pathlib import Path
from unittest import TestCase, mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from pyhiv.loading import discover_fasta_files, read_input_fastas
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
        self.assertEqual(len(fastas), 3, "Expected 5 total sequences from test FASTAs")

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
        self.assertEqual(result[0].annotations["source_file"], "sample.fasta")
        self.assertEqual(result[1].annotations["source_file"], "sample.fasta")

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

    @mock.patch("logging.warning")
    def test_no_supported_fasta_files(self, mock_warning):
        """Ensure line 37 is executed when directory contains files but none are FASTA."""
        # Create some non-FASTA files
        (self.input_path / "file1.txt").write_text("hello")
        (self.input_path / "file2.csv").write_text("1,2,3")

        result = read_input_fastas(self.input_path)
        self.assertEqual(result, [])
        # Check that the warning for no supported FASTA files was triggered
        mock_warning.assert_called_once()
        self.assertIn("No FASTA files with supported extensions found", mock_warning.call_args[0][0])

    def test_discover_fasta_files_finds_nested_files(self):
        """Should recursively find FASTA files inside subdirectories."""
        nested_dir = self.input_path / "subdir"
        nested_dir.mkdir()
        nested_file = nested_dir / "nested.fasta"
        nested_file.write_text(">seq1\nACGT\n")

        result = discover_fasta_files(self.input_path)
        self.assertEqual(result, [nested_file])

    def test_discover_fasta_files_excludes_given_directories(self):
        """Should skip files under any directory listed in exclude_dirs."""
        kept_dir = self.input_path / "keep"
        kept_dir.mkdir()
        kept_file = kept_dir / "sample.fasta"
        kept_file.write_text(">seq1\nACGT\n")

        excluded_dir = self.input_path / "PyHIV_results"
        excluded_dir.mkdir()
        excluded_file = excluded_dir / "best_alignment_sample_seq1.fasta"
        excluded_file.write_text(">Reference X\nACGT\n>seq1\nACGT\n")

        result = discover_fasta_files(self.input_path, exclude_dirs=[excluded_dir])
        self.assertEqual(result, [kept_file])

    def test_read_input_fastas_uses_relative_path_as_source_file(self):
        """source_file should be the path relative to the input folder, not just the basename."""
        nested_dir = self.input_path / "batch1"
        nested_dir.mkdir()
        SeqIO.write([SeqRecord(Seq("ACGT"), id="seq1")], nested_dir / "sample.fasta", "fasta")

        result = read_input_fastas(self.input_path)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].annotations["source_file"], "batch1/sample.fasta")

    def test_read_input_fastas_excludes_output_dir(self):
        """read_input_fastas should not re-read files under an excluded output directory."""
        self._write_fasta("sample.fasta", [SeqRecord(Seq("ACGT"), id="seq1")])

        output_dir = self.input_path / "PyHIV_results"
        output_dir.mkdir()
        SeqIO.write(
            [SeqRecord(Seq("ACGT"), id="old_seq")],
            output_dir / "best_alignment_old_seq.fasta",
            "fasta",
        )

        result = read_input_fastas(self.input_path, exclude_dirs=[output_dir])
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].id, "seq1")

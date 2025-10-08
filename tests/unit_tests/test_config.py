import os
import tempfile
import shutil
import unittest
from pathlib import Path
from unittest.mock import patch

from pyhiv.config import (
    get_reference_base_dir,
    get_reference_paths,
    validate_reference_paths,
)


class TestConfig(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: shutil.rmtree(self.tmp_dir, ignore_errors=True))

    @patch.dict(os.environ, {"REFERENCE_GENOMES_DIR": "/custom/path"})
    def test_get_reference_base_dir_from_env(self):
        """Should use REFERENCE_GENOMES_DIR from environment if set."""
        base = get_reference_base_dir()
        self.assertEqual(base, Path("/custom/path"))

    @patch.dict(os.environ, {}, clear=True)
    def test_get_reference_base_dir_default(self):
        """Should use default path if REFERENCE_GENOMES_DIR not set."""
        base = get_reference_base_dir()
        expected_default = Path(__file__).resolve().parent.parent / "pyhiv" / "reference_genomes"
        # Check that it ends with reference_genomes (don’t hardcode full path)
        self.assertTrue(str(base).endswith("reference_genomes"))

    def test_get_reference_paths_from_base_dir(self):
        """Should construct correct paths given a base directory."""
        paths = get_reference_paths(self.tmp_dir)
        self.assertEqual(paths["REFERENCE_GENOMES_DIR"], self.tmp_dir)
        self.assertEqual(paths["REFERENCE_GENOMES_FASTAS_DIR"], self.tmp_dir / "reference_fastas")
        self.assertEqual(paths["HXB2_GENOME_FASTA_DIR"], self.tmp_dir / "HXB2_fasta")
        self.assertEqual(paths["SEQUENCES_WITH_LOCATION"], self.tmp_dir / "sequences_with_locations.tsv")

    @patch("pyhiv.config.get_reference_base_dir", return_value=Path("/mock/base"))
    def test_get_reference_paths_uses_get_reference_base_dir(self, mock_base):
        """Should fall back to get_reference_base_dir if no base_dir is provided."""
        paths = get_reference_paths()
        self.assertEqual(paths["REFERENCE_GENOMES_DIR"], Path("/mock/base"))
        mock_base.assert_called_once()

    def test_validate_reference_paths_success(self):
        """Should not raise if all required paths exist."""
        (self.tmp_dir / "reference_fastas").mkdir()
        (self.tmp_dir / "HXB2_fasta").mkdir()
        (self.tmp_dir / "sequences_with_locations.tsv").write_text("dummy")

        paths = get_reference_paths(self.tmp_dir)
        validate_reference_paths(paths)

    def test_validate_reference_paths_missing_file(self):
        """Should raise FileNotFoundError if a required path is missing."""
        (self.tmp_dir / "reference_fastas").mkdir()

        paths = get_reference_paths(self.tmp_dir)
        with self.assertRaises(FileNotFoundError):
            validate_reference_paths(paths)

    @patch("pyhiv.config.get_reference_paths")
    def test_validate_reference_paths_calls_get_reference_paths(self, mock_get_paths):
        """Should call get_reference_paths when no argument is passed."""
        mock_get_paths.return_value = {
            "REFERENCE_GENOMES_DIR": self.tmp_dir,
            "REFERENCE_GENOMES_FASTAS_DIR": self.tmp_dir,
            "HXB2_GENOME_FASTA_DIR": self.tmp_dir,
            "SEQUENCES_WITH_LOCATION": self.tmp_dir,
        }
        # Should not raise
        validate_reference_paths()
        mock_get_paths.assert_called_once()

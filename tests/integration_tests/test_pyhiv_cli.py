from unittest import TestCase
import shutil
from unittest.mock import patch

import pytest
from click.testing import CliRunner

from pyhiv import __version__
from pyhiv.cli import cli, validate_n_jobs
from tests import TEST_DIR

DATA_DIR = TEST_DIR / "data" / "fastas"
REFERENCE_BASE = TEST_DIR / "data" / "references"


class TestPyHIVCLI(TestCase):

    def setUp(self):
        self.output_dir = TEST_DIR / "output_cli"
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.runner = CliRunner()

    def tearDown(self):
        shutil.rmtree(self.output_dir, ignore_errors=True)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV")
    def test_run_cli_defaults(self, mock_pyhiv):
        """Test running CLI with default settings."""
        result = self.runner.invoke(
            cli, ["run", str(DATA_DIR)]
        )

        self.assertEqual(result.exit_code, 0)
        mock_pyhiv.assert_called_once()
        self.assertIn("Processing", result.output)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV")
    def test_run_cli_with_options(self, mock_pyhiv):
        """Test CLI with subtyping, splitting, output dir, and n-jobs."""
        result = self.runner.invoke(
            cli, [
                "run", str(DATA_DIR),
                "--no-subtyping",
                "--no-splitting",
                "-o", str(self.output_dir),
                "-j", "2",
                "--verbose"
            ]
        )

        self.assertEqual(result.exit_code, 0)
        mock_pyhiv.assert_called_once_with(
            fastas_dir=str(DATA_DIR),
            subtyping=False,
            splitting=False,
            output_dir=str(self.output_dir),
            n_jobs=2
        )
        self.assertIn(f"PyHIV v{__version__}", result.output)

    def test_validate_cli_no_files(self):
        """Test validate command with empty directory."""
        empty_dir = TEST_DIR / "empty_dir"
        empty_dir.mkdir(exist_ok=True)
        result = self.runner.invoke(cli, ["validate", str(empty_dir)])
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("✗ No FASTA files found", result.output)
        shutil.rmtree(empty_dir)

    def test_validate_cli_with_files(self):
        """Test validate command with FASTA files."""
        result = self.runner.invoke(cli, ["validate", str(DATA_DIR)])
        self.assertEqual(result.exit_code, 0)
        self.assertIn("✓ Found", result.output)
        # Should list at least one file
        self.assertRegex(result.output, r"\s• .+")

    def test_run_cli_verbose_and_quiet_conflict(self):
        """Test that CLI raises UsageError if both --verbose and --quiet are used."""
        result = self.runner.invoke(
            cli, ["run", str(DATA_DIR), "--verbose", "--quiet"]
        )
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Cannot use --verbose and --quiet together", result.output)

    @patch("pyhiv.PyHIV", side_effect=ImportError("mocked import error"))
    def test_run_cli_import_error(self, mock_pyhiv):
        result = self.runner.invoke(cli, ["run", str(DATA_DIR)])
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Error: Could not import PyHIV module", result.output)

    @patch("pyhiv.PyHIV", side_effect=Exception("mocked exception"))
    def test_run_cli_general_exception(self, mock_pyhiv):
        """Test CLI handles general exceptions gracefully."""
        result = self.runner.invoke(cli, ["run", str(DATA_DIR)])
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Error during processing: mocked exception", result.output)

    def test_validate_n_jobs_invalid(self):
        """Test that validate_n_jobs raises BadParameter for n_jobs < 1"""
        ctx = None
        param = None
        with pytest.raises(Exception) as exc:
            validate_n_jobs(ctx, param, 0)
        assert "must be at least 1" in str(exc.value)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV")  # <-- patch pyhiv.PyHIV, NOT pyhiv.cli.PyHIV
    def test_run_cli_output_dir_exists_warning(self, mock_pyhiv):
        existing_dir = TEST_DIR / "existing_dir"
        existing_dir.mkdir(exist_ok=True)
        result = self.runner.invoke(cli, ["run", str(DATA_DIR), "-o", str(existing_dir)])
        self.assertIn("Warning: Output directory", result.output)
        shutil.rmtree(existing_dir)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV")
    def test_run_cli_verbose_output_files(self, mock_pyhiv):
        output_dir = self.output_dir
        (output_dir / "final_table.tsv").touch()
        for i in range(5):
            (output_dir / f"best_alignment_{i}.fasta").touch()

        result = self.runner.invoke(cli, ["run", str(DATA_DIR), "-o", str(output_dir), "--verbose"])
        self.assertIn("final_table.tsv", result.output)
        self.assertIn("best_alignment_0.fasta", result.output)
        self.assertIn("... and 2 more alignment file(s)", result.output)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV")
    def test_run_cli_verbose_no_alignment_files(self, mock_pyhiv):
        output_dir = self.output_dir
        (output_dir / "final_table.tsv").touch()
        result = self.runner.invoke(cli, ["run", str(DATA_DIR), "-o", str(output_dir), "--verbose"])
        self.assertIn("final_table.tsv", result.output)

    def test_run_cli_no_fasta_files(self):
        """Trigger the num_files == 0 branch."""
        empty_dir = TEST_DIR / "empty_dir_cli"
        empty_dir.mkdir(exist_ok=True)
        result = self.runner.invoke(cli, ["run", str(empty_dir)])
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Error: No FASTA files found", result.output)
        shutil.rmtree(empty_dir)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV", side_effect=KeyboardInterrupt)
    def test_run_cli_keyboard_interrupt(self, mock_pyhiv):
        """Test that CLI handles KeyboardInterrupt (Ctrl+C) gracefully."""
        result = self.runner.invoke(cli, ["run", str(DATA_DIR)])
        self.assertEqual(result.exit_code, 130)
        self.assertIn("Processing interrupted by user.", result.output)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.PyHIV", side_effect=Exception("mocked exception"))
    def test_run_cli_general_exception_verbose(self, mock_pyhiv):
        """Test CLI handles general exceptions and prints traceback in verbose mode."""
        import traceback

        result = self.runner.invoke(cli, ["run", str(DATA_DIR), "--verbose"])

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("Error during processing: mocked exception", result.output)

        # The traceback should appear in output
        self.assertIn("Traceback (most recent call last):", result.output)
        self.assertIn("Exception: mocked exception", result.output)


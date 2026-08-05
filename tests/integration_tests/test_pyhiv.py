from pathlib import Path
from unittest import TestCase
import shutil
from unittest.mock import patch

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyhiv import (
    NO_SUBTYPING_LABEL,
    PyHIV,
    SEQUENCE_TOO_LONG_WARNING,
    normalize_reference_groups,
    selected_reference_accessions,
    summarize_closest_subtypes,
)
from tests import TEST_DIR

DATA_DIR = TEST_DIR / "data" / "fastas"
REFERENCE_BASE = TEST_DIR / "data" / "references"


class TestPyHIV(TestCase):

    def setUp(self):
        """Prepare output directory for results."""
        self.output_dir = TEST_DIR / "output"
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.output_dir, ignore_errors=True)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    def test_pyhiv_with_real_fastas_splitting(self):
        """Run PyHIV on real FASTAs with splitting enabled."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=12,
            reporting=True,
            alignment_tool="PyFamsa",
        )

        # Check that best alignment files are created
        alignment_files = list(self.output_dir.glob("best_alignment_*.fasta"))
        self.assertTrue(len(alignment_files) > 0)

        # Check that final table exists and has correct columns
        table_file = self.output_dir / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        expected_cols = [
            'Sequence', 'Reference', 'Group', 'Subtype', 'Closest Subtypes',
            'Most Matching Gene Region', 'Present Gene Regions'
        ]
        self.assertListEqual(list(table.columns), expected_cols)
        self.assertTrue(len(table) > 0)

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    def test_pyhiv_with_real_fastas_no_splitting(self):
        """Run PyHIV on real FASTAs with splitting disabled."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=False,
            output_dir=str(self.output_dir),
            n_jobs=1,
            alignment_tool="PyFamsa",
        )

        table_file = self.output_dir / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')

        # Columns specific to splitting should be dropped
        self.assertListEqual(list(table.columns), ['Sequence', 'Reference', 'Group', 'Subtype', 'Closest Subtypes'])

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    def test_pyhiv_no_subtyping_uses_uniform_table_labels(self):
        """HXB2 alignment should not report reference group/subtype as subtype results."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=False,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=1,
            reporting=False,
            alignment_tool="PyFamsa",
        )

        table = pd.read_csv(self.output_dir / "final_table.tsv", sep='\t')
        self.assertTrue(len(table) > 0)
        self.assertTrue((table["Group"] == NO_SUBTYPING_LABEL).all())
        self.assertTrue((table["Subtype"] == NO_SUBTYPING_LABEL).all())
        self.assertTrue((table["Closest Subtypes"] == NO_SUBTYPING_LABEL).all())

    @patch.dict("os.environ", {"REFERENCE_GENOMES_DIR": str(REFERENCE_BASE)})
    @patch("pyhiv.align_with_references", return_value=None)
    def test_best_alignment_is_none(self, mock_align):
        """Test that PyHIV handles cases when align_with_references returns None."""
        # Run PyHIV on a real FASTA directory, alignment will always return None
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=1
        )

        # The output directory should exist but no best_alignment files should be created
        alignment_files = list(Path(self.output_dir).glob("best_alignment_*.fasta"))
        self.assertEqual(len(alignment_files), 0)

        # The final_table.tsv should still be created but empty (or just headers)
        table_file = Path(self.output_dir) / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        self.assertEqual(len(table), 0)  # No rows because no alignments succeeded

    @patch("pyhiv.align_with_references")
    @patch("pyhiv.read_input_fastas")
    @patch("pyhiv.validate_reference_paths")
    @patch("pyhiv.get_reference_paths")
    @patch("logging.warning")
    def test_sequence_longer_than_hiv_genome_is_skipped(
        self,
        mock_warning,
        mock_paths,
        mock_validate,
        mock_read_fastas,
        mock_align,
    ):
        """Sequences with more than 12000 nt should be skipped before alignment."""
        mock_paths.return_value = {
            "SEQUENCES_WITH_LOCATION": REFERENCE_BASE / "sequences_with_locations.tsv",
            "REFERENCE_GENOMES_FASTAS_DIR": REFERENCE_BASE / "reference_fastas",
            "HXB2_GENOME_FASTA_DIR": REFERENCE_BASE / "HXB2_fasta",
        }
        mock_read_fastas.return_value = [SeqRecord(Seq("A" * 12001), id="too_long")]

        PyHIV(
            fastas_dir=str(DATA_DIR),
            output_dir=str(self.output_dir),
            reporting=False,
        )

        mock_validate.assert_called_once()
        mock_align.assert_not_called()
        warning_calls = [call_args[0] for call_args in mock_warning.call_args_list]
        self.assertIn(
            ("%s Skipping sequence '%s'.", SEQUENCE_TOO_LONG_WARNING, "too_long"),
            warning_calls,
        )

        table_file = Path(self.output_dir) / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        self.assertEqual(len(table), 0)

    def test_normalize_reference_groups_accepts_comma_separated_groups(self):
        self.assertEqual(normalize_reference_groups("M,N,O,P"), ("M", "N", "O", "P"))

    def test_selected_reference_accessions_filters_by_group(self):
        reference_sequences = pd.DataFrame([
            {"accession": "acc_m", "group": "M"},
            {"accession": "acc_o", "group": "O"},
        ])

        self.assertEqual(
            selected_reference_accessions(reference_sequences, ("O",)),
            {"acc_o"},
        )

    def test_summarize_closest_subtypes_returns_top_unique_subtypes(self):
        metadata_by_accession = {
            "acc_b1": {"group": "M", "subtype": "B"},
            "acc_b2": {"group": "M", "subtype": "B"},
            "acc_c": {"group": "M", "subtype": "C"},
            "acc_o": {"group": "O", "subtype": "O"},
        }

        result = summarize_closest_subtypes(
            [
                (100, "acc_b1-B"),
                (98, "acc_b2-B"),
                (95, "acc_c-C"),
                (90, "acc_o-O"),
            ],
            metadata_by_accession,
            top_n=3,
        )

        self.assertEqual(result, "M:B (score=100); M:C (score=95); O:O (score=90)")

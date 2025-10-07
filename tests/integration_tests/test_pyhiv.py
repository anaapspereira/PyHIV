from unittest import TestCase
import shutil
import pandas as pd

from pyhiv import PyHIV
from tests import TEST_DIR

DATA_DIR = TEST_DIR / "data" / "fastas"


class TestPyHIV(TestCase):

    def setUp(self):
        """Prepare output directory for results."""
        self.output_dir = TEST_DIR / "output"
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.output_dir, ignore_errors=True)

    def test_pyhiv_with_real_fastas_splitting(self):
        """Run PyHIV on real FASTAs with splitting enabled."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=True,
            output_dir=str(self.output_dir),
            n_jobs=1
        )

        # Check that best alignment files are created
        alignment_files = list(self.output_dir.glob("best_alignment_*.fasta"))
        self.assertTrue(len(alignment_files) > 0)

        # Check that final table exists and has correct columns
        table_file = self.output_dir / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')
        expected_cols = [
            'Sequence', 'Reference', 'Subtype', 'Most Matching Gene Region', 'Present Gene Regions'
        ]
        self.assertListEqual(list(table.columns), expected_cols)
        self.assertTrue(len(table) > 0)

    def test_pyhiv_with_real_fastas_no_splitting(self):
        """Run PyHIV on real FASTAs with splitting disabled."""
        PyHIV(
            fastas_dir=str(DATA_DIR),
            subtyping=True,
            splitting=False,
            output_dir=str(self.output_dir),
            n_jobs=1
        )

        table_file = self.output_dir / "final_table.tsv"
        self.assertTrue(table_file.exists())
        table = pd.read_csv(table_file, sep='\t')

        # Columns specific to splitting should be dropped
        self.assertListEqual(list(table.columns), ['Sequence', 'Reference', 'Subtype'])

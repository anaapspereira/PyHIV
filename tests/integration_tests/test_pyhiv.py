from unittest import TestCase, skip

from pyhiv import PyHIV
from tests import TEST_DIR


@skip('mafft not installed in actions!')
class TestPyHIV(TestCase):

    def test_pyhiv(self):
        fastas_dir = TEST_DIR / 'data' / 'fastas'

        PyHIV(fastas_dir=fastas_dir, subtyping=True, splitting=True, output_dir='PyHIV_results/', n_jobs=8)
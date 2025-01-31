from unittest import TestCase, skip

from hivseqsplit import HIMAPS
from tests import TEST_DIR


@skip('mafft not installed in actions!')
class TestHIMAPS(TestCase):

    def test_himaps(self):
        fastas_dir = TEST_DIR / 'data' / 'fastas'

        HIMAPS(fastas_dir=fastas_dir, subtyping=True, splitting=True, output_dir='HIMAPS_results/')
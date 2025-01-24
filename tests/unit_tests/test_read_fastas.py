import os
from unittest import TestCase

from hivseqsplit.loading import read_input_fastas
from tests import TEST_DIR


class TestReadFastas(TestCase):

    def test_read_fastas(self):
        fastas_path = os.path.join(TEST_DIR, 'data/fastas')
        fastas = read_input_fastas(fastas_path)
        self.assertEqual(len(fastas), 5)

    def test_read_fastas_invalid_path(self):
        fastas_path = os.path.join(TEST_DIR, 'data/invalid_path')
        with self.assertRaises(FileNotFoundError):
            read_input_fastas(fastas_path)
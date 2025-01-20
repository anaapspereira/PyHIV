import os
import shutil
from unittest import TestCase

from hivseqsplit.loading import ncbi_sequences_download


class TestDownloadNCBISequences(TestCase):

    def tearDown(self):
       if os.path.exists("output/"):
           shutil.rmtree("output/")

    def test_download_ncbi_sequences(self):
        fastas = ncbi_sequences_download()
        self.assertEquals(len(fastas), 10)

        fastas = ncbi_sequences_download('../data/config.yml')
        self.assertEquals(len(fastas), 10)

        fastas = ncbi_sequences_download('../data/config2.yml')
        self.assertEquals(len(fastas), 10)

        with self.assertRaises(ValueError):
            ncbi_sequences_download('../data/invalid_config.yml')

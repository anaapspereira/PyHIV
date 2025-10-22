from unittest import TestCase

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages

from pyhiv.report.pdf_generator import render_sequence_page


class TestReportingPdfGenerator(TestCase):

    def test_render_sequence_page_smoke(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            from pathlib import Path
            pdf_path = Path(d) / "test.pdf"
            with PdfPages(pdf_path) as pdf:
                render_sequence_page(
                    pdf=pdf,
                    sequence="seq1",
                    accession="K03455",
                    subtype="B",
                    mm_region="gag",
                    present_regions=["gag"],
                    features_aln={"gag": (10, 50)},
                    ref_header="Reference K03455-B",
                    ref_seq_aligned="ACGT" * 30,
                    user_header=">seq1",
                    user_seq_aligned="ACGT" * 25,
                    y_positions=None,
                )
            self.assertTrue(pdf_path.exists() and pdf_path.stat().st_size > 0)



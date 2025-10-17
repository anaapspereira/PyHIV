from unittest import TestCase
from pathlib import Path

import pandas as pd

from pyhiv.report import PyHIVReporter


class TestReportingIntegration(TestCase):

    def write_alignment_fasta(self, path: Path, ref_name: str, ref_seq: str, seq_name: str, seq: str):
        path.write_text(
            f">Reference {ref_name}\n{ref_seq}\n>{seq_name}\n{seq}\n",
            encoding="utf-8"
        )

    def test_generate_report_end_to_end(self):
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            dpath = Path(d)
            # Prepare minimal final_table
            final_table = pd.DataFrame(
                [
                    {
                        "Sequence": "seq1",
                        "Reference": "K03455",
                        "Subtype": "B",
                        "Most Matching Gene Region": "gag",
                        "Present Gene Regions": "['gag']",
                    }
                ]
            )
            final_table_path = dpath / "final_table.tsv"
            final_table.to_csv(final_table_path, sep="\t", index=False)

            # sequences_with_locations
            swl = pd.DataFrame(
                [
                    {
                        "accession": "K03455",
                        "features": str({"gag": (1, 50)}),
                    }
                ]
            )
            swl_path = dpath / "sequences_with_locations.tsv"
            swl.to_csv(swl_path, sep="\t", index=False)

            # Alignment fasta in output dir
            output_dir = dpath
            ref_seq = "ACGT" * 60
            user_seq = "-" * 10 + ("ACGT" * 40) + "-" * 10
            self.write_alignment_fasta(output_dir / "best_alignment_seq1.fasta", "K03455-B", ref_seq, "seq1", user_seq)

            reporter = PyHIVReporter(output_dir)
            pdf_path = reporter.generate_report(
                final_table_path=final_table_path,
                sequences_with_locations_path=swl_path,
                output_pdf_name="report.pdf",
            )

            self.assertTrue(pdf_path.exists() and pdf_path.stat().st_size > 0)



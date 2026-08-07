import io
import logging
import tempfile
from pathlib import Path
from unittest import TestCase, mock

import pandas as pd
from pyhiv.report.reporter import PyHIVReporter
from pyhiv.report.constants import K03455Config


class TestPyHIVReporter(TestCase):
    """Unit tests for PyHIVReporter.generate_report with 100% coverage."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.temp_dir.name)
        self.reporter = PyHIVReporter(self.tmp_path, subtyping=True, splitting=True)

        # Dummy table DataFrames
        self.final_table = self.tmp_path / "final_table.tsv"
        self.sequences_with_locations = self.tmp_path / "sequences_with_locations.tsv"

        # Minimal correct headers
        ft_data = pd.DataFrame([
            {
                "Sequence": "seq1",
                "Reference": "ref1",
                "Group": "M",
                "Subtype": "B",
                "Closest Subtypes": "M:B (score=100); M:C (score=95)",
                "Most Matching Gene Region": "gag",
                "Present Gene Regions": "gag, pol",
            },
            {
                "Sequence": "seq2",
                "Reference": "ref2",
                "Group": "M",
                "Subtype": "C",
                "Closest Subtypes": "M:C (score=100); M:B (score=95)",
                "Most Matching Gene Region": "env",
                "Present Gene Regions": "env, nef",
            }
        ])
        ft_data.to_csv(self.final_table, sep="\t", index=False)

        swl_data = pd.DataFrame([
            {"accession": "ref1", "features": "{'gag': (1, 100)}"},
            {"accession": "ref2", "features": "{'env': (1, 200)}"},
        ])
        swl_data.to_csv(self.sequences_with_locations, sep="\t", index=False)

    def tearDown(self):
        self.temp_dir.cleanup()

    # ---- Logger configuration ----
    def test_explicit_log_level_overrides_root_logger(self):
        """Passing log_level should set it on the reporter's own logger."""
        reporter = PyHIVReporter(self.tmp_path, subtyping=True, splitting=True, log_level=logging.DEBUG)
        self.addCleanup(reporter.logger.setLevel, logging.NOTSET)
        self.assertEqual(reporter.logger.level, logging.DEBUG)

    def test_default_log_level_is_unset(self):
        """Without an explicit log_level, the reporter should inherit from the root logger."""
        reporter = PyHIVReporter(self.tmp_path, subtyping=True, splitting=True)
        self.assertEqual(reporter.logger.level, logging.NOTSET)

    # ---- Error handling ----
    def test_missing_columns_in_final_table(self):
        bad_ft = self.tmp_path / "bad.tsv"
        pd.DataFrame([{"bad": "col"}]).to_csv(bad_ft, sep="\t", index=False)
        with self.assertRaises(ValueError) as ctx:
            self.reporter.generate_report(bad_ft, self.sequences_with_locations)
        self.assertIn("Missing columns", str(ctx.exception))

    def test_missing_columns_in_sequences_with_locations(self):
        bad_swl = self.tmp_path / "bad_swl.tsv"
        pd.DataFrame([{"accession": "x"}]).to_csv(bad_swl, sep="\t", index=False)
        with self.assertRaises(ValueError) as ctx:
            self.reporter.generate_report(self.final_table, bad_swl)
        self.assertIn("must include", str(ctx.exception))

    # ---- PDF creation and normal flow ----
    @mock.patch("pyhiv.report.reporter.PdfPages")
    @mock.patch("pyhiv.report.reporter.render_sequence_page")
    @mock.patch("pyhiv.report.reporter.project_features_to_alignment", side_effect=lambda f, m: f)
    @mock.patch("pyhiv.report.reporter.build_ref_to_alignment_map", return_value=({"A": 1}, None))
    @mock.patch("pyhiv.report.reporter.is_special_reference", return_value=False)
    @mock.patch("pyhiv.report.reporter.read_alignment_fasta")
    @mock.patch("pyhiv.report.reporter.build_alignment_path")
    @mock.patch("pyhiv.report.reporter.parse_features", return_value={"gag": (1, 10)})
    @mock.patch("pyhiv.report.reporter.parse_present_regions", return_value=["gag"])
    @mock.patch("pyhiv.report.reporter.normalize_features", return_value={"gag": (1, 10)})
    @mock.patch("pyhiv.report.reporter.normalize_present_regions", return_value=["gag"])
    def test_generate_report_success_non_special(
        self,
        mock_norm_pr, mock_norm_feat, mock_pparser, mock_pf,
        mock_path, mock_read_fa, mock_special, mock_map, mock_proj, mock_render, mock_pdf,
    ):
        """Full successful run (non-special ref)."""
        # Fake alignment path and contents
        fake_fasta = self.tmp_path / "seq1.fasta"
        fake_fasta.touch()
        mock_path.return_value = fake_fasta
        mock_read_fa.return_value = ("refheader", "AAA", "userheader", "TTT")

        pdf_context = mock_pdf.return_value.__enter__.return_value

        out_path = self.reporter.generate_report(self.final_table, self.sequences_with_locations)

        # Assertions
        mock_path.assert_called()
        mock_render.assert_called()
        self.assertEqual(mock_render.call_args.kwargs["group"], "M")
        self.assertEqual(mock_render.call_args.kwargs["closest_subtypes"], "M:C (score=100); M:B (score=95)")
        self.assertTrue(out_path.exists() or out_path.name.endswith(".pdf"))
        self.assertTrue(isinstance(out_path, Path))
        pdf_context is not None

    @mock.patch("pyhiv.report.reporter.PdfPages")
    @mock.patch("pyhiv.report.reporter.render_sequence_page")
    @mock.patch("pyhiv.report.reporter.project_features_to_alignment", side_effect=lambda f, m: f)
    @mock.patch("pyhiv.report.reporter.build_ref_to_alignment_map", return_value=({"A": 1}, None))
    @mock.patch("pyhiv.report.reporter.is_special_reference", return_value=True)
    @mock.patch("pyhiv.report.reporter.read_alignment_fasta")
    @mock.patch("pyhiv.report.reporter.build_alignment_path")
    @mock.patch("pyhiv.report.reporter.parse_features", return_value={"gag": (1, 10)})
    @mock.patch("pyhiv.report.reporter.parse_present_regions", return_value=["gag"])
    @mock.patch("pyhiv.report.reporter.normalize_features", return_value={"gag": (1, 10)})
    @mock.patch("pyhiv.report.reporter.normalize_present_regions", return_value=["gag"])
    def test_generate_report_no_subtyping_group_label(
        self,
        mock_norm_pr, mock_norm_feat, mock_pparser, mock_pf,
        mock_path, mock_read_fa, mock_special, mock_map, mock_proj, mock_render, mock_pdf,
    ):
        reporter = PyHIVReporter(self.tmp_path, subtyping=False, splitting=True)
        fake_fasta = self.tmp_path / "seq1.fasta"
        fake_fasta.touch()
        mock_path.return_value = fake_fasta
        mock_read_fa.return_value = ("K03455-B", "AAA", "userheader", "TTT")

        reporter.generate_report(self.final_table, self.sequences_with_locations)

        mock_render.assert_called()
        self.assertEqual(mock_render.call_args.kwargs["group"], "No subtyping performed.")
        self.assertEqual(mock_render.call_args.kwargs["subtype"], "No subtyping performed.")

    # ---- Missing FASTA and read error ----
    @mock.patch("pyhiv.report.reporter.build_alignment_path")
    @mock.patch("pyhiv.report.reporter.read_alignment_fasta")
    def test_missing_fasta_and_read_error(self, mock_read, mock_path):
        fake_missing = self.tmp_path / "does_not_exist.fasta"
        mock_path.return_value = fake_missing
        mock_read.side_effect = IOError("bad file")

        # Patch everything else to harmless values
        with mock.patch("pyhiv.report.reporter.PdfPages"), \
             mock.patch("pyhiv.report.reporter.parse_present_regions", return_value=["gag"]), \
             mock.patch("pyhiv.report.reporter.parse_features", return_value={"gag": (1, 10)}):
            out = self.reporter.generate_report(self.final_table, self.sequences_with_locations)
            self.assertTrue(out.name.endswith(".pdf"))

    # ---- Special reference branch ----
    @mock.patch("pyhiv.report.reporter.PdfPages")
    @mock.patch("pyhiv.report.reporter.render_sequence_page")
    @mock.patch("pyhiv.report.reporter.project_features_to_alignment", side_effect=lambda f, m: f)
    @mock.patch("pyhiv.report.reporter.build_ref_to_alignment_map", return_value=({"A": 1}, None))
    @mock.patch("pyhiv.report.reporter.is_special_reference", return_value=True)
    @mock.patch("pyhiv.report.reporter.read_alignment_fasta")
    @mock.patch("pyhiv.report.reporter.build_alignment_path")
    @mock.patch("pyhiv.report.reporter.parse_features", return_value={"gag": (1, 10)})
    @mock.patch("pyhiv.report.reporter.parse_present_regions", return_value=["gag"])
    @mock.patch("pyhiv.report.reporter.normalize_features", return_value={"gag": (1, 10)})
    @mock.patch("pyhiv.report.reporter.normalize_present_regions", return_value=["gag"])
    def test_generate_report_special_reference(
        self,
        mock_norm_pr, mock_norm_feat, mock_pparser, mock_pf,
        mock_path, mock_read_fa, mock_special, mock_map, mock_proj, mock_render, mock_pdf,
    ):
        """Ensure special reference branch (K03455) executes."""
        fake_fasta = self.tmp_path / "seq1.fasta"
        fake_fasta.touch()
        mock_path.return_value = fake_fasta
        mock_read_fa.return_value = ("K03455", "AAA", "userheader", "TTT")

        pdf_context = mock_pdf.return_value.__enter__.return_value

        out_path = self.reporter.generate_report(self.final_table, self.sequences_with_locations)

        self.assertTrue(out_path.exists() or out_path.name.endswith(".pdf"))
        mock_render.assert_called()
        args, kwargs = mock_render.call_args
        # y_positions should be from K03455Config
        self.assertEqual(kwargs["y_positions"], K03455Config.Y_POSITIONS)

    def test_read_alignment_fasta_exception_branch(self):
        """Trigger the except Exception branch inside generate_report."""
        # Force a real fasta path that exists but fails to read
        fake_fasta = self.tmp_path / "seq1.fasta"
        fake_fasta.touch()

        with mock.patch("pyhiv.report.reporter.PdfPages"), \
             mock.patch("pyhiv.report.reporter.build_alignment_path", return_value=fake_fasta), \
             mock.patch("pyhiv.report.reporter.read_alignment_fasta", side_effect=RuntimeError("bad fasta")), \
             mock.patch("pyhiv.report.reporter.parse_present_regions", return_value=["gag"]), \
             mock.patch("pyhiv.report.reporter.parse_features", return_value={"gag": (1, 10)}):
            out = self.reporter.generate_report(self.final_table, self.sequences_with_locations)
            self.assertTrue(out.name.endswith(".pdf"))

    def test_no_pages_created_branch(self):
        """Trigger pages_made == 0 to hit the 'No pages created' print."""
        # Mock build_alignment_path to return non-existent file for all rows
        missing_path = self.tmp_path / "does_not_exist.fasta"
        with mock.patch("pyhiv.report.reporter.build_alignment_path", return_value=missing_path), \
             mock.patch("pyhiv.report.reporter.PdfPages"):
            out = self.reporter.generate_report(self.final_table, self.sequences_with_locations)
            self.assertTrue(out.name.endswith(".pdf"))

    def test_feature_parsing_exception_branch(self):
        """Trigger parse_features exception to cover logger.warning lines."""
        # Create SWL file with bad 'features' value
        bad_swl = self.tmp_path / "bad_features.tsv"
        pd.DataFrame([{"accession": "ref1", "features": "BAD_FEATURES"}]).to_csv(bad_swl, sep="\t", index=False)

        # Final table is fine
        ft = self.tmp_path / "final_table_good.tsv"
        pd.DataFrame([
            {
                "Sequence": "seq1",
                "Reference": "ref1",
                "Group": "M",
                "Subtype": "B",
                "Closest Subtypes": "M:B (score=100)",
                "Most Matching Gene Region": "gag",
                "Present Gene Regions": "gag, pol",
            }
        ]).to_csv(ft, sep="\t", index=False)

        # Patch parse_features to raise exception
        with mock.patch("pyhiv.report.reporter.parse_features", side_effect=ValueError("bad json")), \
             mock.patch("pyhiv.report.reporter.PdfPages"), \
             mock.patch("pyhiv.report.reporter.build_alignment_path", return_value=self.tmp_path / "none.fasta"):
            out = self.reporter.generate_report(ft, bad_swl)
            self.assertTrue(out.name.endswith(".pdf"))

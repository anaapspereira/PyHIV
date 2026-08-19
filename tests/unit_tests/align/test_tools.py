import os
import stat
import sys
import tempfile
import types
from pathlib import Path
from unittest import TestCase, mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyhiv.align import tools


class TestAlignmentTools(TestCase):
    def test_normalize_alignment_tool_defaults_to_edlib_hw(self):
        self.assertEqual(tools.normalize_alignment_tool(None), tools.EDLIB_HW)

    def test_normalize_alignment_tool_is_case_insensitive(self):
        self.assertEqual(tools.normalize_alignment_tool("parasail-nw"), tools.PARASAIL_NW)
        self.assertEqual(tools.normalize_alignment_tool("mafft"), tools.MAFFT)
        self.assertEqual(tools.normalize_alignment_tool("pyfamsa"), tools.PYFAMSA)

    def test_normalize_alignment_tool_accepts_parasail_alias(self):
        self.assertEqual(tools.normalize_alignment_tool("parasail"), tools.PARASAIL_NW)

    def test_normalize_alignment_tool_rejects_unknown_tool(self):
        with self.assertRaises(ValueError):
            tools.normalize_alignment_tool("unknown")

    def test_normalize_alignment_tool_rejects_edlib_nw(self):
        with self.assertRaises(ValueError):
            tools.normalize_alignment_tool("edlib-NW")

    def test_align_sequences_dispatches_to_selected_tool(self):
        test_seq = SeqRecord(Seq("AAA"), id="test")
        ref_seq = SeqRecord(Seq("AAA"), id="ref")

        def fake_aligner(received_test_seq, received_ref_seq):
            self.assertIs(received_test_seq, test_seq)
            self.assertIs(received_ref_seq, ref_seq)
            return "AAA", "AAA"

        original = tools._ALIGNERS[tools.EDLIB_HW]
        tools._ALIGNERS[tools.EDLIB_HW] = fake_aligner
        try:
            self.assertEqual(tools.align_sequences(test_seq, ref_seq), ("AAA", "AAA"))
        finally:
            tools._ALIGNERS[tools.EDLIB_HW] = original

    def test_apply_cigar_supports_match_insert_and_delete(self):
        test_aligned, ref_aligned = tools._apply_cigar("ACGT", "AGT", "1=1I2=")
        self.assertEqual(test_aligned, "ACGT")
        self.assertEqual(ref_aligned, "A-GT")

        test_aligned, ref_aligned = tools._apply_cigar("AGT", "ACGT", "1=1D2=")
        self.assertEqual(test_aligned, "A-GT")
        self.assertEqual(ref_aligned, "ACGT")

    @mock.patch("pyhiv.align.tools.logging.debug")
    @mock.patch("pyhiv.align.tools.logging.warning")
    def test_reference_projection_logs_dropped_insertions_as_debug(self, mock_warning, mock_debug):
        result = tools._drop_query_insertions("ACGT", "A-GT")

        self.assertEqual(result, ("AGT", "AGT"))
        mock_warning.assert_not_called()
        mock_debug.assert_called_once_with(
            "Dropped %d query insertion column(s) to keep alignment in reference coordinates.",
            1,
        )

    def test_parse_pair_fasta_alignment(self):
        text = """
>query
ACGTACGT
>reference
ACGT-CGT
"""
        self.assertEqual(
            tools._parse_pair_fasta_alignment(text),
            ("ACGTACGT", "ACGT-CGT"),
        )

    def test_parse_pair_fasta_alignment_rejects_missing_alignment(self):
        with self.assertRaises(ValueError):
            tools._parse_pair_fasta_alignment(">query\nAAA\n")

    def test_get_mafft_binary_uses_environment_override(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            executable = Path(tmp_dir) / "mafft"
            executable.write_text("#!/bin/sh\n")
            executable.chmod(executable.stat().st_mode | stat.S_IXUSR)

            with mock.patch.dict(os.environ, {tools.MAFFT_BIN_ENV: str(executable)}):
                self.assertEqual(tools.get_mafft_binary(), str(executable))

    def test_get_mafft_binary_uses_path(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with mock.patch("pyhiv.align.tools.shutil.which", return_value="/usr/local/bin/mafft"):
                self.assertEqual(tools.get_mafft_binary(), "/usr/local/bin/mafft")

    def test_get_mafft_binary_returns_none_when_not_found(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            with mock.patch("pyhiv.align.tools.shutil.which", return_value=None):
                self.assertIsNone(tools.get_mafft_binary())

    def test_write_pair_fasta_writes_both_records(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "pair.fasta"
            tools._write_pair_fasta(path, "query", "ACGT", "reference", "AC-T")

            self.assertEqual(
                path.read_text(),
                ">query\nACGT\n>reference\nAC-T\n",
            )

    def test_parse_pair_fasta_alignment_rejects_content_before_header(self):
        with self.assertRaises(ValueError):
            tools._parse_pair_fasta_alignment("ACGT\n>reference\nACGT\n")

    def test_parse_pair_fasta_alignment_rejects_mismatched_lengths(self):
        with self.assertRaises(ValueError):
            tools._parse_pair_fasta_alignment(">query\nACGT\n>reference\nAC\n")

    def test_sequence_to_string_uppercases(self):
        seq_record = SeqRecord(Seq("acgt"), id="s")
        self.assertEqual(tools._sequence_to_string(seq_record), "ACGT")

    def test_pick_edlib_location_prefers_closest_length_then_earliest_start(self):
        result = tools._pick_edlib_location(
            locations=[(0, 5), (2, 4)],
            query_length=3,
        )
        self.assertEqual(result, (2, 4))

    def test_pick_edlib_location_rejects_all_invalid_locations(self):
        with self.assertRaises(ValueError):
            tools._pick_edlib_location(
                locations=[(-1, 5), (3, 1)],
                query_length=3,
            )

    def test_drop_query_insertions_rejects_mismatched_lengths(self):
        with self.assertRaises(ValueError):
            tools._drop_query_insertions("AC", "ACG")

    def test_apply_cigar_appends_trailing_unconsumed_test_sequence(self):
        test_aligned, ref_aligned = tools._apply_cigar("ACGTAA", "ACGT", "4=")
        self.assertEqual(test_aligned, "ACGTAA")
        self.assertEqual(ref_aligned, "ACGT--")

    def test_apply_cigar_appends_trailing_unconsumed_reference(self):
        test_aligned, ref_aligned = tools._apply_cigar("ACGT", "ACGTAA", "4=")
        self.assertEqual(test_aligned, "ACGT--")
        self.assertEqual(ref_aligned, "ACGTAA")


class TestEdlibHwAlign(TestCase):
    def test_edlib_hw_align_places_query_inside_reference(self):
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("TTTACGTTTT"), id="ref")

        query_aligned, ref_aligned = tools.edlib_hw_align(test_seq, ref_seq)

        self.assertEqual(ref_aligned, "TTTACGTTTT")
        self.assertEqual(query_aligned, "---ACGT---")

    def test_edlib_hw_align_rejects_empty_test_sequence(self):
        test_seq = SeqRecord(Seq(""), id="test")
        ref_seq = SeqRecord(Seq("ACGT"), id="ref")

        with self.assertRaises(ValueError):
            tools.edlib_hw_align(test_seq, ref_seq)

    def test_edlib_hw_align_rejects_empty_reference_sequence(self):
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq(""), id="ref")

        with self.assertRaises(ValueError):
            tools.edlib_hw_align(test_seq, ref_seq)

    def _fake_edlib_module(self, align_side_effect):
        fake_module = types.SimpleNamespace()
        fake_module.align = mock.Mock(side_effect=align_side_effect)
        fake_module.getNiceAlignment = mock.Mock(
            return_value={"query_aligned": "ACGT", "target_aligned": "ACGT"}
        )
        return fake_module

    def test_edlib_hw_align_raises_when_no_location_found(self):
        fake_module = self._fake_edlib_module(
            [{"editDistance": -1, "locations": None}],
        )
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("TTTTTTTTTT"), id="ref")

        with mock.patch.dict(sys.modules, {"edlib": fake_module}):
            with self.assertRaises(ValueError):
                tools.edlib_hw_align(test_seq, ref_seq)

    def test_edlib_hw_align_raises_for_out_of_range_location(self):
        fake_module = self._fake_edlib_module(
            [{"editDistance": 0, "locations": [(2, 20)]}],
        )
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("TTTTTTTTTT"), id="ref")

        with mock.patch.dict(sys.modules, {"edlib": fake_module}):
            with self.assertRaises(ValueError):
                tools.edlib_hw_align(test_seq, ref_seq)

    def test_edlib_hw_align_raises_when_window_alignment_fails(self):
        fake_module = self._fake_edlib_module(
            [
                {"editDistance": 0, "locations": [(0, 3)]},
                {"editDistance": -1, "locations": []},
            ],
        )
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("TTTTTTTTTT"), id="ref")

        with mock.patch.dict(sys.modules, {"edlib": fake_module}):
            with self.assertRaises(ValueError):
                tools.edlib_hw_align(test_seq, ref_seq)


class TestParasailNwAlign(TestCase):
    def test_parasail_nw_align_uses_cigar_from_traceback(self):
        fake_result = mock.Mock()
        fake_result.cigar = types.SimpleNamespace(decode=b"4=")

        fake_module = types.SimpleNamespace()
        fake_module.matrix_create = mock.Mock(return_value="matrix")
        fake_module.nw_trace_scan_16 = mock.Mock(return_value=fake_result)

        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("ACGT"), id="ref")

        with mock.patch.dict(sys.modules, {"parasail": fake_module}):
            result = tools.parasail_nw_align(test_seq, ref_seq)

        self.assertEqual(result, ("ACGT", "ACGT"))
        fake_module.matrix_create.assert_called_once_with(tools.DNA_ALPHABET, 2, -1)
        fake_module.nw_trace_scan_16.assert_called_once_with(
            "ACGT", "ACGT", 5, 1, "matrix"
        )


class TestMafftAlign(TestCase):
    def test_mafft_align_raises_when_no_binary_available(self):
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("ACGT"), id="ref")

        with mock.patch("pyhiv.align.tools.get_mafft_binary", return_value=None):
            with self.assertRaises(RuntimeError):
                tools.mafft_align(test_seq, ref_seq)

    def test_mafft_align_rejects_empty_test_sequence(self):
        test_seq = SeqRecord(Seq(""), id="test")
        ref_seq = SeqRecord(Seq("ACGT"), id="ref")

        with mock.patch("pyhiv.align.tools.get_mafft_binary", return_value="mafft"):
            with mock.patch("pyhiv.align.tools.subprocess.run") as mock_run:
                with self.assertRaises(ValueError):
                    tools.mafft_align(test_seq, ref_seq)
                mock_run.assert_not_called()

    def test_mafft_align_rejects_empty_reference_sequence(self):
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq(""), id="ref")

        with mock.patch("pyhiv.align.tools.get_mafft_binary", return_value="mafft"):
            with mock.patch("pyhiv.align.tools.subprocess.run") as mock_run:
                with self.assertRaises(ValueError):
                    tools.mafft_align(test_seq, ref_seq)
                mock_run.assert_not_called()

    def test_mafft_align_parses_successful_output(self):
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("AC-T"), id="ref")

        fake_completed = mock.Mock(
            returncode=0,
            stdout=">query\nACGT\n>reference\nAC-T\n",
            stderr="",
        )

        with mock.patch("pyhiv.align.tools.get_mafft_binary", return_value="mafft"):
            with mock.patch(
                "pyhiv.align.tools.subprocess.run", return_value=fake_completed
            ) as mock_run:
                result = tools.mafft_align(test_seq, ref_seq)

        self.assertEqual(result, ("ACGT", "AC-T"))
        called_command = mock_run.call_args.args[0]
        self.assertEqual(called_command[0], "mafft")

    def test_mafft_align_raises_runtime_error_on_failure(self):
        test_seq = SeqRecord(Seq("ACGT"), id="test")
        ref_seq = SeqRecord(Seq("ACGT"), id="ref")

        fake_completed = mock.Mock(returncode=1, stdout="", stderr="boom")

        with mock.patch("pyhiv.align.tools.get_mafft_binary", return_value="mafft"):
            with mock.patch(
                "pyhiv.align.tools.subprocess.run", return_value=fake_completed
            ):
                with self.assertRaisesRegex(RuntimeError, "boom"):
                    tools.mafft_align(test_seq, ref_seq)

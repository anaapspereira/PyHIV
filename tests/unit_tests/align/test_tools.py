import os
import stat
import tempfile
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

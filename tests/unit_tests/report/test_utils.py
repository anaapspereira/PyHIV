import tempfile
from pathlib import Path
from unittest import TestCase

from pyhiv.report.utils import (
    ungap,
    first_last_nongap_idx,
    parse_present_regions,
    parse_features,
    build_ref_to_alignment_map,
    project_features_to_alignment,
    canon_label,
    is_special_reference,
    get_numeric_offsets_non_special,
    read_alignment_fasta,
    build_alignment_path,
    normalize_features,
    normalize_present_regions,
)
from pyhiv.report.constants import NumericOffsets, K03455Config


class TestReportingUtils(TestCase):

    def setUp(self):
        self.tmp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        for f in self.tmp_dir.glob("*"):
            try:
                f.unlink()
            except Exception:
                pass

    def test_ungap_and_first_last(self):
        self.assertEqual(ungap("A-C.G-"), "ACG")
        self.assertEqual(first_last_nongap_idx("--A.C-"), (2, 4))
        self.assertEqual(first_last_nongap_idx("----"), (0, 0))
        self.assertEqual(first_last_nongap_idx(""), (0, 0))

    def test_parse_present_regions(self):
        cases = [
            ("['gag','pol']", ["gag", "pol"]),
            ("gag, pol", ["gag", "pol"]),
            ("-", []),
            (None, []),
        ]
        for cell, expected in cases:
            self.assertEqual(parse_present_regions(cell), expected)

        # Fallback branch: invalid Python literal
        self.assertEqual(parse_present_regions("[gag,pol"), ["[gag", "pol"])

    def test_parse_features_and_projection(self):
        cell = {"gag": (1, 10), "pol": (20, 30)}
        features = parse_features(cell)
        self.assertEqual(features, {"gag": (1, 10), "pol": (20, 30)})

        ref = "ACGT" * 10  # no gaps
        mapping, _ = build_ref_to_alignment_map(ref)
        self.assertEqual(mapping[1], 0)
        self.assertEqual(mapping[10], 9)

        projected = project_features_to_alignment(features, mapping)
        self.assertEqual(projected["gag"], (0, 9))

    def test_parse_features_none_and_nan(self):
        self.assertEqual(parse_features(None), {})
        self.assertEqual(parse_features(float("nan")), {})
        self.assertEqual(parse_features("None"), {})

    def test_read_alignment_fasta_file_not_found(self):
        missing = self.tmp_dir / "missing.fasta"
        with self.assertRaises(FileNotFoundError):
            read_alignment_fasta(missing)

    def test_read_alignment_fasta_invalid_or_incomplete(self):
        f = self.tmp_dir / "bad.fasta"
        # Empty file
        f.write_text("")
        with self.assertRaises(ValueError):
            read_alignment_fasta(f)
        # Only one sequence
        f.write_text(">reference\nACGT\n")
        with self.assertRaises(ValueError):
            read_alignment_fasta(f)

    def test_read_alignment_fasta_reference_in_second_header(self):
        f = self.tmp_dir / "ok.fasta"
        f.write_text(">user\nTTTT\n>reference genome\nAAAA\n")
        ref_header, ref_seq, user_header, user_seq = read_alignment_fasta(f)
        self.assertIn("reference", ref_header.lower())
        self.assertEqual(user_seq, "TTTT")

    def test_canon_label_and_special_ref(self):
        self.assertEqual(canon_label("5 ltr"), "5' LTR")
        self.assertEqual(canon_label("tat i"), "tat 1")
        self.assertEqual(canon_label("REV II"), "rev 2")
        self.assertEqual(canon_label("env"), "env")

        self.assertTrue(is_special_reference("K03455", "ref header"))
        self.assertTrue(is_special_reference("X", "K03455-B something"))
        self.assertFalse(is_special_reference("X", "Y"))
        # Empty args branch
        self.assertFalse(is_special_reference("", ""))

    def test_normalize_features_special_and_non_special(self):
        feats = {"tat I": (1, 5), "unknown": (10, 20)}
        non_special = normalize_features(feats, special=False)
        self.assertIn("unknown", non_special)

        special = normalize_features(feats, special=True)
        self.assertIn("tat 1", special)
        self.assertNotIn("unknown", special)

    def test_normalize_present_regions_special_and_non_special(self):
        regions = ["tat i", "pol", "other"]
        non_special = normalize_present_regions(regions, special=False)
        self.assertEqual(non_special, regions)

        special = normalize_present_regions(regions, special=True)
        self.assertIn("tat 1", special)
        self.assertNotIn("other", special)

    def test_get_numeric_offsets_non_special(self):
        self.assertEqual(
            get_numeric_offsets_non_special("tat 1"),
            NumericOffsets.GENE_OFFSET_MAP["tat 1"],
        )
        self.assertEqual(
            get_numeric_offsets_non_special("rev 1"),
            NumericOffsets.GENE_OFFSET_MAP["rev 1"],
        )
        self.assertEqual(
            get_numeric_offsets_non_special("unknown"),
            NumericOffsets.DEFAULT_OFFSETS,
        )
        self.assertEqual(
            get_numeric_offsets_non_special(None),
            NumericOffsets.DEFAULT_OFFSETS,
        )

    def test_build_alignment_path_fallback(self):
        seq = "SEQX"
        base = self.tmp_dir
        best = base / f"best_alignment_{seq}.fasta"
        alt = base / f"{seq}.fasta"
        alt.touch()  # fallback file exists
        result = build_alignment_path(seq, base)
        self.assertEqual(result, alt)

    def test_get_k03455_offsets(self):
        """Ensure K03455Config.get_k03455_offsets returns correct offsets or default."""
        # Known gene
        self.assertEqual( K03455Config.get_k03455_offsets("tat 1"), K03455Config.K03455_NUMERIC_OFFSETS["tat 1"], )
        # Another known gene
        self.assertEqual( K03455Config.get_k03455_offsets("vpu"), K03455Config.K03455_NUMERIC_OFFSETS["vpu"], )
        # Unknown gene → default
        self.assertEqual( K03455Config.get_k03455_offsets("foo"), K03455Config.DEFAULT_K03455_OFFSETS, )


    def test_get_numeric_offsets_non_special_monkeypatched(self):
        from pyhiv.report import utils

        called = {}
        original_get_offsets = utils.NumericOffsets.get_offsets
        try:
            def fake_get_offsets(gene):
                called["gene"] = gene
                return (1.1, 2.2)

            utils.NumericOffsets.get_offsets = fake_get_offsets
            result = utils.get_numeric_offsets_non_special("abc")
            self.assertEqual(result, (1.1, 2.2))
            self.assertEqual(called["gene"], "abc")
        finally:
            utils.NumericOffsets.get_offsets = original_get_offsets

    def test_read_alignment_fasta_with_blank_lines(self):
        # Ensure lines with blanks are skipped properly
        f = self.tmp_dir / "weird.fasta"
        f.write_text("\n>ref reference\nAAAA\n\n>user\nTTTT\n\n")
        from pyhiv.report.utils import read_alignment_fasta
        ref_h, ref_s, user_h, user_s = read_alignment_fasta(f)
        self.assertIn("ref", ref_h.lower())
        self.assertEqual(ref_s, "AAAA")
        self.assertEqual(user_s, "TTTT")

    def test_canon_label_case_insensitive_and_no_match(self):
        from pyhiv.report.utils import canon_label
        self.assertEqual(canon_label("ENV"), "env")  # case-insensitive match
        self.assertIsNone(canon_label("madeup"))  # unmatched returns None


    def test_build_alignment_path_best_exists(self):
        from pyhiv.report.utils import build_alignment_path
        best = self.tmp_dir / "best_alignment_SEQY.fasta"
        best.touch()
        result = build_alignment_path("SEQY", self.tmp_dir)
        self.assertEqual(result, best)

    def test_read_alignment_fasta_no_reference_in_headers_defaults_first(self):
        """Covers branch where neither header contains 'reference' -> idx_ref=0."""
        f = self.tmp_dir / "no_ref.fasta"
        f.write_text(">alpha\nAAAA\n>beta\nTTTT\n")
        ref_h, ref_s, user_h, user_s = read_alignment_fasta(f)
        self.assertEqual(ref_h, "alpha")
        self.assertEqual(ref_s, "AAAA")
        self.assertEqual(user_h, "beta")
        self.assertEqual(user_s, "TTTT")

    def test_parse_features_from_string_literal(self):
        """Covers ast.literal_eval path in parse_features."""
        cell = "{'gag': (1, 10), 'pol': (20, 30)}"
        out = parse_features(cell)
        self.assertEqual(out, {"gag": (1, 10), "pol": (20, 30)})

    def test_project_features_to_alignment_exclusions(self):
        """Covers branches where (a) endpoints missing and (b) aend < astart."""
        ref = "ACGT" * 5  # positions 1..20 map to indices 0..19
        ref_map, _ = build_ref_to_alignment_map(ref)

        feats = {
            "ok": (3, 7),        # included
            "missing": (100, 101),  # start/end not in map -> excluded
            "reversed": (10, 1), # both present but aend < astart -> excluded
        }
        proj = project_features_to_alignment(feats, ref_map)
        self.assertIn("ok", proj)
        self.assertNotIn("missing", proj)
        self.assertNotIn("reversed", proj)
        self.assertEqual(proj["ok"], (ref_map[3], ref_map[7]))

    def test_canon_label_targets_without_regex_via_monkeypatch(self):
        """
        Force canon_label to skip regex patterns to exercise:
          - direct TARGET_REGIONS hit
          - case-insensitive TARGET_REGIONS hit
          - final None return
          - early empty/whitespace return
        """
        from pyhiv.report import utils
        original_patterns = utils._CANON_PATTERNS
        try:
            utils._CANON_PATTERNS = []  # bypass regex block entirely

            # Pick a real target from config
            some_target = next(iter(K03455Config.TARGET_REGIONS))

            # Direct exact match branch
            self.assertEqual(utils.canon_label(some_target), some_target)

            # Case-insensitive match branch
            self.assertEqual(utils.canon_label(some_target.upper()), some_target)

            # Empty/whitespace -> early None
            self.assertIsNone(utils.canon_label("   "))

            # Completely unknown -> final None
            self.assertIsNone(utils.canon_label("not-a-real-target"))
        finally:
            utils._CANON_PATTERNS = original_patterns

    def test_parse_present_regions_non_list_non_str_hits_final_return(self):
        self.assertEqual(parse_present_regions("123"), [])

    def test_canon_label_case_insensitive_target_return_branch(self):
        from pyhiv.report import utils
        original_patterns = utils._CANON_PATTERNS
        original_targets = K03455Config.TARGET_REGIONS
        try:
            # Disable regex so we don't exit early via regex matches
            utils._CANON_PATTERNS = []

            # Use a synthetic target that won't match any regex but *will* match case-insensitively
            K03455Config.TARGET_REGIONS = ("FooBar",)

            # Lowercase input ensures exact-match check fails, then loop hits case-insensitive equality -> return target
            self.assertEqual(utils.canon_label("foobar"), "FooBar")
        finally:
            utils._CANON_PATTERNS = original_patterns
            K03455Config.TARGET_REGIONS = original_targets

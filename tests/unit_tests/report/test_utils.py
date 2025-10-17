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
)
from pyhiv.report.constants import NumericOffsets


class TestReportingUtils(TestCase):

    def test_ungap_and_first_last(self):
        self.assertEqual(ungap("A-C.G-"), "ACG")
        self.assertEqual(first_last_nongap_idx("--A.C-"), (2, 4))
        self.assertEqual(first_last_nongap_idx("----"), (0, 0))

    def test_parse_present_regions(self):
        cases = [
            ("['gag','pol']", ["gag", "pol"]),
            ("gag, pol", ["gag", "pol"]),
            ("-", []),
            (None, []),
        ]
        for cell, expected in cases:
            self.assertEqual(parse_present_regions(cell), expected)

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

    def test_canon_label_and_special_ref(self):
        self.assertEqual(canon_label("5 ltr"), "5' LTR")
        self.assertEqual(canon_label("tat i"), "tat 1")
        self.assertEqual(canon_label("REV II"), "rev 2")
        self.assertEqual(canon_label("env"), "env")

        self.assertTrue(is_special_reference("K03455", "ref header"))
        self.assertTrue(is_special_reference("X", "K03455-B something"))
        self.assertFalse(is_special_reference("X", "Y"))

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

from unittest import TestCase, mock

from pyhiv.split.split import get_gene_region, get_present_gene_regions


class TestGeneRegions(TestCase):

    @mock.patch("pyhiv.align.calculate_alignment_score")
    def test_get_gene_region_empty_ranges(self, mock_score):
        """Should return empty list if gene_ranges is empty."""
        result = get_gene_region("AAA", "AAA", {})
        self.assertEqual(result, [])
        mock_score.assert_not_called()

    @mock.patch("pyhiv.split.split.calculate_alignment_score")
    def test_get_gene_region_single_max(self, mock_score):
        mock_score.side_effect = [5, 10, 3]  # scores for three genes
        gene_ranges = {"g1": (1, 3), "g2": (1, 3), "g3": (1, 3)}
        result = get_gene_region("AAA", "AAA", gene_ranges)
        self.assertEqual(result, ["g2"])
        self.assertEqual(mock_score.call_count, 3)

    @mock.patch("pyhiv.split.split.calculate_alignment_score")
    def test_get_gene_region_multiple_max(self, mock_score):
        mock_score.side_effect = [10, 10, 5]
        gene_ranges = {"g1": (1, 3), "g2": (1, 3), "g3": (1, 3)}
        result = get_gene_region("AAA", "AAA", gene_ranges)
        self.assertCountEqual(result, ["g1", "g2"])

    @mock.patch("pyhiv.align.calculate_alignment_score")
    def test_get_gene_region_max_score_none(self, mock_score):
        """Should handle None max_score gracefully."""
        mock_score.side_effect = []
        gene_ranges = {}
        result = get_gene_region("AAA", "AAA", gene_ranges)
        self.assertEqual(result, [])

    def test_get_present_gene_regions_all_gaps(self):
        """Should return empty list if all genes contain only '-'."""
        test_aligned = "------"
        gene_ranges = {"g1": (1, 3), "g2": (4, 6)}
        result = get_present_gene_regions(test_aligned, gene_ranges)
        self.assertEqual(result, [])

    def test_get_present_gene_regions_partial_bases(self):
        """Should return genes that have at least one base."""
        test_aligned = "A--G--"
        gene_ranges = {"g1": (1, 3), "g2": (4, 6)}
        result = get_present_gene_regions(test_aligned, gene_ranges)
        self.assertEqual(result, ["g1", "g2"])

    def test_get_present_gene_regions_empty_ranges(self):
        """Should return empty list if gene_ranges is empty."""
        test_aligned = "AAA"
        result = get_present_gene_regions(test_aligned, {})
        self.assertEqual(result, [])

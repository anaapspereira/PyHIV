from unittest import TestCase, mock

from pyhiv.split.split import get_gene_region, get_present_gene_regions, map_ref_coords_to_alignment


class TestSplitting(TestCase):

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
        gene_ranges = {"g1": (0, 2), "g2": (3, 5)}
        result = get_present_gene_regions(test_aligned, gene_ranges)
        self.assertEqual(result, ["g1", "g2"])

    def test_get_present_gene_regions_empty_ranges(self):
        """Should return empty list if gene_ranges is empty."""
        test_aligned = "AAA"
        result = get_present_gene_regions(test_aligned, {})
        self.assertEqual(result, [])


class TestMapRefCoordsToAlignment(TestCase):

    def test_basic_mapping(self):
        """Should correctly map reference coordinates ignoring gaps."""
        ref_aligned = "A-CGT--A"
        # Positions without gaps: A(1), C(2), G(3), T(4), A(5)
        # Corresponding alignment indices (0-based): 0, 2, 3, 4, 7
        expected = {1: 0, 2: 2, 3: 3, 4: 4, 5: 7}
        result = map_ref_coords_to_alignment(ref_aligned)
        self.assertEqual(result, expected)

    def test_all_gaps(self):
        """Should return empty dict when all positions are gaps."""
        ref_aligned = "-----"
        result = map_ref_coords_to_alignment(ref_aligned)
        self.assertEqual(result, {})

    def test_empty_input(self):
        """Should handle empty input string gracefully."""
        result = map_ref_coords_to_alignment("")
        self.assertEqual(result, {})

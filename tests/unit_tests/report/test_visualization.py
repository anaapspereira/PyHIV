from unittest import TestCase

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pyhiv.report.visualization import plot_gene_axes


class TestReportingVisualization(TestCase):

    def test_plot_gene_axes_basic(self):
        fig, ax = plt.subplots(figsize=(4, 3))
        genes = {"gag": (10, 50), "pol": (60, 120)}
        plot_gene_axes(ax, genes_ranges=genes, alignment_start=15, alignment_end=100, y_positions=None)
        # should have one patch per gene
        self.assertGreaterEqual(len(ax.patches), 2)
        plt.close(fig)

    def test_plot_gene_axes_k03455_connectors(self):
        # y_positions provided => K03455 mode
        fig, ax = plt.subplots(figsize=(4, 3))
        genes = {
            "tat 1": (10, 20),
            "tat 2": (25, 35),
            "rev 1": (40, 50),
            "rev 2": (55, 65),
            "gag": (70, 100),
            "pol": (110, 140),
            "vif": (150, 180),
            "vpu": (190, 220),
            "vpr": (230, 260),
            "env": (270, 300),
            "nef": (310, 340),
            "5' LTR": (350, 360),
            "3' LTR": (370, 380),
            "unknown": (390, 400),
        }
        y_positions = {g: i * 0.1 for i, g in enumerate(genes)}
        plot_gene_axes(ax, genes, alignment_start=10, alignment_end=380, y_positions=y_positions)
        # connectors for tat and rev should be drawn
        self.assertGreater(len(ax.lines), 0)
        self.assertGreater(len(ax.patches), 0)
        plt.close(fig)

    def test_plot_gene_axes_single_connector(self):
        # only one tat gene => single connector path
        fig, ax = plt.subplots(figsize=(3, 2))
        genes = {"tat 1": (5, 15)}
        y_positions = {"tat 1": 0.1}
        plot_gene_axes(ax, genes, alignment_start=5, alignment_end=15, y_positions=y_positions)
        # one patch should exist
        self.assertEqual(len(ax.patches), 1)
        plt.close(fig)

    def test_plot_gene_axes_alignment_point(self):
        # alignment_start == alignment_end
        fig, ax = plt.subplots(figsize=(3, 2))
        genes = {"gag": (10, 50)}
        plot_gene_axes(ax, genes_ranges=genes, alignment_start=50, alignment_end=50, y_positions=None)
        self.assertGreaterEqual(len(ax.patches), 1)
        plt.close(fig)

    def test_plot_gene_axes_missing_gene_in_y_positions(self):
        fig, ax = plt.subplots(figsize=(3, 2))
        genes = {"gag": (10, 50), "pol": (60, 120)}
        # y_positions missing one gene → triggers "continue"
        y_positions = {"gag": 0.1}
        plot_gene_axes(ax, genes_ranges=genes, alignment_start=10, alignment_end=100, y_positions=y_positions)
        plt.close(fig)

    def test_plot_gene_axes_empty_genes(self):
        fig, ax = plt.subplots(figsize=(3, 2))
        # No genes → triggers x_min_genes, x_max_genes = 0, 0
        plot_gene_axes(ax, genes_ranges={}, alignment_start=0, alignment_end=0, y_positions=None)
        plt.close(fig)

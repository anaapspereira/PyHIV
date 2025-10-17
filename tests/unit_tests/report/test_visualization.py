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
        self.assertGreaterEqual(len(ax.patches), 2)
        plt.close(fig)



"""
PDF report generation for PyHIV results.
"""

import textwrap
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

from .constants import PageLayout, MetadataConfig
from .utils import ungap, first_last_nongap_idx
from .visualization import plot_gene_axes


def render_sequence_page(
    pdf: PdfPages,
    sequence: str,
    accession: str,
    subtype: str,
    mm_region: str,
    present_regions: List[str],
    features_aln: Dict[str, Tuple[int, int]],
    ref_seq_aligned: str,
    user_seq_aligned: str,
    y_positions: Optional[Dict[str, float]] = None
):
    """
    Render a single sequence page in the PDF report.

    Parameters
    ----------
    pdf : PdfPages
        The PdfPages object to save the figure into.
    sequence : str
        The name or identifier of the sequence.
    accession : str
        The accession number of the sequence.
    subtype : str
        The subtype of the sequence.
    mm_region : str
        The most matching region of the sequence.
    present_regions : List[str]
        List of present regions in the sequence.
    features_aln : Dict[str, Tuple[int, int]]
        Dictionary of gene features with their alignment coordinate ranges.
    ref_seq_aligned : str
        The reference sequence aligned (with gaps).
    user_seq_aligned : str
        The user's sequence aligned (with gaps).
    y_positions : Optional[Dict[str, float]], optional
        Fixed y-positions for gene lanes, by default None (auto lanes).
    """
    fig = plt.figure(figsize=PageLayout.FIGSIZE)
    gs = GridSpec(
        2, 1,  # one column
        height_ratios=PageLayout.GRID_HEIGHT_RATIOS,
        figure=fig,
        hspace=PageLayout.HSPACE
    )

    # Top: metadata (wrapped) — no bar chart
    ax_meta = fig.add_subplot(gs[0, 0])
    ax_meta.axis("off")
    ax_meta.text(0.5, MetadataConfig.TITLE_Y, f"PyHIV Report — {sequence}", ha="center", va="top",
                 fontsize=16, fontweight="bold", transform=ax_meta.transAxes)

    ref_len_nt = len(ungap(ref_seq_aligned))
    usr_len_nt = len(ungap(user_seq_aligned))

    meta_lines = [
        f"Sequence: {sequence}",
        f"Subtype: {subtype}",
        f"Reference Accession: {accession}",
        f"Most matching region: {mm_region or '-'}",
        f"Present regions ({len(present_regions)}): {', '.join(present_regions) if present_regions else '-'}",
        f"Lengths — nt (no gaps) Ref|Seq: {ref_len_nt} | {usr_len_nt}",
    ]
    wrapped = "\n".join(textwrap.fill(l, width=MetadataConfig.WRAP, subsequent_indent='   ') for l in meta_lines)

    ax_meta.text(0.0, MetadataConfig.INFO_TOP_Y, wrapped, ha="left", va="top",
                 family="monospace", fontsize=MetadataConfig.FONTSIZE, transform=ax_meta.transAxes, wrap=True)

    # Bottom: gene panel (full width)
    ax_map = fig.add_subplot(gs[1, 0])

    a_start, a_end = first_last_nongap_idx(user_seq_aligned)

    plot_gene_axes(
        ax=ax_map,
        genes_ranges=features_aln,
        alignment_start=a_start,
        alignment_end=a_end,
        y_positions=y_positions,
    )

    pdf.savefig(fig)
    plt.close(fig)

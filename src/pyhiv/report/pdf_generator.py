"""
PDF report generation for PyHIV results.
"""

import textwrap
from typing import Dict, List, Optional

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
    features_aln: Dict[str, tuple],
    ref_header: str,
    ref_seq_aligned: str,
    user_header: str,
    user_seq_aligned: str,
    y_positions: Optional[Dict[str, float]] = None
):
    """Render a single sequence page in the PDF report."""
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
        f"Accession: {accession}",
        f"Most matching region: {mm_region or '-'}",
        f"Present regions ({len(present_regions)}): {', '.join(present_regions) if present_regions else '-'}",
        f"Lengths — nt (no gaps) Ref|Seq: {ref_len_nt} | {usr_len_nt}",
    ]
    wrapped = "\n".join(textwrap.fill(l, width=MetadataConfig.WRAP, subsequent_indent='   ') for l in meta_lines)

    ax_meta.text(0.0, MetadataConfig.INFO_TOP_Y, wrapped, ha="left", va="top",
                 family="monospace", fontsize=MetadataConfig.FONTSIZE, transform=ax_meta.transAxes, wrap=True)

    # Bottom: gene panel (full width)
    ax_map = fig.add_subplot(gs[1, 0])

    a_start, a_end = first_last_nongap_idx(user_seq_aligned)  # alignment columns of the user's non-gap span

    plot_gene_axes(
        ax=ax_map,
        genes_ranges=features_aln,     # already in ALIGNMENT coordinates
        alignment_start=a_start,
        alignment_end=a_end,
        y_positions=y_positions,       # fixed dict (K03455) or None for auto lanes
    )

    pdf.savefig(fig)
    plt.close(fig)

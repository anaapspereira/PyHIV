"""
Gene visualization and plotting functions for PyHIV reporting module.
"""

from typing import Dict, List, Tuple, Optional
import matplotlib.patches as mpatches

from .constants import GenePanelConfig
from .utils import get_numeric_offsets_non_special


def _assign_lanes_nonoverlap(items: List[Tuple[str, Tuple[int, int]]], min_gap: int = 50) -> Dict[str, float]:
    """
    Greedy interval graph coloring to avoid overlaps (only for non-K03455).
    items: list[(gene, (start, end))] sorted by start
    min_gap: minimal separation (alignment columns) to share a lane
    Returns: dict[gene] -> y_base in [0..1]
    """
    lanes_last_end = []   # last end per lane
    y_lane = {}           # gene -> lane idx

    for gene, (s, e) in items:
        placed = False
        for li, last_end in enumerate(lanes_last_end):
            if s >= last_end + min_gap:
                lanes_last_end[li] = e
                y_lane[gene] = li
                placed = True
                break
        if not placed:
            lanes_last_end.append(e)
            y_lane[gene] = len(lanes_last_end) - 1

    base, step = 0.2, 0.28
    return {g: base + lane * step for g, lane in y_lane.items()}


def plot_gene_axes(
    ax,
    genes_ranges: Dict[str, Tuple[int, int]],
    alignment_start: int,
    alignment_end: int,
    y_positions: Optional[Dict[str, float]] = None,
):
    """Plot gene visualization with alignment information."""
    items = sorted(genes_ranges.items(), key=lambda x: x[1][0])

    # K03455: keep fixed map; others: compute non-overlapping lanes
    non_special = y_positions is None
    if non_special:
        y_positions = _assign_lanes_nonoverlap(items, min_gap=50)

    y_positions_scaled = {g: y * GenePanelConfig.Y_SCALE for g, y in y_positions.items()}

    tat_centers, rev_centers = [], []

    for gene, (start, end) in items:
        if gene not in y_positions_scaled:
            continue
        y = y_positions_scaled[gene]
        height = 0.15
        rect = mpatches.Rectangle(
            (start, y - height / 2),
            max(1, end - start),
            height,
            color="#778da9",
            ec="#6c757d",
        )
        ax.add_patch(rect)

        cx = (start + end) / 2
        gl = gene.lower()
        is_tat = gl.startswith("tat")
        is_rev = gl.startswith("rev")

        # Label inside rectangle:
        # - non-K03455: all genes (including tat/rev)
        # - K03455: all except tat/rev (they use connectors)
        if non_special or not (is_tat or is_rev):
            ax.text(cx, y, gene, ha="center", va="center",
                    fontsize=8, fontweight="bold", color="white")
        else:
            # K03455: collect centers for connectors
            if is_tat:
                tat_centers.append((cx, y))
            if is_rev:
                rev_centers.append((cx, y))

        # --- numeric label offsets (start/end) ---
        if non_special:
            # NON-K03455 → configurable map
            offset_down, offset_up = get_numeric_offsets_non_special(gene)
        else:
            # K03455 → original behavior (unchanged)
            if gene == "tat 1":       offset_down, offset_up = -0.15, 0.15
            elif gene == "tat 2":     offset_down, offset_up = -0.15, 0.20
            elif gene == "rev 1":     offset_down, offset_up =  0.20, 0.25
            elif gene == "rev 2":     offset_down, offset_up =  0.20, 0.25
            elif gene == "gag":       offset_down, offset_up =  0.30, 0.20
            elif gene == "vif":       offset_down, offset_up = -0.15, -0.20
            elif gene == "5' LTR":    offset_down, offset_up = -0.15, -0.20
            elif gene == "pol":       offset_down, offset_up = -0.15, -0.20
            elif gene == "vpu":       offset_down, offset_up =  0.35, 0.15
            elif gene == "env":       offset_down, offset_up = -0.15, -0.20
            elif gene == "vpr":       offset_down, offset_up =  0.20, -0.20
            elif gene == "nef":       offset_down, offset_up = -0.15, -0.20
            elif gene == "3' LTR":    offset_down, offset_up =  0.25, -0.20
            else:                      offset_down, offset_up = -0.15, 0.35

        ax.annotate(str(start), xy=(start, y), xycoords="data",
                    xytext=(start, y + offset_down), textcoords="data",
                    ha="center", va="top", fontsize=7, color="#1b263b", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#6c757d", lw=0.5))
        ax.annotate(str(end), xy=(end, y), xycoords="data",
                    xytext=(end, y + offset_up), textcoords="data",
                    ha="center", va="bottom", fontsize=7, color="#1b263b", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#6c757d", lw=0.5))

    # ----- tat/rev connectors only for K03455 -----
    def draw_connector(centers, extra_height, label):
        if len(centers) < 1:
            return None
        if len(centers) >= 2:
            (x1, y1), (x2, y2) = centers[0], centers[-1]
            y_line = max(y1, y2) + extra_height
            ax.plot([x1, x2], [y_line, y_line], color="#778da9", lw=1.2)
            ax.plot([x1, x1], [y1, y_line], color="#778da9", lw=1.2)
            ax.plot([x2, x2], [y2, y_line], color="#778da9", lw=1.2)
            ax.text((x1 + x2) / 2, y_line + 0.0, label, ha="center", va="bottom",
                    fontsize=9, fontweight="bold", color="#1b263b")
            return y_line
        else:
            (x, y) = centers[0]
            y_line = y + extra_height
            ax.text(x, y_line, label, ha="center", va="bottom",
                    fontsize=9, fontweight="bold", color="#1b263b")
            return y_line

    y_line_tat = draw_connector(tat_centers, GenePanelConfig.TAT_CONNECTOR, "tat") if not non_special else None
    y_line_rev = draw_connector(rev_centers, GenePanelConfig.REV_CONNECTOR, "rev") if not non_special else None

    # ----- red alignment span -----
    tops = list(y_positions_scaled.values()) if y_positions_scaled else [0.0]
    if y_line_tat is not None: tops.append(y_line_tat)
    if y_line_rev is not None: tops.append(y_line_rev)
    y_alignment = max(tops) + GenePanelConfig.ALIGNMENT_CLEARANCE

    ax.plot([alignment_start, alignment_end], [y_alignment, y_alignment], color="#eb5e28", linewidth=4)

    if alignment_start == alignment_end:
        # single point
        ax.annotate(str(alignment_end), xy=(alignment_end, y_alignment), xycoords="data",
                    xytext=(alignment_end, y_alignment + 0.12), textcoords="data",
                    ha="center", va="bottom", fontsize=7, color="#eb5e28", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#eb5e28", lw=0.5))
    else:
        # start below, end above (reduces collision)
        ax.annotate(str(alignment_start), xy=(alignment_start, y_alignment), xycoords="data",
                    xytext=(alignment_start, y_alignment - 0.15), textcoords="data",
                    ha="center", va="bottom", fontsize=7, color="#eb5e28", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#eb5e28", lw=0.5))
        ax.annotate(str(alignment_end), xy=(alignment_end, y_alignment), xycoords="data",
                    xytext=(alignment_end, y_alignment + 0.12), textcoords="data",
                    ha="center", va="bottom", fontsize=7, color="#eb5e28", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#eb5e28", lw=0.5))

    # ----- baseline & limits -----
    if items:
        xs = [start for _, (start, _) in items] + [end for _, (_, end) in items]
        x_min_genes, x_max_genes = min(xs), max(xs)
    else:
        x_min_genes, x_max_genes = 0, 0

    if non_special:
        # Force the visible domain to reach at least 10,000 or alignment_end
        hard_max = max(GenePanelConfig.NON_K03455_X_MAX_DEFAULT, alignment_end, x_max_genes)
        pad_x = max(GenePanelConfig.X_PAD_MIN, int(0.02 * max(1, hard_max)))  # keep some air to the right
        ax.set_xlim(0, hard_max + pad_x)
        baseline_y = (min(y_positions_scaled.values()) if y_positions_scaled else 0.0) - 0.32
        ax.plot([0, hard_max], [baseline_y, baseline_y], color="#778da9")
    else:
        # K03455: keep auto-fit
        x_min = x_min_genes
        x_max = x_max_genes
        width = max(1, x_max - x_min)
        pad_x = max(GenePanelConfig.X_PAD_MIN, int(0.02 * width))
        ax.set_xlim(x_min - pad_x, x_max + pad_x)
        baseline_y = (min(y_positions_scaled.values()) if y_positions_scaled else 0.0) - 0.32
        ax.plot([x_min, x_max], [baseline_y, baseline_y], color="#778da9")

    y_min = (min(y_positions_scaled.values()) if y_positions_scaled else 0.0) - GenePanelConfig.BOTTOM_MARGIN
    y_max = y_alignment + GenePanelConfig.TOP_MARGIN
    ax.set_ylim(y_min, y_max)
    ax.axis("off")

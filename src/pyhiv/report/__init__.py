"""
PyHIV reporting module.

This module provides functionality to generate PDF reports from PyHIV analysis results.
The reports include sequence metadata and gene visualization plots.
"""

from .reporter import PyHIVReporter
from .constants import (
    PageLayout, MetadataConfig, GenePanelConfig,
    NumericOffsets, K03455Config,
)
from .utils import (
    ungap, first_last_nongap_idx, read_alignment_fasta,
    parse_present_regions, parse_features, is_special_reference,
    canon_label, normalize_features, normalize_present_regions,
    build_ref_to_alignment_map, project_features_to_alignment,
    get_numeric_offsets_non_special, build_alignment_path
)
from .visualization import plot_gene_axes, _assign_lanes_nonoverlap
from .pdf_generator import render_sequence_page

__all__ = [
    "PyHIVReporter",
    # Config classes
    "PageLayout", "MetadataConfig", "GenePanelConfig",
    "NumericOffsets", "K03455Config",
    # Utilities
    "ungap", "first_last_nongap_idx", "read_alignment_fasta",
    "parse_present_regions", "parse_features", "is_special_reference",
    "canon_label", "normalize_features", "normalize_present_regions",
    "build_ref_to_alignment_map", "project_features_to_alignment",
    "get_numeric_offsets_non_special", "build_alignment_path",
    # Viz + rendering
    "plot_gene_axes", "_assign_lanes_nonoverlap",
    "render_sequence_page",
]

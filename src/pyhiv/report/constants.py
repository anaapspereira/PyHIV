"""
Constants and configuration for PyHIV reporting module.
"""

# ================================
# PAGE LAYOUT CONFIGURATION
# ================================
class PageLayout:
    """Page layout and spacing configuration."""
    FIGSIZE = (11.69, 9.2)
    GRID_HEIGHT_RATIOS = [0.9, 1.9]
    HSPACE = 0.42

# ================================
# METADATA BLOCK CONFIGURATION
# ================================
class MetadataConfig:
    """Metadata block display configuration."""
    WRAP = 75
    FONTSIZE = 9.5
    TITLE_Y = 1.06
    INFO_TOP_Y = 0.78

# ================================
# GENE PANEL CONFIGURATION
# ================================
class GenePanelConfig:
    """Gene panel spacing and positioning configuration."""
    Y_SCALE = 1.35
    TAT_CONNECTOR = 0.07 * Y_SCALE
    REV_CONNECTOR = 0.45 * Y_SCALE
    ALIGNMENT_CLEARANCE = 0.18 * Y_SCALE
    BOTTOM_MARGIN = 0.55
    TOP_MARGIN = 0.30
    X_PAD_MIN = 60
    NON_K03455_X_MAX_DEFAULT = 10_000


# ================================
# NUMERIC LABEL OFFSETS
# ================================
class NumericOffsets:
    """Vertical offsets for numeric labels (non-K03455 only)."""

    # Format: gene_name: (down_offset, up_offset)
    GENE_OFFSET_MAP = {
        "tat 1": (0.15, 0.15),
        "tat 2": (0.15, 0.20),
        "rev 1": (0.20, 0.25),
        "rev 2": (0.20, 0.25),
        "gag": (-0.15, -0.15),
        "vif": (-0.15, -0.20),
        "5' ltr": (-0.15, -0.20),
        "pol": (-0.15, -0.20),
        "vpu": (0.15, 0.15),
        "env": (0.15, 0.15),
        "vpr": (0.15, -0.20),
        "nef": (-0.15, -0.20),
        "gag-pol": (-0.15, -0.20),
        "3' ltr": (0.25, -0.20),
    }
    DEFAULT_OFFSETS = (-0.15, 0.35)

    @classmethod
    def get_offsets(cls, gene: str) -> tuple[float, float]:
        """
        Get offsets for a gene (case-insensitive).

        Parameters
        ----------
        gene : str
            Gene name.
        """
        if not gene:
            return cls.DEFAULT_OFFSETS

        key = gene.strip().lower()
        if key in cls.GENE_OFFSET_MAP:
            return cls.GENE_OFFSET_MAP[key]
        if key.startswith("tat"):
            return cls.GENE_OFFSET_MAP.get("tat 1", cls.DEFAULT_OFFSETS)
        if key.startswith("rev"):
            return cls.GENE_OFFSET_MAP.get("rev 1", cls.DEFAULT_OFFSETS)
        return cls.DEFAULT_OFFSETS


# ================================
# K03455 SPECIAL REFERENCE CONFIG
# ================================
class K03455Config:
    """Configuration for K03455 special reference handling."""

    TARGET_REGIONS = [
        "5' LTR", "gag", "pol", "vif", "vpr", "vpu",
        "tat 1", "tat 2", "rev 1", "rev 2", "env", "nef", "3' LTR",
    ]

    Y_POSITIONS = {
        "5' LTR": 0.0,
        "gag": 0.4,
        "pol": 0.0,
        "vif": 0.4,
        "vpr": 0.8,
        "vpu": 0.6,
        "tat 1": 0.2,
        "tat 2": 0.4,
        "rev 1": 1.0,
        "rev 2": 1.0,
        "env": 0.0,
        "nef": 0.8,
        "3' LTR": 0.0,
    }

    # K03455-specific numeric label offsets
    K03455_NUMERIC_OFFSETS = {
        "tat 1": (-0.15, 0.15),
        "tat 2": (-0.15, 0.20),
        "rev 1": (0.20, 0.25),
        "rev 2": (0.20, 0.25),
        "gag": (0.30, 0.20),
        "vif": (-0.15, -0.20),
        "5' LTR": (-0.15, -0.20),
        "pol": (-0.15, -0.20),
        "vpu": (0.35, 0.15),
        "env": (-0.15, -0.20),
        "vpr": (0.20, -0.20),
        "nef": (-0.15, -0.20),
        "3' LTR": (0.25, -0.20),
    }
    DEFAULT_K03455_OFFSETS = (-0.15, 0.35)

    @classmethod
    def get_k03455_offsets(cls, gene: str) -> tuple[float, float]:
        """
        Get K03455-specific offsets for a gene.

        Parameters
        ----------
        gene : str
            Gene name.
        """
        return cls.K03455_NUMERIC_OFFSETS.get(gene, cls.DEFAULT_K03455_OFFSETS)

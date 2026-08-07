import re
from pathlib import Path

from pyhiv.align.align_with_reference import calculate_alignment_score


FEATURE_OUTPUT_PATHS = {
    "5' ltr": ("ltr", "5"),
    "5' ltr r": ("ltr", "5", "r"),
    "5' ltr u3": ("ltr", "5", "u3"),
    "5' ltr u5": ("ltr", "5", "u5"),
    "3' ltr": ("ltr", "3"),
    "3' ltr r": ("ltr", "3", "r"),
    "3' ltr u3": ("ltr", "3", "u3"),
    "3' ltr u5": ("ltr", "3", "u5"),
    "tar": ("ltr", "5", "tar"),
    "gag-pol": ("gag-pol",),
    "gag-pol fusion": ("gag-pol",),
    "gag_pol precursor": ("gag-pol",),
    "p17 (matrix)": ("gag", "p17-matrix"),
    "p24 (capsid)": ("gag", "p24-capsid"),
    "p7 (nucleocapsid)": ("gag", "p7-nucleocapsid"),
    "p6": ("gag", "p6"),
    "pol cds": ("pol",),
    "protease": ("pol", "protease"),
    "p51 (rt)": ("pol", "p51-rt"),
    "p15 (rnase h)": ("pol", "p15-rnase-h"),
    "p31 (integrase)": ("pol", "p31-integrase"),
    "vif cds": ("vif",),
    "vpr cds": ("vpr",),
    "tat cds (plus intron)": ("tat",),
    "tat exon 1": ("tat", "exon-1"),
    "tat exon 2": ("tat", "exon-2"),
    "rev cds (plus intron)": ("rev",),
    "rev exon 1": ("rev", "exon-1"),
    "rev exon 2": ("rev", "exon-2"),
    "vpu cds": ("vpu",),
    "vpu (defective)": ("vpu", "defective"),
    "env cds": ("env",),
    "gp160": ("env", "gp160"),
    "gp120": ("env", "gp120"),
    "gp41": ("env", "gp41"),
    "v1": ("env", "v1"),
    "v2": ("env", "v2"),
    "v3": ("env", "v3"),
    "v4": ("env", "v4"),
    "v5": ("env", "v5"),
    "rre": ("env", "rre"),
    "nef cds": ("nef",),
}


def safe_output_label(value: str) -> str:
    """Sanitize an arbitrary string into a stable, filesystem-friendly label."""
    label = str(value or "").strip()
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", label)
    label = label.strip("._-")
    return label or "unknown"


def sequence_output_label(file_name: str | None, sequence_name: str) -> str:
    """
    Build the output filename label for a sequence.

    Combines the source FASTA file stem with the sequence id so that
    same-named records from different input files don't collide. Used both
    when PyHIV writes alignment/split files and when PyHIVReporter looks
    them back up, so the two must stay in sync.
    """
    sequence_label = safe_output_label(sequence_name)
    if not file_name or file_name == "-":
        return sequence_label

    file_label = safe_output_label(Path(file_name).stem)
    if not file_label:
        return sequence_label
    return f"{file_label}_{sequence_label}"


def slugify_feature_name(feature_name: str) -> str:
    """Convert a feature name into a stable filesystem-friendly label."""
    slug = feature_name.strip().lower().replace("_", "-")
    slug = slug.replace("'", "")
    slug = re.sub(r"[^a-z0-9]+", "-", slug)
    return slug.strip("-") or "unknown"


def feature_output_location(feature_name: str, hierarchical: bool) -> tuple[Path, str]:
    """
    Return the output directory and filename suffix for a feature.

    HXB2-only splitting keeps the original feature names so fine-grained regions
    such as v3 remain directly visible. Subtyping output uses a stable hierarchy
    to group equivalent features under their parent genes.
    """
    if not hierarchical:
        return Path(feature_name), feature_name

    key = feature_name.strip().lower()
    parts = FEATURE_OUTPUT_PATHS.get(key)
    if parts is None:
        parts = (slugify_feature_name(feature_name),)

    return Path(*parts), parts[-1]


def map_ref_coords_to_alignment(ref_aligned: str) -> dict:
    """
    Build a mapping from reference coordinates without gaps (GenBank) to alignment columns with gaps.

    Parameters
    ----------
    ref_aligned : str
        The aligned reference sequence (may contain '-' characters representing gaps).

    Returns
    -------
    dict
        A dictionary mapping 1-based reference positions (without gaps) to 0-based alignment positions (with gaps).
    """
    mapping = {}
    ref_pos = 0
    for aln_pos, base in enumerate(ref_aligned):
        if base != "-":
            ref_pos += 1
            mapping[ref_pos] = aln_pos
    return mapping


def get_gene_region(test_aligned: str, ref_aligned: str, aligned_gene_ranges: dict) -> list:
    """
    Identify the gene region(s) with the highest alignment score.

    Parameters
    ----------
    test_aligned : str
        The aligned test sequence (with gaps).
    ref_aligned : str
        The aligned reference sequence (with gaps).
    aligned_gene_ranges : dict
        Dictionary mapping gene names to (start, end) positions in the alignment coordinates (0-based).

    Returns
    -------
    list
        A list of gene names corresponding to the region(s) with the highest alignment score.
        If multiple genes share the same maximum score, all of them are returned.
    """
    if not aligned_gene_ranges:
        return []

    gene_scores = {
        gene: calculate_alignment_score(
            test_aligned[start:end+1], ref_aligned[start:end+1]
        )
        for gene, (start, end) in aligned_gene_ranges.items()
    }

    max_score = max(gene_scores.values(), default=None)
    return [gene for gene, score in gene_scores.items() if score == max_score] if max_score is not None else []

def get_present_gene_regions(test_aligned: str, aligned_gene_ranges: dict) -> list:
    """
    Identify gene regions that contain at least one base (non-gap) in the aligned test sequence.

    Parameters
    ----------
    test_aligned : str
        The aligned test sequence (with gaps).
    aligned_gene_ranges : dict
        Dictionary mapping gene names to (start, end) positions in the alignment coordinates (0-based).

    Returns
    -------
    list
        A list of gene names where the test sequence contains non-gap characters within the region.
    """
    return [
        gene for gene, (start, end) in aligned_gene_ranges.items()
        if any(base != '-' for base in test_aligned[start:end+1])
    ]

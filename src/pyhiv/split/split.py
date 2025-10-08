from pyhiv.align.align_with_reference import calculate_alignment_score

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
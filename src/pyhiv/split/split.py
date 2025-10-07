from pyhiv.align.align_with_reference import calculate_alignment_score


def get_gene_region(test_aligned: str, ref_aligned: str, gene_ranges: dict) -> list:
    """
    Identify the gene region with the highest alignment score.
    If multiple regions have the same highest score, return all of them.

    Parameters
    ----------
    test_aligned : str
        The aligned test sequence.
    ref_aligned : str
        The aligned reference sequence.
    gene_ranges : dict
        Dictionary mapping gene names to their (start, end) positions.

    Returns
    -------
    list
        List of gene regions with the highest alignment score.
    """
    if not gene_ranges:
        return []
    gene_scores = {
        gene: calculate_alignment_score(test_aligned[start - 1:end], ref_aligned[start - 1:end])
        for gene, (start, end) in gene_ranges.items()
    }

    max_score = max(gene_scores.values(), default=None)
    return [gene for gene, score in gene_scores.items() if score == max_score] if max_score is not None else []


def get_present_gene_regions(test_aligned: str, gene_ranges: dict) -> list:
    """
    Identify gene regions that contain actual base pair letters instead of gaps ('-').

    Parameters
    ----------
    test_aligned : str
        The aligned test sequence.
    gene_ranges : dict
        Dictionary mapping gene names to their (start, end) positions.

    Returns
    -------
    list
        List of gene regions that contain base pair letters.
    """
    return [
        gene for gene, (start, end) in gene_ranges.items()
        if any(base != '-' for base in test_aligned[start - 1:end])
    ]

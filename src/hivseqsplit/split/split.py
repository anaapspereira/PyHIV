from hivseqsplit.align.align_with_reference import calculate_alignment_score


# Check what is the gene region with most matches and return the gene region. If more than two gene regions have the same score, return all of them.
def get_gene_region(test_aligned, ref_aligned, gene_ranges):
    gene_scores = {}
    for gene, (start, end) in gene_ranges.items():
        gene_seq = test_aligned[start - 1:end]
        score = calculate_alignment_score(gene_seq, ref_aligned[start - 1:end])
        gene_scores[gene] = score

    max_score = max(gene_scores.values())
    gene_region = [gene for gene, score in gene_scores.items() if score == max_score]
    return gene_region

# Check what gene regions have base pair letters instead of '-' in the test sequence and return the gene regions.
def get_present_gene_regions(test_aligned, gene_ranges):
    present_gene_regions = []
    for gene, (start, end) in gene_ranges.items():
        gene_seq = test_aligned[start - 1:end]
        #check if there is any other value other than '-' in that region
        if any([base != '-' for base in gene_seq]):
            present_gene_regions.append(gene)
    #remove duplicates
    present_gene_regions = list(set(present_gene_regions))
    return present_gene_regions

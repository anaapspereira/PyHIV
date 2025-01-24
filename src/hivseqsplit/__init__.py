__version__ = "0.0.1"

import ast
import os

import pandas as pd

from hivseqsplit.align import align_with_references
from hivseqsplit.loading import read_input_fastas, REFERENCE_GENOMES_DIR, REFERENCE_GENOMES_FASTAS_DIR, \
    HXB2_GENOME_FASTA_DIR
from hivseqsplit.split import get_gene_region, get_present_gene_regions


def HIMAPS(fastas_dir: str, subtyping: bool = True, splitting: bool = True, output_dir: str = None):
    user_fastas = read_input_fastas(fastas_dir)
    reference_sequences = pd.read_csv(os.path.join(REFERENCE_GENOMES_DIR, 'sequences_with_locations.tsv'), sep='\t')
    final_table = pd.DataFrame(
        columns=['Sequence', 'Reference', 'Subtype', 'Most Matching Gene Region', 'Present Gene Regions'])
    for fasta in user_fastas:
        if subtyping:
            best_alignment = align_with_references(fasta, references_dir=REFERENCE_GENOMES_FASTAS_DIR)
        else:
            best_alignment = align_with_references(fasta, references_dir=HXB2_GENOME_FASTA_DIR)
        sequence_name = fasta.id
        if best_alignment is not None:
            test_aligned, ref_aligned, ref_file = best_alignment
            # get gene ranges from reference file based on ref_file name
            gene_ranges = ast.literal_eval(
                reference_sequences[reference_sequences['accession'] == ref_file.split('-')[0]][
                    'features'].values[0])
            if not os.path.exists('data/final_results/'):
                os.makedirs('data/final_results/', exist_ok=True)

            # save a fasta file with the best alignment
            final_alignment_file = f'data/final_results/best_alignment_{sequence_name}.fasta'
            with open(final_alignment_file, 'w') as output_file:
                output_file.write(f">Reference {ref_file.split('.')[0]}\n{ref_aligned}\n")
                output_file.write(f'>{sequence_name}\n{test_aligned}\n')

        if splitting:
            # get gene region with most matches
            region = get_gene_region(test_aligned, ref_aligned, gene_ranges)
            # get gene regions with base pair letters
            present_regions = get_present_gene_regions(test_aligned, gene_ranges)

            # save gene regions fasta in each present regions specific gene region folder
            for gene in present_regions:
                if not os.path.exists(f'data/final_results/{gene}/'):
                    os.makedirs(f'data/final_results/{gene}/', exist_ok=True)
                with open(f'data/final_results/{gene}/{sequence_name}_{gene}.fasta', 'w') as output_file:
                    output_file.write(
                        f'>{sequence_name}\n{test_aligned[gene_ranges[gene][0] - 1:gene_ranges[gene][1]]}\n')
                    # save the results in a final global table
            row = pd.Series([sequence_name, ref_file.split('.')[0], ref_file.split('-')[1].split('.')[0],
                             str(region).strip("[]"), str(present_regions).strip("[]")], index=final_table.columns)
            final_table.loc[len(final_table)] = row
        else:
            row = pd.Series([sequence_name, ref_file.split('.')[0], ref_file.split('-')[1].split('.')[0],
                             "-", "-"], index=final_table.columns)
            final_table.loc[len(final_table)] = row
    if not splitting:
        final_table.drop(columns=['Most Matching Gene Region', 'Present Gene Regions'], inplace=True)
    final_table.to_csv('data/final_results/final_table.tsv', sep='\t', index=False)

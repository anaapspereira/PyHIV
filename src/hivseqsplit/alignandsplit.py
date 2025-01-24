import os
import pandas as pd
from Bio import SeqIO
import subprocess
from Bio import AlignIO
import ast
from Bio.SeqRecord import SeqRecord


def main():
    #check if sequence folder does not exist to convert the data
    reference_sequences = pd.read_csv('data/sequenceswithlocations.tsv', sep='\t')
    try:
        if not os.path.exists('data/reference_sequences/') or not os.listdir('data/reference_sequences/'):
            print('Converting data')
            convert_to_fasta_files(reference_sequences)
            print('Data converted and saved in data/reference_sequences/')
    except:
        print('Error converting data. Check sequencesewithlocations.tsv')
        return None

    # Align all sequences from user_sequences folder with the reference sequences
    try:
        final_table = pd.DataFrame(columns=['Sequence', 'Reference', 'Subtype', 'Most Matching Gene Region', 'Present Gene Regions'])

        for file in os.listdir('data/user_sequences/'):
            test_sequence = f'data/user_sequences/{file}'
            sequence_name = file.split('.')[0]
            best_alignment = align_with_references(test_sequence)
            if best_alignment is not None:
                test_aligned, ref_aligned, ref_file = best_alignment
                print(f"Best alignment for {file} found with {ref_file}")
                # get gene ranges from reference file based on ref_file name
                gene_ranges = ast.literal_eval(
                    reference_sequences[reference_sequences['accession'] == ref_file.split('-')[0]][
                        'features'].values[0])
                if not os.path.exists('data/final_results/'):
                    os.mkdir('data/final_results')

                # save a fasta file with the best alignment
                final_alignment_file = f'data/final_results/best_alignment_{sequence_name}.fasta'
                with open(final_alignment_file, 'w') as output_file:
                    output_file.write(f">Reference {ref_file.split('.')[0]}\n{ref_aligned}\n")
                    output_file.write(f'>{sequence_name}\n{test_aligned}\n')

                # final_gene_file = f'data/final_results/gene_sequences_{sequence_name}.tsv'
                # with open(final_gene_file, 'w') as output_file:
                #     output_file.write('Gene\tSequence\n')
                #     for gene, (start, end) in gene_ranges.items():
                #         gene_seq = test_aligned[start - 1:end]
                #         output_file.write(f'{gene}\t{gene_seq}\n')

                # get gene region with most matches
                region = get_gene_region(test_aligned, ref_aligned, gene_ranges)
                # get gene regions with base pair letters
                present_regions = get_present_gene_regions(test_aligned, gene_ranges)

                # save gene regions fasta in each present regions specific gene region folder
                for gene in present_regions:
                    if not os.path.exists(f'data/final_results/{gene}/'):
                        os.mkdir(f'data/final_results/{gene}/')
                    with open(f'data/final_results/{gene}/{sequence_name}_{gene}.fasta', 'w') as output_file:
                        output_file.write(f'>{sequence_name}\n{test_aligned[gene_ranges[gene][0] - 1:gene_ranges[gene][1]]}\n')
                # save the results in a final global table
                row = pd.Series([file, ref_file.split('.')[0], ref_file.split('-')[1].split('.')[0],
                                 str(region).strip("[]"), str(present_regions).strip("[]")], index=final_table.columns)
                final_table = final_table.append(row, ignore_index=True)
            else:
                print(f"No best alignment found for file {file}.")

        final_table.to_csv('data/final_results/final_table.tsv', sep='\t', index=False)

    except:
        print('Error aligning sequences. Check user_sequences folder')
        return None


if __name__ == "__main__":
    main()
import os
import pandas as pd
from Bio import Entrez, SeqIO
import subprocess
from Bio import AlignIO
import ast
from Bio.SeqRecord import SeqRecord


# Convert the sequences from sequenceswithlocations.tsv to individual fasta files
def convert_to_fasta_files(input_file):
    try:
        if not os.path.exists('data/reference_sequences/'):
            os.mkdir("data/reference_sequences")
        sequences = input_file
        for index, row in sequences.iterrows():
            sequence = row['sequence']
            subtype = row['subtype']
            accession = row['accession']
            with open(f'data/reference_sequences/{subtype}_{accession}.fasta', 'w') as output_file:
                output_file.write(f'>{subtype}_{accession}\n')
                output_file.write(f'{sequence}\n')
    except Exception as e:
        print(f"Error converting data: {e}")

# Align the test sequence with the reference sequence using MAFFT
def mafft_align(test_seq, ref_seq):
    temp_input_file = "temp_input.fasta"
    temp_output_file = "temp_output.fasta"
    try:
        with open(temp_input_file, 'w') as temp_file:
            SeqIO.write([SeqRecord(ref_seq, id="ref_seq"), SeqRecord(test_seq, id="test_seq")], temp_file, "fasta")

        mafft_cmd = ["mafft", "--localpair", "--maxiterate", "1000", temp_input_file]
        with open(temp_output_file, 'w') as output_file:
            subprocess.run(mafft_cmd, stdout=output_file, stderr=subprocess.PIPE, check=True)

        alignment = AlignIO.read(temp_output_file, "fasta")
        ref_aligned = None
        test_aligned = None
        for record in alignment:
            if record.id == "ref_seq":
                ref_aligned = str(record.seq)
            elif record.id == "test_seq":
                test_aligned = str(record.seq)

        if ref_aligned is None or test_aligned is None:
            print("Error: Reference sequence or test sequence not found in the alignment.")
            return None, None

        return test_aligned, ref_aligned
    except subprocess.CalledProcessError as e:
        print(f"Error in MAFFT alignment: {e}")
        return None, None
    finally:
        if os.path.exists(temp_input_file):
            os.remove(temp_input_file)
        if os.path.exists(temp_output_file):
            os.remove(temp_output_file)

# Calculate the alignment score between two sequences
def calculate_alignment_score(seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1
    return score


#Align testing sequence with sequences from references folder to find best alignment with MAFFT
def align_with_references(test_sequence, references_dir='data/reference_sequences/'):
    # Read the test sequence
    try:
        test_seq = SeqIO.read(test_sequence, "fasta")
    # Initialize the best alignment score to a high value
        best_score = -1
        best_alignment = None

        # Iterate over each reference sequence
        for ref_file in os.listdir(references_dir):
            ref_seq = SeqIO.read(f'{references_dir}{ref_file}', "fasta")
            # Align the test sequence with the reference sequence
            test_aligned, ref_aligned = mafft_align(test_seq.seq, ref_seq.seq)
            # Calculate the alignment score
            score = calculate_alignment_score(test_aligned, ref_aligned)
            # Update the best alignment if the current score is better
            if score > best_score:
                best_score = score
                best_alignment = (test_aligned, ref_aligned, ref_file)
        return best_alignment
    except Exception as e:
        print(f"Error reading the test sequence file: {e}")
        return None

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
                    reference_sequences[reference_sequences['accession'] == ref_file.split('_')[1].split('.')[0]][
                        'features'].values[0])
                if not os.path.exists('data/final_results/'):
                    os.mkdir('data/final_results')

                # save a fasta file with the best alignment
                final_alignment_file = f'data/final_results/best_alignment_{sequence_name}.fasta'
                with open(final_alignment_file, 'w') as output_file:
                    output_file.write(f">Reference {ref_file.split('.')[0]}\n{ref_aligned}\n")
                    output_file.write(f'>{sequence_name}\n{test_aligned}\n')
                final_gene_file = f'data/final_results/gene_sequences_{sequence_name}.tsv'
                with open(final_gene_file, 'w') as output_file:
                    output_file.write('Gene\tSequence\n')
                    for gene, (start, end) in gene_ranges.items():
                        gene_seq = test_aligned[start - 1:end]
                        output_file.write(f'{gene}\t{gene_seq}\n')

                # get gene region with most matches
                region = get_gene_region(test_aligned, ref_aligned, gene_ranges)
                # get gene regions with base pair letters
                present_regions = get_present_gene_regions(test_aligned, gene_ranges)

                row = pd.Series([file, ref_file.split('.')[0], ref_file.split('_')[0],
                                 str(region), str(present_regions)], index=final_table.columns)
                final_table = final_table.append(row, ignore_index=True)
            else:
                print(f"No best alignment found for file {file}.")

        final_table.to_csv('data/final_results/final_table.tsv', sep='\t', index=False)

    except:
        print('Error aligning sequences. Check user_sequences folder')
        return None




if __name__ == "__main__":
    main()
# Convert the sequences from sequenceswithlocations.tsv to individual fasta files
import os


def convert_to_fasta_files(input_file):
    try:
        if not os.path.exists('reference_genomes/reference_fastas/'):
            os.mkdir("reference_genomes/reference_fastas/")
        sequences = input_file
        for index, row in sequences.iterrows():
            sequence = row['sequence']
            subtype = row['subtype']
            accession = row['accession']
            with open(f'reference_genomes/reference_fastas/{accession}-{subtype}.fasta', 'w') as output_file:
                output_file.write(f'>{accession}-{subtype}\n')
                output_file.write(f'{sequence}\n')
    except Exception as e:
        print(f"Error converting data: {e}")

if __name__ == '__main__':
    import pandas as pd
    file = pd.read_csv('reference_genomes/sequences_with_locations.tsv', sep='\t')
    convert_to_fasta_files(file)
import os
from Bio import SeqIO

#Open hiv-db.fasta file and separate the sequences in the file into individual fasta files
def separate_sequences(input_file):
    try:
        if not os.path.exists('data/user_sequences/'):
            os.mkdir("data/user_sequences")
        sequences = SeqIO.parse(input_file, "fasta")
        for sequence in sequences:
            accession = sequence.id.split(".")[-1]
            with open(f'data/user_sequences/{accession}.fasta', 'w') as output_file:
                SeqIO.write(sequence, output_file, "fasta")
    except Exception as e:
        print(f"Error separating sequences: {e}")

def main():
    input_file = "data/hiv-db.fasta"
    separate_sequences(input_file)

if __name__ == "__main__":
    main()
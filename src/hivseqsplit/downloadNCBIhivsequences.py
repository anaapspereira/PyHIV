from Bio import Entrez, SeqIO

#download all complete HIV genome sequences from NCBI
def download_hiv_sequences(search_term, output_folder, output_file):
    # Set the email address
    Entrez.email = "youremail@example.com"
    # Define the search term
    search_term = search_term
    # Search the database
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10)
    record = Entrez.read(handle)
    # Get the list of IDs
    id_list = record["IdList"]
    # Fetch the records
    records = []
    # save the records in individual fasta files
    for i in id_list:
        handle = Entrez.efetch(db="nucleotide", id=i, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        SeqIO.write(record, f"{output_folder}{i}.fasta", "fasta")
        records.append(record)


def main():
    # Define the search term
    search_term = "HIV-1[Organism] AND complete genome[Title]"
    # Define the output folder
    output_folder = "data/user_sequences/"
    # Define the output file
    output_file = "data/user_sequences/hiv_sequences.fasta"
    # Download the sequences
    download_hiv_sequences(search_term, output_folder, output_file)

if __name__ == "__main__":
    main()


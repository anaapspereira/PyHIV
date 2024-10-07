import streamlit as st
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

#Get accession numbers from the Los Alamos reference genome file
def get_accession_numbers():
    try:
        references = pd.read_csv('data/HIV1_REF_2023_genome_DNA.csv', header=None)
        references[['subtype', 'country', 'year', 'isolate']] = references[0].str.split('.', expand=True).iloc[:, 1:5]
        references['accession'] = references[0].str.strip().str.split('.').str[-1]
        references['align fasta'] = references[1]
        references.drop(columns=[0,1], inplace=True)
        return references
    except:
        print('Error reading data/HIV1_REF_2023_genome_DNA.csv')
        return None

#Download the sequence features for each accession number from NCBI
def fetch_features(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        locations = get_gene_locations(record.features)
        return locations
    except Exception as e:
        print(f"Error fetching {accession_id}: {e}")
        return None

#download the sequence fasta for each accession number from NCBI
def fetch_sequence(accession_id):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record.seq
    except Exception as e:
        print(f"Error fetching {accession_id}: {e}")
        return None

#Get gene locations for each subtype from feature data
def get_gene_locations(features):
    gene_ranges = {}
    for feature in features:
        if feature.type == "gene":
            gene_name = feature.qualifiers.get("gene", ["locations"])[0]
            start = feature.location.start.position
            end = feature.location.end.position
            gene_ranges[gene_name] = (start, end)
    return gene_ranges

# Check what is the gene region with most matches and return the gene region
def get_gene_region(test_seq, ref_seq, gene_ranges):
    best_score = -1
    best_gene = None
    for gene, (start, end) in gene_ranges.items():
        score = calculate_alignment_score(test_seq[start:end], ref_seq[start:end])
        if score > best_score:
            best_score = score
            best_gene = gene
    return best_gene


def download_data():
    # Get accession numbers from the Los Alamos reference genome file
    references = get_accession_numbers()
    # Download the sequence data (features) for each accession number from NCBI
    Entrez.email = "your_email@example.com"

    references['features'] = references['accession'].apply(lambda x: fetch_features(x))
    references['sequence'] = references['accession'].apply(lambda x: fetch_sequence(x))

    #remove rows with missing features or sequences or empty dictionaries in features or sequences

    references = references[references['features'].notnull() & references['features'].apply(lambda x: bool(x))]
    references = references[references['sequence'].notnull() & references['sequence'].apply(lambda x: bool(x))]

    return references

# Define the directory to save uploaded files
uploadfolder = 'data/user_sequences'

# Create the upload folder if it does not exist
if not os.path.exists(uploadfolder):
    os.makedirs(uploadfolder)

# Function to save the uploaded FASTA files
def save_fasta_files(uploaded_files):
    #clean folder first
    for file in os.listdir(uploadfolder):
        os.remove(os.path.join(uploadfolder, file))
    saved_files = []
    for uploaded_file in uploaded_files:
        with open(os.path.join(uploadfolder, uploaded_file.name), "wb") as f:
            f.write(uploaded_file.getbuffer())
        saved_files.append(os.path.join(uploadfolder, uploaded_file.name))
    return saved_files

st.title("HIV-1 FASTA Sequence Alignment and Splitting")

# File uploader allows multiple files
uploaded_files = st.file_uploader("Upload FASTA files", type=["fasta"], accept_multiple_files=True)

if uploaded_files:
    st.write(f"Uploaded {len(uploaded_files)} file(s)")
    saved_files = save_fasta_files(uploaded_files)

    # Display button to start alignment
    if st.button("Start Alignment"):
        st.write("Alignment started...")


        #check if sequencewithlocations.tsv does not exist to download the data

        if not os.path.exists('data/sequenceswithlocations.tsv'):
            st.write('Downloading reference sequences')
            data = download_data()
            data.to_csv('data/sequenceswithlocations.tsv', sep='\t', index=False)
            st.write('Data downloaded, processing sequences...')
        else:
            data = pd.read_csv('data/sequenceswithlocations.tsv', sep='\t')

        reference_sequences = data
        try:
            if not os.path.exists('data/reference_sequences/') or not os.listdir('data/reference_sequences/'):
                st.write('Converting sequences')
                convert_to_fasta_files(reference_sequences)
                st.write('Sequences converted and saved, aligning sequences...')
        except:
            st.write('Error converting data. Check sequencesewithlocations.tsv')

        try:
            final_table = pd.DataFrame(columns=['Sequence', 'Reference', 'Subtype', 'Region'])
            # Align all sequences from user_sequences folder with the reference sequences

            for i, file in enumerate(os.listdir('data/user_sequences/')):
                test_sequence = f'data/user_sequences/{file}'
                sequence_name = file.split('.')[0]
                #final_table['Sequence'] = sequence_name
                best_alignment = align_with_references(test_sequence)
                if best_alignment is not None:
                    test_aligned, ref_aligned, ref_file = best_alignment
                    st.write(f"Best alignment for {file} found with {ref_file}")
                    # get gene ranges from reference file based on ref_file name
                    gene_ranges = ast.literal_eval(
                        reference_sequences[reference_sequences['accession'] == ref_file.split('_')[1].split('.')[0]][
                            'features'].values[0])
                    #row_data = {
                    #    'Sequence': file,
                    #    'Reference': ref_file.split('.')[0],
                    #    'Subtype': ref_file.split('_')[0],
                    #    'Region': get_gene_region(test_aligned, ref_aligned, gene_ranges)
                    #}

                    #final_table = final_table.append(row_data, ignore_index=True)
                    row = pd.Series([file, ref_file.split('.')[0], ref_file.split('_')[0], get_gene_region(test_aligned, ref_aligned, gene_ranges)], index=final_table.columns)
                    final_table = final_table.append(row, ignore_index=True)

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
                else:
                    st.write(f"No best alignment found for file {file}.")
            final_table.to_csv('data/final_results/final_table.tsv', sep='\t', index=False)
        except:
            st.write('Error aligning sequences. Check user_sequences folder')

        st.write("Alignment completed.")
        st.session_state['alignment_completed'] = True

    # Show Results button logic with session state management
    if 'alignment_completed' not in st.session_state:
        st.session_state['alignment_completed'] = False

    # Button to show results, and update session state when clicked
    if st.session_state['alignment_completed']:
        if st.button("Show Results"):
            try:
                # Load the final results file
                final_table = pd.read_csv('data/final_results/final_table.tsv', sep='\t')
                st.write("Results:")
                st.write(final_table)
            except Exception as e:
                st.write(f"Error loading results: {e}")



    #Clean up the screen

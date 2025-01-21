import pandas as pd
from Bio import Entrez, SeqIO
import os

#Get accession numbers from the Los Alamos reference genome file
def get_accession_numbers():
    try:
        references = pd.read_csv('data/HIV1_REF_2023_genome_DNA.csv', header=None)
        references[['subtype', 'country', 'year', 'isolate']] = references[0].str.split('.', expand=True).iloc[:, 1:5]
        references['accession'] = references[0].str.strip().str.split('.').str[-1]
        references['aligned fasta'] = references[1]
        references.drop(columns=[0,1], inplace=True)
        return references
    except:
        print('Error reading data/HIV1_REF_2023_genome_DNA.csv')
        return None

#get gene locations for each subtype sequence reference from NCBI
def fetch_features(accession_id):
    try:
        if accession_id =="K03455":
            locations = get_hxb2_locations()
            return (locations)
        else:
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            locations = get_gene_locations(record.features)
            return locations
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

def get_hxb2_locations():
    #get hxb2 location from HXB2_coordinates.csv
    hxb2_locations = pd.read_csv('data/HXB2_coordinates.csv', sep='\t')
    #set start location to coordinate value before -
    hxb2_locations[['start', 'end']] = hxb2_locations['HXB2 coordinates'].str.extract(r'(\d+)\s*-\s*(\d+)').astype(int)
    hxb2_locations['gene'] = hxb2_locations['Fragment'].str.replace('"', '\'')
    hxb2_locations = hxb2_locations[['gene','start','end']]
    locations = hxb2_locations.set_index('gene')[['start', 'end']].apply(tuple, axis=1).to_dict()
    return locations

def download_data():
    # Get accession numbers from the Los Alamos reference genome file
    references = get_accession_numbers()
    # Download the sequence data (features) for each accession number from NCBI
    Entrez.email = "your_email@example.com"

    references['features'] = references['accession'].apply(lambda x: fetch_features(x))
    references['sequence'] = references['aligned fasta'].str.replace('-','')

    #remove rows with missing features or empty
    references = references[references['features'].notnull() & references['features'].apply(lambda x: bool(x))]

    #remove rows with subtype equal to 'N','P', 'CPZ', 'GOR' # or 'O'
    references = references[~references['subtype'].isin(['N','P', 'O', 'CPZ', 'GOR'])]

    return references

def main():

    #check if sequencewithlocations.tsv does not exist to download the data

    if not os.path.exists('data/sequenceswithlocations.tsv'):
        print('Downloading data')
        data = download_data()
        print('Data downloaded, saving to data/sequenceswithlocations.tsv')
        data.to_csv('data/sequenceswithlocations.tsv', sep='\t', index=False)
    else:
        print('Data already downloaded')


if __name__ == "__main__":
    main()
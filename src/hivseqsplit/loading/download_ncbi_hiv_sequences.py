import logging
import os
from typing import List

import yaml
from Bio import Entrez, SeqIO
from pydantic import BaseModel, EmailStr, Field, field_validator

DEFAULT_NCBI_CONFIG = os.path.join(os.path.dirname(__file__), "data/default_ncbi_config.yml")
DEFAULT_SEQUENCES_PATH = os.path.join(os.path.dirname(__file__), "data/ncbi_hiv_sequences/")
DEFAULT_SEARCH_TERM = "HIV-1[Organism] AND complete genome[Title]"


class ConfigNCBI(BaseModel):
    """
    Configuration model for downloading sequences from NCBI.

    This class defines and validates the configuration required for querying
    the NCBI Entrez API to fetch sequences. It includes the email for API access,
    the search term to query the NCBI database, and the output folder where the
    fetched sequences will be saved.

    Attributes
    ----------
    email : EmailStr
        The email address to use with the Entrez API for identification.
    search_term : str
        The search query to be used for the NCBI Entrez search.
    output_folder : str
        The directory path where the resulting sequences will be saved.
    """
    EMAIL: str = Field("youremail@example.com", description="Email address for NCBI queries")
    SEARCH_TERM: str = Field("HIV-1[Organism] AND complete genome[Title]",
                             description="Search term for the NCBI database")
    OUTPUT_FOLDER: str = Field("ncbi_hiv_sequences/", description="Path to the folder where outputs will be saved")

    @classmethod
    @field_validator("output_folder", mode="before")
    def ensure_trailing_slash(cls, value: str) -> str:
        """Ensure the output folder path ends with a slash."""
        return value if value.endswith("/") else f"{value}/"

    class Config:
        schema_extra = {
            "example": {
                "EMAIL": "youremail@example.com",
                "SEARCH_TERM": "HIV-1[Organism] AND complete genome[Title]",
                "OUTPUT_FOLDER": "ncbi_hiv_sequences/"
            }
        }

def load_config_from_yaml(file_path: str) -> ConfigNCBI:
    """
    Load configuration from a YAML file and validate it using ConfigNCBI.

    Parameters
    ----------
    file_path: str
        Path to the YAML configuration file.

    Returns
    -------
        ConfigNCBI:
            Validated configuration object.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Configuration file not found at {file_path}")

    with open(file_path, 'r') as file:
        config_data = yaml.safe_load(file)

    if not isinstance(config_data, dict):
        raise ValueError("Invalid configuration format. Expected a dictionary.")

    return ConfigNCBI(**config_data)


def download_hiv_sequences(config: ConfigNCBI) -> List[str]:
    """
    Downloads HIV sequences from the NCBI nucleotide database and saves them as FASTA files.

    Parameters
    ----------
    config: ConfigNCBI
        Configuration settings for the download.

    Returns
    -------
    list
        List of paths to the saved FASTA files.

    Raises
    ------
    RuntimeError
        If the search query fails.
    ValueError
        If no sequences are found for the given search term.
    """
    # Set the email address for Entrez
    Entrez.email = config.EMAIL
    logging.info("Entrez email set for NCBI queries.")

    # Ensure the output folder exists
    os.makedirs(config.OUTPUT_FOLDER, exist_ok=True)
    logging.info(f"Output folder verified: {config.OUTPUT_FOLDER}")

    # Search the database
    try:
        with Entrez.esearch(db="nucleotide", term=config.SEARCH_TERM, retmax=10) as handle:
            record = Entrez.read(handle)
        logging.info(f"Search completed for term: {config.SEARCH_TERM}")
    except Exception as e:
        logging.error(f"Error during Entrez search: {e}")
        raise RuntimeError("Failed to search NCBI database.") from e

    # Get the list of IDs
    id_list = record.get("IdList", [])
    if not id_list:
        logging.warning("No sequences found for the given search term.")
        raise ValueError("No sequences found for the given search term.")

    # Fetch and save records
    saved_files = []
    for seq_id in id_list:
        try:
            with Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text") as handle:
                record = SeqIO.read(handle, "genbank")

            output_file = os.path.join(config.OUTPUT_FOLDER, f"{seq_id}.fasta")
            SeqIO.write(record, output_file, "fasta")
            saved_files.append(output_file)

            logging.info(f"Sequence {seq_id} saved to {output_file}")
        except Exception as e:
            logging.error(f"Error fetching or saving sequence {seq_id}: {e}")

    logging.info(f"Total sequences saved: {len(saved_files)}")
    return saved_files

def ncbi_sequences_download(config_path: str = None) -> List[str]:
    """
    Download HIV sequences from NCBI using the provided configuration file.

    Parameters
    ----------
    config_path: str
        Path to the YAML configuration file.

    Returns
    -------
    list
        List of paths to the saved FASTA files.

    Raises
    ------
    FileNotFoundError
        If the configuration file is not found.
    """
    config = load_config_from_yaml(config_path or DEFAULT_NCBI_CONFIG)
    # Download the sequences only if search term is different from default, otherwise use the default
    if config.SEARCH_TERM == DEFAULT_SEARCH_TERM:
        logging.info("Using default search term for NCBI sequences.")
        # list all files in the default folder
        return [os.path.join(DEFAULT_SEQUENCES_PATH, file) for file in os.listdir(DEFAULT_SEQUENCES_PATH)]
    return download_hiv_sequences(config)

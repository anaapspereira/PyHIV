import os

from .read_fastas import read_input_fastas

REFERENCE_GENOMES_DIR = os.path.join(os.path.dirname(__file__), 'reference_genomes')
REFERENCE_GENOMES_FASTAS_DIR = os.path.join(REFERENCE_GENOMES_DIR, 'reference_fastas')

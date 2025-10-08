__version__ = "0.0.3"

import ast
from pathlib import Path
import pandas as pd

from pyhiv.align import align_with_references
from pyhiv.config import get_reference_paths, validate_reference_paths
from pyhiv.loading import read_input_fastas
from pyhiv.split import get_gene_region, get_present_gene_regions, map_ref_coords_to_alignment

FINAL_TABLE_COLUMNS = ['Sequence', 'Reference', 'Subtype', 'Most Matching Gene Region', 'Present Gene Regions']


def PyHIV(fastas_dir: str, subtyping: bool = True, splitting: bool = True,
          output_dir: str = None, n_jobs: int = None):
    """
    Main function to run the PyHIV pipeline.
    It aligns the user sequences with the reference sequences and saves the
    best alignment in a fasta file. If subtyping is True, it aligns the user
    sequences with the reference sequences from the HIV-1 subtyping tool.
    If splitting is True, it splits the user sequences into gene regions
    and saves them in specific folders. It also saves a final table with the results.
    """
    paths = get_reference_paths()
    validate_reference_paths(paths)

    fastas_dir = Path(fastas_dir)
    output_dir = Path(output_dir) if output_dir else Path('PyHIV_results')
    output_dir.mkdir(parents=True, exist_ok=True)

    user_fastas = read_input_fastas(fastas_dir)
    reference_sequences = pd.read_csv(paths["SEQUENCES_WITH_LOCATION"], sep='\t')

    final_table = pd.DataFrame(columns=FINAL_TABLE_COLUMNS)

    for fasta in user_fastas:
        reference_dir = (
            paths["REFERENCE_GENOMES_FASTAS_DIR"]
            if subtyping else
            paths["HXB2_GENOME_FASTA_DIR"]
        )
        best_alignment = align_with_references(fasta, references_dir=reference_dir, n_jobs=n_jobs)

        if best_alignment is None:
            continue

        sequence_name = fasta.id
        test_aligned, ref_aligned, ref_file = best_alignment



        # Extract reference information
        ref_file_parts = Path(ref_file).stem.split('-')
        accession = ref_file_parts[0]
        subtype = ref_file_parts[1] if len(ref_file_parts) > 1 else "Unknown"

        # Retrieve gene ranges
        gene_ranges = ast.literal_eval(
            reference_sequences.loc[
                reference_sequences['accession'] == accession, 'features'
            ].values[0]
        )

        # save a fasta file with the best alignment
        final_alignment_file = output_dir / f"best_alignment_{sequence_name}.fasta"
        with open(final_alignment_file, 'w') as output_file:
            output_file.write(
                f">Reference {Path(ref_file).stem}\n{ref_aligned}\n>{sequence_name}\n{test_aligned}\n"
            )

        if splitting:

            mapping = map_ref_coords_to_alignment(ref_aligned)

            aligned_gene_ranges = {
                gene: (mapping[start], mapping[end])
                for gene, (start, end) in gene_ranges.items()
                if start in mapping and end in mapping
            }

            # get gene region with most matches
            region = get_gene_region(test_aligned, ref_aligned, aligned_gene_ranges)
            # get gene regions with base pair letters
            present_regions = get_present_gene_regions(test_aligned, aligned_gene_ranges)

            # Save gene regions fasta in each region-specific folder
            for gene in present_regions:
                gene_path = output_dir / gene
                gene_path.mkdir(parents=True, exist_ok=True)
                gene_file = gene_path / f"{sequence_name}_{gene}.fasta"
                with open(gene_file, 'w') as output_file:
                    aln_start, aln_end = aligned_gene_ranges[gene]
                    seq_fragment = test_aligned[aln_start:aln_end+1]
                    output_file.write(f'>{sequence_name}\n{seq_fragment}\n')

            # Save the results in the final global table
            row_data = [sequence_name, accession, subtype, str(region).strip("[]"), str(present_regions).strip("[]")]
        else:
            row_data = [sequence_name, accession, subtype, "-", "-"]

        final_table = pd.concat(
            [final_table, pd.DataFrame([row_data], columns=final_table.columns)],
            ignore_index=True
        )
    if not splitting:
        final_table.drop(columns=['Most Matching Gene Region', 'Present Gene Regions'], inplace=True)

    final_table.to_csv(output_dir / 'final_table.tsv', sep='\t', index=False)

__version__ = "0.1.0"

import ast
from pathlib import Path
import pandas as pd

from pyhiv.align import align_with_references, reference_accession_from_name
from pyhiv.config import get_reference_paths, validate_reference_paths
from pyhiv.loading import read_input_fastas
from pyhiv.split import (
    feature_output_location,
    get_gene_region,
    get_present_gene_regions,
    map_ref_coords_to_alignment,
    slugify_feature_name,
)
from pyhiv.report import PyHIVReporter
import logging

from pyhiv.align import DEFAULT_ALIGNMENT_TOOL, DEFAULT_KMER_SIZE, DEFAULT_REFERENCE_TOP_K

FINAL_TABLE_COLUMNS = [
    'Sequence',
    'Reference',
    'Group',
    'Subtype',
    'Closest Subtypes',
    'Most Matching Gene Region',
    'Present Gene Regions',
]
MAX_HIV1_SEQUENCE_LENGTH = 12000
SEQUENCE_TOO_LONG_WARNING = "The submitted sequence is longer than the HIV-1 genome."
DEFAULT_REFERENCE_GROUPS = ("M",)
SUPPORTED_REFERENCE_GROUPS = ("M", "N", "O", "P")
DEFAULT_CLOSEST_SUBTYPES_COUNT = 3
NO_SUBTYPING_LABEL = "No subtyping performed."

def PyHIV(fastas_dir: str, subtyping: bool = True, splitting: bool = True,
          output_dir: str = None, n_jobs: int = None, reporting: bool = True,
          alignment_tool: str = DEFAULT_ALIGNMENT_TOOL,
          kmer_size: int = DEFAULT_KMER_SIZE,
          reference_top_k: int = DEFAULT_REFERENCE_TOP_K,
          reference_groups=None):
    """
    Main function to run the PyHIV pipeline.
    It aligns the user sequences with the reference sequences and saves the
    best alignment in a fasta file. If subtyping is True, it aligns the user
    sequences with the reference sequences from the HIV-1 subtyping tool.
    If splitting is True, it splits the user sequences into gene regions
    and saves them in specific folders. It also saves a final table with the results.
    If reporting is True, it generates a PDF report with visualizations.
    alignment_tool selects the alignment backend. Defaults to edlib-HW.
    kmer_size and reference_top_k control the k-mer candidate reference
    prefilter used before alignment.
    reference_groups controls which HIV-1 reference groups are eligible for
    subtyping. When splitting is enabled, omitted groups default to group M.
    When splitting is disabled, omitted groups allow all FASTA references.
    Use ("M", "N", "O", "P") to include all known groups.
    """
    paths = get_reference_paths()
    validate_reference_paths(paths)

    fastas_dir = Path(fastas_dir)
    output_dir = Path(output_dir) if output_dir else Path('PyHIV_results')
    output_dir.mkdir(parents=True, exist_ok=True)

    user_fastas = read_input_fastas(fastas_dir)
    reference_sequences = pd.read_csv(paths["SEQUENCES_WITH_LOCATION"], sep='\t')
    metadata_by_accession = build_reference_metadata(reference_sequences)
    reference_group_filter_requested = reference_groups is not None
    selected_reference_groups = normalize_reference_groups(reference_groups)
    allowed_reference_accessions = selected_reference_accessions(
        reference_sequences,
        selected_reference_groups,
    ) if subtyping and (splitting or reference_group_filter_requested) else None

    final_table = pd.DataFrame(columns=FINAL_TABLE_COLUMNS)

    for fasta in user_fastas:
        sequence_name = fasta.id
        sequence_length = len(str(fasta.seq).replace("-", ""))
        if sequence_length > MAX_HIV1_SEQUENCE_LENGTH:
            logging.warning("%s Skipping sequence '%s'.", SEQUENCE_TOO_LONG_WARNING, sequence_name)
            continue

        reference_dir = (
            paths["REFERENCE_GENOMES_FASTAS_DIR"]
            if subtyping else
            paths["HXB2_GENOME_FASTA_DIR"]
        )
        best_alignment = align_with_references(
            fasta,
            references_dir=reference_dir,
            n_jobs=n_jobs,
            alignment_tool=alignment_tool,
            kmer_size=kmer_size,
            reference_top_k=reference_top_k,
            allowed_reference_accessions=allowed_reference_accessions,
            include_alignment_scores=True,
        )

        if best_alignment is None:
            continue

        test_aligned, ref_aligned, ref_file, alignment_scores = best_alignment

        # Extract reference information
        ref_file_parts = Path(ref_file).stem.split('-')
        accession = ref_file_parts[0]
        if subtyping:
            reference_metadata = metadata_by_accession.get(accession, {})
            group = reference_metadata.get("group", "Unknown")
            subtype = reference_metadata.get(
                "subtype",
                ref_file_parts[1] if len(ref_file_parts) > 1 else "Unknown",
            )
            closest_subtypes = summarize_closest_subtypes(
                alignment_scores,
                metadata_by_accession,
                top_n=DEFAULT_CLOSEST_SUBTYPES_COUNT,
            )
        else:
            group = NO_SUBTYPING_LABEL
            subtype = NO_SUBTYPING_LABEL
            closest_subtypes = NO_SUBTYPING_LABEL

        # save a fasta file with the best alignment
        final_alignment_file = output_dir / f"best_alignment_{sequence_name}.fasta"
        with open(final_alignment_file, 'w') as output_file:
            output_file.write(
                f">Reference {Path(ref_file).stem}\n{ref_aligned}\n>{sequence_name}\n{test_aligned}\n"
            )

        if splitting:
            # Retrieve gene ranges
            gene_ranges = ast.literal_eval(
                reference_sequences.loc[
                    reference_sequences['accession'] == accession, 'features'
                ].values[0]
            )

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
            written_region_files: set[Path] = set()
            for gene in present_regions:
                relative_gene_path, file_suffix = feature_output_location(
                    gene,
                    hierarchical=subtyping,
                )
                gene_path = output_dir / relative_gene_path
                gene_path.mkdir(parents=True, exist_ok=True)
                gene_file = gene_path / f"{sequence_name}_{file_suffix}.fasta"
                if gene_file in written_region_files:
                    original_suffix = slugify_feature_name(gene)
                    gene_file = gene_path / f"{sequence_name}_{original_suffix}.fasta"
                written_region_files.add(gene_file)
                with open(gene_file, 'w') as output_file:
                    aln_start, aln_end = aligned_gene_ranges[gene]
                    seq_fragment = test_aligned[aln_start:aln_end+1]
                    output_file.write(f'>{sequence_name}\n{seq_fragment}\n')

            # Save the results in the final global table
            row_data = [
                sequence_name,
                accession,
                group,
                subtype,
                closest_subtypes,
                str(region).strip("[]"),
                str(present_regions).strip("[]"),
            ]
        else:
            row_data = [sequence_name, accession, group, subtype, closest_subtypes, "-", "-"]

        final_table = pd.concat(
            [final_table, pd.DataFrame([row_data], columns=final_table.columns)],
            ignore_index=True
        )
    if not splitting:
        final_table.drop(columns=['Most Matching Gene Region', 'Present Gene Regions'], inplace=True)

    final_table.to_csv(output_dir / 'final_table.tsv', sep='\t', index=False)
    
    # Generate PDF report if requested
    if reporting:
        try:
            reporter = PyHIVReporter(output_dir, subtyping=subtyping, splitting=splitting)
            sequences_with_locations_path = paths["SEQUENCES_WITH_LOCATION"]
            pdf_path = reporter.generate_report(
                final_table_path=output_dir / 'final_table.tsv',
                sequences_with_locations_path=sequences_with_locations_path,
                output_pdf_name="PyHIV_report_all_sequences.pdf"
            )
            logging.info(f"PDF report generated: {pdf_path}")
        except Exception as e: # pragma: no cover
            logging.exception("Error generating PDF report — continuing without it.")


def normalize_reference_groups(reference_groups) -> tuple[str, ...]:
    if reference_groups is None:
        reference_groups = DEFAULT_REFERENCE_GROUPS

    if isinstance(reference_groups, str):
        raw_groups = reference_groups.replace(";", ",").split(",")
    else:
        raw_groups = []
        for group in reference_groups:
            raw_groups.extend(str(group).replace(";", ",").split(","))

    normalized = tuple(dict.fromkeys(group.strip().upper() for group in raw_groups if group.strip()))
    if not normalized:
        raise ValueError("reference_groups must include at least one group.")

    if normalized == ("ALL",):
        return SUPPORTED_REFERENCE_GROUPS

    unsupported = [group for group in normalized if group not in SUPPORTED_REFERENCE_GROUPS]
    if unsupported:
        supported = ", ".join(SUPPORTED_REFERENCE_GROUPS)
        raise ValueError(
            f"Unsupported reference group(s): {', '.join(unsupported)}. "
            f"Supported groups: {supported}."
        )

    return normalized


def build_reference_metadata(reference_sequences: pd.DataFrame) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    has_group = "group" in reference_sequences.columns
    has_subtype = "subtype" in reference_sequences.columns

    for _, row in reference_sequences.iterrows():
        accession = str(row["accession"])
        metadata[accession] = {
            "group": str(row["group"]) if has_group else "Unknown",
            "subtype": str(row["subtype"]) if has_subtype else "Unknown",
        }

    return metadata


def summarize_closest_subtypes(
    alignment_scores: list[tuple[int, str]],
    metadata_by_accession: dict[str, dict[str, str]],
    top_n: int = DEFAULT_CLOSEST_SUBTYPES_COUNT,
) -> str:
    closest: list[str] = []
    seen: set[tuple[str, str]] = set()

    for score, reference_name in sorted(alignment_scores, key=lambda item: (-item[0], item[1])):
        accession = reference_accession_from_name(reference_name)
        metadata = metadata_by_accession.get(accession, {})
        group = metadata.get("group", "Unknown")
        subtype = metadata.get("subtype", "Unknown")
        key = (group, subtype)
        if key in seen:
            continue

        seen.add(key)
        closest.append(f"{group}:{subtype} (score={score})")
        if len(closest) == top_n:
            break

    return "; ".join(closest) if closest else "-"


def selected_reference_accessions(
    reference_sequences: pd.DataFrame,
    reference_groups: tuple[str, ...],
) -> set[str] | None:
    if "group" not in reference_sequences.columns:
        logging.warning(
            "sequences_with_locations.tsv has no 'group' column; using all listed reference accessions."
        )
        return set(reference_sequences["accession"].astype(str))

    selected = set(reference_sequences.loc[
        reference_sequences["group"].astype(str).str.upper().isin(reference_groups),
        "accession",
    ].astype(str))

    if not selected:
        raise ValueError(
            f"No reference accessions found for group(s): {', '.join(reference_groups)}."
        )

    return selected

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
    sequence_output_label,
    slugify_feature_name,
)
from pyhiv.report import PyHIVReporter
import logging

from pyhiv.align import DEFAULT_ALIGNMENT_TOOL, DEFAULT_KMER_SIZE, DEFAULT_REFERENCE_TOP_K

FINAL_TABLE_BASE_COLUMNS = [
    'File Name',
    'Sequence',
    'Reference',
    'Group',
    'Subtype',
    'Closest Subtypes',
]
FINAL_TABLE_SPLITTING_COLUMNS = [
    'Splitting Reference',
    'Most Matching Gene Region',
    'Present Gene Regions',
]
MAX_HIV1_SEQUENCE_LENGTH = 12000
SEQUENCE_TOO_LONG_WARNING = "The submitted sequence is longer than the HIV-1 genome."
DEFAULT_REFERENCE_GROUPS = ("M",)
SUPPORTED_REFERENCE_GROUPS = ("M", "N", "O", "P")
DEFAULT_CLOSEST_SUBTYPES_COUNT = 3
NO_SUBTYPING_LABEL = "No subtyping performed."
SPLITTING_MODE_NONE = "none"
SPLITTING_MODE_SUBTYPE = "subtype"
SPLITTING_MODE_HXB2 = "hxb2"
SPLITTING_ALIGNMENT_PREFIX = "splitting_alignment"
SPLITTING_MODE_ALIASES = {
    "0": SPLITTING_MODE_NONE,
    "false": SPLITTING_MODE_NONE,
    "no": SPLITTING_MODE_NONE,
    "none": SPLITTING_MODE_NONE,
    "off": SPLITTING_MODE_NONE,
    "1": SPLITTING_MODE_SUBTYPE,
    "true": SPLITTING_MODE_SUBTYPE,
    "yes": SPLITTING_MODE_SUBTYPE,
    "on": SPLITTING_MODE_SUBTYPE,
    "subtype": SPLITTING_MODE_SUBTYPE,
    "subtyping": SPLITTING_MODE_SUBTYPE,
    "reference": SPLITTING_MODE_HXB2,
    "ref": SPLITTING_MODE_HXB2,
    "hxb2": SPLITTING_MODE_HXB2,
}


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
    If splitting is enabled, it splits the user sequences into gene regions
    and saves them in specific folders. splitting accepts booleans for backwards
    compatibility, or "subtype", "hxb2"/"reference", and "none".
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
    splitting_mode = normalize_splitting_mode(splitting, subtyping=subtyping)
    should_split = splitting_mode != SPLITTING_MODE_NONE

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
        require_features=False,
    ) if subtyping and (should_split or reference_group_filter_requested) else None

    final_table_columns = final_table_columns_for_splitting(should_split)
    final_table = pd.DataFrame(columns=final_table_columns)

    for fasta in user_fastas:
        sequence_name = fasta.id
        file_name = fasta.annotations.get("source_file", "-")
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
        splitting_test_aligned = test_aligned
        splitting_ref_aligned = ref_aligned
        splitting_ref_file = ref_file

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
        output_label = sequence_output_label(file_name, sequence_name)
        final_alignment_file = output_dir / f"best_alignment_{output_label}.fasta"
        write_alignment_file(final_alignment_file, Path(ref_file).stem, ref_aligned, sequence_name, test_aligned)

        if should_split:
            split_against_hxb2 = (
                splitting_mode == SPLITTING_MODE_HXB2 and subtyping
            ) or (
                splitting_mode == SPLITTING_MODE_SUBTYPE
                and subtyping
                and not has_reference_features(reference_sequences, accession)
            )

            if split_against_hxb2:
                hxb2_alignment = align_with_references(
                    fasta,
                    references_dir=paths["HXB2_GENOME_FASTA_DIR"],
                    n_jobs=n_jobs,
                    alignment_tool=alignment_tool,
                    kmer_size=kmer_size,
                    reference_top_k=reference_top_k,
                )
                if hxb2_alignment is None:
                    row_data = [
                        file_name,
                        sequence_name,
                        accession,
                        group,
                        subtype,
                        closest_subtypes,
                        "-",
                        "-",
                        "-",
                    ]
                    final_table = pd.concat(
                        [final_table, pd.DataFrame([row_data], columns=final_table.columns)],
                        ignore_index=True
                    )
                    continue

                splitting_test_aligned, splitting_ref_aligned, splitting_ref_file = hxb2_alignment
                splitting_alignment_file = output_dir / f"{SPLITTING_ALIGNMENT_PREFIX}_{output_label}.fasta"
                write_alignment_file(
                    splitting_alignment_file,
                    Path(splitting_ref_file).stem,
                    splitting_ref_aligned,
                    sequence_name,
                    splitting_test_aligned,
                )

            splitting_accession = Path(splitting_ref_file).stem.split('-')[0]
            # Retrieve gene ranges
            gene_ranges = reference_features(reference_sequences, splitting_accession)

            mapping = map_ref_coords_to_alignment(splitting_ref_aligned)

            aligned_gene_ranges = {
                gene: (mapping[start], mapping[end])
                for gene, (start, end) in gene_ranges.items()
                if start in mapping and end in mapping
            }

            # get gene region with most matches
            region = get_gene_region(splitting_test_aligned, splitting_ref_aligned, aligned_gene_ranges)
            # get gene regions with base pair letters
            present_regions = get_present_gene_regions(splitting_test_aligned, aligned_gene_ranges)

            # Save gene regions fasta in each region-specific folder
            written_region_files: set[Path] = set()
            for gene in present_regions:
                relative_gene_path, file_suffix = feature_output_location(
                    gene,
                    hierarchical=splitting_mode == SPLITTING_MODE_SUBTYPE,
                )
                gene_path = output_dir / relative_gene_path
                gene_path.mkdir(parents=True, exist_ok=True)
                gene_file = gene_path / f"{output_label}_{file_suffix}.fasta"
                if gene_file in written_region_files:
                    original_suffix = slugify_feature_name(gene)
                    gene_file = gene_path / f"{output_label}_{original_suffix}.fasta"
                written_region_files.add(gene_file)
                with open(gene_file, 'w') as output_file:
                    aln_start, aln_end = aligned_gene_ranges[gene]
                    seq_fragment = splitting_test_aligned[aln_start:aln_end+1]
                    output_file.write(f'>{sequence_name}\n{seq_fragment}\n')

            # Save the results in the final global table
            row_data = [
                file_name,
                sequence_name,
                accession,
                group,
                subtype,
                closest_subtypes,
                splitting_accession,
                str(region).strip("[]"),
                str(present_regions).strip("[]"),
            ]
        else:
            row_data = [file_name, sequence_name, accession, group, subtype, closest_subtypes]

        final_table = pd.concat(
            [final_table, pd.DataFrame([row_data], columns=final_table.columns)],
            ignore_index=True
        )

    final_table.to_csv(output_dir / 'final_table.tsv', sep='\t', index=False)
    
    # Generate PDF report if requested
    if reporting:
        try:
            reporter = PyHIVReporter(output_dir, subtyping=subtyping, splitting=should_split)
            sequences_with_locations_path = paths["SEQUENCES_WITH_LOCATION"]
            pdf_path = reporter.generate_report(
                final_table_path=output_dir / 'final_table.tsv',
                sequences_with_locations_path=sequences_with_locations_path,
                output_pdf_name="PyHIV_report_all_sequences.pdf"
            )
            logging.info(f"PDF report generated: {pdf_path}")
        except Exception as e: # pragma: no cover
            logging.exception("Error generating PDF report — continuing without it.")


def final_table_columns_for_splitting(should_split: bool) -> list[str]:
    if should_split:
        return FINAL_TABLE_BASE_COLUMNS + FINAL_TABLE_SPLITTING_COLUMNS
    return FINAL_TABLE_BASE_COLUMNS.copy()


def normalize_splitting_mode(splitting, subtyping: bool = True) -> str:
    if isinstance(splitting, bool):
        if not splitting:
            return SPLITTING_MODE_NONE
        return SPLITTING_MODE_SUBTYPE if subtyping else SPLITTING_MODE_HXB2

    mode = str(splitting).strip().lower().replace("-", "_")
    mode = mode.replace("_", "")
    try:
        normalized = SPLITTING_MODE_ALIASES[mode]
    except KeyError as exc:
        supported = "subtype, hxb2/reference, none"
        raise ValueError(f"Unsupported splitting mode: {splitting!r}. Supported modes: {supported}.") from exc

    if normalized == SPLITTING_MODE_SUBTYPE and not subtyping:
        return SPLITTING_MODE_HXB2
    return normalized


def write_alignment_file(
    path: Path,
    reference_name: str,
    ref_aligned: str,
    sequence_name: str,
    test_aligned: str,
) -> None:
    with open(path, 'w') as output_file:
        output_file.write(
            f">Reference {reference_name}\n{ref_aligned}\n>{sequence_name}\n{test_aligned}\n"
        )


def reference_features(reference_sequences: pd.DataFrame, accession: str) -> dict:
    matches = reference_sequences.loc[
        reference_sequences['accession'].astype(str) == str(accession),
        'features'
    ]
    if matches.empty:
        raise ValueError(f"No annotated features found for splitting reference accession: {accession}.")
    return ast.literal_eval(matches.values[0])


def has_reference_features(reference_sequences: pd.DataFrame, accession: str) -> bool:
    if "features" not in reference_sequences.columns:
        return False

    matches = reference_sequences.loc[
        reference_sequences["accession"].astype(str) == str(accession),
        "features",
    ]
    if matches.empty:
        return False

    return reference_has_features(matches.values[0])


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
    require_features: bool = False,
) -> set[str] | None:
    selected_rows = reference_sequences

    if "group" not in reference_sequences.columns:
        logging.warning(
            "sequences_with_locations.tsv has no 'group' column; using all listed reference accessions."
        )
    else:
        selected_rows = selected_rows.loc[
            selected_rows["group"].astype(str).str.upper().isin(reference_groups)
        ]

    if require_features:
        if "features" not in reference_sequences.columns:
            raise ValueError(
                "sequences_with_locations.tsv must include a 'features' column when splitting is enabled."
            )
        selected_rows = selected_rows.loc[
            selected_rows["features"].apply(reference_has_features)
        ]

    selected = set(selected_rows["accession"].astype(str))

    if not selected:
        detail = " with annotated features" if require_features else ""
        raise ValueError(
            f"No reference accessions{detail} found for group(s): {', '.join(reference_groups)}."
        )

    return selected


def reference_has_features(features) -> bool:
    if features is None or (isinstance(features, float) and pd.isna(features)):
        return False

    if isinstance(features, dict):
        return bool(features)

    text = str(features).strip()
    if not text or text.lower() in {"none", "nan", "null"}:
        return False

    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return True

    if parsed is None:
        return False
    if isinstance(parsed, dict):
        return bool(parsed)
    return True

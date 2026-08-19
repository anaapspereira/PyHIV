"""
Main reporting class for PyHIV results.
"""

from pathlib import Path
import logging

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .pdf_generator import render_sequence_page
from pyhiv.report.constants import K03455Config
from pyhiv.report.utils import (
    read_alignment_fasta, parse_present_regions, parse_features,
    is_special_reference, normalize_features, normalize_present_regions,
    build_ref_to_alignment_map, project_features_to_alignment,
    build_alignment_path
)


class PyHIVReporter:
    """Main class for generating PyHIV PDF reports."""
    
    def __init__(self, output_dir: Path, subtyping: bool, splitting: bool, log_level: int | None = None):
        """
        Initialize the reporter with output directory and logger.

        Parameters
        ----------
        output_dir : Path
            Directory to save the generated PDF report.
        subtyping : bool
            Whether subtyping was performed.
        splitting : bool
            Whether splitting was performed.
        log_level : int, optional
            If given, overrides this logger's level explicitly. By default
            (None), no handler or level is set here, so this logger
            propagates to and inherits its effective level from the root
            logger — e.g. the level the CLI configures from -v/-q — instead
            of always emitting INFO-level messages regardless of those
            flags.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.subtyping = subtyping
        self.splitting = splitting

        self.logger = logging.getLogger(self.__class__.__name__)
        if log_level is not None:
            self.logger.setLevel(log_level)
    
    def generate_report(
        self,
        final_table_path: Path,
        sequences_with_locations_path: Path,
        output_pdf_name: str = "PyHIV_report_all_sequences.pdf"
    ) -> Path:
        """
        Generate PDF report from PyHIV results.

        Parameters
        ----------
        final_table_path : Path
            Path to final_table.tsv file.
        sequences_with_locations_path : Path
            Path to sequences_with_locations.tsv file.
        output_pdf_name : str, optional
            Name of the output PDF file, by default "PyHIV_report_all_sequences.pdf"

        Returns
        -------
        Path
            Path to the generated PDF report.
        """
        # Read input data
        self.logger.info(f"Reading {final_table_path} and {sequences_with_locations_path}...")
        ft = pd.read_csv(final_table_path, sep="\t")
        if self.splitting:
            required = ["Sequence", "Reference", "Subtype", "Most Matching Gene Region", "Present Gene Regions"]
        else:
            required = ["Sequence", "Reference", "Subtype"]
        missing = [c for c in required if c not in ft.columns]
        if missing:
            raise ValueError(f"Missing columns in final_table: {missing}")

        swl = pd.read_csv(sequences_with_locations_path, sep="\t")
        if "accession" not in swl.columns or "features" not in swl.columns:
            raise ValueError("sequences_with_locations.tsv must include 'accession' and 'features' columns")
        
        # Parse features by accession
        features_by_acc = {}
        for _, row in swl.iterrows():
            acc = str(row["accession"])
            try:
                features_by_acc[acc] = parse_features(row["features"])
            except Exception as e:
                self.logger.warning(f"Failed to parse features for accession {acc}: {e}")
                features_by_acc[acc] = {}

        # Generate PDF
        output_pdf_path = self.output_dir / output_pdf_name
        self.logger.info(f"Generating PDF report: {output_pdf_path}")
        with PdfPages(output_pdf_path) as pdf:
            pages_made = 0

            for _, r in ft.iterrows():
                file_name = str(r["File Name"]) if "File Name" in r else "-"
                sequence = str(r["Sequence"])
                accession = str(r["Reference"])
                splitting_accession = str(r.get("Splitting Reference", accession))
                group = str(r["Group"]) if self.subtyping and "Group" in r else "No subtyping performed."
                subtype = str(r["Subtype"]) if self.subtyping else "No subtyping performed."
                closest_subtypes = str(r["Closest Subtypes"]) if "Closest Subtypes" in r else "-"
                subtype_score_warning_value = r.get("Subtype Score Warning", "")
                subtype_score_warning = (
                    ""
                    if pd.isna(subtype_score_warning_value)
                    else str(subtype_score_warning_value)
                )
                mm_region = str(r["Most Matching Gene Region"]) if "Most Matching Gene Region" in r else "-"
                present_regions_raw = parse_present_regions(r.get("Present Gene Regions", "")) if self.splitting else []

                # Find alignment file
                alignment_prefix = (
                    "splitting_alignment"
                    if self.splitting and splitting_accession not in {"", "-", accession}
                    else "best_alignment"
                )
                fasta_path = build_alignment_path(
                    sequence,
                    self.output_dir,
                    prefix=alignment_prefix,
                    file_name=file_name,
                )
                if not fasta_path.exists():
                    self.logger.warning(f"Alignment FASTA not found for {sequence}: {fasta_path}")
                    continue

                try:
                    ref_header, ref_seq_aln, user_header, user_seq_aln = read_alignment_fasta(fasta_path)
                except Exception as e:
                    self.logger.error(f"Error reading {fasta_path}: {e}")
                    continue

                special = is_special_reference(splitting_accession, ref_header)

                ref_map, _ = build_ref_to_alignment_map(ref_seq_aln)

                raw_features = features_by_acc.get(splitting_accession, {})
                features_genomic = normalize_features(raw_features, special)
                present_regions = normalize_present_regions(present_regions_raw, special)

                features_aln = project_features_to_alignment(features_genomic, ref_map)

                if special:
                    features_aln = {k: v for k, v in features_aln.items() if k in K03455Config.TARGET_REGIONS}
                    present_regions = [g for g in present_regions if g in features_aln]
                    y_pos = K03455Config.Y_POSITIONS
                else:
                    y_pos = None  # non-K03455 → auto lanes + configurable numeric offsets + fixed x >= 10,000

                render_sequence_page(
                    pdf=pdf,
                    sequence=sequence,
                    file_name=file_name,
                    accession=accession,
                    group=group,
                    subtype=subtype,
                    closest_subtypes=closest_subtypes,
                    subtype_score_warning=subtype_score_warning,
                    mm_region=mm_region if mm_region != "-" else "",
                    present_regions=present_regions,
                    features_aln=features_aln,
                    ref_seq_aligned=ref_seq_aln,
                    user_seq_aligned=user_seq_aln,
                    y_positions=y_pos,
                    splitting_accession=splitting_accession if self.splitting else None,
                )

                pages_made += 1
                self.logger.info(
                    f"Added page for {sequence} (special={special}, features={list(features_aln.keys())})"
                )

        if pages_made == 0:
            self.logger.warning("No pages created. Check your paths and file formats.")
        else:
            self.logger.info(f"PDF created: {output_pdf_path} (pages: {pages_made})")
        
        return output_pdf_path

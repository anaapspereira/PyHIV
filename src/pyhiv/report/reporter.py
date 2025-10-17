"""
Main reporting class for PyHIV results.
"""

from pathlib import Path

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .constants import K03455Config
from .utils import (
    read_alignment_fasta, parse_present_regions, parse_features,
    is_special_reference, normalize_features, normalize_present_regions,
    build_ref_to_alignment_map, project_features_to_alignment,
    build_alignment_path
)
from .pdf_generator import render_sequence_page


class PyHIVReporter:
    """Main class for generating PyHIV PDF reports."""
    
    def __init__(self, output_dir: Path):
        """Initialize the reporter with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_report(
        self,
        final_table_path: Path,
        sequences_with_locations_path: Path,
        output_pdf_name: str = "PyHIV_report_all_sequences.pdf"
    ) -> Path:
        """
        Generate PDF report from PyHIV results.
        
        Args:
            final_table_path: Path to final_table.tsv
            sequences_with_locations_path: Path to sequences_with_locations.tsv
            output_pdf_name: Name of output PDF file
            
        Returns:
            Path to generated PDF file
        """
        # Read input data
        ft = pd.read_csv(final_table_path, sep="\t")
        required = ["Sequence", "Reference", "Subtype", "Most Matching Gene Region", "Present Gene Regions"]
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
            except Exception:
                features_by_acc[acc] = {}

        # Generate PDF
        output_pdf_path = self.output_dir / output_pdf_name
        with PdfPages(output_pdf_path) as pdf:
            pages_made = 0

            for _, r in ft.iterrows():
                sequence = str(r["Sequence"])
                accession = str(r["Reference"])
                subtype = str(r["Subtype"])
                mm_region = str(r["Most Matching Gene Region"]) if "Most Matching Gene Region" in r else "-"
                present_regions_raw = parse_present_regions(r.get("Present Gene Regions", ""))

                # Find alignment file
                fasta_path = build_alignment_path(sequence, self.output_dir)
                if not fasta_path.exists():
                    print(f"[!] Alignment FASTA not found for {sequence}: {fasta_path}")
                    continue

                try:
                    ref_header, ref_seq_aln, user_header, user_seq_aln = read_alignment_fasta(fasta_path)
                except Exception as e:
                    print(f"[!] Error reading {fasta_path}: {e}")
                    continue

                special = is_special_reference(accession, ref_header)

                ref_map, _ = build_ref_to_alignment_map(ref_seq_aln)

                raw_features = features_by_acc.get(accession, {})
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
                    accession=accession,
                    subtype=subtype,
                    mm_region=mm_region if mm_region != "-" else "",
                    present_regions=present_regions,
                    features_aln=features_aln,
                    ref_header=ref_header,
                    ref_seq_aligned=ref_seq_aln,
                    user_header=user_header,
                    user_seq_aligned=user_seq_aln,
                    y_positions=y_pos,
                )

                pages_made += 1
                print(f"[OK] Added page for: {sequence} (special={special}, features={list(features_aln.keys())})")

        if pages_made == 0:
            print("[!] No pages created. Check your paths and file formats.")
        else:
            print(f"[OK] Single PDF created: {output_pdf_path} (pages: {pages_made})")
        
        return output_pdf_path

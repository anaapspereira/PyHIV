"""
Utility functions for PyHIV reporting module.
"""

import ast
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd


def ungap(seq: str) -> str:
    """Remove gaps from sequence."""
    return seq.replace("-", "").replace(".", "")


def first_last_nongap_idx(seq: str) -> Tuple[int, int]:
    """Find first and last non-gap positions in sequence."""
    if not seq or set(seq).issubset({"-", "."}):
        return 0, 0
    try:
        first = next((i for i, c in enumerate(seq) if c not in "-."), 0)
    except StopIteration:
        return 0, 0
    try:
        last = len(seq) - 1 - next((i for i, c in enumerate(reversed(seq)) if c not in "-."), 0)
    except StopIteration:
        last = first
    return first, last


def read_alignment_fasta(fpath: Path) -> Tuple[str, str, str, str]:
    """Read alignment FASTA file and return headers and sequences."""
    if not fpath.exists():
        raise FileNotFoundError(f"Alignment FASTA not found: {fpath}")

    headers, seqs, cur, cur_header = [], [], [], None
    with open(fpath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_header is not None:
                    seqs.append("".join(cur))
                    cur = []
                cur_header = line[1:].strip()
                headers.append(cur_header)
            else:
                cur.append(line)
        if cur_header is not None:
            seqs.append("".join(cur))

    if len(seqs) != 2:
        raise ValueError(f"Expected 2 sequences in {fpath}, found {len(seqs)}.")

    idx_ref = 0 if "reference" in headers[0].lower() else (1 if "reference" in headers[1].lower() else 0)
    idx_usr = 1 - idx_ref
    return headers[idx_ref], seqs[idx_ref], headers[idx_usr], seqs[idx_usr]


def parse_present_regions(cell) -> List[str]:
    """Parse present regions from table cell."""
    if cell is None:
        return []
    cell = str(cell).strip()
    if not cell or cell == "-":
        return []
    try:
        val = ast.literal_eval(cell)
        if isinstance(val, (list, tuple)):
            return [str(x).strip().strip("'").strip('"') for x in val]
        if isinstance(val, str):
            cell = val
    except Exception:
        pass
    return [p.strip().strip("'").strip('"') for p in cell.split(",") if p.strip()]


def parse_features(cell) -> Dict[str, Tuple[int, int]]:
    """Parse features from table cell."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return {}
    if isinstance(cell, dict):
        return {str(k): (int(v[0]), int(v[1])) for k, v in cell.items()}
    d = ast.literal_eval(str(cell))
    return {str(k): (int(v[0]), int(v[1])) for k, v in d.items()}


def is_special_reference(accession: str, ref_header: str) -> bool:
    """Check if reference is special (K03455)."""
    return (accession or "").strip() == "K03455" or "K03455-B" in (ref_header or "")


def canon_label(label: str) -> Optional[str]:
    """Canonicalize gene label for K03455."""
    from .constants import K03455Config
    
    # Patterns for canonicalization
    patterns = [
        (re.compile(r"^\s*5\s*'? *ltr\s*$", re.I), "5' LTR"),
        (re.compile(r"^\s*gag\s*$", re.I), "gag"),
        (re.compile(r"^\s*pol(\s*cds)?\s*$", re.I), "pol"),
        (re.compile(r"^\s*vif(\s*cds)?\s*$", re.I), "vif"),
        (re.compile(r"^\s*vpr(\s*cds)?\s*$", re.I), "vpr"),
        (re.compile(r"^\s*vpu(\s*cds)?\s*$", re.I), "vpu"),
        # Accept either numeric or roman numerals (i/ii) for exon indices
        (re.compile(r"^\s*tat(\s*exon)?\s*(?:1|i)\s*$", re.I), "tat 1"),
        (re.compile(r"^\s*tat(\s*exon)?\s*(?:2|ii)\s*$", re.I), "tat 2"),
        (re.compile(r"^\s*rev(\s*exon)?\s*(?:1|i)\s*$", re.I), "rev 1"),
        (re.compile(r"^\s*rev(\s*exon)?\s*(?:2|ii)\s*$", re.I), "rev 2"),
        (re.compile(r"^\s*env(\s*cds)?\s*$", re.I), "env"),
        (re.compile(r"^\s*nef(\s*cds)?\s*$", re.I), "nef"),
        (re.compile(r"^\s*3\s*'? *ltr\s*$", re.I), "3' LTR"),
    ]
    
    s = (label or "").strip()
    for pat, out in patterns:
        if pat.match(s):
            return out
    if s in K03455Config.TARGET_REGIONS:
        return s
    s_lower = s.lower()
    for x in K03455Config.TARGET_REGIONS:
        if x.lower() == s_lower:
            return x
    return None


def normalize_features(raw_features: Dict[str, Tuple[int, int]], special: bool) -> Dict[str, Tuple[int, int]]:
    """Normalize features based on reference type."""
    from .constants import K03455Config
    
    raw_features = raw_features or {}
    if not special:
        return {str(k): (int(v[0]), int(v[1])) for k, v in raw_features.items()}
    
    out = {}
    for k, (s, e) in raw_features.items():
        canon = canon_label(k)
        if canon in K03455Config.TARGET_REGIONS:
            out[canon] = (int(s), int(e))
    return {k: v for k, v in out.items() if k in K03455Config.TARGET_REGIONS}


def normalize_present_regions(regions: List[str], special: bool) -> List[str]:
    """Normalize present regions based on reference type."""
    from .constants import K03455Config
    
    regions = regions or []
    if not special:
        return regions
    
    out = []
    for r in regions:
        canon = canon_label(r)
        if canon in K03455Config.TARGET_REGIONS:
            out.append(canon)
    return out


def build_ref_to_alignment_map(ref_aligned: str) -> Tuple[Dict[int, int], int]:
    """Build mapping from reference coordinates to alignment coordinates."""
    mapping = {}
    ref_pos = 0
    for aln_idx, ch in enumerate(ref_aligned):
        if ch not in "-.":
            ref_pos += 1
            mapping[ref_pos] = aln_idx
    return mapping, len(ref_aligned)


def project_features_to_alignment(features_genomic: Dict[str, Tuple[int, int]], ref_map: Dict[int, int]) -> Dict[str, Tuple[int, int]]:
    """Project genomic features to alignment coordinates."""
    projected = {}
    for gene, (gstart, gend) in features_genomic.items():
        if gstart in ref_map and gend in ref_map:
            astart, aend = ref_map[gstart], ref_map[gend]
            if aend >= astart:
                projected[gene] = (astart, aend)
    return projected


def get_numeric_offsets_non_special(gene: str) -> Tuple[float, float]:
    """Get numeric offsets for non-K03455 references."""
    from .constants import NumericOffsets
    
    if not gene:
        return NumericOffsets.DEFAULT_OFFSETS
    key = gene.strip().lower()
    if key in NumericOffsets.GENE_OFFSET_MAP:
        return NumericOffsets.GENE_OFFSET_MAP[key]
    if key.startswith("tat"):
        return NumericOffsets.GENE_OFFSET_MAP.get("tat 1", NumericOffsets.DEFAULT_OFFSETS)
    if key.startswith("rev"):
        return NumericOffsets.GENE_OFFSET_MAP.get("rev 1", NumericOffsets.DEFAULT_OFFSETS)
    return NumericOffsets.DEFAULT_OFFSETS


def build_alignment_path(sequence: str, alignments_dir: Path) -> Path:
    """Build path to alignment FASTA file."""
    p = alignments_dir / f"best_alignment_{sequence}.fasta"
    return p if p.exists() else (alignments_dir / f"{sequence}.fasta")

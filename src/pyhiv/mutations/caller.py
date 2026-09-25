"""Mutation calling against HIVDB Consensus B protein references."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyhiv import __version__
from pyhiv.loading import discover_fasta_files
from pyhiv.mutations.drm import (
    DRMComponent,
    UNRESOLVED_CODON_STATUSES,
    annotate_drm,
    drm_evaluation_status,
    drm_status,
    load_drm_catalog,
)
from pyhiv.mutations.references import (
    REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    normalize_gene,
    reference_for_gene,
)


MUTATION_TSV_COLUMNS = [
    "sequence_id",
    "gene",
    "position",
    "ref_aa",
    "alt_aa",
    "mutation",
    "mutation_type",
    "mixture",
    "insertion",
    "deletion",
    "stop",
    "hxb2_ref_codon",
    "query_codon",
    "possible_alt_aas",
    "codon_status",
    "inserted_nts",
    "deleted_nt_count",
    "phase_at_start",
    "phase_at_end",
    "reference_system",
    "coordinate_system",
    "is_drm",
    "drm_status",
    "drm_evaluation_status",
    "drm_source",
    "drm_version",
    "drm_catalog_sha256",
    "drm_components",
    "mutation_types",
    "is_accessory",
    "drm_class",
    "drug_class",
    "caller",
    "caller_version",
    "algorithm_version",
    "qc_status",
]
DRM_SCREENING_TSV_COLUMNS = [
    "sequence_id",
    "gene",
    "position",
    "ref_aa",
    "alt_aa",
    "codon_status",
    "status",
]
DRM_SCREENING_SUMMARY_TSV_COLUMNS = [
    "sequence_id",
    "gene",
    "drm_evaluation_status",
    "resolved_wt",
    "resolved_mutant",
    "unresolved",
    "not_covered",
    "positions",
]
MUTATION_POSITION_QC_TSV_COLUMNS = [
    "sequence_id",
    "gene",
    "HXB2_position",
    "Consensus_B_ref_aa",
    "observed_aa",
    "mutation",
    "query_codon",
    "codon_status",
    "coverage_status",
    "ambiguity",
    "insertion",
    "deletion",
    "stop",
    "is_drm",
    "drm_class",
    "drug_class",
    "qc_status",
]

AA_ALPHABET = set("ABCDEFGHIKLMNPQRSTVWXYZ*-")
NT_ALPHABET = set("ACGTRYSWKMBDHVN-")
HXB2_GENE_COORDINATES = {
    "PR": (2253, 2549),
    "RT": (2550, 4229),
    "IN": (4230, 5093),
    "CA": (1186, 1878),
}
IUPAC_DNA = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

DEFAULT_NT_ALIGNMENT_BACKEND = "semiglobal_mismatch2_open20_trimmed"

ALIGNMENT_BACKENDS = {
    "builtin",
    "postalign",
    "semiglobal_open20",
    "semiglobal_mismatch2_open20",
    "semiglobal_mismatch2_open20_trimmed",
    "semiglobal_codonaware",
    "semiglobal_codonaware_local",
    "semiglobal_segmented",
    "semiglobal_anchorblocks",
    "semiglobal_coverage_confidence",
}
POSTALIGN_PROGRAM_ENV = "POSTALIGN_PROGRAM"
DRM_SCREENING_RESOLVED_WT = "RESOLVED_WT"
DRM_SCREENING_RESOLVED_MUTANT = "RESOLVED_MUTANT"
DRM_SCREENING_UNRESOLVED = "UNRESOLVED"
DRM_SCREENING_NOT_COVERED = "NOT_COVERED"
SEMIGLOBAL_MATCH_SCORE = 2
SEMIGLOBAL_GAP_OPEN_SCORE = -20
SEMIGLOBAL_GAP_EXTEND_SCORE = -1
CODONAWARE_MISMATCH_SCORE = -2
CODONAWARE_FRAMESHIFT_GAP_SCORE = -12



@dataclass(frozen=True)
class MutationCall:
    sequence_id: str
    gene: str
    position: int
    ref_aa: str
    alt_aa: str
    mutation: str
    mutation_type: str
    mixture: bool
    insertion: bool
    deletion: bool
    stop: bool
    hxb2_ref_codon: str
    query_codon: str
    possible_alt_aas: str
    codon_status: str
    inserted_nts: str
    deleted_nt_count: int
    phase_at_start: int
    phase_at_end: int
    reference_system: str
    coordinate_system: str
    is_drm: bool | None
    drm_source: str
    drm_version: str
    drm_catalog_sha256: str
    drm_components: tuple[DRMComponent, ...]
    mutation_types: str
    is_accessory: bool | None
    drm_class: str
    drug_class: str
    caller: str
    caller_version: str
    algorithm_version: str
    qc_status: str

    @property
    def drm_status(self) -> str:
        return drm_status(self.is_drm)

    @property
    def drm_evaluation_status(self) -> str:
        if not self.drm_components:
            return "UNRESOLVED" if self.is_drm is None else "COMPLETE"
        return drm_evaluation_status(self.drm_components)

    def to_tsv_row(self) -> dict[str, object]:
        row = asdict(self)
        row["drm_status"] = self.drm_status
        row["drm_evaluation_status"] = self.drm_evaluation_status
        row["drm_components"] = json.dumps([
            {**asdict(component), "drm_status": component.drm_status}
            for component in self.drm_components
        ], sort_keys=True, separators=(",", ":"))
        for field in ("is_drm", "is_accessory"):
            if row[field] is None:
                row[field] = "UNRESOLVED"
        return {
            key: str(value).lower() if isinstance(value, bool) else value
            for key, value in row.items()
        }


@dataclass(frozen=True)
class DRMScreeningPosition:
    sequence_id: str
    gene: str
    position: int
    ref_aa: str
    alt_aa: str
    codon_status: str
    status: str


@dataclass(frozen=True)
class DRMScreeningSummary:
    sequence_id: str
    gene: str
    status: str
    positions: tuple[DRMScreeningPosition, ...]


@dataclass(frozen=True)
class MutationPositionQC:
    sequence_id: str
    gene: str
    HXB2_position: int
    Consensus_B_ref_aa: str
    observed_aa: str
    mutation: str
    query_codon: str
    codon_status: str
    coverage_status: str
    ambiguity: str
    insertion: bool
    deletion: bool
    stop: bool
    is_drm: bool | None
    drm_class: str
    drug_class: str
    qc_status: str

    def to_tsv_row(self) -> dict[str, object]:
        row = asdict(self)
        if row["is_drm"] is None:
            row["is_drm"] = "UNRESOLVED"
        return {
            key: str(value).lower() if isinstance(value, bool) else value
            for key, value in row.items()
        }


def call_mutations_for_record(
    record: SeqRecord,
    gene: str,
    sequence_type: str = "auto",
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[MutationCall]:
    """Call amino-acid mutations for one sequence record."""
    return call_mutations(
        sequence=str(record.seq),
        sequence_id=record.id,
        gene=gene,
        sequence_type=sequence_type,
        reference_system=reference_system,
        alignment_backend=alignment_backend,
    )


def call_mutations_for_record_with_drm_screening(
    record: SeqRecord,
    gene: str,
    sequence_type: str = "auto",
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], DRMScreeningSummary]:
    return call_mutations_with_drm_screening(
        sequence=str(record.seq),
        sequence_id=record.id,
        gene=gene,
        sequence_type=sequence_type,
        reference_system=reference_system,
        alignment_backend=alignment_backend,
    )


def call_mutations_for_record_with_qc(
    record: SeqRecord,
    gene: str,
    sequence_type: str = "auto",
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], DRMScreeningSummary, list[MutationPositionQC]]:
    return call_mutations_with_qc(
        sequence=str(record.seq),
        sequence_id=record.id,
        gene=gene,
        sequence_type=sequence_type,
        reference_system=reference_system,
        alignment_backend=alignment_backend,
    )


def call_mutations(
    sequence: str,
    sequence_id: str,
    gene: str,
    sequence_type: str = "auto",
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[MutationCall]:
    """Call amino-acid substitutions, mixtures, insertions, deletions and stops."""
    canonical_gene = normalize_gene(gene)
    alignment_backend = normalize_alignment_backend(alignment_backend)
    if reference_system != REFERENCE_SYSTEM_HIVDB_CONSENSUS_B:
        raise ValueError(f"Unsupported mutation reference system: {reference_system!r}.")

    reference = reference_for_gene(canonical_gene)
    cleaned = "".join(str(sequence).upper().split())
    effective_sequence_type = sequence_type.lower()
    if effective_sequence_type == "auto":
        effective_sequence_type = infer_sequence_type(cleaned)
    if effective_sequence_type == "nt":
        return mutations_from_unaligned_dna(
            sequence_id=sequence_id,
            gene=canonical_gene,
            sequence=cleaned,
            reference_protein=reference,
            reference_system=reference_system,
            alignment_backend=alignment_backend,
        )

    protein, qc_status = protein_sequence_for_calling(
        cleaned,
        reference,
        sequence_type=effective_sequence_type,
    )
    query_aligned, ref_aligned = align_proteins(protein, reference)
    return mutations_from_alignment(
        sequence_id=sequence_id,
        gene=canonical_gene,
        query_aligned=query_aligned,
        ref_aligned=ref_aligned,
        qc_status=qc_status,
        reference_system=reference_system,
    )


def call_mutations_with_drm_screening(
    sequence: str,
    sequence_id: str,
    gene: str,
    sequence_type: str = "auto",
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], DRMScreeningSummary]:
    """Call mutations and summarize DRM screening over all catalogue positions."""
    canonical_gene = normalize_gene(gene)
    alignment_backend = normalize_alignment_backend(alignment_backend)
    if reference_system != REFERENCE_SYSTEM_HIVDB_CONSENSUS_B:
        raise ValueError(f"Unsupported mutation reference system: {reference_system!r}.")

    reference = reference_for_gene(canonical_gene)
    cleaned = "".join(str(sequence).upper().split())
    effective_sequence_type = sequence_type.lower()
    if effective_sequence_type == "auto":
        effective_sequence_type = infer_sequence_type(cleaned)
    if effective_sequence_type == "nt":
        return mutations_and_drm_screening_from_unaligned_dna(
            sequence_id=sequence_id,
            gene=canonical_gene,
            sequence=cleaned,
            reference_protein=reference,
            reference_system=reference_system,
            alignment_backend=alignment_backend,
        )

    protein, qc_status = protein_sequence_for_calling(
        cleaned,
        reference,
        sequence_type=effective_sequence_type,
    )
    query_aligned, ref_aligned = align_proteins(protein, reference)
    calls = mutations_from_alignment(
        sequence_id=sequence_id,
        gene=canonical_gene,
        query_aligned=query_aligned,
        ref_aligned=ref_aligned,
        qc_status=qc_status,
        reference_system=reference_system,
    )
    screening = summarize_drm_screening(
        sequence_id=sequence_id,
        gene=canonical_gene,
        codon_calls=codon_calls_from_protein_alignment(query_aligned, ref_aligned, reference),
        mutation_calls=calls,
        reference_protein=reference,
    )
    return calls, screening


def call_mutations_with_qc(
    sequence: str,
    sequence_id: str,
    gene: str,
    sequence_type: str = "auto",
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], DRMScreeningSummary, list[MutationPositionQC]]:
    """Call mutations, DRM screening, and per-position coverage/QC rows."""
    canonical_gene = normalize_gene(gene)
    alignment_backend = normalize_alignment_backend(alignment_backend)
    if reference_system != REFERENCE_SYSTEM_HIVDB_CONSENSUS_B:
        raise ValueError(f"Unsupported mutation reference system: {reference_system!r}.")

    reference = reference_for_gene(canonical_gene)
    cleaned = "".join(str(sequence).upper().split())
    effective_sequence_type = sequence_type.lower()
    if effective_sequence_type == "auto":
        effective_sequence_type = infer_sequence_type(cleaned)
    if effective_sequence_type == "nt":
        return mutations_screening_and_position_qc_from_unaligned_dna(
            sequence_id=sequence_id,
            gene=canonical_gene,
            sequence=cleaned,
            reference_protein=reference,
            reference_system=reference_system,
            alignment_backend=alignment_backend,
        )

    protein, qc_status = protein_sequence_for_calling(
        cleaned,
        reference,
        sequence_type=effective_sequence_type,
    )
    query_aligned, ref_aligned = align_proteins(protein, reference)
    calls = mutations_from_alignment(
        sequence_id=sequence_id,
        gene=canonical_gene,
        query_aligned=query_aligned,
        ref_aligned=ref_aligned,
        qc_status=qc_status,
        reference_system=reference_system,
    )
    codon_calls = codon_calls_from_protein_alignment(query_aligned, ref_aligned, reference)
    screening = summarize_drm_screening(
        sequence_id=sequence_id,
        gene=canonical_gene,
        codon_calls=codon_calls,
        mutation_calls=calls,
        reference_protein=reference,
    )
    position_qc = mutation_position_qc_rows(
        sequence_id=sequence_id,
        gene=canonical_gene,
        codon_calls=codon_calls,
        mutation_calls=calls,
        reference_protein=reference,
        qc_status=qc_status,
    )
    return calls, screening, position_qc


def protein_sequence_for_calling(
    sequence: str,
    reference_protein: str,
    sequence_type: str = "auto",
) -> tuple[str, str]:
    """Return the protein sequence and a coarse QC status."""
    sequence_type = sequence_type.lower()
    if sequence_type not in {"auto", "nt", "aa"}:
        raise ValueError("sequence_type must be one of: auto, nt, aa.")

    cleaned = "".join(str(sequence).upper().split())
    if sequence_type == "auto":
        sequence_type = infer_sequence_type(cleaned)

    if sequence_type == "aa":
        return cleaned, "PASS"

    nucleotide_sequence = cleaned.replace("-", "")
    translated_frames: list[tuple[float, str, str]] = []
    for frame in range(3):
        framed = nucleotide_sequence[frame:]
        trimmed_length = len(framed) - (len(framed) % 3)
        coding = framed[:trimmed_length]
        protein = translate_ambiguous_dna(coding)
        qc_status = "PASS" if len(framed) % 3 == 0 else "WARN"
        translated_frames.append((alignment_score(protein, reference_protein), protein, qc_status))

    translated_frames.sort(key=lambda item: item[0], reverse=True)
    return translated_frames[0][1], translated_frames[0][2]


def infer_sequence_type(sequence: str) -> str:
    """Infer nucleotide versus amino-acid sequence input."""
    letters = set(sequence)
    if letters <= NT_ALPHABET and len(sequence.replace("-", "")) >= 3:
        return "nt"
    return "aa"


def translate_ambiguous_dna(sequence: str) -> str:
    """Translate DNA, preserving ambiguous codons as X."""
    sequence = sequence.replace("-", "")
    amino_acids = []
    for index in range(0, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        translated = translate_codon(codon)
        amino_acids.append(translated if len(translated) == 1 else "X")
    return "".join(amino_acids)


def translate_reference_aligned_dna(sequence: str) -> str:
    """Translate a gapped coding-region alignment in reference codon columns."""
    amino_acids = []
    for index in range(0, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        if codon == "---":
            amino_acids.append("-")
        elif "-" in codon:
            amino_acids.append("X")
        else:
            amino_acids.append(translate_codon(codon))
    return "".join(amino_acids)


def qc_status_for_protein(protein: str) -> str:
    return "WARN" if "*" in protein or "X" in protein else "PASS"


def translate_codon(codon: str) -> str:
    """Translate a concrete or IUPAC-ambiguous codon to one or more amino acids."""
    codon = codon.upper()
    if len(codon) != 3 or any(base not in IUPAC_DNA for base in codon):
        return "X"

    possibilities = [""]
    for base in codon:
        possibilities = [
            prefix + resolved
            for prefix in possibilities
            for resolved in IUPAC_DNA[base]
        ]

    amino_acids = {str(Seq(possible).translate()) for possible in possibilities}
    if len(amino_acids) > 4:
        return "X"
    if len(amino_acids) == 1:
        return next(iter(amino_acids))
    return "".join(sorted(amino_acids))


def mutations_from_unaligned_dna(
    sequence_id: str,
    gene: str,
    sequence: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[MutationCall]:
    """Call mutations from raw/unaligned nucleotide sequence.

    Raw HIV fragments are first mapped to HXB2 nucleotide coordinates, but
    mutation calls are then made independently per HXB2 codon against the
    HIVDB Consensus B amino acid.  Importantly, a 1- or 2-nt gap introduced by
    the nucleotide alignment is treated as a *local uncertain codon* here; it
    is not propagated as a frameshift through every downstream position.

    Long-range frame state is only interpreted when the caller is given an
    explicit pre-existing HXB2 alignment via ``call_mutations_from_hxb2_alignment``.
    This keeps raw-fragment calling conservative and avoids hundreds of false X
    calls caused by alternative gap placement in divergent HIV sequences.
    """
    reference_nt = hxb2_reference_nt_for_gene(gene)
    query_aligned, ref_aligned = align_dna_to_hxb2_columns(
        sequence,
        reference_nt,
        alignment_backend=alignment_backend,
    )
    force_not_covered_codons = (
        coverage_confidence_not_covered_codons(query_aligned, ref_aligned, reference_protein)
        if normalize_alignment_backend(alignment_backend) == "semiglobal_coverage_confidence"
        else None
    )
    return mutations_from_raw_hxb2_alignment(
        sequence_id=sequence_id,
        gene=gene,
        query_aligned=query_aligned,
        hxb2_ref_aligned=ref_aligned,
        reference_protein=reference_protein,
        reference_system=reference_system,
        force_not_covered_codons=force_not_covered_codons,
    )


def mutations_and_drm_screening_from_unaligned_dna(
    sequence_id: str,
    gene: str,
    sequence: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], DRMScreeningSummary]:
    reference_nt = hxb2_reference_nt_for_gene(gene)
    query_aligned, ref_aligned = align_dna_to_hxb2_columns(
        sequence,
        reference_nt,
        alignment_backend=alignment_backend,
    )
    force_not_covered_codons = (
        coverage_confidence_not_covered_codons(query_aligned, ref_aligned, reference_protein)
        if normalize_alignment_backend(alignment_backend) == "semiglobal_coverage_confidence"
        else None
    )
    return mutations_and_drm_screening_from_raw_hxb2_alignment(
        sequence_id=sequence_id,
        gene=gene,
        query_aligned=query_aligned,
        hxb2_ref_aligned=ref_aligned,
        reference_protein=reference_protein,
        reference_system=reference_system,
        force_not_covered_codons=force_not_covered_codons,
    )


def mutations_screening_and_position_qc_from_unaligned_dna(
    sequence_id: str,
    gene: str,
    sequence: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], DRMScreeningSummary, list[MutationPositionQC]]:
    reference_nt = hxb2_reference_nt_for_gene(gene)
    query_aligned, ref_aligned = align_dna_to_hxb2_columns(
        sequence,
        reference_nt,
        alignment_backend=alignment_backend,
    )
    force_not_covered_codons = (
        coverage_confidence_not_covered_codons(query_aligned, ref_aligned, reference_protein)
        if normalize_alignment_backend(alignment_backend) == "semiglobal_coverage_confidence"
        else None
    )
    return mutations_screening_and_position_qc_from_raw_hxb2_alignment(
        sequence_id=sequence_id,
        gene=gene,
        query_aligned=query_aligned,
        hxb2_ref_aligned=ref_aligned,
        reference_protein=reference_protein,
        reference_system=reference_system,
        force_not_covered_codons=force_not_covered_codons,
    )


def mutations_from_raw_hxb2_alignment(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    hxb2_ref_aligned: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    force_not_covered_codons: set[int] | None = None,
) -> list[MutationCall]:
    """Call mutations from a newly-created raw-query to HXB2 alignment.

    This is deliberately codon-local.  HXB2 determines the coordinate columns;
    Consensus B determines whether the observed amino acid is a mutation.
    """
    codon_calls = codon_calls_from_raw_nucleotide_alignment(
        query_aligned=query_aligned,
        ref_aligned=hxb2_ref_aligned,
        reference_protein=reference_protein,
        force_not_covered_codons=force_not_covered_codons,
    )
    rows = mutation_calls_from_raw_codon_calls(
        sequence_id=sequence_id,
        gene=gene,
        codon_calls=codon_calls,
        reference_protein=reference_protein,
        reference_system=reference_system,
    )
    rows.extend(
        insertion_calls_from_nucleotide_alignment(
            sequence_id=sequence_id,
            gene=gene,
            query_aligned=query_aligned,
            ref_aligned=hxb2_ref_aligned,
            qc_status=qc_status_for_codon_calls(codon_calls),
            reference_system=reference_system,
        )
    )
    return rows


def mutations_and_drm_screening_from_raw_hxb2_alignment(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    hxb2_ref_aligned: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    force_not_covered_codons: set[int] | None = None,
) -> tuple[list[MutationCall], DRMScreeningSummary]:
    codon_calls = codon_calls_from_raw_nucleotide_alignment(
        query_aligned=query_aligned,
        ref_aligned=hxb2_ref_aligned,
        reference_protein=reference_protein,
        force_not_covered_codons=force_not_covered_codons,
    )
    qc_status = qc_status_for_codon_calls(codon_calls)
    rows = mutation_calls_from_raw_codon_calls(
        sequence_id=sequence_id,
        gene=gene,
        codon_calls=codon_calls,
        reference_protein=reference_protein,
        reference_system=reference_system,
    )
    rows.extend(
        insertion_calls_from_nucleotide_alignment(
            sequence_id=sequence_id,
            gene=gene,
            query_aligned=query_aligned,
            ref_aligned=hxb2_ref_aligned,
            qc_status=qc_status,
            reference_system=reference_system,
        )
    )
    screening = summarize_drm_screening(
        sequence_id=sequence_id,
        gene=gene,
        codon_calls=codon_calls,
        mutation_calls=rows,
        reference_protein=reference_protein,
    )
    return rows, screening


def mutations_screening_and_position_qc_from_raw_hxb2_alignment(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    hxb2_ref_aligned: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
    force_not_covered_codons: set[int] | None = None,
) -> tuple[list[MutationCall], DRMScreeningSummary, list[MutationPositionQC]]:
    codon_calls = codon_calls_from_raw_nucleotide_alignment(
        query_aligned=query_aligned,
        ref_aligned=hxb2_ref_aligned,
        reference_protein=reference_protein,
        force_not_covered_codons=force_not_covered_codons,
    )
    qc_status = qc_status_for_codon_calls(codon_calls)
    rows = mutation_calls_from_raw_codon_calls(
        sequence_id=sequence_id,
        gene=gene,
        codon_calls=codon_calls,
        reference_protein=reference_protein,
        reference_system=reference_system,
    )
    rows.extend(
        insertion_calls_from_nucleotide_alignment(
            sequence_id=sequence_id,
            gene=gene,
            query_aligned=query_aligned,
            ref_aligned=hxb2_ref_aligned,
            qc_status=qc_status,
            reference_system=reference_system,
        )
    )
    screening = summarize_drm_screening(
        sequence_id=sequence_id,
        gene=gene,
        codon_calls=codon_calls,
        mutation_calls=rows,
        reference_protein=reference_protein,
    )
    position_qc = mutation_position_qc_rows(
        sequence_id=sequence_id,
        gene=gene,
        codon_calls=codon_calls,
        mutation_calls=rows,
        reference_protein=reference_protein,
        qc_status=qc_status,
    )
    return rows, screening, position_qc


def qc_status_for_codon_calls(codon_calls: list[dict[str, object]]) -> str:
    protein_tokens = [call["alt_aa"] for call in codon_calls]
    return "WARN" if any("*" in token or "X" in token for token in protein_tokens) else "PASS"


def mutation_calls_from_raw_codon_calls(
    sequence_id: str,
    gene: str,
    codon_calls: list[dict[str, object]],
    reference_protein: str,
    reference_system: str,
) -> list[MutationCall]:
    qc_status = qc_status_for_codon_calls(codon_calls)
    rows: list[MutationCall] = []
    for index, codon_call in enumerate(codon_calls[:len(reference_protein)]):
        status = codon_call["codon_status"]
        if status == "NOT_COVERED":
            continue

        alt_aa = codon_call["alt_aa"]
        ref_aa = reference_protein[index]
        position = index + 1
        if status == "DELETION":
            mutation_type = "deletion"
        elif "PARTIAL_CODON" in status:
            mutation_type = "partial_codon"
        else:
            mutation_type = "substitution"

        if alt_aa == "-" or should_report_substitution(ref_aa, normalize_alt_aa(ref_aa, alt_aa)):
            rows.append(
                build_mutation_call(
                    sequence_id=sequence_id,
                    gene=gene,
                    position=position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    mutation_type=mutation_type,
                    qc_status=qc_status,
                    reference_system=reference_system,
                    hxb2_ref_codon=codon_call["hxb2_ref_codon"],
                    query_codon=codon_call["query_codon"],
                    possible_alt_aas=codon_call["possible_alt_aas"],
                    codon_status=status,
                    inserted_nts="",
                    deleted_nt_count=codon_call["deleted_nt_count"],
                    phase_at_start=0,
                    phase_at_end=0,
                )
            )

    return rows


def mutations_from_reference_aligned_dna(
    sequence_id: str,
    gene: str,
    sequence: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
) -> list[MutationCall]:
    """Call mutations from nucleotide data already in HXB2 reference columns.

    Gaps in ``sequence`` are meaningful alignment columns and are preserved
    exactly. This helper is for reference-column alignments without insertion
    columns on the HXB2 side. For a full pairwise alignment containing HXB2
    gap columns (insertions in the query), use ``call_mutations_from_hxb2_alignment``.
    """
    reference_nt = hxb2_reference_nt_for_gene(gene)
    query_aligned = "".join(str(sequence).upper().split())
    if len(query_aligned) != len(reference_nt):
        raise ValueError(
            "Reference-aligned DNA must have exactly one query column per HXB2 "
            "reference nucleotide. For alignments with insertion columns, use "
            "call_mutations_from_hxb2_alignment()."
        )
    return mutations_from_hxb2_aligned_dna(
        sequence_id=sequence_id,
        gene=gene,
        query_aligned=query_aligned,
        hxb2_ref_aligned=reference_nt,
        reference_protein=reference_protein,
        reference_system=reference_system,
    )


def call_mutations_from_hxb2_alignment(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    hxb2_ref_aligned: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
) -> list[MutationCall]:
    """Call mutations from a pre-existing query-to-HXB2 nucleotide alignment.

    HXB2 supplies only nucleotide coordinates. Amino-acid mutation calls are
    always defined against the HIVDB Consensus B protein for the requested gene.
    """
    canonical_gene = normalize_gene(gene)
    if reference_system != REFERENCE_SYSTEM_HIVDB_CONSENSUS_B:
        raise ValueError(f"Unsupported mutation reference system: {reference_system!r}.")
    return mutations_from_hxb2_aligned_dna(
        sequence_id=sequence_id,
        gene=canonical_gene,
        query_aligned=query_aligned,
        hxb2_ref_aligned=hxb2_ref_aligned,
        reference_protein=reference_for_gene(canonical_gene),
        reference_system=reference_system,
    )


def mutations_from_hxb2_aligned_dna(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    hxb2_ref_aligned: str,
    reference_protein: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
) -> list[MutationCall]:
    """Call mutations from an existing query-to-HXB2 nucleotide alignment.

    The ungapped reference side must be exactly the HXB2 nucleotide sequence for
    the gene. This prevents an arbitrary or offset alignment from silently being
    interpreted as HXB2 coordinates.
    """
    query_aligned = "".join(str(query_aligned).upper().split())
    hxb2_ref_aligned = "".join(str(hxb2_ref_aligned).upper().split())
    if len(query_aligned) != len(hxb2_ref_aligned):
        raise ValueError("Aligned query and HXB2 reference must have the same length.")

    expected_hxb2 = hxb2_reference_nt_for_gene(gene)
    if hxb2_ref_aligned.replace("-", "") != expected_hxb2:
        raise ValueError(
            "The ungapped reference alignment does not match the HXB2 gene "
            f"reference for {normalize_gene(gene)}."
        )

    codon_calls = codon_calls_from_nucleotide_alignment(
        query_aligned=query_aligned,
        ref_aligned=hxb2_ref_aligned,
        reference_protein=reference_protein,
    )
    protein_tokens = [call["alt_aa"] for call in codon_calls]
    qc_status = "WARN" if any("*" in token or "X" in token for token in protein_tokens) else "PASS"
    rows: list[MutationCall] = []
    for index, codon_call in enumerate(codon_calls[:len(reference_protein)]):
        if codon_call["codon_status"] == "NOT_COVERED":
            continue

        alt_aa = codon_call["alt_aa"]
        ref_aa = reference_protein[index]
        position = index + 1
        if alt_aa == "-":
            mutation_type = "deletion"
        elif "FRAMESHIFT" in codon_call["codon_status"]:
            mutation_type = "frameshift"
        elif "PARTIAL_CODON" in codon_call["codon_status"]:
            mutation_type = "partial_codon"
        else:
            mutation_type = "substitution"
        if alt_aa == "-" or should_report_substitution(ref_aa, normalize_alt_aa(ref_aa, alt_aa)):
            rows.append(
                build_mutation_call(
                    sequence_id=sequence_id,
                    gene=gene,
                    position=position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    mutation_type=mutation_type,
                    qc_status=qc_status,
                    reference_system=reference_system,
                    hxb2_ref_codon=codon_call["hxb2_ref_codon"],
                    query_codon=codon_call["query_codon"],
                    possible_alt_aas=codon_call["possible_alt_aas"],
                    codon_status=codon_call["codon_status"],
                    inserted_nts=codon_call["inserted_nts"],
                    deleted_nt_count=codon_call["deleted_nt_count"],
                    phase_at_start=codon_call["phase_at_start"],
                    phase_at_end=codon_call["phase_at_end"],
                )
            )
    rows.extend(
        insertion_calls_from_nucleotide_alignment(
            sequence_id=sequence_id,
            gene=gene,
            query_aligned=query_aligned,
            ref_aligned=hxb2_ref_aligned,
            qc_status=qc_status,
            reference_system=reference_system,
        )
    )
    return rows


@lru_cache(maxsize=1)
def hxb2_genome_sequence() -> str:
    path = Path(__file__).resolve().parents[1] / "loading" / "reference_genomes" / "HXB2_fasta" / "K03455-B.fasta"
    return "".join(
        line.strip().upper()
        for line in path.read_text().splitlines()
        if line and not line.startswith(">")
    )


def hxb2_reference_nt_for_gene(gene: str) -> str:
    canonical_gene = normalize_gene(gene)
    start, end = HXB2_GENE_COORDINATES[canonical_gene]
    return hxb2_genome_sequence()[start - 1:end]


def normalize_alignment_backend(alignment_backend: str) -> str:
    backend = str(alignment_backend or DEFAULT_NT_ALIGNMENT_BACKEND).strip().lower()
    if backend not in ALIGNMENT_BACKENDS:
        raise ValueError(
            f"Unsupported alignment_backend {alignment_backend!r}; "
            f"expected one of {sorted(ALIGNMENT_BACKENDS)}."
        )
    return backend


def postalign_program() -> str:
    """Return the PostAlign executable path or raise an actionable error.

    The Stanford Sierra stack also supports POSTALIGN_PROGRAM as the path to
    the PostAlign command, so PyHIV uses the same environment variable.
    """
    configured = os.environ.get(POSTALIGN_PROGRAM_ENV, "").strip()
    if configured:
        candidate = Path(configured).expanduser()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
        raise RuntimeError(
            f"{POSTALIGN_PROGRAM_ENV} is set to {configured!r}, but that file "
            "does not exist or is not executable."
        )

    discovered = shutil.which("postalign")
    if discovered:
        return discovered

    raise RuntimeError(
        "alignment_backend='postalign' requires the PostAlign executable. "
        "Install the 'post-align' package and ensure 'postalign' is on PATH, "
        f"or set {POSTALIGN_PROGRAM_ENV} to the executable path."
    )


def refine_hxb2_alignment_with_postalign(
    query_aligned: str,
    ref_aligned: str,
    reference_nt: str,
) -> tuple[str, str]:
    """Refine an existing pairwise HXB2 alignment with PostAlign.

    PostAlign is used only for codon-aware gap placement.  HXB2 remains the
    coordinate reference and mutation calls remain defined against the HIVDB
    Consensus B protein elsewhere in this module.  No position-specific gap
    scores are supplied.
    """
    query_aligned = "".join(str(query_aligned).upper().split())
    ref_aligned = "".join(str(ref_aligned).upper().split())
    reference_nt = "".join(str(reference_nt).upper().split())

    if len(query_aligned) != len(ref_aligned):
        raise ValueError("Pairwise query/reference alignment lengths differ.")
    if ref_aligned.replace("-", "") != reference_nt:
        raise ValueError("The initial reference alignment is not the supplied HXB2 reference.")

    executable = postalign_program()
    with tempfile.TemporaryDirectory(prefix="pyhiv-postalign-") as tmpdir:
        tmp = Path(tmpdir)
        input_path = tmp / "input_alignment.fasta"
        output_path = tmp / "postaligned.fasta"
        input_path.write_text(
            f">HXB2\n{ref_aligned}\n>QUERY\n{query_aligned}\n"
        )

        command = [
            executable,
            "-i", str(input_path),
            "-o", str(output_path),
            "-f", "MSA",
            # PostAlign reads the first FASTA record as the reference header.
            # Supplying a file also works with versions that reject bare headers.
            "-r", str(input_path),
            "-q",
            "codon-alignment",
            "1",
            str(len(reference_nt)),
            "save-fasta",
            "--pairwise",
            "--no-modifiers",
        ]
        try:
            completed = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            stderr = (exc.stderr or "").strip()
            stdout = (exc.stdout or "").strip()
            detail = stderr or stdout or f"exit status {exc.returncode}"
            raise RuntimeError(f"PostAlign codon alignment failed: {detail}") from exc

        if not output_path.exists():
            detail = (completed.stderr or completed.stdout or "").strip()
            raise RuntimeError(
                "PostAlign completed without creating the expected FASTA output"
                + (f": {detail}" if detail else ".")
            )

        with output_path.open() as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        by_id = {record.id: str(record.seq).upper() for record in records}
        try:
            refined_ref = by_id["HXB2"]
            refined_query = by_id["QUERY"]
        except KeyError as exc:
            raise RuntimeError(
                "PostAlign FASTA output did not contain both HXB2 and QUERY sequences."
            ) from exc

    if len(refined_query) != len(refined_ref):
        raise RuntimeError("PostAlign returned unequal query/reference alignment lengths.")
    if refined_ref.replace("-", "") != reference_nt:
        raise RuntimeError("PostAlign changed the ungapped HXB2 reference sequence.")
    if refined_query.replace("-", "") != query_aligned.replace("-", ""):
        raise RuntimeError("PostAlign changed the ungapped query sequence.")
    return refined_query, refined_ref


def align_dna_to_hxb2_columns(
    sequence: str,
    reference_nt: str,
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[str, str]:
    """Align raw nucleotide input to HXB2 reference columns.

    ``builtin`` returns the ordinary pairwise nucleotide alignment.
    ``postalign`` first creates that alignment, then asks Stanford PostAlign to
    refine gap placement using its codon-aware algorithm.
    ``semiglobal_*`` backends are experimental mutation-only aligners with
    terminal gaps unpenalized.
    """
    backend = normalize_alignment_backend(alignment_backend)
    cleaned = "".join(str(sequence).upper().split())
    query = cleaned.replace("-", "")
    if backend == "semiglobal_open20":
        return align_nucleotides_semiglobal(
            query,
            reference_nt,
            mismatch_score=-1,
            gap_open_score=SEMIGLOBAL_GAP_OPEN_SCORE,
            gap_extend_score=SEMIGLOBAL_GAP_EXTEND_SCORE,
        )
    if backend == "semiglobal_mismatch2_open20":
        return align_nucleotides_semiglobal(
            query,
            reference_nt,
            mismatch_score=CODONAWARE_MISMATCH_SCORE,
            gap_open_score=SEMIGLOBAL_GAP_OPEN_SCORE,
            gap_extend_score=SEMIGLOBAL_GAP_EXTEND_SCORE,
        )
    if backend == "semiglobal_mismatch2_open20_trimmed":
        return align_nucleotides_semiglobal_trimmed(
            query,
            reference_nt,
            mismatch_score=CODONAWARE_MISMATCH_SCORE,
            gap_open_score=SEMIGLOBAL_GAP_OPEN_SCORE,
            gap_extend_score=SEMIGLOBAL_GAP_EXTEND_SCORE,
        )
    if backend == "semiglobal_codonaware":
        return align_nucleotides_semiglobal_codonaware(query, reference_nt)
    if backend == "semiglobal_codonaware_local":
        return align_nucleotides_semiglobal_codonaware_local(query, reference_nt)
    if backend == "semiglobal_segmented":
        return align_nucleotides_semiglobal_segmented(query, reference_nt)
    if backend == "semiglobal_anchorblocks":
        return align_nucleotides_semiglobal_anchorblocks(query, reference_nt)
    if backend == "semiglobal_coverage_confidence":
        return align_nucleotides_semiglobal_coverage_confidence(query, reference_nt)

    query_aligned, ref_aligned = align_nucleotides(query, reference_nt)
    if backend == "postalign":
        query_aligned, ref_aligned = refine_hxb2_alignment_with_postalign(
            query_aligned=query_aligned,
            ref_aligned=ref_aligned,
            reference_nt=reference_nt,
        )
    return query_aligned, ref_aligned

def align_nucleotides(query: str, reference: str) -> tuple[str, str]:
    """Globally align raw coding DNA to HXB2 with affine gap penalties.

    A unit-cost edit-distance alignment is a poor default for divergent HIV
    coding regions because a mismatch and a 1-nt indel can have the same cost.
    That can create artificial single-nucleotide gaps, which then look like
    long frameshifts.  Affine gap penalties make isolated substitutions much
    cheaper than opening a gap while still allowing genuine indels when the
    surrounding sequence supports them.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -8
    aligner.extend_gap_score = -1

    alignments = aligner.align(reference, query)
    try:
        alignment = alignments[0]
    except IndexError:
        raise ValueError("Could not align nucleotide sequence to HXB2 reference.")

    ref_aligned, query_aligned = format_pairwise_alignment(alignment, reference, query)
    return query_aligned, ref_aligned


def align_nucleotides_semiglobal(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
) -> tuple[str, str]:
    """Semiglobally align raw coding DNA to HXB2 coordinates.

    Terminal gaps are free so complete pol sequences can be aligned directly to
    gene-sized HXB2 coordinate references. Internal gap scoring is ordinary
    affine scoring without codon awareness.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = SEMIGLOBAL_MATCH_SCORE
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_open_score
    aligner.extend_gap_score = gap_extend_score
    aligner.query_end_gap_score = 0
    aligner.target_end_gap_score = 0

    alignments = aligner.align(reference, query)
    try:
        alignment = alignments[0]
    except IndexError:
        raise ValueError("Could not align nucleotide sequence to HXB2 reference.")

    ref_aligned, query_aligned = format_pairwise_alignment(alignment, reference, query)
    return query_aligned, ref_aligned


def trim_terminal_query_overhangs(
    query_aligned: str,
    ref_aligned: str,
) -> tuple[str, str]:
    """Drop terminal columns where the query extends outside HXB2 coordinates."""
    ref_columns = [index for index, base in enumerate(ref_aligned) if base != "-"]
    if not ref_columns:
        return query_aligned, ref_aligned
    window_start = ref_columns[0]
    window_end = ref_columns[-1] + 1
    return query_aligned[window_start:window_end], ref_aligned[window_start:window_end]


def align_nucleotides_semiglobal_trimmed(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
) -> tuple[str, str]:
    """Semiglobal mismatch2/open20 alignment with terminal query overhangs removed."""
    query_aligned, ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    return trim_terminal_query_overhangs(query_aligned, ref_aligned)


def align_nucleotides_semiglobal_codonaware(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    frameshift_gap_score: int = CODONAWARE_FRAMESHIFT_GAP_SCORE,
) -> tuple[str, str]:
    """Semiglobal nucleotide alignment with codon-aware tie-breaking.

    HXB2 is used only as the nucleotide coordinate reference. The primary score
    is the ordinary semiglobal nucleotide score with terminal gaps unpenalized.
    Codon awareness is used only as a secondary criterion among alignments with
    the same primary score: internal gaps whose length is not a multiple of
    three add a frame cost, but frameshift-length gaps remain allowed.
    """
    query = "".join(str(query).upper().split()).replace("-", "")
    reference = "".join(str(reference).upper().split()).replace("-", "")
    if not reference:
        return query, "-" * len(query)
    if not query:
        return "-" * len(reference), reference

    baseline_query_aligned, baseline_ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )

    def primary_alignment_score(query_aligned: str, ref_aligned: str) -> int:
        score = 0
        index = 0
        alignment_length = len(query_aligned)
        while index < alignment_length:
            query_nt = query_aligned[index]
            ref_nt = ref_aligned[index]
            if query_nt != "-" and ref_nt != "-":
                score += SEMIGLOBAL_MATCH_SCORE if query_nt == ref_nt else mismatch_score
                index += 1
                continue
            gap_start = index
            while index < alignment_length and (
                query_aligned[index] == "-" or ref_aligned[index] == "-"
            ):
                index += 1
            if gap_start == 0 or index == alignment_length:
                continue
            gap_length = index - gap_start
            score += gap_open_score + gap_extend_score * (gap_length - 1)
        return score

    frameshift_cost = max(0, -frameshift_gap_score)

    def frame_alignment_cost(query_aligned: str, ref_aligned: str) -> int:
        cost = 0
        index = 0
        alignment_length = len(query_aligned)
        while index < alignment_length:
            query_nt = query_aligned[index]
            ref_nt = ref_aligned[index]
            if query_nt != "-" and ref_nt != "-":
                index += 1
                continue
            gap_start = index
            gap_side = "query" if query_nt == "-" else "reference"
            while index < alignment_length:
                same_query_gap = gap_side == "query" and query_aligned[index] == "-"
                same_ref_gap = gap_side == "reference" and ref_aligned[index] == "-"
                if not (same_query_gap or same_ref_gap):
                    break
                index += 1
            if gap_start == 0 or index == alignment_length:
                continue
            if (index - gap_start) % 3 != 0:
                cost += frameshift_cost
        return cost

    ref_columns = [
        index for index, base in enumerate(baseline_ref_aligned) if base != "-"
    ]
    if ref_columns:
        window_start = ref_columns[0]
        window_end = ref_columns[-1] + 1
        baseline_window_query = baseline_query_aligned[window_start:window_end]
        baseline_window_ref = baseline_ref_aligned[window_start:window_end]
    else:
        baseline_window_query = baseline_query_aligned
        baseline_window_ref = baseline_ref_aligned

    baseline_primary_score = primary_alignment_score(
        baseline_window_query,
        baseline_window_ref,
    )
    if frame_alignment_cost(baseline_window_query, baseline_window_ref) == 0:
        return baseline_window_query, baseline_window_ref

    query = baseline_window_query.replace("-", "") or query

    # Banded DP with explicit modulo-3 gap states. State 0 is an ordinary
    # aligned-column state. States 1..3 are active gaps in the query
    # (reference bases deleted from the query) with length modulo 3 equal to
    # 1, 2, 0. States 4..6 are active gaps in the reference (query insertion)
    # with the same phase convention. The DP optimizes lexicographically:
    # first the semiglobal nucleotide score, then the cumulative frame cost.
    # Terminal gaps can finish at the matrix edge without a closing frame cost.
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise ImportError("NumPy is required for alignment_backend='semiglobal_codonaware'.") from exc

    ref_len = len(reference)
    query_len = len(query)
    band = max(32, abs(ref_len - query_len) + 12, int(max(ref_len, query_len) * 0.025))
    width = 2 * band + 1
    neg = -10**9
    inf = 10**9
    states = 7
    nonzero_phase_states = {1, 2, 4, 5}

    scores = np.full((ref_len + 1, width, states), neg, dtype=np.int32)
    frame_costs = np.full((ref_len + 1, width, states), inf, dtype=np.int32)
    traces = np.full((ref_len + 1, width, states), 255, dtype=np.uint8)

    def band_index(i: int, j: int) -> int:
        return j - i + band

    def in_band(i: int, j: int) -> bool:
        k = band_index(i, j)
        return 0 <= k < width

    def score_at(i: int, j: int, state: int) -> int:
        if i < 0 or j < 0 or not in_band(i, j):
            return neg
        return int(scores[i, band_index(i, j), state])

    def frame_cost_at(i: int, j: int, state: int) -> int:
        if i < 0 or j < 0 or not in_band(i, j):
            return inf
        return int(frame_costs[i, band_index(i, j), state])

    def better(
        candidate_score: int,
        candidate_frame_cost: int,
        current_score: int,
        current_frame_cost: int,
    ) -> bool:
        return candidate_score > current_score or (
            candidate_score == current_score and candidate_frame_cost < current_frame_cost
        )

    def set_state(
        i: int,
        j: int,
        state: int,
        candidate_score: int,
        candidate_frame_cost: int,
        previous_state: int,
    ) -> None:
        if candidate_score <= neg // 2 or candidate_frame_cost >= inf:
            return
        k = band_index(i, j)
        if better(
            candidate_score,
            candidate_frame_cost,
            int(scores[i, k, state]),
            int(frame_costs[i, k, state]),
        ):
            scores[i, k, state] = candidate_score
            frame_costs[i, k, state] = candidate_frame_cost
            traces[i, k, state] = previous_state

    for i in range(ref_len + 1):
        if in_band(i, 0):
            k = band_index(i, 0)
            scores[i, k, 0] = 0
            frame_costs[i, k, 0] = 0
    for j in range(query_len + 1):
        if in_band(0, j):
            k = band_index(0, j)
            scores[0, k, 0] = 0
            frame_costs[0, k, 0] = 0

    def best_closed(i: int, j: int) -> tuple[int, int, int]:
        best_score = score_at(i, j, 0)
        best_frame_cost = frame_cost_at(i, j, 0)
        best_state = 0
        for state in range(1, states):
            value = score_at(i, j, state)
            if value <= neg // 2:
                continue
            frame_cost = frame_cost_at(i, j, state)
            if state in nonzero_phase_states:
                frame_cost += frameshift_cost
            if better(value, frame_cost, best_score, best_frame_cost):
                best_score = value
                best_frame_cost = frame_cost
                best_state = state
        return best_score, best_frame_cost, best_state

    for i in range(1, ref_len + 1):
        j_start = max(1, i - band)
        j_end = min(query_len, i + band)
        for j in range(j_start, j_end + 1):
            diagonal_score, diagonal_frame_cost, diagonal_state = best_closed(i - 1, j - 1)
            if diagonal_score > neg // 2:
                base_score = SEMIGLOBAL_MATCH_SCORE if reference[i - 1] == query[j - 1] else mismatch_score
                set_state(
                    i,
                    j,
                    0,
                    diagonal_score + base_score,
                    diagonal_frame_cost,
                    diagonal_state,
                )

            set_state(
                i,
                j,
                1,
                score_at(i - 1, j, 0) + gap_open_score,
                frame_cost_at(i - 1, j, 0),
                0,
            )
            set_state(
                i,
                j,
                1,
                score_at(i - 1, j, 3) + gap_extend_score,
                frame_cost_at(i - 1, j, 3),
                3,
            )
            set_state(
                i,
                j,
                2,
                score_at(i - 1, j, 1) + gap_extend_score,
                frame_cost_at(i - 1, j, 1),
                1,
            )
            set_state(
                i,
                j,
                3,
                score_at(i - 1, j, 2) + gap_extend_score,
                frame_cost_at(i - 1, j, 2),
                2,
            )

            set_state(
                i,
                j,
                4,
                score_at(i, j - 1, 0) + gap_open_score,
                frame_cost_at(i, j - 1, 0),
                0,
            )
            set_state(
                i,
                j,
                4,
                score_at(i, j - 1, 6) + gap_extend_score,
                frame_cost_at(i, j - 1, 6),
                6,
            )
            set_state(
                i,
                j,
                5,
                score_at(i, j - 1, 4) + gap_extend_score,
                frame_cost_at(i, j - 1, 4),
                4,
            )
            set_state(
                i,
                j,
                6,
                score_at(i, j - 1, 5) + gap_extend_score,
                frame_cost_at(i, j - 1, 5),
                5,
            )

    best_score = neg
    best_frame_cost = inf
    best_i = ref_len
    best_j = query_len
    best_state = 0
    edge_cells = [(ref_len, j) for j in range(query_len + 1)] + [
        (i, query_len) for i in range(ref_len + 1)
    ]
    for i, j in edge_cells:
        if not in_band(i, j):
            continue
        k = band_index(i, j)
        for state in range(states):
            value = int(scores[i, k, state])
            frame_cost = int(frame_costs[i, k, state])
            if better(value, frame_cost, best_score, best_frame_cost):
                best_score = value
                best_frame_cost = frame_cost
                best_i = i
                best_j = j
                best_state = state

    if best_score <= neg // 2:
        raise ValueError("Could not align nucleotide sequence to HXB2 reference.")

    ref_parts: list[str] = []
    query_parts: list[str] = []

    for i in range(ref_len, best_i, -1):
        ref_parts.append(reference[i - 1])
        query_parts.append("-")
    for j in range(query_len, best_j, -1):
        ref_parts.append("-")
        query_parts.append(query[j - 1])

    i, j, state = best_i, best_j, best_state
    while i > 0 and j > 0:
        k = band_index(i, j)
        previous_state = int(traces[i, k, state])
        if state == 0:
            ref_parts.append(reference[i - 1])
            query_parts.append(query[j - 1])
            i -= 1
            j -= 1
        elif state in {1, 2, 3}:
            ref_parts.append(reference[i - 1])
            query_parts.append("-")
            i -= 1
        else:
            ref_parts.append("-")
            query_parts.append(query[j - 1])
            j -= 1
        state = previous_state if previous_state != 255 else 0

    while i > 0:
        ref_parts.append(reference[i - 1])
        query_parts.append("-")
        i -= 1
    while j > 0:
        ref_parts.append("-")
        query_parts.append(query[j - 1])
        j -= 1

    candidate_query_aligned = "".join(reversed(query_parts))
    candidate_ref_aligned = "".join(reversed(ref_parts))
    candidate_primary_score = primary_alignment_score(
        candidate_query_aligned,
        candidate_ref_aligned,
    )
    if candidate_primary_score < baseline_primary_score:
        return baseline_window_query, baseline_window_ref
    return candidate_query_aligned, candidate_ref_aligned


def _alignment_score_and_frame_cost(
    query_aligned: str,
    ref_aligned: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    frameshift_gap_score: int = CODONAWARE_FRAMESHIFT_GAP_SCORE,
    terminal_gaps_free: bool = True,
) -> tuple[int, int]:
    primary_score = 0
    frame_cost = 0
    frameshift_cost = max(0, -frameshift_gap_score)
    index = 0
    alignment_length = len(query_aligned)
    while index < alignment_length:
        query_nt = query_aligned[index]
        ref_nt = ref_aligned[index]
        if query_nt != "-" and ref_nt != "-":
            primary_score += SEMIGLOBAL_MATCH_SCORE if query_nt == ref_nt else mismatch_score
            index += 1
            continue

        gap_start = index
        gap_side = "query" if query_nt == "-" else "reference"
        while index < alignment_length:
            same_query_gap = gap_side == "query" and query_aligned[index] == "-"
            same_ref_gap = gap_side == "reference" and ref_aligned[index] == "-"
            if not (same_query_gap or same_ref_gap):
                break
            index += 1
        gap_length = index - gap_start
        if terminal_gaps_free and (gap_start == 0 or index == alignment_length):
            continue
        primary_score += gap_open_score + gap_extend_score * (gap_length - 1)
        if gap_length % 3 != 0:
            frame_cost += frameshift_cost
    return primary_score, frame_cost


def _internal_gap_runs(query_aligned: str, ref_aligned: str) -> list[tuple[int, int, str]]:
    runs: list[tuple[int, int, str]] = []
    index = 0
    alignment_length = len(query_aligned)
    while index < alignment_length:
        query_gap = query_aligned[index] == "-"
        ref_gap = ref_aligned[index] == "-"
        if not query_gap and not ref_gap:
            index += 1
            continue
        gap_side = "query" if query_gap else "reference"
        start = index
        while index < alignment_length:
            same_query_gap = gap_side == "query" and query_aligned[index] == "-"
            same_ref_gap = gap_side == "reference" and ref_aligned[index] == "-"
            if not (same_query_gap or same_ref_gap):
                break
            index += 1
        if start != 0 and index != alignment_length:
            runs.append((start, index, gap_side))
    return runs


def _align_nucleotides_global_codonaware_tiebreak(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    frameshift_gap_score: int = CODONAWARE_FRAMESHIFT_GAP_SCORE,
) -> tuple[str, str]:
    """Globally align a short local window, using frame cost only as tie-break."""
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise ImportError("NumPy is required for local codon-aware alignment.") from exc

    ref_len = len(reference)
    query_len = len(query)
    if not reference:
        return query, "-" * len(query)
    if not query:
        return "-" * len(reference), reference

    band = max(8, abs(ref_len - query_len) + 6, max(ref_len, query_len))
    width = 2 * band + 1
    neg = -10**9
    inf = 10**9
    states = 7
    nonzero_phase_states = {1, 2, 4, 5}
    frameshift_cost = max(0, -frameshift_gap_score)

    scores = np.full((ref_len + 1, width, states), neg, dtype=np.int32)
    frame_costs = np.full((ref_len + 1, width, states), inf, dtype=np.int32)
    traces = np.full((ref_len + 1, width, states), 255, dtype=np.uint8)

    def band_index(i: int, j: int) -> int:
        return j - i + band

    def in_band(i: int, j: int) -> bool:
        k = band_index(i, j)
        return 0 <= k < width

    def score_at(i: int, j: int, state: int) -> int:
        if i < 0 or j < 0 or not in_band(i, j):
            return neg
        return int(scores[i, band_index(i, j), state])

    def frame_cost_at(i: int, j: int, state: int) -> int:
        if i < 0 or j < 0 or not in_band(i, j):
            return inf
        return int(frame_costs[i, band_index(i, j), state])

    def better(
        candidate_score: int,
        candidate_frame_cost: int,
        current_score: int,
        current_frame_cost: int,
    ) -> bool:
        return candidate_score > current_score or (
            candidate_score == current_score and candidate_frame_cost < current_frame_cost
        )

    def set_state(
        i: int,
        j: int,
        state: int,
        candidate_score: int,
        candidate_frame_cost: int,
        previous_state: int,
    ) -> None:
        if candidate_score <= neg // 2 or candidate_frame_cost >= inf:
            return
        k = band_index(i, j)
        if better(
            candidate_score,
            candidate_frame_cost,
            int(scores[i, k, state]),
            int(frame_costs[i, k, state]),
        ):
            scores[i, k, state] = candidate_score
            frame_costs[i, k, state] = candidate_frame_cost
            traces[i, k, state] = previous_state

    def best_closed(i: int, j: int) -> tuple[int, int, int]:
        best_score = score_at(i, j, 0)
        best_frame_cost = frame_cost_at(i, j, 0)
        best_state = 0
        for state in range(1, states):
            value = score_at(i, j, state)
            if value <= neg // 2:
                continue
            frame_cost = frame_cost_at(i, j, state)
            if state in nonzero_phase_states:
                frame_cost += frameshift_cost
            if better(value, frame_cost, best_score, best_frame_cost):
                best_score = value
                best_frame_cost = frame_cost
                best_state = state
        return best_score, best_frame_cost, best_state

    if in_band(0, 0):
        scores[0, band_index(0, 0), 0] = 0
        frame_costs[0, band_index(0, 0), 0] = 0

    for i in range(ref_len + 1):
        j_start = max(0, i - band)
        j_end = min(query_len, i + band)
        for j in range(j_start, j_end + 1):
            if i == 0 and j == 0:
                continue
            if i > 0 and j > 0:
                diagonal_score, diagonal_frame_cost, diagonal_state = best_closed(i - 1, j - 1)
                if diagonal_score > neg // 2:
                    base_score = SEMIGLOBAL_MATCH_SCORE if reference[i - 1] == query[j - 1] else mismatch_score
                    set_state(i, j, 0, diagonal_score + base_score, diagonal_frame_cost, diagonal_state)
            if i > 0:
                set_state(
                    i, j, 1,
                    score_at(i - 1, j, 0) + gap_open_score,
                    frame_cost_at(i - 1, j, 0),
                    0,
                )
                set_state(
                    i, j, 1,
                    score_at(i - 1, j, 3) + gap_extend_score,
                    frame_cost_at(i - 1, j, 3),
                    3,
                )
                set_state(
                    i, j, 2,
                    score_at(i - 1, j, 1) + gap_extend_score,
                    frame_cost_at(i - 1, j, 1),
                    1,
                )
                set_state(
                    i, j, 3,
                    score_at(i - 1, j, 2) + gap_extend_score,
                    frame_cost_at(i - 1, j, 2),
                    2,
                )
            if j > 0:
                set_state(
                    i, j, 4,
                    score_at(i, j - 1, 0) + gap_open_score,
                    frame_cost_at(i, j - 1, 0),
                    0,
                )
                set_state(
                    i, j, 4,
                    score_at(i, j - 1, 6) + gap_extend_score,
                    frame_cost_at(i, j - 1, 6),
                    6,
                )
                set_state(
                    i, j, 5,
                    score_at(i, j - 1, 4) + gap_extend_score,
                    frame_cost_at(i, j - 1, 4),
                    4,
                )
                set_state(
                    i, j, 6,
                    score_at(i, j - 1, 5) + gap_extend_score,
                    frame_cost_at(i, j - 1, 5),
                    5,
                )

    best_score, best_frame_cost, best_state = best_closed(ref_len, query_len)
    if best_score <= neg // 2:
        raise ValueError("Could not locally align nucleotide window to HXB2 reference.")

    ref_parts: list[str] = []
    query_parts: list[str] = []
    i, j, state = ref_len, query_len, best_state
    while i > 0 or j > 0:
        k = band_index(i, j)
        previous_state = int(traces[i, k, state]) if in_band(i, j) else 255
        if state == 0:
            ref_parts.append(reference[i - 1])
            query_parts.append(query[j - 1])
            i -= 1
            j -= 1
        elif state in {1, 2, 3}:
            ref_parts.append(reference[i - 1])
            query_parts.append("-")
            i -= 1
        else:
            ref_parts.append("-")
            query_parts.append(query[j - 1])
            j -= 1
        state = previous_state if previous_state != 255 else 0
    return "".join(reversed(query_parts)), "".join(reversed(ref_parts))


def align_nucleotides_semiglobal_codonaware_local(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    frameshift_gap_score: int = CODONAWARE_FRAMESHIFT_GAP_SCORE,
    flank_columns: int = 36,
) -> tuple[str, str]:
    """Semiglobal baseline with local codon-aware refinements near 1/2 nt gaps."""
    query = "".join(str(query).upper().split()).replace("-", "")
    reference = "".join(str(reference).upper().split()).replace("-", "")
    if not reference:
        return query, "-" * len(query)
    if not query:
        return "-" * len(reference), reference

    baseline_query_aligned, baseline_ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    ref_columns = [index for index, base in enumerate(baseline_ref_aligned) if base != "-"]
    if ref_columns:
        window_start = ref_columns[0]
        window_end = ref_columns[-1] + 1
        query_aligned = baseline_query_aligned[window_start:window_end]
        ref_aligned = baseline_ref_aligned[window_start:window_end]
    else:
        query_aligned = baseline_query_aligned
        ref_aligned = baseline_ref_aligned

    baseline_primary, _ = _alignment_score_and_frame_cost(
        query_aligned,
        ref_aligned,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
        frameshift_gap_score=frameshift_gap_score,
        terminal_gaps_free=True,
    )

    max_passes = 6
    for _ in range(max_passes):
        changed = False
        for gap_start, gap_end, _gap_side in _internal_gap_runs(query_aligned, ref_aligned):
            gap_length = gap_end - gap_start
            if gap_length not in {1, 2}:
                continue
            segment_start = max(0, gap_start - flank_columns)
            segment_end = min(len(query_aligned), gap_end + flank_columns)
            original_query_segment = query_aligned[segment_start:segment_end]
            original_ref_segment = ref_aligned[segment_start:segment_end]
            if not original_query_segment.replace("-", "") or not original_ref_segment.replace("-", ""):
                continue
            original_primary, original_frame = _alignment_score_and_frame_cost(
                original_query_segment,
                original_ref_segment,
                mismatch_score=mismatch_score,
                gap_open_score=gap_open_score,
                gap_extend_score=gap_extend_score,
                frameshift_gap_score=frameshift_gap_score,
                terminal_gaps_free=False,
            )
            if original_frame == 0:
                continue
            candidate_query_segment, candidate_ref_segment = _align_nucleotides_global_codonaware_tiebreak(
                original_query_segment.replace("-", ""),
                original_ref_segment.replace("-", ""),
                mismatch_score=mismatch_score,
                gap_open_score=gap_open_score,
                gap_extend_score=gap_extend_score,
                frameshift_gap_score=frameshift_gap_score,
            )
            if (
                candidate_query_segment[0] == "-"
                or candidate_ref_segment[0] == "-"
                or candidate_query_segment[-1] == "-"
                or candidate_ref_segment[-1] == "-"
            ):
                continue
            candidate_primary, candidate_frame = _alignment_score_and_frame_cost(
                candidate_query_segment,
                candidate_ref_segment,
                mismatch_score=mismatch_score,
                gap_open_score=gap_open_score,
                gap_extend_score=gap_extend_score,
                frameshift_gap_score=frameshift_gap_score,
                terminal_gaps_free=False,
            )
            if candidate_primary < original_primary or candidate_frame >= original_frame:
                continue
            proposed_query = query_aligned[:segment_start] + candidate_query_segment + query_aligned[segment_end:]
            proposed_ref = ref_aligned[:segment_start] + candidate_ref_segment + ref_aligned[segment_end:]
            proposed_primary, _ = _alignment_score_and_frame_cost(
                proposed_query,
                proposed_ref,
                mismatch_score=mismatch_score,
                gap_open_score=gap_open_score,
                gap_extend_score=gap_extend_score,
                frameshift_gap_score=frameshift_gap_score,
                terminal_gaps_free=True,
            )
            if proposed_primary < baseline_primary:
                continue
            query_aligned = proposed_query
            ref_aligned = proposed_ref
            changed = True
            break
        if not changed:
            break

    return query_aligned, ref_aligned



def _align_nucleotides_global_affine(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
) -> tuple[str, str]:
    if not reference:
        return query, "-" * len(query)
    if not query:
        return "-" * len(reference), reference
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = SEMIGLOBAL_MATCH_SCORE
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = gap_open_score
    aligner.extend_gap_score = gap_extend_score
    alignments = aligner.align(reference, query)
    try:
        alignment = alignments[0]
    except IndexError:
        raise ValueError("Could not align nucleotide segment to HXB2 reference.")
    ref_aligned, query_aligned = format_pairwise_alignment(alignment, reference, query)
    return query_aligned, ref_aligned


def _unique_kmer_positions(sequence: str, k: int) -> dict[str, int]:
    counts: dict[str, int] = {}
    positions: dict[str, int] = {}
    for index in range(0, len(sequence) - k + 1):
        kmer = sequence[index:index + k]
        if set(kmer) - {"A", "C", "G", "T"}:
            continue
        counts[kmer] = counts.get(kmer, 0) + 1
        positions[kmer] = index
    return {kmer: positions[kmer] for kmer, count in counts.items() if count == 1}


def _monotonic_exact_anchors(
    query: str,
    reference: str,
    k: int = 18,
    min_anchors: int = 3,
) -> list[tuple[int, int, int]]:
    query_positions = _unique_kmer_positions(query, k)
    reference_positions = _unique_kmer_positions(reference, k)
    pairs = sorted(
        (ref_pos, query_positions[kmer])
        for kmer, ref_pos in reference_positions.items()
        if kmer in query_positions
    )
    if len(pairs) < min_anchors:
        return []

    # Longest increasing subsequence on query coordinates after sorting by HXB2.
    tails: list[int] = []
    tails_indices: list[int] = []
    predecessors = [-1] * len(pairs)
    for index, (_ref_pos, query_pos) in enumerate(pairs):
        left, right = 0, len(tails)
        while left < right:
            middle = (left + right) // 2
            if tails[middle] < query_pos:
                left = middle + 1
            else:
                right = middle
        if left:
            predecessors[index] = tails_indices[left - 1]
        if left == len(tails):
            tails.append(query_pos)
            tails_indices.append(index)
        else:
            tails[left] = query_pos
            tails_indices[left] = index

    if not tails_indices:
        return []
    chain_indices: list[int] = []
    cursor = tails_indices[-1]
    while cursor != -1:
        chain_indices.append(cursor)
        cursor = predecessors[cursor]
    chain = [pairs[index] for index in reversed(chain_indices)]
    if len(chain) < min_anchors:
        return []

    merged: list[tuple[int, int, int]] = []
    for ref_pos, query_pos in chain:
        if not merged:
            merged.append((ref_pos, query_pos, k))
            continue
        last_ref, last_query, last_length = merged[-1]
        ref_delta = ref_pos - last_ref
        query_delta = query_pos - last_query
        if ref_delta == query_delta and 0 < ref_delta <= last_length:
            merged[-1] = (last_ref, last_query, max(last_length, ref_delta + k))
            continue
        if ref_pos >= last_ref + last_length and query_pos >= last_query + last_length:
            merged.append((ref_pos, query_pos, k))

    anchors = [anchor for anchor in merged if anchor[2] >= k]
    return anchors if len(anchors) >= min_anchors else []


def align_nucleotides_semiglobal_segmented(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    anchor_k: int = 18,
    min_anchors: int = 3,
) -> tuple[str, str]:
    """Semiglobal baseline split into independently aligned anchored segments.

    The baseline is `semiglobal_mismatch2_open20` with terminal query overhangs
    trimmed. Exact unique nucleotide k-mers provide monotonic anchors between
    the trimmed query and HXB2. Anchors are fixed, and the intervening segments
    are aligned independently with the same affine nucleotide scoring.
    """
    query = "".join(str(query).upper().split()).replace("-", "")
    reference = "".join(str(reference).upper().split()).replace("-", "")
    if not reference:
        return query, "-" * len(query)
    if not query:
        return "-" * len(reference), reference

    baseline_query_aligned, baseline_ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    ref_columns = [index for index, base in enumerate(baseline_ref_aligned) if base != "-"]
    if ref_columns:
        window_start = ref_columns[0]
        window_end = ref_columns[-1] + 1
        baseline_query_aligned = baseline_query_aligned[window_start:window_end]
        baseline_ref_aligned = baseline_ref_aligned[window_start:window_end]
    trimmed_query = baseline_query_aligned.replace("-", "")

    anchors = _monotonic_exact_anchors(
        trimmed_query,
        reference,
        k=anchor_k,
        min_anchors=min_anchors,
    )
    if not anchors:
        return baseline_query_aligned, baseline_ref_aligned

    query_parts: list[str] = []
    ref_parts: list[str] = []
    query_cursor = 0
    ref_cursor = 0
    for ref_start, query_start, anchor_length in anchors:
        if ref_start < ref_cursor or query_start < query_cursor:
            continue
        segment_query = trimmed_query[query_cursor:query_start]
        segment_ref = reference[ref_cursor:ref_start]
        segment_query_aligned, segment_ref_aligned = _align_nucleotides_global_affine(
            segment_query,
            segment_ref,
            mismatch_score=mismatch_score,
            gap_open_score=gap_open_score,
            gap_extend_score=gap_extend_score,
        )
        query_parts.append(segment_query_aligned)
        ref_parts.append(segment_ref_aligned)
        query_anchor = trimmed_query[query_start:query_start + anchor_length]
        ref_anchor = reference[ref_start:ref_start + anchor_length]
        if query_anchor != ref_anchor or len(query_anchor) != anchor_length:
            return baseline_query_aligned, baseline_ref_aligned
        query_parts.append(query_anchor)
        ref_parts.append(ref_anchor)
        query_cursor = query_start + anchor_length
        ref_cursor = ref_start + anchor_length

    tail_query_aligned, tail_ref_aligned = _align_nucleotides_global_affine(
        trimmed_query[query_cursor:],
        reference[ref_cursor:],
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    query_parts.append(tail_query_aligned)
    ref_parts.append(tail_ref_aligned)

    segmented_query_aligned = "".join(query_parts)
    segmented_ref_aligned = "".join(ref_parts)
    if segmented_query_aligned.replace("-", "") != trimmed_query:
        return baseline_query_aligned, baseline_ref_aligned
    if segmented_ref_aligned.replace("-", "") != reference:
        return baseline_query_aligned, baseline_ref_aligned
    return segmented_query_aligned, segmented_ref_aligned



def _exact_run_anchors_from_alignment(
    query_aligned: str,
    ref_aligned: str,
    min_length: int = 18,
) -> list[tuple[int, int, int]]:
    anchors: list[tuple[int, int, int]] = []
    ref_pos = 0
    query_pos = 0
    run_ref_start = 0
    run_query_start = 0
    run_length = 0

    def flush_run() -> None:
        nonlocal run_length
        if run_length >= min_length:
            anchors.append((run_ref_start, run_query_start, run_length))
        run_length = 0

    for query_nt, ref_nt in zip(query_aligned, ref_aligned):
        current_ref_pos = ref_pos
        current_query_pos = query_pos
        if ref_nt != "-":
            ref_pos += 1
        if query_nt != "-":
            query_pos += 1

        if (
            query_nt == ref_nt
            and query_nt in {"A", "C", "G", "T"}
        ):
            if run_length == 0:
                run_ref_start = current_ref_pos
                run_query_start = current_query_pos
            run_length += 1
        else:
            flush_run()
    flush_run()
    return anchors



def align_nucleotides_semiglobal_anchorblocks(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    anchor_k: int = 18,
    min_anchors: int = 2,
    min_unsupported_ref_nt: int = 45,
) -> tuple[str, str]:
    """Anchor-block alignment that does not force unsupported deletions.

    The baseline is `semiglobal_mismatch2_open20` with terminal trimming. Exact
    unique nucleotide k-mers define monotonic anchors. Intervals between anchors
    are aligned independently when they contain query sequence. A long reference
    interval with no intervening query nucleotides is represented as uninformative
    `N` coverage rather than as a forced internal deletion. This uses only
    nucleotide evidence and never consults Sierra, DRMs, subtype, Consensus B, or
    position-specific rules.
    """
    query = "".join(str(query).upper().split()).replace("-", "")
    reference = "".join(str(reference).upper().split()).replace("-", "")
    if not reference:
        return query, "-" * len(query)
    if not query:
        return "-" * len(reference), reference

    baseline_query_aligned, baseline_ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    ref_columns = [index for index, base in enumerate(baseline_ref_aligned) if base != "-"]
    if ref_columns:
        window_start = ref_columns[0]
        window_end = ref_columns[-1] + 1
        baseline_query_aligned = baseline_query_aligned[window_start:window_end]
        baseline_ref_aligned = baseline_ref_aligned[window_start:window_end]
    trimmed_query = baseline_query_aligned.replace("-", "")

    anchors = _exact_run_anchors_from_alignment(
        baseline_query_aligned,
        baseline_ref_aligned,
        min_length=anchor_k,
    )
    if len(anchors) < min_anchors:
        anchors = _monotonic_exact_anchors(
            trimmed_query,
            reference,
            k=anchor_k,
            min_anchors=min_anchors,
        )
    if len(anchors) < min_anchors:
        return baseline_query_aligned, baseline_ref_aligned

    query_parts: list[str] = []
    ref_parts: list[str] = []
    query_cursor = 0
    ref_cursor = 0
    used_unsupported_interval = False
    for ref_start, query_start, anchor_length in anchors:
        if ref_start < ref_cursor or query_start < query_cursor:
            continue
        segment_query = trimmed_query[query_cursor:query_start]
        segment_ref = reference[ref_cursor:ref_start]
        if (
            not segment_query
            and len(segment_ref) >= min_unsupported_ref_nt
            and query_cursor > 0
            and ref_cursor > 0
        ):
            query_parts.append("N" * len(segment_ref))
            ref_parts.append(segment_ref)
            used_unsupported_interval = True
        else:
            segment_query_aligned, segment_ref_aligned = _align_nucleotides_global_affine(
                segment_query,
                segment_ref,
                mismatch_score=mismatch_score,
                gap_open_score=gap_open_score,
                gap_extend_score=gap_extend_score,
            )
            query_parts.append(segment_query_aligned)
            ref_parts.append(segment_ref_aligned)
        query_anchor = trimmed_query[query_start:query_start + anchor_length]
        ref_anchor = reference[ref_start:ref_start + anchor_length]
        if query_anchor != ref_anchor or len(query_anchor) != anchor_length:
            return baseline_query_aligned, baseline_ref_aligned
        query_parts.append(query_anchor)
        ref_parts.append(ref_anchor)
        query_cursor = query_start + anchor_length
        ref_cursor = ref_start + anchor_length

    tail_query_aligned, tail_ref_aligned = _align_nucleotides_global_affine(
        trimmed_query[query_cursor:],
        reference[ref_cursor:],
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    query_parts.append(tail_query_aligned)
    ref_parts.append(tail_ref_aligned)

    block_query_aligned = "".join(query_parts)
    block_ref_aligned = "".join(ref_parts)
    if not used_unsupported_interval:
        return baseline_query_aligned, baseline_ref_aligned
    if block_query_aligned.replace("-", "").replace("N", "") != trimmed_query:
        return baseline_query_aligned, baseline_ref_aligned
    if block_ref_aligned.replace("-", "") != reference:
        return baseline_query_aligned, baseline_ref_aligned
    return block_query_aligned, block_ref_aligned



def align_nucleotides_semiglobal_coverage_confidence(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
) -> tuple[str, str]:
    """Baseline semiglobal alignment plus terminal trimming.

    This backend intentionally does not realign or move internal gaps. It uses
    `semiglobal_mismatch2_open20` as the coordinate alignment, trims terminal
    query overhang columns, and applies coverage-confidence masking later at
    codon interpretation time.
    """
    query_aligned, ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    return trim_terminal_query_overhangs(query_aligned, ref_aligned)


def anchorless_codon_runs(anchor_supported: list[bool]) -> list[tuple[int, int, int]]:
    runs: list[tuple[int, int, int]] = []
    start: int | None = None
    for index, supported in enumerate(anchor_supported):
        if not supported and start is None:
            start = index
        elif supported and start is not None:
            runs.append((start, index - 1, index - start))
            start = None
    if start is not None:
        runs.append((start, len(anchor_supported) - 1, len(anchor_supported) - start))
    return runs


def exact_anchor_supported_codons(
    query_aligned: str,
    ref_aligned: str,
    reference_protein: str,
    min_anchor_nt: int = 12,
) -> list[bool]:
    supported = [False] * len(reference_protein)
    ref_nt_position = 0
    run_ref_start = 0
    run_length = 0

    def flush_run() -> None:
        nonlocal run_length
        if run_length >= min_anchor_nt:
            for nt_position in range(run_ref_start, run_ref_start + run_length):
                codon_index = nt_position // 3
                if 0 <= codon_index < len(supported):
                    supported[codon_index] = True
        run_length = 0

    for query_nt, ref_nt in zip(query_aligned, ref_aligned):
        current_ref_position = ref_nt_position
        if ref_nt != "-":
            ref_nt_position += 1
        if query_nt == ref_nt and query_nt in {"A", "C", "G", "T"}:
            if run_length == 0:
                run_ref_start = current_ref_position
            run_length += 1
        else:
            flush_run()
    flush_run()
    return supported


def coverage_confidence_not_covered_codons(
    query_aligned: str,
    ref_aligned: str,
    reference_protein: str,
    min_anchor_nt: int = 12,
    min_anchorless_fraction: float = 0.40,
    min_anchorless_codons_floor: int = 30,
) -> set[int]:
    """Return low-confidence codons using only nucleotide anchor support.

    A codon is low confidence when it belongs to a long HXB2 interval with no
    exact nucleotide anchor support. The minimum interval length scales with the
    gene length, so the rule is generic across PR, RT and IN rather than tied to
    one benchmark-derived coordinate.
    """
    anchor_supported = exact_anchor_supported_codons(
        query_aligned,
        ref_aligned,
        reference_protein,
        min_anchor_nt=min_anchor_nt,
    )
    minimum_run = max(
        min_anchorless_codons_floor,
        int(round(len(reference_protein) * min_anchorless_fraction)),
    )
    not_covered: set[int] = set()
    for start, end, length in anchorless_codon_runs(anchor_supported):
        if length >= minimum_run:
            not_covered.update(range(start, end + 1))
    return not_covered



def semiglobal_hxb2_query_window(
    query: str,
    reference: str,
    mismatch_score: int = CODONAWARE_MISMATCH_SCORE,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
) -> str:
    """Return the query span aligned to HXB2 by a semiglobal coordinate pass."""
    query_aligned, ref_aligned = align_nucleotides_semiglobal(
        query,
        reference,
        mismatch_score=mismatch_score,
        gap_open_score=gap_open_score,
        gap_extend_score=gap_extend_score,
    )
    ref_columns = [index for index, base in enumerate(ref_aligned) if base != "-"]
    if not ref_columns:
        return query
    start = ref_columns[0]
    end = ref_columns[-1] + 1
    window = "".join(base for base in query_aligned[start:end] if base != "-")
    return window or query


def internal_gap_score(
    start: int,
    length: int,
    sequence_length: int,
    gap_open_score: int = SEMIGLOBAL_GAP_OPEN_SCORE,
    gap_extend_score: int = SEMIGLOBAL_GAP_EXTEND_SCORE,
    frameshift_gap_score: int = CODONAWARE_FRAMESHIFT_GAP_SCORE,
) -> int:
    if length <= 0:
        return 0
    if start == 0 or start == sequence_length:
        return 0
    score = gap_open_score + gap_extend_score * (length - 1)
    if length % 3 != 0:
        score += frameshift_gap_score
    return score


def codon_calls_from_nucleotide_alignment(
    query_aligned: str,
    ref_aligned: str,
    reference_protein: str,
) -> list[dict[str, object]]:
    """Map an HXB2 nucleotide alignment into codon-level observations.

    Coverage, ambiguity and frame state are inferred generically from alignment
    columns. Terminal missing sequence does not alter frame; only internal indels
    inside the covered span do.
    """
    if len(query_aligned) != len(ref_aligned):
        raise ValueError("Aligned query and reference must have the same length.")

    codons = [
        {
            "hxb2_ref_codon": [],
            "query_codon": [],
            "inserted_nts": [],
            "internal_deleted_nt_count": 0,
            "missing_nt_count": 0,
            "phase_values": [],
            "phase_at_start": 0,
            "phase_at_end": 0,
        }
        for _ in range(len(reference_protein))
    ]

    # Find the informative covered span on HXB2. N-only terminal padding is not
    # considered informative, but internal N codons remain part of coverage.
    informative_ref_positions: list[int] = []
    ref_nt_position = 0
    for ref_nt, query_nt in zip(ref_aligned, query_aligned):
        if ref_nt == "-":
            continue
        if query_nt not in {"-", "N"}:
            informative_ref_positions.append(ref_nt_position)
        ref_nt_position += 1

    if informative_ref_positions:
        first_covered_nt = min(informative_ref_positions)
        last_covered_nt = max(informative_ref_positions)
    else:
        first_covered_nt = len(reference_protein) * 3
        last_covered_nt = -1

    phase = 0
    ref_nt_position = 0
    previous_codon_index = 0

    for ref_nt, query_nt in zip(ref_aligned, query_aligned):
        if ref_nt == "-":
            # An insertion is anchored to the preceding HXB2 codon. Only an
            # insertion within the covered region is allowed to change frame.
            if query_nt != "-" and codons:
                insertion_index = min(max(previous_codon_index, 0), len(codons) - 1)
                within_covered_span = first_covered_nt <= ref_nt_position <= last_covered_nt + 1
                if within_covered_span:
                    codons[insertion_index]["inserted_nts"].append(query_nt)
                    phase = (phase + 1) % 3
                    codons[insertion_index]["phase_values"].append(phase)
                    codons[insertion_index]["phase_at_end"] = phase
            continue

        codon_index = ref_nt_position // 3
        if codon_index < len(codons):
            codon = codons[codon_index]
            if not codon["hxb2_ref_codon"]:
                codon["phase_at_start"] = phase

            codon["hxb2_ref_codon"].append(ref_nt)
            codon["query_codon"].append(query_nt)

            within_covered_span = first_covered_nt <= ref_nt_position <= last_covered_nt
            if query_nt == "-":
                if within_covered_span:
                    codon["internal_deleted_nt_count"] += 1
                    phase = (phase - 1) % 3
                else:
                    codon["missing_nt_count"] += 1

            codon["phase_values"].append(phase)
            codon["phase_at_end"] = phase
            previous_codon_index = codon_index

        ref_nt_position += 1

    first_covered_codon = first_covered_nt // 3 if last_covered_nt >= 0 else len(codons)
    last_covered_codon = last_covered_nt // 3 if last_covered_nt >= 0 else -1

    calls: list[dict[str, object]] = []
    for index, codon in enumerate(codons):
        ref_codon = "".join(codon["hxb2_ref_codon"])
        query_codon = "".join(codon["query_codon"])
        inserted_nts = "".join(codon["inserted_nts"])
        deleted_nt_count = int(codon["internal_deleted_nt_count"])
        missing_nt_count = int(codon["missing_nt_count"])
        phase_at_start = int(codon["phase_at_start"])
        phase_at_end = int(codon["phase_at_end"])

        possible_alt_aas = possible_amino_acids_for_query_codon(query_codon)
        status = codon_status(
            index=index,
            query_codon=query_codon,
            inserted_nts=inserted_nts,
            deleted_nt_count=deleted_nt_count,
            missing_nt_count=missing_nt_count,
            phase_at_start=phase_at_start,
            phase_at_end=phase_at_end,
            first_covered_codon=first_covered_codon,
            last_covered_codon=last_covered_codon,
        )
        alt_aa = amino_acid_token_from_possible(possible_alt_aas, status)
        calls.append(
            {
                "hxb2_ref_codon": ref_codon,
                "query_codon": query_codon,
                "possible_alt_aas": "".join(sorted(possible_alt_aas)),
                "inserted_nts": inserted_nts,
                "deleted_nt_count": deleted_nt_count,
                "phase_at_start": phase_at_start,
                "phase_at_end": phase_at_end,
                "codon_status": status,
                "alt_aa": alt_aa,
            }
        )
    return calls


def codon_calls_from_raw_nucleotide_alignment(
    query_aligned: str,
    ref_aligned: str,
    reference_protein: str,
    force_not_covered_codons: set[int] | None = None,
) -> list[dict[str, object]]:
    """Create conservative codon calls from a raw-query HXB2 alignment.

    Each HXB2 codon is interpreted independently. Partial nucleotide gaps are
    local uncertainty and are never allowed to turn every downstream codon into
    a frameshift call. This mirrors the stable pre-frameshift-propagation caller
    while retaining the richer QC/debug fields.
    """
    if len(query_aligned) != len(ref_aligned):
        raise ValueError("Aligned query and reference must have the same length.")

    codons = [
        {"hxb2_ref_codon": [], "query_codon": []}
        for _ in range(len(reference_protein))
    ]
    ref_nt_position = 0
    for ref_nt, query_nt in zip(ref_aligned, query_aligned):
        if ref_nt == "-":
            continue
        codon_index = ref_nt_position // 3
        if codon_index < len(codons):
            codons[codon_index]["hxb2_ref_codon"].append(ref_nt)
            codons[codon_index]["query_codon"].append(query_nt)
        ref_nt_position += 1

    rendered = ["".join(codon["query_codon"]) for codon in codons]
    informative = [
        i for i, codon in enumerate(rendered)
        if any(base not in {"-", "N"} for base in codon)
    ]
    first_covered = min(informative) if informative else len(codons)
    last_covered = max(informative) if informative else -1
    forced_not_covered = force_not_covered_codons or set()

    calls: list[dict[str, object]] = []
    for index, codon in enumerate(codons):
        hxb2_ref_codon = "".join(codon["hxb2_ref_codon"])
        query_codon = "".join(codon["query_codon"])
        deleted_nt_count = query_codon.count("-")
        possible_alt_aas = possible_amino_acids_for_query_codon(query_codon)

        if index in forced_not_covered or index < first_covered or index > last_covered:
            status = "NOT_COVERED"
        elif query_codon == "---":
            status = "DELETION"
        else:
            statuses: list[str] = []
            if "-" in query_codon:
                statuses.append("PARTIAL_CODON")
            if any(base not in {"A", "C", "G", "T", "-"} for base in query_codon):
                statuses.append("AMBIGUOUS")
            # A codon that reaches a terminal N-only tail is partially covered,
            # not a biological mixture. Keep its possibilities for debugging but
            # render the compact amino-acid call as X.
            if "N" in query_codon and index == last_covered:
                trailing = rendered[index + 1:]
                if all((not item) or set(item) <= {"N", "-"} for item in trailing):
                    statuses.append("PARTIAL_COVERAGE")
            status = "|".join(statuses) if statuses else "OBSERVED"

        alt_aa = amino_acid_token_from_possible(possible_alt_aas, status)
        calls.append(
            {
                "hxb2_ref_codon": hxb2_ref_codon,
                "query_codon": query_codon,
                "possible_alt_aas": "".join(sorted(possible_alt_aas)),
                "inserted_nts": "",
                "deleted_nt_count": deleted_nt_count,
                "phase_at_start": 0,
                "phase_at_end": 0,
                "codon_status": status,
                "alt_aa": alt_aa,
            }
        )
    return calls


def possible_amino_acids_for_query_codon(codon: str) -> set[str]:
    """Return every amino acid compatible with an IUPAC query codon.

    Gaps are expanded as unknown nucleotides only for diagnostic purposes. The
    codon status determines whether those possibilities are biologically usable.
    """
    if len(codon) != 3:
        return set()
    if codon == "---":
        return {"-"}
    if any(base not in IUPAC_DNA and base != "-" for base in codon):
        return set()

    possibilities = [""]
    for base in codon:
        choices = IUPAC_DNA["N"] if base == "-" else IUPAC_DNA[base]
        possibilities = [
            prefix + resolved
            for prefix in possibilities
            for resolved in choices
        ]
    return {str(Seq(possible).translate()) for possible in possibilities}


def amino_acid_token_from_possible(possible_alt_aas: set[str], codon_status: str = "OBSERVED") -> str:
    if possible_alt_aas == {"-"}:
        return "-"
    if any(status in codon_status for status in ("NOT_COVERED", "PARTIAL_CODON", "PARTIAL_COVERAGE", "FRAMESHIFT")):
        return "X"
    if not possible_alt_aas:
        return "X"
    if len(possible_alt_aas) == 1:
        return next(iter(possible_alt_aas))
    if len(possible_alt_aas) <= 4:
        return "".join(sorted(possible_alt_aas))
    return "X"


def codon_status(
    index: int,
    query_codon: str,
    inserted_nts: str,
    deleted_nt_count: int,
    missing_nt_count: int,
    phase_at_start: int,
    phase_at_end: int,
    first_covered_codon: int,
    last_covered_codon: int,
) -> str:
    """Classify a codon using position-independent alignment rules."""
    if index < first_covered_codon or index > last_covered_codon:
        return "NOT_COVERED"

    statuses: list[str] = []
    query_nt_count = sum(1 for nt in query_codon if nt != "-")

    # Missing terminal bases mean partial coverage, not a biological deletion.
    if missing_nt_count:
        statuses.append("PARTIAL_CODON")

    if query_nt_count == 0 and deleted_nt_count == 3:
        statuses.append("DELETION")
    elif deleted_nt_count:
        statuses.append("PARTIAL_CODON")

    if inserted_nts:
        statuses.append("INSERTION")

    # IUPAC ambiguity is preserved independently of whether it becomes X in the
    # compact amino-acid representation. Internal NNN is therefore AMBIGUOUS,
    # while terminal N-only padding is already classified as NOT_COVERED above.
    observed_bases = [base for base in query_codon if base != "-"]
    if any(base not in {"A", "C", "G", "T"} for base in observed_bases):
        statuses.append("AMBIGUOUS")

    # Whole-codon indels are in-frame. Non-triplet indels create a frameshift; a
    # compensating downstream indel can restore phase without erasing the event.
    nontriplet_indel_here = (len(inserted_nts) % 3 != 0) or (deleted_nt_count % 3 != 0)
    frame_affected = phase_at_start != 0 or phase_at_end != 0 or nontriplet_indel_here
    pure_inframe_deletion = deleted_nt_count == 3 and not inserted_nts and phase_at_start == phase_at_end == 0
    if frame_affected and not pure_inframe_deletion:
        statuses.append("FRAMESHIFT_RESTORED" if phase_at_end == 0 else "FRAMESHIFT")

    # Deduplicate while preserving a stable, meaningful order.
    ordered = []
    for status in statuses:
        if status not in ordered:
            ordered.append(status)
    return "|".join(ordered) if ordered else "OBSERVED"


def insertion_calls_from_nucleotide_alignment(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    ref_aligned: str,
    qc_status: str,
    reference_system: str,
) -> list[MutationCall]:
    rows: list[MutationCall] = []
    ref_nt_position = 0
    index = 0
    while index < len(ref_aligned):
        if ref_aligned[index] != "-":
            ref_nt_position += 1
            index += 1
            continue

        inserted = []
        while index < len(ref_aligned) and ref_aligned[index] == "-":
            if query_aligned[index] != "-":
                inserted.append(query_aligned[index])
            index += 1

        if inserted and len(inserted) % 3 == 0:
            inserted_nt = "".join(inserted)
            inserted_aa = translate_ambiguous_dna(inserted_nt)
            if inserted_aa:
                position = max(ref_nt_position // 3, 1)
                rows.append(
                    build_mutation_call(
                        sequence_id=sequence_id,
                        gene=gene,
                        position=position,
                        ref_aa="-",
                        alt_aa=inserted_aa,
                        mutation_type="insertion",
                        qc_status=qc_status,
                        reference_system=reference_system,
                        hxb2_ref_codon="",
                        query_codon=inserted_nt,
                        possible_alt_aas=inserted_aa,
                        codon_status="INSERTION",
                        inserted_nts=inserted_nt,
                        deleted_nt_count=0,
                        phase_at_start=0,
                        phase_at_end=0,
                    )
                )
    return rows


def mask_terminal_unsequenced_tokens(protein_tokens: list[str], unsequenced_tokens: list[bool]) -> list[str]:
    """Do not report leading/trailing fully unsequenced codons as mutations."""
    masked = list(protein_tokens)
    start = 0
    while start < len(masked) and unsequenced_tokens[start]:
        masked[start] = "."
        start += 1
    end = len(masked) - 1
    while end >= 0 and unsequenced_tokens[end]:
        masked[end] = "."
        end -= 1
    return masked


def align_proteins(query: str, reference: str) -> tuple[str, str]:
    """Globally align a called protein to its reference."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1
    alignment = aligner.align(reference, query)[0]
    ref_aligned, query_aligned = format_pairwise_alignment(alignment, reference, query)
    return query_aligned, ref_aligned


def format_pairwise_alignment(alignment, reference: str, query: str) -> tuple[str, str]:
    """Format a Bio.Align alignment as two gapped strings."""
    ref_parts: list[str] = []
    query_parts: list[str] = []
    ref_cursor = 0
    query_cursor = 0

    for (ref_start, ref_end), (query_start, query_end) in zip(*alignment.aligned):
        if ref_start > ref_cursor:
            ref_parts.append(reference[ref_cursor:ref_start])
            query_parts.append("-" * (ref_start - ref_cursor))
        if query_start > query_cursor:
            ref_parts.append("-" * (query_start - query_cursor))
            query_parts.append(query[query_cursor:query_start])

        ref_parts.append(reference[ref_start:ref_end])
        query_parts.append(query[query_start:query_end])
        ref_cursor = ref_end
        query_cursor = query_end

    if ref_cursor < len(reference):
        ref_parts.append(reference[ref_cursor:])
        query_parts.append("-" * (len(reference) - ref_cursor))
    if query_cursor < len(query):
        ref_parts.append("-" * (len(query) - query_cursor))
        query_parts.append(query[query_cursor:])

    return "".join(ref_parts), "".join(query_parts)


def alignment_score(query: str, reference: str) -> float:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -6
    aligner.extend_gap_score = -1
    return float(aligner.score(reference, query))


def mutations_from_alignment(
    sequence_id: str,
    gene: str,
    query_aligned: str,
    ref_aligned: str,
    qc_status: str,
    reference_system: str = REFERENCE_SYSTEM_HIVDB_CONSENSUS_B,
) -> list[MutationCall]:
    """Convert aligned protein strings into mutation rows."""
    if len(query_aligned) != len(ref_aligned):
        raise ValueError("Aligned query and reference must have the same length.")

    rows: list[MutationCall] = []
    ref_position = 0
    index = 0
    while index < len(ref_aligned):
        ref_aa = ref_aligned[index]
        alt_aa = query_aligned[index]

        if ref_aa == "-":
            inserted = []
            while index < len(ref_aligned) and ref_aligned[index] == "-":
                if query_aligned[index] != "-":
                    inserted.append(query_aligned[index])
                index += 1
            if inserted:
                position = max(ref_position, 1)
                rows.append(
                    build_mutation_call(
                        sequence_id=sequence_id,
                        gene=gene,
                        position=position,
                        ref_aa="-",
                        alt_aa="".join(inserted),
                        mutation_type="insertion",
                        qc_status=qc_status,
                        reference_system=reference_system,
                    )
                )
            continue

        ref_position += 1
        if alt_aa == "-":
            rows.append(
                build_mutation_call(
                    sequence_id=sequence_id,
                    gene=gene,
                    position=ref_position,
                    ref_aa=ref_aa,
                    alt_aa="-",
                    mutation_type="deletion",
                    qc_status=qc_status,
                    reference_system=reference_system,
                )
            )
        elif should_report_substitution(ref_aa, alt_aa):
            rows.append(
                build_mutation_call(
                    sequence_id=sequence_id,
                    gene=gene,
                    position=ref_position,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    mutation_type="substitution",
                    qc_status=qc_status,
                    reference_system=reference_system,
                )
            )
        index += 1

    return rows


def codon_calls_from_protein_alignment(
    query_aligned: str,
    ref_aligned: str,
    reference_protein: str,
) -> list[dict[str, object]]:
    codon_calls = [
        {
            "alt_aa": "X",
            "codon_status": "NOT_COVERED",
            "hxb2_ref_codon": "",
            "query_codon": "",
            "possible_alt_aas": "",
            "inserted_nts": "",
            "deleted_nt_count": 0,
            "phase_at_start": 0,
            "phase_at_end": 0,
        }
        for _ in reference_protein
    ]
    ref_position = 0
    for query_aa, ref_aa in zip(query_aligned, ref_aligned):
        if ref_aa == "-":
            continue
        ref_position += 1
        if ref_position > len(codon_calls):
            continue
        if query_aa == "-":
            alt_aa, status = "-", "DELETION"
        elif query_aa == "X":
            alt_aa, status = "X", "PARTIAL_CODON"
        else:
            alt_aa, status = query_aa, "OBSERVED"
        codon_calls[ref_position - 1] = {
            **codon_calls[ref_position - 1],
            "alt_aa": alt_aa,
            "codon_status": status,
            "possible_alt_aas": alt_aa,
        }
    return codon_calls


def coverage_status_for_codon_status(codon_status: str) -> str:
    tokens = {token for token in str(codon_status).split("|") if token}
    if "NOT_COVERED" in tokens:
        return "NOT_COVERED"
    if tokens.intersection({"PARTIAL_CODON", "PARTIAL_COVERAGE"}):
        return "PARTIAL"
    return "COVERED"


def ambiguity_for_codon_status(codon_status: str) -> str:
    tokens = [token for token in str(codon_status).split("|") if token]
    uncertain = [
        token for token in tokens
        if token in {
            "NOT_COVERED", "PARTIAL_CODON", "PARTIAL_COVERAGE",
            "AMBIGUOUS", "FRAMESHIFT", "FRAMESHIFT_RESTORED",
        }
    ]
    return "|".join(uncertain) if uncertain else "NONE"


def combined_tri_state(values: Iterable[bool | None]) -> bool | None:
    values = list(values)
    if any(value is True for value in values):
        return True
    if any(value is None for value in values):
        return None
    return False


def mutation_position_qc_rows(
    sequence_id: str,
    gene: str,
    codon_calls: list[dict[str, object]],
    mutation_calls: Iterable[MutationCall],
    reference_protein: str,
    qc_status: str,
) -> list[MutationPositionQC]:
    canonical_gene = normalize_gene(gene)
    calls_by_position: dict[int, list[MutationCall]] = {}
    for call in mutation_calls:
        if call.sequence_id == sequence_id and call.gene == canonical_gene:
            calls_by_position.setdefault(call.position, []).append(call)

    rows: list[MutationPositionQC] = []
    for index, ref_aa in enumerate(reference_protein, start=1):
        codon_call = codon_calls[index - 1] if index <= len(codon_calls) else {}
        base_codon_status = str(codon_call.get("codon_status", "NOT_COVERED"))
        position_calls = sorted(
            calls_by_position.get(index, []),
            key=lambda call: (call.insertion, call.deletion, call.mutation),
        )
        status_tokens = [token for token in base_codon_status.split("|") if token]
        if any(call.insertion for call in position_calls) and "INSERTION" not in status_tokens:
            status_tokens.append("INSERTION")
        codon_status = "|".join(status_tokens) if status_tokens else "OBSERVED"
        observed_aa = normalize_alt_aa(ref_aa, str(codon_call.get("alt_aa", "X")))
        mutation = ", ".join(call.mutation for call in position_calls)
        if not mutation:
            if coverage_status_for_codon_status(codon_status) == "NOT_COVERED":
                mutation = "NOT_COVERED"
            elif should_report_substitution(ref_aa, observed_aa):
                mutation = format_mutation(canonical_gene, index, ref_aa, observed_aa, "substitution")
            else:
                mutation = "-"

        insertion = any(call.insertion for call in position_calls) or "INSERTION" in status_tokens
        deletion = any(call.deletion for call in position_calls) or "DELETION" in status_tokens
        stop = any(call.stop for call in position_calls) or "*" in observed_aa
        if position_calls:
            is_drm = combined_tri_state(call.is_drm for call in position_calls)
            drm_class = ";".join(sorted({call.drm_class for call in position_calls if call.drm_class}))
            drug_class = ";".join(sorted({call.drug_class for call in position_calls if call.drug_class}))
        elif mutation == "-":
            is_drm = False
            drm_class = ""
            drug_class = ""
        else:
            drm = annotate_drm(
                canonical_gene, index, observed_aa,
                codon_status=codon_status,
                mutation_type="deletion" if deletion else "substitution",
            )
            is_drm = drm.is_drm
            drm_class = drm.drm_class
            drug_class = drm.drug_class

        rows.append(
            MutationPositionQC(
                sequence_id=sequence_id,
                gene=canonical_gene,
                HXB2_position=index,
                Consensus_B_ref_aa=ref_aa,
                observed_aa=observed_aa,
                mutation=mutation,
                query_codon=str(codon_call.get("query_codon", "")),
                codon_status=codon_status,
                coverage_status=coverage_status_for_codon_status(codon_status),
                ambiguity=ambiguity_for_codon_status(codon_status),
                insertion=insertion,
                deletion=deletion,
                stop=stop,
                is_drm=is_drm,
                drm_class=drm_class,
                drug_class=drug_class,
                qc_status=qc_status,
            )
        )
    return rows


@lru_cache(maxsize=None)
def drm_catalog_positions(gene: str) -> tuple[int, ...]:
    canonical_gene = normalize_gene(gene)
    return tuple(sorted({position for catalog_gene, position, _ in load_drm_catalog()
                         if catalog_gene == canonical_gene}))


def summarize_drm_screening(
    sequence_id: str,
    gene: str,
    codon_calls: list[dict[str, object]],
    mutation_calls: Iterable[MutationCall],
    reference_protein: str | None = None,
) -> DRMScreeningSummary:
    canonical_gene = normalize_gene(gene)
    reference = reference_protein or reference_for_gene(canonical_gene)
    calls_by_position: dict[int, list[MutationCall]] = {}
    for call in mutation_calls:
        if call.sequence_id == sequence_id and call.gene == canonical_gene:
            calls_by_position.setdefault(call.position, []).append(call)

    positions = tuple(
        summarize_drm_screening_position(
            sequence_id=sequence_id,
            gene=canonical_gene,
            position=position,
            reference_protein=reference,
            codon_calls=codon_calls,
            mutation_calls=calls_by_position.get(position, []),
        )
        for position in drm_catalog_positions(canonical_gene)
    )
    return DRMScreeningSummary(
        sequence_id=sequence_id,
        gene=canonical_gene,
        status=drm_screening_status_for_positions(positions),
        positions=positions,
    )


def summarize_drm_screening_position(
    sequence_id: str,
    gene: str,
    position: int,
    reference_protein: str,
    codon_calls: list[dict[str, object]],
    mutation_calls: Iterable[MutationCall],
) -> DRMScreeningPosition:
    ref_aa = reference_protein[position - 1] if position <= len(reference_protein) else ""
    if position > len(codon_calls) or not ref_aa:
        return DRMScreeningPosition(
            sequence_id, gene, position, ref_aa, "X", "NOT_COVERED",
            DRM_SCREENING_NOT_COVERED,
        )

    codon_call = codon_calls[position - 1]
    codon_status = str(codon_call.get("codon_status", ""))
    alt_aa = normalize_alt_aa(ref_aa, str(codon_call.get("alt_aa", "X")))
    status_tokens = {token for token in codon_status.split("|") if token}
    if "NOT_COVERED" in status_tokens:
        screening_status = DRM_SCREENING_NOT_COVERED
    elif UNRESOLVED_CODON_STATUSES.intersection(status_tokens):
        screening_status = DRM_SCREENING_UNRESOLVED
    elif any(call.insertion or call.deletion for call in mutation_calls):
        screening_status = DRM_SCREENING_RESOLVED_MUTANT
    elif annotate_drm(gene, position, alt_aa, codon_status=codon_status).drm_evaluation_status != "COMPLETE":
        screening_status = DRM_SCREENING_UNRESOLVED
    elif should_report_substitution(ref_aa, alt_aa):
        screening_status = DRM_SCREENING_RESOLVED_MUTANT
    else:
        screening_status = DRM_SCREENING_RESOLVED_WT
    return DRMScreeningPosition(
        sequence_id, gene, position, ref_aa, alt_aa, codon_status, screening_status,
    )


def drm_screening_status_for_positions(positions: Iterable[DRMScreeningPosition]) -> str:
    statuses = [position.status for position in positions]
    if not statuses:
        return "COMPLETE"
    resolved = any(status in {DRM_SCREENING_RESOLVED_WT, DRM_SCREENING_RESOLVED_MUTANT}
                   for status in statuses)
    unresolved = any(status in {DRM_SCREENING_UNRESOLVED, DRM_SCREENING_NOT_COVERED}
                     for status in statuses)
    if resolved and unresolved:
        return "PARTIAL"
    if unresolved:
        return "UNRESOLVED"
    return "COMPLETE"


def drm_screening_status_by_sequence(
    summaries: Iterable[DRMScreeningSummary],
) -> dict[str, str]:
    positions_by_sequence: dict[str, list[DRMScreeningPosition]] = {}
    for summary in summaries:
        positions_by_sequence.setdefault(summary.sequence_id, []).extend(summary.positions)
    return {
        sequence_id: drm_screening_status_for_positions(positions)
        for sequence_id, positions in positions_by_sequence.items()
    }


def write_drm_screening_tsv(
    summaries: Iterable[DRMScreeningSummary],
    positions_output_path: Path,
    summary_output_path: Path,
) -> None:
    summaries = list(summaries)
    positions_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(positions_output_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=DRM_SCREENING_TSV_COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        for summary in summaries:
            for position in summary.positions:
                writer.writerow({
                    "sequence_id": position.sequence_id,
                    "gene": position.gene,
                    "position": position.position,
                    "ref_aa": position.ref_aa,
                    "alt_aa": position.alt_aa,
                    "codon_status": position.codon_status,
                    "status": position.status,
                })

    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_output_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=DRM_SCREENING_SUMMARY_TSV_COLUMNS, delimiter="\t"
        )
        writer.writeheader()
        for summary in summaries:
            counts = {
                status: sum(position.status == status for position in summary.positions)
                for status in (
                    DRM_SCREENING_RESOLVED_WT,
                    DRM_SCREENING_RESOLVED_MUTANT,
                    DRM_SCREENING_UNRESOLVED,
                    DRM_SCREENING_NOT_COVERED,
                )
            }
            writer.writerow({
                "sequence_id": summary.sequence_id,
                "gene": summary.gene,
                "drm_evaluation_status": summary.status,
                "resolved_wt": counts[DRM_SCREENING_RESOLVED_WT],
                "resolved_mutant": counts[DRM_SCREENING_RESOLVED_MUTANT],
                "unresolved": counts[DRM_SCREENING_UNRESOLVED],
                "not_covered": counts[DRM_SCREENING_NOT_COVERED],
                "positions": len(summary.positions),
            })


def should_report_substitution(ref_aa: str, alt_aa: str) -> bool:
    if alt_aa == ".":
        return False
    if len(alt_aa) > 1:
        return set(alt_aa) != {ref_aa}
    return alt_aa != ref_aa


def build_mutation_call(
    sequence_id: str,
    gene: str,
    position: int,
    ref_aa: str,
    alt_aa: str,
    mutation_type: str,
    qc_status: str,
    reference_system: str,
    hxb2_ref_codon: str = "",
    query_codon: str = "",
    possible_alt_aas: str = "",
    codon_status: str = "",
    inserted_nts: str = "",
    deleted_nt_count: int = 0,
    phase_at_start: int = 0,
    phase_at_end: int = 0,
) -> MutationCall:
    alt_aa = normalize_alt_aa(ref_aa, alt_aa)
    status_tokens = {token for token in codon_status.split("|") if token}
    is_unreliable_codon = bool({
        "NOT_COVERED", "PARTIAL_COVERAGE", "PARTIAL_CODON",
        "FRAMESHIFT", "FRAMESHIFT_RESTORED",
    } & status_tokens)
    mixture = (
        mutation_type == "substitution"
        and not is_unreliable_codon
        and alt_aa not in {"", "X", "-"}
        and len(set(alt_aa)) > 1
    )
    insertion = mutation_type == "insertion" or "INSERTION" in status_tokens or bool(inserted_nts)
    deletion = mutation_type == "deletion" or "DELETION" in status_tokens
    stop = "*" in alt_aa
    mutation = format_mutation(gene, position, ref_aa, alt_aa, mutation_type)
    drm = annotate_drm(gene, position, alt_aa,
                       codon_status=codon_status, mutation_type=mutation_type)
    return MutationCall(
        sequence_id=sequence_id,
        gene=gene,
        position=position,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        mutation=mutation,
        mutation_type=mutation_type,
        mixture=mixture,
        insertion=insertion,
        deletion=deletion,
        stop=stop,
        hxb2_ref_codon=hxb2_ref_codon,
        query_codon=query_codon,
        possible_alt_aas=possible_alt_aas,
        codon_status=codon_status,
        inserted_nts=inserted_nts,
        deleted_nt_count=deleted_nt_count,
        phase_at_start=phase_at_start,
        phase_at_end=phase_at_end,
        reference_system=reference_system,
        coordinate_system="HXB2",
        is_drm=drm.is_drm,
        drm_source=drm.drm_source,
        drm_version=drm.drm_version,
        drm_catalog_sha256=drm.drm_catalog_sha256,
        drm_components=drm.components,
        mutation_types=drm.mutation_types,
        is_accessory=drm.is_accessory,
        drm_class=drm.drm_class,
        drug_class=drm.drug_class,
        caller="PyHIV",
        caller_version=__version__,
        algorithm_version="hxb2-codon-local-v2",
        qc_status=qc_status,
    )


def normalize_alt_aa(ref_aa: str, alt_aa: str) -> str:
    if len(alt_aa) <= 1:
        return alt_aa
    unique = sorted(set(alt_aa))
    if ref_aa in unique:
        unique.remove(ref_aa)
        return ref_aa + "".join(unique)
    return "".join(unique)


def format_mutation(gene: str, position: int, ref_aa: str, alt_aa: str, mutation_type: str) -> str:
    if mutation_type == "deletion":
        return f"{gene}:{ref_aa}{position}del"
    if mutation_type == "insertion":
        return f"{gene}:{position}ins{alt_aa}"
    return f"{gene}:{ref_aa}{position}{alt_aa}"


def call_mutations_for_fasta_files(
    paths: Iterable[Path],
    gene: str | None = None,
    sequence_type: str = "auto",
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[MutationCall]:
    """Call mutations for records from FASTA files."""
    calls: list[MutationCall] = []
    for path in paths:
        inferred_gene = normalize_gene(gene) if gene else infer_gene_from_path(path)
        if inferred_gene is None:
            continue
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                calls.extend(
                    call_mutations_for_record(
                        record,
                        inferred_gene,
                        sequence_type=sequence_type,
                        alignment_backend=alignment_backend,
                    )
                )
    return calls


def call_mutations_for_fasta_files_with_drm_screening(
    paths: Iterable[Path],
    gene: str | None = None,
    sequence_type: str = "auto",
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], list[DRMScreeningSummary]]:
    """Call mutations and DRM screening summaries for records from FASTA files."""
    calls: list[MutationCall] = []
    summaries: list[DRMScreeningSummary] = []
    for path in paths:
        inferred_gene = normalize_gene(gene) if gene else infer_gene_from_path(path)
        if inferred_gene is None:
            continue
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                record_calls, screening = call_mutations_for_record_with_drm_screening(
                    record,
                    inferred_gene,
                    sequence_type=sequence_type,
                    alignment_backend=alignment_backend,
                )
                calls.extend(record_calls)
                summaries.append(screening)
    return calls, summaries


def call_mutations_for_fasta_files_with_qc(
    paths: Iterable[Path],
    gene: str | None = None,
    sequence_type: str = "auto",
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> tuple[list[MutationCall], list[DRMScreeningSummary], list[MutationPositionQC]]:
    """Call mutations, DRM screening summaries, and position-level QC rows."""
    calls: list[MutationCall] = []
    summaries: list[DRMScreeningSummary] = []
    position_qc_rows: list[MutationPositionQC] = []
    for path in paths:
        inferred_gene = normalize_gene(gene) if gene else infer_gene_from_path(path)
        if inferred_gene is None:
            continue
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                record_calls, screening, record_position_qc = call_mutations_for_record_with_qc(
                    record,
                    inferred_gene,
                    sequence_type=sequence_type,
                    alignment_backend=alignment_backend,
                )
                calls.extend(record_calls)
                summaries.append(screening)
                position_qc_rows.extend(record_position_qc)
    return calls, summaries, position_qc_rows


def call_mutations_for_directory(
    input_dir: Path,
    gene: str | None = None,
    sequence_type: str = "auto",
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[MutationCall]:
    """Call mutations for all supported FASTA files under a directory."""
    fasta_files = discover_fasta_files(input_dir)
    return call_mutations_for_fasta_files(
        fasta_files,
        gene=gene,
        sequence_type=sequence_type,
        alignment_backend=alignment_backend,
    )


def infer_gene_from_path(path: Path) -> str | None:
    """Infer a supported gene from a PyHIV split-output path."""
    candidates = [path.stem, *[part for part in path.parts[-4:-1]]]
    for candidate in candidates:
        cleaned = candidate.replace("_", "-")
        for token in (cleaned, cleaned.split("-")[-1]):
            try:
                return normalize_gene(token)
            except ValueError:
                continue
    return None


def sequence_gene_pairs_from_fasta_files(
    paths: Iterable[Path],
    gene: str | None = None,
) -> list[tuple[str, str]]:
    """Return sequence/gene pairs represented by split FASTA inputs."""
    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for path in paths:
        inferred_gene = normalize_gene(gene) if gene else infer_gene_from_path(path)
        if inferred_gene is None:
            continue
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                key = (record.id, inferred_gene)
                if key not in seen:
                    pairs.append(key)
                    seen.add(key)
    return pairs


def mutation_matrix_cell(call: MutationCall) -> str:
    """Return a compact mutation label for one matrix cell."""
    prefix = f"{call.gene}:"
    if call.mutation.startswith(prefix):
        return call.mutation[len(prefix):]
    return call.mutation


def write_mutation_matrices_tsv(
    calls: Iterable[MutationCall],
    output_dir: Path,
    sequence_gene_pairs: Iterable[tuple[str, str]] | None = None,
) -> list[Path]:
    """Write one reference-position mutation matrix per gene.

    The first row is the amino-acid reference. Sequence rows contain compact
    mutation labels at mutated positions and "-" at positions without a reported
    mutation.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    calls = list(calls)
    ordered_pairs: list[tuple[str, str]] = []
    seen_pairs: set[tuple[str, str]] = set()
    for sequence_id, gene in sequence_gene_pairs or ():
        canonical_gene = normalize_gene(gene)
        key = (str(sequence_id), canonical_gene)
        if key not in seen_pairs:
            ordered_pairs.append(key)
            seen_pairs.add(key)
    for call in calls:
        key = (call.sequence_id, call.gene)
        if key not in seen_pairs:
            ordered_pairs.append(key)
            seen_pairs.add(key)

    genes = [
        gene for gene in supported_gene_names()
        if any(pair_gene == gene for _, pair_gene in ordered_pairs)
        or any(call.gene == gene for call in calls)
    ]
    calls_by_gene_sequence_position: dict[tuple[str, str, int], list[MutationCall]] = {}
    for call in calls:
        calls_by_gene_sequence_position.setdefault(
            (call.gene, call.sequence_id, call.position), []
        ).append(call)

    paths: list[Path] = []
    for gene in genes:
        reference = reference_for_gene(gene)
        fields = ["Sequence"] + [
            f"{gene}:{ref_aa}{position}"
            for position, ref_aa in enumerate(reference, start=1)
        ]
        rows: list[dict[str, object]] = [
            {
                "Sequence": "REFERENCE",
                **{
                    fields[position]: reference[position - 1]
                    for position in range(1, len(reference) + 1)
                },
            }
        ]
        sequence_ids = [
            sequence_id for sequence_id, pair_gene in ordered_pairs
            if pair_gene == gene
        ]
        for sequence_id in sequence_ids:
            row = {"Sequence": sequence_id}
            for position in range(1, len(reference) + 1):
                position_calls = sorted(
                    calls_by_gene_sequence_position.get((gene, sequence_id, position), []),
                    key=lambda call: (call.position, call.mutation),
                )
                row[fields[position]] = (
                    ", ".join(mutation_matrix_cell(call) for call in position_calls)
                    if position_calls else "-"
                )
            rows.append(row)

        output_path = output_dir / f"mutation_matrix_{gene}.tsv"
        with open(output_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        paths.append(output_path)
    return paths



def write_mutation_position_qc_tsv(
    rows: Iterable[MutationPositionQC],
    output_path: Path,
) -> None:
    """Write one row per sequence/gene/reference position for KG ingestion."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=MUTATION_POSITION_QC_TSV_COLUMNS, delimiter="	"
        )
        writer.writeheader()
        for row in rows:
            tsv_row = row.to_tsv_row()
            writer.writerow({
                column: tsv_row[column]
                for column in MUTATION_POSITION_QC_TSV_COLUMNS
            })


def write_mutations_tsv(calls: Iterable[MutationCall], output_path: Path) -> None:
    """Write mutation calls to a TSV file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MUTATION_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for call in calls:
            row = call.to_tsv_row()
            writer.writerow({column: row[column] for column in MUTATION_TSV_COLUMNS})


def supported_gene_names() -> tuple[str, ...]:
    return ("PR", "RT", "IN", "CA")

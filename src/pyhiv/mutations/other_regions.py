"""HXB2 mutation calling for regions outside PR, RT, IN and CA.

This module deliberately keeps non-Stanford/non-HIVfacts resistance markers
separate from the existing DRM caller.  It reports amino-acid mutations against
HXB2-derived protein references for additional coding regions, and annotates
only markers present in the local ``other_region_resistance.tsv`` catalogue.
"""

from __future__ import annotations

import csv
import logging
import os
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Iterable

from Bio import SeqIO

from pyhiv import __version__
from pyhiv.mutations.caller import (
    DEFAULT_NT_ALIGNMENT_BACKEND,
    align_dna_to_hxb2_columns,
    amino_acid_token_from_possible,
    codon_calls_from_nucleotide_alignment,
    hxb2_genome_sequence,
    infer_sequence_type,
    normalize_alt_aa,
    should_report_substitution,
    translate_ambiguous_dna,
)


REFERENCE_SYSTEM_HXB2 = "HXB2"
CATALOG_FILENAME = "other_region_resistance.tsv"
OTHER_REGION_CATALOG_ENV_VAR = "PYHIV_OTHER_REGION_RESISTANCE_CATALOG"
RESISTANCE_STATUS_MARKER_MATCH = "CATALOG_MARKER_MATCH"
RESISTANCE_STATUS_EVALUATED_NO_MATCH = "EVALUATED_NO_MATCH"
RESISTANCE_STATUS_NOT_EVALUATED_NO_CATALOG = "NOT_EVALUATED_NO_CATALOG_MARKERS"
NO_REGION_CATALOG_MARKERS_MESSAGE = "No resistance markers are described for this region in the PyHIV other-region catalogue; mutations are reported against HXB2 only."
NO_MATCHING_CATALOG_MARKER_MESSAGE = "No matching resistance marker was found in the PyHIV other-region catalogue."
MATCHING_CATALOG_MARKER_MESSAGE = "Mutation matched a marker in the PyHIV other-region resistance catalogue."


@dataclass(frozen=True)
class OtherRegionReference:
    region: str
    display_name: str
    start: int
    end: int
    aliases: tuple[str, ...]


BROAD_OTHER_REGION_TARGETS = {
    "GAG": ("MA", "NC", "P6"),
    "GAG-CDS": ("MA", "NC", "P6"),
    "GAG-POL": ("MA", "NC", "P6"),
    "GAG-POL-FUSION": ("MA", "NC", "P6"),
    "GAG-POL-PRECURSOR": ("MA", "NC", "P6"),
}


OTHER_REGION_REFERENCES: dict[str, OtherRegionReference] = {
    "MA": OtherRegionReference("MA", "p17 matrix", 790, 1185, ("p17", "p17-matrix", "matrix")),
    "NC": OtherRegionReference("NC", "p7 nucleocapsid", 1921, 2085, ("p7", "p7-nucleocapsid", "nucleocapsid")),
    "P6": OtherRegionReference("P6", "p6", 2134, 2292, ("p6",)),
    "RH": OtherRegionReference("RH", "p15 RNase H", 3870, 4229, ("p15", "p15-rnase-h", "rnase-h", "rnaseh")),
    "VIF": OtherRegionReference("VIF", "Vif", 5041, 5619, ("vif", "vif-cds")),
    "VPR": OtherRegionReference("VPR", "Vpr", 5559, 5850, ("vpr", "vpr-cds")),
    "VPU": OtherRegionReference("VPU", "Vpu", 6062, 6310, ("vpu", "vpu-cds")),
    "ENV": OtherRegionReference("ENV", "Env gp160", 6225, 8795, ("env", "env-cds", "gp160")),
    "GP120": OtherRegionReference("GP120", "gp120", 6225, 7757, ("gp120",)),
    "GP41": OtherRegionReference("GP41", "gp41", 7758, 8795, ("gp41",)),
    "NEF": OtherRegionReference("NEF", "Nef", 8797, 9417, ("nef", "nef-cds")),
}

SUPPORTED_OTHER_REGIONS = tuple(OTHER_REGION_REFERENCES)


@dataclass(frozen=True)
class OtherRegionResistanceMarker:
    region: str
    position: int
    ref_aa: str
    alt_aa: str
    drug_class: str
    drug: str
    resistance_scope: str
    evidence_level: str
    source: str
    source_url: str
    note: str


@dataclass(frozen=True)
class OtherRegionMutationCall:
    sequence_id: str
    region: str
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
    reference_system: str
    coordinate_system: str
    resistance_marker: bool
    resistance_evaluation_status: str
    resistance_evaluation_note: str
    resistance_scope: str
    drug_class: str
    drug: str
    evidence_level: str
    resistance_source: str
    resistance_note: str
    caller: str
    caller_version: str
    algorithm_version: str
    qc_status: str

    def to_tsv_row(self) -> dict[str, object]:
        return self.__dict__.copy()


OTHER_REGION_MUTATION_TSV_COLUMNS = [
    "sequence_id",
    "region",
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
    "reference_system",
    "coordinate_system",
    "resistance_marker",
    "resistance_evaluation_status",
    "resistance_evaluation_note",
    "resistance_scope",
    "drug_class",
    "drug",
    "evidence_level",
    "resistance_source",
    "resistance_note",
    "caller",
    "caller_version",
    "algorithm_version",
    "qc_status",
]

OTHER_REGION_RESISTANCE_TSV_COLUMNS = [
    "sequence_id",
    "region",
    "position",
    "ref_aa",
    "alt_aa",
    "mutation",
    "drug_class",
    "drug",
    "resistance_scope",
    "evidence_level",
    "resistance_source",
    "resistance_note",
]

OTHER_REGION_WARNING_TSV_COLUMNS = [
    "sequence_id",
    "region",
    "warning_type",
    "message",
    "reference_system",
    "catalog",
]


def normalize_other_region(region: str) -> str:
    key = _clean_region_token(region)
    for canonical, reference in OTHER_REGION_REFERENCES.items():
        if key == _clean_region_token(canonical):
            return canonical
        if key in {_clean_region_token(alias) for alias in reference.aliases}:
            return canonical
    raise ValueError(f"Unsupported other-region mutation target: {region!r}.")


def infer_other_region_from_path(path: Path) -> str | None:
    targets = infer_other_region_targets_from_path(path, include_broad=False)
    return targets[0] if targets else None


def infer_other_region_targets_from_path(path: Path, include_broad: bool = True) -> tuple[str, ...]:
    candidates: list[str] = []
    candidates.append(path.stem)
    candidates.extend(path.parts[-5:-1])
    for candidate in reversed(candidates):
        for token in _candidate_tokens(candidate):
            try:
                return (normalize_other_region(token),)
            except ValueError:
                continue
    if include_broad:
        broad_candidates = [path.parent.name, path.stem]
        for candidate in broad_candidates:
            for token in _broad_region_path_tokens(candidate):
                broad_targets = BROAD_OTHER_REGION_TARGETS.get(_clean_region_token(token).upper())
                if broad_targets:
                    return broad_targets
    return ()


def hxb2_reference_nt_for_other_region(region: str) -> str:
    reference = OTHER_REGION_REFERENCES[normalize_other_region(region)]
    return hxb2_genome_sequence()[reference.start - 1:reference.end]


def hxb2_reference_protein_for_other_region(region: str) -> str:
    nt = hxb2_reference_nt_for_other_region(region)
    coding_length = len(nt) - (len(nt) % 3)
    return translate_ambiguous_dna(nt[:coding_length]).rstrip("*")


def call_other_region_mutations(
    sequence: str,
    sequence_id: str,
    region: str,
    sequence_type: str = "auto",
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[OtherRegionMutationCall]:
    canonical_region = normalize_other_region(region)
    reference_nt = hxb2_reference_nt_for_other_region(canonical_region)
    reference_protein = hxb2_reference_protein_for_other_region(canonical_region)
    cleaned = "".join(str(sequence).upper().split())
    effective_sequence_type = sequence_type.lower()
    if effective_sequence_type == "auto":
        effective_sequence_type = infer_sequence_type(cleaned)
    if effective_sequence_type != "nt":
        raise ValueError("Other-region mutation calling currently requires nucleotide input.")

    query_aligned, ref_aligned = align_dna_to_hxb2_columns(
        cleaned,
        reference_nt,
        alignment_backend=alignment_backend,
    )
    codon_calls = codon_calls_from_nucleotide_alignment(
        query_aligned=query_aligned,
        ref_aligned=ref_aligned,
        reference_protein=reference_protein,
    )
    protein_tokens = [call["alt_aa"] for call in codon_calls]
    qc_status = "WARN" if any("*" in token or "X" in token for token in protein_tokens) else "PASS"
    calls: list[OtherRegionMutationCall] = []
    for index, codon_call in enumerate(codon_calls[:len(reference_protein)]):
        if codon_call["codon_status"] == "NOT_COVERED":
            continue
        ref_aa = reference_protein[index]
        alt_aa = str(codon_call["alt_aa"])
        if alt_aa == "-":
            mutation_type = "deletion"
        elif "FRAMESHIFT" in str(codon_call["codon_status"]):
            mutation_type = "frameshift"
        elif "PARTIAL_CODON" in str(codon_call["codon_status"]):
            mutation_type = "partial_codon"
        else:
            mutation_type = "substitution"
        rendered_alt = normalize_alt_aa(ref_aa, alt_aa)
        if rendered_alt == "-" or should_report_substitution(ref_aa, rendered_alt):
            calls.append(
                build_other_region_mutation_call(
                    sequence_id=sequence_id,
                    region=canonical_region,
                    position=index + 1,
                    ref_aa=ref_aa,
                    alt_aa=rendered_alt,
                    mutation_type=mutation_type,
                    qc_status=qc_status,
                    hxb2_ref_codon=str(codon_call["hxb2_ref_codon"]),
                    query_codon=str(codon_call["query_codon"]),
                    possible_alt_aas=str(codon_call["possible_alt_aas"]),
                    codon_status=str(codon_call["codon_status"]),
                    inserted_nts=str(codon_call["inserted_nts"]),
                    deleted_nt_count=int(codon_call["deleted_nt_count"]),
                )
            )
    calls.extend(
        insertion_calls_from_other_region_alignment(
            sequence_id=sequence_id,
            region=canonical_region,
            query_aligned=query_aligned,
            ref_aligned=ref_aligned,
            qc_status=qc_status,
        )
    )
    return calls


def call_other_region_mutations_for_fasta_files(
    paths: Iterable[Path],
    region: str | None = None,
    sequence_type: str = "auto",
    alignment_backend: str = DEFAULT_NT_ALIGNMENT_BACKEND,
) -> list[OtherRegionMutationCall]:
    calls: list[OtherRegionMutationCall] = []
    for path in paths:
        inferred_regions = (normalize_other_region(region),) if region else infer_other_region_targets_from_path(path)
        if not inferred_regions:
            continue
        with open(path, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
        for inferred_region in inferred_regions:
            for record in records:
                try:
                    calls.extend(
                        call_other_region_mutations(
                            sequence=str(record.seq),
                            sequence_id=record.id,
                            region=inferred_region,
                            sequence_type=sequence_type,
                            alignment_backend=alignment_backend,
                        )
                    )
                except ValueError as exc:
                    logging.warning(
                        "Skipping other-region mutation calling for %s in %s: %s",
                        record.id,
                        path,
                        exc,
                    )
    return calls


def sequence_region_pairs_from_fasta_files(
    paths: Iterable[Path],
    region: str | None = None,
) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for path in paths:
        inferred_regions = (normalize_other_region(region),) if region else infer_other_region_targets_from_path(path)
        if not inferred_regions:
            continue
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                for inferred_region in inferred_regions:
                    key = (record.id, inferred_region)
                    if key not in seen:
                        seen.add(key)
                        pairs.append(key)
    return pairs


def build_other_region_mutation_call(
    sequence_id: str,
    region: str,
    position: int,
    ref_aa: str,
    alt_aa: str,
    mutation_type: str,
    qc_status: str,
    hxb2_ref_codon: str = "",
    query_codon: str = "",
    possible_alt_aas: str = "",
    codon_status: str = "",
    inserted_nts: str = "",
    deleted_nt_count: int = 0,
) -> OtherRegionMutationCall:
    alt_aa = normalize_alt_aa(ref_aa, alt_aa)
    mutation = format_other_region_mutation(region, position, ref_aa, alt_aa, mutation_type)
    markers = annotate_other_region_resistance(region, position, ref_aa, alt_aa, mutation_type)
    marker = markers[0] if markers else None
    if marker is not None:
        resistance_evaluation_status = RESISTANCE_STATUS_MARKER_MATCH
        resistance_evaluation_note = MATCHING_CATALOG_MARKER_MESSAGE
    elif has_other_region_resistance_markers(region):
        resistance_evaluation_status = RESISTANCE_STATUS_EVALUATED_NO_MATCH
        resistance_evaluation_note = NO_MATCHING_CATALOG_MARKER_MESSAGE
    else:
        resistance_evaluation_status = RESISTANCE_STATUS_NOT_EVALUATED_NO_CATALOG
        resistance_evaluation_note = NO_REGION_CATALOG_MARKERS_MESSAGE
    is_unreliable = any(
        token in {"NOT_COVERED", "PARTIAL_COVERAGE", "PARTIAL_CODON", "FRAMESHIFT", "FRAMESHIFT_RESTORED"}
        for token in codon_status.split("|")
    )
    mixture = (
        mutation_type == "substitution"
        and not is_unreliable
        and alt_aa not in {"", "X", "-"}
        and len(set(alt_aa)) > 1
    )
    return OtherRegionMutationCall(
        sequence_id=sequence_id,
        region=region,
        position=position,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        mutation=mutation,
        mutation_type=mutation_type,
        mixture=mixture,
        insertion=mutation_type == "insertion" or bool(inserted_nts),
        deletion=mutation_type == "deletion",
        stop="*" in alt_aa,
        hxb2_ref_codon=hxb2_ref_codon,
        query_codon=query_codon,
        possible_alt_aas=possible_alt_aas,
        codon_status=codon_status,
        inserted_nts=inserted_nts,
        deleted_nt_count=deleted_nt_count,
        reference_system=REFERENCE_SYSTEM_HXB2,
        coordinate_system="HXB2",
        resistance_marker=marker is not None,
        resistance_evaluation_status=resistance_evaluation_status,
        resistance_evaluation_note=resistance_evaluation_note,
        resistance_scope=marker.resistance_scope if marker else "",
        drug_class=marker.drug_class if marker else "",
        drug=marker.drug if marker else "",
        evidence_level=marker.evidence_level if marker else "",
        resistance_source=marker.source if marker else "",
        resistance_note=marker.note if marker else "",
        caller="PyHIV",
        caller_version=__version__,
        algorithm_version="hxb2-other-region-codon-v1",
        qc_status=qc_status,
    )


def insertion_calls_from_other_region_alignment(
    sequence_id: str,
    region: str,
    query_aligned: str,
    ref_aligned: str,
    qc_status: str,
) -> list[OtherRegionMutationCall]:
    calls: list[OtherRegionMutationCall] = []
    ref_nt_position = 0
    index = 0
    while index < len(ref_aligned):
        if ref_aligned[index] != "-":
            ref_nt_position += 1
            index += 1
            continue
        inserted: list[str] = []
        while index < len(ref_aligned) and ref_aligned[index] == "-":
            if query_aligned[index] != "-":
                inserted.append(query_aligned[index])
            index += 1
        if inserted and len(inserted) % 3 == 0:
            inserted_nts = "".join(inserted)
            inserted_aa = translate_ambiguous_dna(inserted_nts)
            if inserted_aa:
                calls.append(
                    build_other_region_mutation_call(
                        sequence_id=sequence_id,
                        region=region,
                        position=max(ref_nt_position // 3, 1),
                        ref_aa="-",
                        alt_aa=inserted_aa,
                        mutation_type="insertion",
                        qc_status=qc_status,
                        query_codon=inserted_nts,
                        possible_alt_aas=inserted_aa,
                        codon_status="INSERTION",
                        inserted_nts=inserted_nts,
                    )
                )
    return calls


def annotate_other_region_resistance(
    region: str,
    position: int,
    ref_aa: str,
    alt_aa: str,
    mutation_type: str,
) -> tuple[OtherRegionResistanceMarker, ...]:
    catalog = load_other_region_resistance_catalog()
    lookup_alt = "del" if mutation_type == "deletion" else ("ins" if mutation_type == "insertion" else alt_aa)
    candidates: list[str]
    if lookup_alt in {"del", "ins"}:
        candidates = [lookup_alt]
    else:
        candidates = [aa for aa in lookup_alt if aa != ref_aa]
    markers: list[OtherRegionResistanceMarker] = []
    for aa in candidates:
        markers.extend(catalog.get((normalize_other_region(region), position, aa), ()))
    return tuple(markers)


def has_other_region_resistance_markers(region: str) -> bool:
    canonical_region = normalize_other_region(region)
    return any(key_region == canonical_region for key_region, _, _ in load_other_region_resistance_catalog())


@lru_cache(maxsize=1)
def load_other_region_resistance_catalog(
    catalog_path: str | Path | None = None,
) -> dict[tuple[str, int, str], tuple[OtherRegionResistanceMarker, ...]]:
    rows_by_key: dict[tuple[str, int, str], list[OtherRegionResistanceMarker]] = {}
    with _open_other_region_catalog(catalog_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row.get("region") or not row.get("position") or not row.get("alt_aa"):
                continue
            marker = OtherRegionResistanceMarker(
                region=normalize_other_region(row["region"]),
                position=int(row["position"]),
                ref_aa=row.get("ref_aa", ""),
                alt_aa=row["alt_aa"],
                drug_class=row.get("drug_class", ""),
                drug=row.get("drug", ""),
                resistance_scope=row.get("resistance_scope", ""),
                evidence_level=row.get("evidence_level", ""),
                source=row.get("source", ""),
                source_url=row.get("source_url", ""),
                note=row.get("note", ""),
            )
            rows_by_key.setdefault((marker.region, marker.position, marker.alt_aa), []).append(marker)
    return {key: tuple(value) for key, value in rows_by_key.items()}


def _open_other_region_catalog(catalog_path: str | Path | None = None):
    configured_path = catalog_path or os.environ.get(OTHER_REGION_CATALOG_ENV_VAR)
    if configured_path:
        return open(Path(configured_path), "r", newline="")
    resource = files("pyhiv.mutations.data").joinpath(CATALOG_FILENAME)
    return resource.open("r", newline="")


def write_other_region_mutations_tsv(
    calls: Iterable[OtherRegionMutationCall],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OTHER_REGION_MUTATION_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for call in calls:
            row = call.to_tsv_row()
            writer.writerow({column: row[column] for column in OTHER_REGION_MUTATION_TSV_COLUMNS})


def write_other_region_warnings_tsv(
    sequence_region_pairs: Iterable[tuple[str, str]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    seen: set[tuple[str, str]] = set()
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OTHER_REGION_WARNING_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for sequence_id, region in sequence_region_pairs:
            canonical_region = normalize_other_region(region)
            key = (str(sequence_id), canonical_region)
            if key in seen or has_other_region_resistance_markers(canonical_region):
                continue
            seen.add(key)
            writer.writerow({
                "sequence_id": sequence_id,
                "region": canonical_region,
                "warning_type": RESISTANCE_STATUS_NOT_EVALUATED_NO_CATALOG,
                "message": NO_REGION_CATALOG_MARKERS_MESSAGE,
                "reference_system": REFERENCE_SYSTEM_HXB2,
                "catalog": CATALOG_FILENAME,
            })


def write_other_region_resistance_tsv(
    calls: Iterable[OtherRegionMutationCall],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OTHER_REGION_RESISTANCE_TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for call in calls:
            if not call.resistance_marker:
                continue
            writer.writerow({
                "sequence_id": call.sequence_id,
                "region": call.region,
                "position": call.position,
                "ref_aa": call.ref_aa,
                "alt_aa": call.alt_aa,
                "mutation": call.mutation,
                "drug_class": call.drug_class,
                "drug": call.drug,
                "resistance_scope": call.resistance_scope,
                "evidence_level": call.evidence_level,
                "resistance_source": call.resistance_source,
                "resistance_note": call.resistance_note,
            })


def write_other_region_mutation_matrices_tsv(
    calls: Iterable[OtherRegionMutationCall],
    output_dir: Path,
    sequence_region_pairs: Iterable[tuple[str, str]] | None = None,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    calls = list(calls)
    ordered_pairs: list[tuple[str, str]] = []
    seen_pairs: set[tuple[str, str]] = set()
    for sequence_id, region in sequence_region_pairs or ():
        canonical_region = normalize_other_region(region)
        key = (str(sequence_id), canonical_region)
        if key not in seen_pairs:
            ordered_pairs.append(key)
            seen_pairs.add(key)
    for call in calls:
        key = (call.sequence_id, call.region)
        if key not in seen_pairs:
            ordered_pairs.append(key)
            seen_pairs.add(key)

    paths: list[Path] = []
    calls_by_region_sequence_position: dict[tuple[str, str, int], list[OtherRegionMutationCall]] = {}
    for call in calls:
        calls_by_region_sequence_position.setdefault(
            (call.region, call.sequence_id, call.position), []
        ).append(call)

    for region in SUPPORTED_OTHER_REGIONS:
        if not any(pair_region == region for _, pair_region in ordered_pairs):
            continue
        reference = hxb2_reference_protein_for_other_region(region)
        fields = ["Sequence"] + [
            f"{region}:{ref_aa}{position}"
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
        for sequence_id, pair_region in ordered_pairs:
            if pair_region != region:
                continue
            row = {"Sequence": sequence_id}
            for position in range(1, len(reference) + 1):
                position_calls = sorted(
                    calls_by_region_sequence_position.get((region, sequence_id, position), []),
                    key=lambda call: (call.position, call.mutation),
                )
                row[fields[position]] = (
                    ", ".join(mutation_matrix_cell(call) for call in position_calls)
                    if position_calls else "-"
                )
            rows.append(row)
        output_path = output_dir / f"other_region_mutation_matrix_{region}.tsv"
        with open(output_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        paths.append(output_path)
    return paths


def mutation_matrix_cell(call: OtherRegionMutationCall) -> str:
    prefix = f"{call.region}:"
    if call.mutation.startswith(prefix):
        return call.mutation[len(prefix):]
    return call.mutation


def format_other_region_mutation(region: str, position: int, ref_aa: str, alt_aa: str, mutation_type: str) -> str:
    if mutation_type == "deletion":
        return f"{region}:{ref_aa}{position}del"
    if mutation_type == "insertion":
        return f"{region}:{position}ins{alt_aa}"
    return f"{region}:{ref_aa}{position}{alt_aa}"


def _candidate_tokens(value: str) -> tuple[str, ...]:
    cleaned = str(value or "")
    stem = Path(cleaned).stem
    tokens = {cleaned, stem}
    tokens.update(part for part in cleaned.replace("_", "-").split("-") if part)
    return tuple(tokens)


def _broad_region_path_tokens(value: str) -> tuple[str, ...]:
    cleaned = _clean_region_token(value).strip("-")
    parts = [part for part in cleaned.split("-") if part]
    tokens = [cleaned]
    if len(parts) >= 2:
        tokens.extend("-".join(parts[index:]) for index in range(1, len(parts)))
    tokens.extend(parts)
    seen: set[str] = set()
    ordered: list[str] = []
    for token in tokens:
        key = token.upper()
        if key not in seen:
            seen.add(key)
            ordered.append(token)
    return tuple(ordered)


def _clean_region_token(value: str) -> str:
    return str(value or "").strip().lower().replace("_", "-").replace(" ", "-")

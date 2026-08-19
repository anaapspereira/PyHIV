from __future__ import annotations

from dataclasses import dataclass
import csv
from collections import Counter
from collections.abc import Iterable, Mapping
from html import unescape
from io import StringIO
from pathlib import Path
import re
import shutil
import tempfile
import time
from typing import Callable
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import warnings

import pandas as pd

from pyhiv.config import get_reference_paths

try:
    from Bio import Entrez, SeqIO
except ImportError:  # pragma: no cover
    raise ImportError("BioPython is required. Please install it via 'pip install biopython'.")


LANL_CRF_URL = "https://www.hiv.lanl.gov/components/sequence/HIV/crfdb/crfs.comp"
LANL_ALIGNMENT_LIST_URL = "https://www.hiv.lanl.gov/content/sequence/NEWALIGN/list.comp"
LANL_ALIGNMENT_DOWNLOAD_URL = "https://www.hiv.lanl.gov/cgi-bin/NEWALIGN/align.cgi"
LANL_ALIGNMENT_PREFIX = "HIV1_REF"
LANL_ALIGNMENT_SUFFIX = "genome_DNA.csv"
NCBI_TOOL_NAME = "pyhiv_reference_dataset_update"
USER_AGENT = "PyHIV-reference-dataset-updater/1.0"
HIV1_NEAR_COMPLETE_MIN_LENGTH = 8000
REFERENCE_TABLE_COLUMNS = [
    "subtype",
    "country",
    "year",
    "isolate",
    "accession",
    "align sequence",
    "sequence",
    "features",
    "group",
]
ACC_PATTERN = (
    r"(?:[A-Z]\d{5}|[A-Z]{2}\d{6}|[A-Z]{2}\d{8}|"
    r"[A-Z]{4}\d{8,10}|[A-Z]{2}_\d{6,9})(?:\.\d+)?"
)
ACC_RE = re.compile(rf"\b({ACC_PATTERN})\b")
STRAIN_ACCESSION_RE = re.compile(
    rf"(?P<strain>[^\s()]+)\s*\((?P<accession>{ACC_PATTERN})\)"
)
LANL_ALIGNMENT_RE = re.compile(
    rf"^{LANL_ALIGNMENT_PREFIX}_(\d{{4}})_{LANL_ALIGNMENT_SUFFIX}$"
)


@dataclass(frozen=True)
class CrfReference:
    label: str
    crf: str
    strain: str
    accession: str


CrfReferenceCollection = (
    Mapping[str, CrfReference | Iterable[CrfReference]]
    | Iterable[CrfReference]
)


class GenBankRecordUnavailable(RuntimeError):
    """Raised when NCBI has no GenBank record for a requested accession."""


class GenBankRecordPartial(RuntimeError):
    """Raised when a GenBank record is not complete or nearly complete."""


class ReconciliationError(RuntimeError):
    """Raised when reference FASTAs cannot be reconciled deterministically."""


class LocalReferenceSequenceConflict(RuntimeError):
    """Raised when local FASTAs for one accession contain different sequences."""


@dataclass(frozen=True)
class LanlAlignment:
    year: int
    content: str

    @property
    def filename(self) -> str:
        return f"{LANL_ALIGNMENT_PREFIX}_{self.year}_{LANL_ALIGNMENT_SUFFIX}"


@dataclass(frozen=True)
class LanlAlignmentUpdateResult:
    year: int | None
    path: Path | None
    updated: bool
    fasta_records: int
    fasta_files_added: int
    fasta_files_updated: int
    fasta_files_removed: int
    fasta_files_reclassified: int = 0


@dataclass(frozen=True)
class UpdateResult:
    added_sequences: int
    added_rows: int
    updated_rows: int
    updated_feature_rows: int
    total_references: int
    sequences_with_locations: Path
    reference_fastas_dir: Path
    lanl_alignment_year: int | None = None
    lanl_alignment_updated: bool = False
    lanl_fasta_records: int = 0
    lanl_fasta_files_added: int = 0
    lanl_fasta_files_updated: int = 0
    lanl_fasta_files_removed: int = 0
    reclassified_references: int = 0
    partial_sequences_rejected: int = 0
    duplicate_rows_removed: int = 0
    crfs_with_multiple_references: int = 0
    pending_ncbi_references: int = 0
    pending_ncbi_accessions: tuple[str, ...] = ()
    failed_sequences: int = 0
    failed_accessions: tuple[str, ...] = ()

    @property
    def up_to_date(self) -> bool:
        return (
            self.added_sequences == 0
            and self.added_rows == 0
            and self.updated_rows == 0
            and self.updated_feature_rows == 0
            and not self.lanl_alignment_updated
            and self.lanl_fasta_files_added == 0
            and self.lanl_fasta_files_updated == 0
            and self.lanl_fasta_files_removed == 0
            and self.reclassified_references == 0
            and self.partial_sequences_rejected == 0
            and self.duplicate_rows_removed == 0
            and self.pending_ncbi_references == 0
            and self.failed_sequences == 0
        )


def update_reference_dataset(
    base_dir: Path | None = None,
    email: str | None = None,
    ncbi_api_key: str | None = None,
    dry_run: bool = False,
    refresh_features: bool = False,
    sleep_seconds: float = 0.34,
    update_lanl_alignment: bool = True,
    lanl_alignment_loader: Callable[[], LanlAlignment] | None = None,
    crf_loader: Callable[[], CrfReferenceCollection] | None = None,
    genbank_fetcher: Callable[[str], tuple[str, dict[str, tuple[int, int]] | None]] | None = None,
) -> UpdateResult:
    paths = get_reference_paths(base_dir)
    reference_base_dir = paths["REFERENCE_GENOMES_DIR"]
    references_dir = paths["REFERENCE_GENOMES_FASTAS_DIR"]
    sequences_path = paths["SEQUENCES_WITH_LOCATION"]

    references_dir.mkdir(parents=True, exist_ok=True)

    lanl_result = LanlAlignmentUpdateResult(
        year=None,
        path=None,
        updated=False,
        fasta_records=0,
        fasta_files_added=0,
        fasta_files_updated=0,
        fasta_files_removed=0,
    )
    alignment_records: list[LanlReferenceRecord] = []
    if update_lanl_alignment:
        alignment = (lanl_alignment_loader or download_latest_lanl_reference_alignment)()
        alignment_records = parse_lanl_alignment_records(alignment.content)
        lanl_alignment_updated = save_lanl_alignment_file(
            reference_base_dir,
            alignment,
            dry_run=dry_run,
        )
        lanl_result = LanlAlignmentUpdateResult(
            year=alignment.year,
            path=reference_base_dir / alignment.filename,
            updated=lanl_alignment_updated,
            fasta_records=len(alignment_records),
            fasta_files_added=0,
            fasta_files_updated=0,
            fasta_files_removed=0,
        )

    raw_crf_references = (crf_loader or load_lanl_crf_references)()
    crf_references = normalize_crf_references(raw_crf_references)
    desired_references = build_desired_reference_map(alignment_records, crf_references)
    local_accessions = reference_accessions_in_dir(references_dir)
    missing_crfs = [
        ref for ref in desired_references.values()
        if ref.sequence is None and ref.accession not in local_accessions
    ]

    added_sequences = 0
    if missing_crfs and not email and genbank_fetcher is None:
        raise ValueError(
            "NCBI email is required to download missing reference sequences. "
            "Pass --email or set NCBI_EMAIL."
        )

    if email:
        Entrez.email = email
    if ncbi_api_key:
        Entrez.api_key = ncbi_api_key

    fetcher = genbank_fetcher or fetch_genbank_sequence_and_features
    fetched_features: dict[str, dict[str, tuple[int, int]] | None] = {}
    partial_sequences_rejected = 0

    pending_ncbi_accessions: list[str] = []
    failed_accessions: list[str] = []
    for crf_ref in missing_crfs:
        try:
            sequence, features = fetcher(crf_ref.accession)
        except GenBankRecordPartial:
            partial_sequences_rejected += 1
            failed_accessions.append(
                f"{crf_ref.subtype}:{crf_ref.accession}:not_complete_or_near_complete"
            )
            continue
        except GenBankRecordUnavailable:
            pending_ncbi_accessions.append(
                f"{crf_ref.subtype}:{crf_ref.accession}:authoritative_reference_pending_ncbi"
            )
            continue
        except HTTPError as exc:
            if exc.code not in {400, 404}:
                raise
            pending_ncbi_accessions.append(
                f"{crf_ref.subtype}:{crf_ref.accession}:authoritative_reference_pending_ncbi"
            )
            continue

        fetched_features[crf_ref.accession] = features
        desired_references[crf_ref.accession] = DesiredReference(
            accession=crf_ref.accession,
            subtype=crf_ref.subtype,
            sequence=sequence,
            source=crf_ref.source,
        )
        added_sequences += 1
        time.sleep(sleep_seconds)

    with tempfile.TemporaryDirectory(
        prefix="reference_update_",
        dir=reference_base_dir,
    ) as tmp_dir_name:
        staging_dir = Path(tmp_dir_name)
        staged_references_dir = staging_dir / references_dir.name
        staged_sequences_path = staging_dir / sequences_path.name
        fasta_result = reconcile_reference_fastas(
            desired_references.values(),
            references_dir,
            dry_run=True,
            output_dir=staged_references_dir,
            fetcher=fetcher,
            fetched_features=fetched_features,
        )
        lanl_result = LanlAlignmentUpdateResult(
            year=lanl_result.year,
            path=lanl_result.path,
            updated=lanl_result.updated,
            fasta_records=lanl_result.fasta_records,
            fasta_files_added=fasta_result["files_added"],
            fasta_files_updated=fasta_result["files_updated"],
            fasta_files_removed=fasta_result["files_removed"],
            fasta_files_reclassified=fasta_result["files_reclassified"],
        )

        table_result = rebuild_sequences_with_locations(
            staged_references_dir,
            sequences_path,
            fetched_features=fetched_features,
            refresh_features=refresh_features,
            dry_run=False,
            fetcher=fetcher,
            sleep_seconds=sleep_seconds,
            output_path=staged_sequences_path,
        )

        if not dry_run:
            commit_staged_reference_dataset(
                staged_references_dir,
                references_dir,
                staged_sequences_path,
                sequences_path,
            )

    return UpdateResult(
        added_sequences=added_sequences,
        added_rows=table_result["added_rows"],
        updated_rows=table_result["updated_rows"],
        updated_feature_rows=table_result["updated_feature_rows"],
        total_references=table_result["total_references"],
        sequences_with_locations=sequences_path,
        reference_fastas_dir=references_dir,
        lanl_alignment_year=lanl_result.year,
        lanl_alignment_updated=lanl_result.updated,
        lanl_fasta_records=lanl_result.fasta_records,
        lanl_fasta_files_added=lanl_result.fasta_files_added,
        lanl_fasta_files_updated=lanl_result.fasta_files_updated,
        lanl_fasta_files_removed=lanl_result.fasta_files_removed,
        reclassified_references=lanl_result.fasta_files_reclassified,
        partial_sequences_rejected=partial_sequences_rejected,
        duplicate_rows_removed=table_result["duplicate_rows_removed"],
        crfs_with_multiple_references=count_crfs_with_multiple_references(raw_crf_references),
        pending_ncbi_references=len(pending_ncbi_accessions),
        pending_ncbi_accessions=tuple(pending_ncbi_accessions),
        failed_sequences=len(failed_accessions),
        failed_accessions=tuple(failed_accessions),
    )


def update_lanl_reference_alignment(
    reference_base_dir: Path,
    references_dir: Path,
    dry_run: bool = False,
    alignment_loader: Callable[[], LanlAlignment] | None = None,
    crf_loader: Callable[[], CrfReferenceCollection] | None = None,
    genbank_fetcher: Callable[[str], tuple[str, dict[str, tuple[int, int]] | None]] | None = None,
) -> LanlAlignmentUpdateResult:
    reference_base_dir.mkdir(parents=True, exist_ok=True)
    references_dir.mkdir(parents=True, exist_ok=True)

    alignment = (alignment_loader or download_latest_lanl_reference_alignment)()
    alignment_path = reference_base_dir / alignment.filename
    updated = save_lanl_alignment_file(reference_base_dir, alignment, dry_run=dry_run)

    alignment_records = parse_lanl_alignment_records(alignment.content)
    raw_crf_references = (crf_loader or load_lanl_crf_references)()
    crf_references = normalize_crf_references(raw_crf_references)
    desired_references = build_desired_reference_map(alignment_records, crf_references)
    fasta_result = reconcile_reference_fastas(
        desired_references.values(),
        references_dir,
        dry_run=dry_run,
        fetcher=genbank_fetcher,
    )

    return LanlAlignmentUpdateResult(
        year=alignment.year,
        path=alignment_path,
        updated=updated,
        fasta_records=len(alignment_records),
        fasta_files_added=fasta_result["files_added"],
        fasta_files_updated=fasta_result["files_updated"],
        fasta_files_removed=fasta_result["files_removed"],
        fasta_files_reclassified=fasta_result.get("files_reclassified", 0),
    )


def save_lanl_alignment_file(
    reference_base_dir: Path,
    alignment: LanlAlignment,
    dry_run: bool = False,
) -> bool:
    reference_base_dir.mkdir(parents=True, exist_ok=True)
    alignment_path = reference_base_dir / alignment.filename
    local_alignment = latest_local_lanl_alignment(reference_base_dir)

    current_content = ""
    if alignment_path.exists():
        current_content = alignment_path.read_text(encoding="utf-8")
    elif local_alignment is not None:
        current_content = local_alignment.read_text(encoding="utf-8")

    updated = current_content != alignment.content or not alignment_path.exists()
    if updated and not dry_run:
        if alignment_path.exists():
            backup_existing_file(alignment_path)
        elif local_alignment is not None:
            backup_existing_file(local_alignment)
        alignment_path.write_text(alignment.content, encoding="utf-8")

    return updated


def latest_local_lanl_alignment(reference_base_dir: Path) -> Path | None:
    candidates: list[tuple[int, Path]] = []
    for path in reference_base_dir.glob(f"{LANL_ALIGNMENT_PREFIX}_*_{LANL_ALIGNMENT_SUFFIX}"):
        match = LANL_ALIGNMENT_RE.match(path.name)
        if match:
            candidates.append((int(match.group(1)), path))
    if not candidates:
        return None
    return max(candidates, key=lambda item: item[0])[1]


def download_latest_lanl_reference_alignment() -> LanlAlignment:
    years = load_lanl_reference_alignment_years()
    if not years:
        raise RuntimeError("Could not find LANL HIV1 REF genome DNA alignment years.")

    year = max(years)
    content = download_lanl_reference_alignment(year)
    return LanlAlignment(year=year, content=content)


def load_lanl_reference_alignment_years() -> list[int]:
    query = urlencode(
        {
            "submit": "Retrieve",
            "server": "HIV",
            "ebolaIF": "",
            "col": "al_date",
            "where": "al_align_type:REF,al_organism:HIV1,al_region:GENOME,al_base_type:DNA",
        }
    )
    request = Request(
        f"{LANL_ALIGNMENT_LIST_URL}?{query}",
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/html,text/plain,*/*",
        },
    )
    with urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", errors="replace")
    return [
        int(value)
        for value in re.findall(r'value\s*=\s*["\']?(\d{4})["\']?', html, flags=re.I)
    ]


def download_lanl_reference_alignment(year: int) -> str:
    payload = urlencode(
        {
            "ALIGN_TYPE": "REF",
            "SUBORGANISM": "HIV1",
            "PRE_USER": "predefined",
            "REGION": "GENOME",
            "GENO_SUB": "ALLM",
            "BASETYPE": "DNA",
            "YEAR": str(year),
            "ORGANISM": "HIV",
            "FORMAT": "csv",
            "down_acc": "1",
            "submit": "Get Alignment",
        }
    ).encode("utf-8")
    request = Request(
        LANL_ALIGNMENT_DOWNLOAD_URL,
        data=payload,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/html,text/plain,*/*",
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    with urlopen(request, timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")
    return extract_lanl_alignment_content(html)


def extract_lanl_alignment_content(html: str) -> str:
    match = re.search(r"<pre[^>]*>(.*?)</pre>", html, flags=re.I | re.S)
    if not match:
        raise RuntimeError("Could not extract LANL alignment content from response.")

    content = unescape(match.group(1)).strip()
    if not content or "," not in content:
        raise RuntimeError("LANL alignment response did not contain CSV sequence rows.")
    return content + "\n"


def build_desired_reference_map(
    alignment_records: list[LanlReferenceRecord],
    crf_references: tuple[CrfReference, ...],
) -> dict[str, DesiredReference]:
    desired: dict[str, DesiredReference] = {}

    for reference in crf_references:
        accession = normalize_accession(reference.accession)
        add_desired_reference(
            desired,
            DesiredReference(
                accession=accession,
                subtype=reference.label,
                sequence=None,
                source="crf",
            ),
        )

    for record in alignment_records:
        accession = normalize_accession(record.accession)
        add_desired_reference(
            desired,
            DesiredReference(
                accession=accession,
                subtype=record.subtype,
                sequence=record.sequence,
                source="alignment",
            ),
        )

    return desired


def add_desired_reference(
    desired: dict[str, DesiredReference],
    reference: DesiredReference,
) -> None:
    accession = normalize_accession(reference.accession)
    existing = desired.get(accession)
    if existing is None:
        desired[accession] = DesiredReference(
            accession=accession,
            subtype=reference.subtype,
            sequence=reference.sequence,
            source=reference.source,
        )
        return

    if existing.source == reference.source:
        if existing.subtype != reference.subtype:
            if existing.source == "crf" and are_same_crf_number_aliases(
                existing.subtype,
                reference.subtype,
            ):
                subtype = preferred_crf_alias_label(existing.subtype, reference.subtype)
                warnings.warn(
                    f"CRF database aliases {existing.subtype} and {reference.subtype} "
                    f"refer to accession {accession}; using {subtype}.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                desired[accession] = DesiredReference(
                    accession=accession,
                    subtype=subtype,
                    sequence=existing.sequence or reference.sequence,
                    source=existing.source,
                )
                return
            raise ReconciliationError(
                f"Conflicting {reference.source} labels for accession {accession}: "
                f"{existing.subtype} and {reference.subtype}."
            )
        if existing.sequence is None and reference.sequence is not None:
            desired[accession] = DesiredReference(
                accession=accession,
                subtype=existing.subtype,
                sequence=reference.sequence,
                source=existing.source,
            )
        return

    if {existing.source, reference.source} == {"alignment", "crf"}:
        alignment_ref = existing if existing.source == "alignment" else reference
        crf_ref = existing if existing.source == "crf" else reference
        if alignment_ref.subtype != crf_ref.subtype:
            warnings.warn(
                f"CRF database label {crf_ref.subtype} overrides alignment "
                f"label {alignment_ref.subtype} for accession {accession}.",
                RuntimeWarning,
                stacklevel=2,
            )
        desired[accession] = DesiredReference(
            accession=accession,
            subtype=crf_ref.subtype,
            sequence=crf_ref.sequence or alignment_ref.sequence,
            source="crf",
        )
        return

    desired[accession] = preferred_desired_reference(existing, reference)


def preferred_desired_reference(
    first: DesiredReference,
    second: DesiredReference,
) -> DesiredReference:
    priority = {"crf": 3, "alignment": 2, "local": 1}
    if priority.get(second.source, 0) > priority.get(first.source, 0):
        return second
    return first


def are_same_crf_number_aliases(first: str, second: str) -> bool:
    first_number = crf_label_number(first)
    second_number = crf_label_number(second)
    return first_number is not None and first_number == second_number


def crf_label_number(label: str) -> str | None:
    match = re.match(r"^(\d{2,3})_[A-Za-z0-9]+$", str(label))
    if match is None:
        return None
    return match.group(1)


def preferred_crf_alias_label(first: str, second: str) -> str:
    return max(
        (str(first), str(second)),
        key=lambda label: (
            len(label),
            label.count("_"),
            label,
        ),
    )


def reconcile_reference_fastas(
    desired_references,
    references_dir: Path,
    dry_run: bool = False,
    output_dir: Path | None = None,
    fetcher: Callable[[str], tuple[str, dict[str, tuple[int, int]] | None]] | None = None,
    fetched_features: dict[str, dict[str, tuple[int, int]] | None] | None = None,
) -> dict[str, int]:
    desired_by_accession = {
        normalize_accession(reference.accession): reference
        for reference in desired_references
    }
    existing_by_accession = reference_fasta_paths_by_accession(references_dir)
    temp_dir = None
    if output_dir is None:
        temp_dir = tempfile.TemporaryDirectory(
            prefix="reference_fastas_",
            dir=references_dir.parent,
        )
        staged_references_dir = Path(temp_dir.name) / references_dir.name
    else:
        staged_references_dir = output_dir
    staged_references_dir.mkdir(parents=True, exist_ok=True)

    files_added = 0
    files_updated = 0
    files_removed = 0
    files_reclassified = 0

    for accession in sorted(desired_by_accession):
        reference = desired_by_accession[accession]
        fasta_path = references_dir / reference_fasta_name(accession, reference.subtype)
        staged_fasta_path = staged_references_dir / fasta_path.name
        local_paths = existing_by_accession.get(accession, [])
        sequence = reference.sequence
        if sequence is None:
            try:
                sequence = local_reference_sequence(fasta_path, local_paths)
            except LocalReferenceSequenceConflict:
                if fetcher is None:
                    raise ReconciliationError(
                        f"Conflicting local FASTA sequences for accession {accession}; "
                        "NCBI fetcher is required to choose an authoritative sequence."
                    )
                try:
                    sequence, features = fetcher(accession)
                except (GenBankRecordPartial, GenBankRecordUnavailable, HTTPError) as exc:
                    raise ReconciliationError(
                        f"Conflicting local FASTA sequences for accession {accession}; "
                        "could not fetch authoritative NCBI record."
                    ) from exc
                if fetched_features is not None:
                    fetched_features[accession] = features
        if sequence is None:
            continue

        stale_paths = [
            path
            for path in local_paths
            if path != fasta_path
        ]
        if stale_paths:
            files_removed += len(stale_paths)
            files_reclassified += 1

        fasta_exists = fasta_path.exists()
        fasta_matches = fasta_exists and reference_fasta_matches(
            fasta_path,
            accession,
            reference.subtype,
            sequence,
        )
        if not fasta_exists:
            if local_paths:
                files_updated += 1
            else:
                files_added += 1
        elif not fasta_matches:
            files_updated += 1

        write_reference_fasta(
            staged_fasta_path,
            accession,
            reference.subtype,
            sequence,
        )

    extra_paths = [
        path
        for accession, paths in existing_by_accession.items()
        if accession not in desired_by_accession
        for path in paths
    ]
    files_removed += len(extra_paths)
    validate_reference_fasta_dir(staged_references_dir)

    if output_dir is None and not dry_run:
        commit_staged_reference_fastas(staged_references_dir, references_dir)

    if temp_dir is not None:
        temp_dir.cleanup()

    return {
        "records": len(desired_by_accession),
        "files_added": files_added,
        "files_updated": files_updated,
        "files_removed": files_removed,
        "files_reclassified": files_reclassified,
    }


def reference_fasta_paths_by_accession(references_dir: Path) -> dict[str, list[Path]]:
    paths_by_accession: dict[str, list[Path]] = {}
    for path in sorted(references_dir.glob("*.fasta")):
        accession = parse_reference_name(path.stem)["accession"]
        if accession not in paths_by_accession:
            paths_by_accession[accession] = []
        paths_by_accession[accession].append(path)
    return paths_by_accession


def local_reference_sequence(final_path: Path, local_paths: list[Path]) -> str | None:
    candidates = [final_path] if final_path in local_paths else []
    candidates.extend(path for path in local_paths if path not in candidates)
    sequences: list[str] = []
    for path in candidates:
        try:
            sequences.append(parse_reference_fasta(path)["sequence"])
        except ValueError:
            continue
    if len(set(sequences)) > 1:
        raise LocalReferenceSequenceConflict(str(final_path))
    if sequences:
        return sequences[0]
    return None


@dataclass(frozen=True)
class LanlReferenceRecord:
    accession: str
    subtype: str
    sequence: str


@dataclass(frozen=True)
class DesiredReference:
    accession: str
    subtype: str
    sequence: str | None
    source: str


def parse_lanl_alignment_records(alignment_csv: str) -> list[LanlReferenceRecord]:
    records: list[LanlReferenceRecord] = []
    reader = csv.reader(StringIO(alignment_csv))
    for row in reader:
        if len(row) < 2:
            continue

        label = row[0].strip()
        sequence = row[1].strip()
        if not label or not sequence or label.lower() in {"name", "sequence"}:
            continue

        parsed = parse_reference_name(label)
        records.append(
            LanlReferenceRecord(
                accession=parsed["accession"],
                subtype=parsed["subtype"],
                sequence=sequence,
            )
        )

    return records


def load_lanl_crf_references() -> dict[str, tuple[CrfReference, ...]]:
    request = Request(
        LANL_CRF_URL,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "text/html,text/plain,*/*",
        },
    )
    with urlopen(request, timeout=60) as response:
        html = response.read().decode("utf-8", errors="replace")

    rows = re.findall(r"<tr\b[^>]*>(.*?)</tr>", html, flags=re.I | re.S)
    if not rows:
        rows = [html]

    found: dict[str, list[CrfReference]] = {}
    rejected_by_label: dict[str, set[str]] = {}
    current_label: str | None = None
    current_crf: str | None = None
    for row_html in rows:
        row_text = strip_html(row_html)
        if is_hiv2_text(row_text):
            current_label = None
            current_crf = None
            continue
        cells = [
            strip_html(cell)
            for cell in re.findall(r"<t[dh]\b[^>]*>(.*?)</t[dh]>", row_html, flags=re.I | re.S)
        ]
        match = re.search(r"\b(CRF(\d{2,3}_[A-Za-z0-9]+))\b", row_text)
        is_summary_row = (
            match is not None
            and any(cell.strip() in {match.group(1), match.group(2)} for cell in cells)
        )
        if is_summary_row:
            label = match.group(2)
            crf = match.group(1)
            current_label = label
            current_crf = crf
            rejected_by_label.setdefault(label, set()).update(
                subgenomic_accessions_from_text(row_text)
            )
            row_references = parse_lanl_crf_row(row_html)
        elif current_label is not None and current_crf is not None:
            label = current_label
            crf = current_crf
            rejected_by_label.setdefault(label, set()).update(
                subgenomic_accessions_from_text(row_text)
            )
            row_references = parse_lanl_crf_detail_text(label, crf, row_text)
        else:
            row_references = ()

        add_crf_references(found, row_references, rejected_by_label)

        if row_references:
            label = row_references[0].label
            rejected = rejected_by_label.get(label, set())
            if rejected and label in found:
                found[label] = [
                    reference
                    for reference in found[label]
                    if reference.accession not in rejected
                ]

    if not found:
        raise RuntimeError(
            "Could not extract CRF reference strains from the LANL CRF table."
        )
    return {
        label: tuple(references)
        for label, references in found.items()
    }


def add_crf_references(
    found: dict[str, list[CrfReference]],
    row_references: tuple[CrfReference, ...],
    rejected_by_label: dict[str, set[str]],
) -> None:
    if not row_references:
        return

    label = row_references[0].label
    rejected = rejected_by_label.get(label, set())
    incoming = [
        reference
        for reference in row_references
        if reference.accession not in rejected
    ]
    if not incoming:
        return

    references = found.setdefault(label, [])
    incoming_accessions = {reference.accession for reference in incoming}
    if len(incoming) > 1 and any(
        existing.accession in incoming_accessions for existing in references
    ):
        references[:] = dedupe_crf_references(
            incoming + [
                existing
                for existing in references
                if existing.accession not in incoming_accessions
            ]
        )
        return

    references[:] = dedupe_crf_references(references + incoming)


def parse_lanl_crf_row(row_html: str) -> tuple[CrfReference, ...]:
    cells = [
        strip_html(cell)
        for cell in re.findall(r"<t[dh]\b[^>]*>(.*?)</t[dh]>", row_html, flags=re.I | re.S)
    ]
    text = strip_html(row_html)
    match = re.search(r"\b(CRF(\d{2,3}_[A-Za-z0-9]+))\b", text)
    if not match:
        return ()

    crf = match.group(1)
    label = match.group(2)
    reference_text = reference_strain_text(cells, text, match)
    references = list(references_from_text(label, crf, reference_text))
    references.extend(genome_context_references_from_text(label, crf, text[match.end():]))
    references.extend(classification_context_references_from_text(label, crf, text[match.end():]))
    rejected_accessions = subgenomic_accessions_from_text(text)
    return tuple(
        reference
        for reference in dedupe_crf_references(references)
        if reference.accession not in rejected_accessions
    )


def parse_lanl_crf_detail_text(label: str, crf: str, text: str) -> tuple[CrfReference, ...]:
    if is_hiv2_text(text):
        return ()
    references = list(genome_context_references_from_text(label, crf, text))
    references.extend(classification_context_references_from_text(label, crf, text))
    rejected_accessions = subgenomic_accessions_from_text(text)
    return tuple(
        reference
        for reference in dedupe_crf_references(references)
        if reference.accession not in rejected_accessions
    )


def is_hiv2_text(text: str) -> bool:
    return re.search(r"\bHIV[- ]?2\b", text, flags=re.I) is not None


def reference_strain_text(cells: list[str], text: str, match: re.Match) -> str:
    for index, cell in enumerate(cells):
        if cell.strip() in {match.group(1), match.group(2)}:
            if index + 1 < len(cells):
                return cells[index + 1]
            break
    return ""


def references_from_text(label: str, crf: str, text: str) -> tuple[CrfReference, ...]:
    references: list[CrfReference] = []
    for reference_match in STRAIN_ACCESSION_RE.finditer(text):
        accession = normalize_accession(reference_match.group("accession"))
        references.append(
            CrfReference(
                label=label,
                crf=crf,
                strain=reference_match.group("strain"),
                accession=accession,
            )
        )

    if references:
        return tuple(references)

    accession_match = ACC_RE.search(text)
    if accession_match is None:
        return ()

    accession = normalize_accession(accession_match.group(1))
    return (
        CrfReference(
            label=label,
            crf=crf,
            strain=extract_reference_strain(text, accession),
            accession=accession,
        ),
    )


def genome_context_references_from_text(label: str, crf: str, text: str) -> tuple[CrfReference, ...]:
    references: list[CrfReference] = []
    for context in genome_reference_contexts(text):
        for accession in accessions_from_text(context):
            references.append(
                CrfReference(
                    label=label,
                    crf=crf,
                    strain="NA",
                    accession=accession,
                )
            )
    return tuple(dedupe_crf_references(references))


def classification_context_references_from_text(
    label: str,
    crf: str,
    text: str,
) -> tuple[CrfReference, ...]:
    references: list[CrfReference] = []
    for context in classification_reference_contexts(label, text):
        for accession in accessions_from_text(context):
            references.append(
                CrfReference(
                    label=label,
                    crf=crf,
                    strain="NA",
                    accession=accession,
                )
            )
    return tuple(dedupe_crf_references(references))


def genome_reference_contexts(text: str) -> list[str]:
    contexts: list[str] = []
    for segment in re.split(r"(?<=\.)\s+|\n+", text):
        normalized = segment.strip()
        if not normalized:
            continue
        for context in non_partial_genome_contexts(normalized):
            contexts.append(context)
    return contexts


def classification_reference_contexts(label: str, text: str) -> list[str]:
    contexts: list[str] = []
    label_pattern = re.escape(label)
    classification_re = re.compile(
        rf"\b(?:classif(?:y|ied)|reclassif(?:y|ied)|assigned|designated)\b"
        rf"[^.;]*\b(?:as|to)\s+(?:CRF)?{label_pattern}\b",
        flags=re.I,
    )
    for segment in re.split(r"(?<=[.;])\s+|\n+", text):
        normalized = segment.strip()
        if not normalized or not ACC_RE.search(normalized):
            continue
        if classification_re.search(normalized):
            contexts.append(normalized)
    return contexts


def subgenomic_accessions_from_text(text: str) -> set[str]:
    rejected: set[str] = set()
    for segment in re.split(r"(?<=[.;])\s+|\n+", text):
        normalized = segment.strip()
        if not normalized:
            continue
        if has_subgenomic_sequence_annotation(normalized.lower()):
            rejected.update(accessions_from_text(normalized))
    return rejected


def non_partial_genome_contexts(segment: str) -> list[str]:
    contexts: list[str] = []
    last_end = 0
    for match in partial_genome_annotation_re().finditer(segment.lower()):
        candidate = segment[last_end:match.start()].strip(" ,;:.")
        if candidate and has_valid_genome_context(candidate.lower()):
            contexts.append(candidate)
        last_end = match.end()

    candidate = segment[last_end:].strip(" ,;:.")
    if candidate and has_valid_genome_context(candidate.lower()):
        contexts.append(candidate)

    return contexts


def has_valid_genome_context(text: str) -> bool:
    valid_patterns = (
        r"\bcomplete genomes?\b",
        r"\bcomplete genome sequences?\b",
        r"\bcomplete sequences?\b",
        r"\bfull[- ]length genomes?\b",
        r"\bfull[- ]length genome sequences?\b",
        r"\bfull[- ]length sequences?\b",
        r"\bnear[- ]complete genomes?\b",
        r"\bnear[- ]complete genome sequences?\b",
        r"\bnear[- ]complete sequences?\b",
        r"\bnear[- ]full[- ]length genomes?\b",
        r"\bnear[- ]full[- ]length genome sequences?\b",
        r"\bnear[- ]full[- ]length sequences?\b",
        r"\bnear complete genomes?\b",
        r"\bnear complete genome sequences?\b",
        r"\bnear complete sequences?\b",
        r"\bnear full length genomes?\b",
        r"\bnear full length genome sequences?\b",
        r"\bnear full length sequences?\b",
        r"\bnflgs?\b",
        r"\bdefining genomes?\b",
        r"\bgenomes? used to (?:define|describe)\b",
        r"\bgenome sequences? used to (?:define|describe)\b",
        r"\bused to (?:define|describe) (?:the )?crf\b",
        r"\bbased on (?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve) (?:near[- ]complete |near[- ]full[- ]length |complete |full[- ]length )?(?:genomes?|genome sequences?|sequences?)\b",
        r"\bdescribed (?:using|based on) (?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve) (?:near[- ]complete |near[- ]full[- ]length |complete |full[- ]length )?(?:genomes?|genome sequences?|sequences?)\b",
        r"\bdefined (?:using|by|based on) (?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve) (?:near[- ]complete |near[- ]full[- ]length |complete |full[- ]length )?(?:genomes?|genome sequences?|sequences?)\b",
    )
    return any(re.search(pattern, text) for pattern in valid_patterns)


def accessions_from_text(text: str) -> tuple[str, ...]:
    accessions: list[str] = []
    for start, end in accession_ranges_from_text(text):
        accessions.extend(expand_accession_range(start, end))

    accessions.extend(
        normalize_accession(match.group(1))
        for match in ACC_RE.finditer(text)
        if not accession_is_parenthesized_reference_strain(text, match)
        and not accession_is_nonterminal_lanl_name_component(text, match)
    )
    return tuple(dedupe_accessions(accessions))


def accession_is_parenthesized_reference_strain(text: str, match: re.Match) -> bool:
    suffix = text[match.end():].lstrip()
    return re.match(rf"^\(\s*{ACC_PATTERN}\s*\)", suffix) is not None


def accession_is_nonterminal_lanl_name_component(text: str, match: re.Match) -> bool:
    start, end = match.span(1)
    token_start = start
    while token_start > 0 and re.match(r"[A-Za-z0-9_.]", text[token_start - 1]):
        token_start -= 1
    token_end = end
    while token_end < len(text) and re.match(r"[A-Za-z0-9_.]", text[token_end]):
        token_end += 1

    token = text[token_start:token_end]
    local_end = end - token_start
    if "." not in token or local_end >= len(token) or token[local_end] != ".":
        return False

    return ACC_RE.search(token[local_end + 1:]) is not None


def accession_ranges_from_text(text: str) -> tuple[tuple[str, str], ...]:
    range_re = re.compile(rf"\b({ACC_PATTERN})\s*(?:-|to)\s*({ACC_PATTERN})\b", flags=re.I)
    return tuple(
        (match.group(1), match.group(2))
        for match in range_re.finditer(text)
    )


def expand_accession_range(start: str, end: str) -> tuple[str, ...]:
    start = normalize_accession(start)
    end = normalize_accession(end)
    start_match = re.match(r"^([A-Z]{1,4}_?)(\d+)$", start)
    end_match = re.match(r"^([A-Z]{1,4}_?)(\d+)$", end)
    if not start_match or not end_match:
        return (start, end)
    if start_match.group(1) != end_match.group(1):
        return (start, end)

    start_number = int(start_match.group(2))
    end_number = int(end_match.group(2))
    if end_number < start_number or end_number - start_number > 500:
        return (start, end)

    prefix = start_match.group(1)
    width = len(start_match.group(2))
    return tuple(
        f"{prefix}{number:0{width}d}"
        for number in range(start_number, end_number + 1)
    )


def dedupe_crf_references(references: list[CrfReference]) -> list[CrfReference]:
    deduped: list[CrfReference] = []
    seen: set[str] = set()
    for reference in references:
        accession = normalize_accession(reference.accession)
        if accession in seen:
            continue
        seen.add(accession)
        deduped.append(
            CrfReference(
                label=reference.label,
                crf=reference.crf,
                strain=reference.strain,
                accession=accession,
            )
        )
    return deduped


def dedupe_accessions(accessions: list[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for accession in accessions:
        accession = normalize_accession(accession)
        if accession in seen:
            continue
        seen.add(accession)
        deduped.append(accession)
    return deduped


def strip_html(value: str) -> str:
    text = re.sub(r"<[^>]+>", " ", value)
    return " ".join(unescape(text).split())


def extract_reference_strain(text: str, accession: str) -> str:
    match = re.search(r"([^\s()]+)\s*\(" + re.escape(accession) + r"(?:\.\d+)?\)", text)
    if match:
        return match.group(1)
    return "NA"


def reference_accessions_in_dir(references_dir: Path) -> set[str]:
    return {
        parse_reference_name(path.stem)["accession"]
        for path in references_dir.glob("*.fasta")
    }


def normalize_crf_references(crf_references: CrfReferenceCollection) -> tuple[CrfReference, ...]:
    values = crf_references.values() if isinstance(crf_references, dict) else crf_references

    normalized: list[CrfReference] = []
    seen: set[str] = set()
    for value in values:
        if isinstance(value, CrfReference):
            references = (value,)
        else:
            references = tuple(value)

        for reference in references:
            key = reference_key(reference.accession, reference.label)
            if key in seen:
                continue
            seen.add(key)
            normalized.append(
                CrfReference(
                    label=reference.label,
                    crf=reference.crf,
                    strain=reference.strain,
                    accession=normalize_accession(reference.accession),
                )
            )

    return tuple(normalized)


def count_crfs_with_multiple_references(crf_references: CrfReferenceCollection) -> int:
    values = crf_references.values() if isinstance(crf_references, dict) else crf_references
    count = 0
    for value in values:
        if isinstance(value, CrfReference):
            continue
        if len(tuple(value)) > 1:
            count += 1
    return count


def reference_key(accession: str, subtype: str) -> str:
    return f"{normalize_accession(accession)}|{str(subtype)}"


def normalize_accession(accession: str) -> str:
    return str(accession).split(".")[0]


def reference_fasta_name(accession: str, subtype: str) -> str:
    safe_accession = sanitize_filename_part(accession)
    safe_subtype = sanitize_filename_part(subtype)
    return f"{safe_accession}-{safe_subtype}.fasta"


def sanitize_filename_part(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value).strip())
    return text.strip("._-") or "NA"


def write_reference_fasta(path: Path, accession: str, subtype: str, sequence: str) -> None:
    path.write_text(reference_fasta_content(accession, subtype, sequence), encoding="utf-8")


def reference_fasta_matches(path: Path, accession: str, subtype: str, sequence: str) -> bool:
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        return False

    record = records[0]
    expected_id = f"{accession}-{subtype}"
    return record.id == expected_id and normalize_sequence(str(record.seq)) == normalize_sequence(sequence)


def reference_fasta_content(accession: str, subtype: str, sequence: str) -> str:
    sequence = normalize_sequence(sequence)
    lines = [f">{accession}-{subtype}"]
    lines.extend(sequence[start:start + 80] for start in range(0, len(sequence), 80))
    return "\n".join(lines) + "\n"


def fetch_genbank_sequence_and_features(accession: str) -> tuple[str, dict[str, tuple[int, int]] | None]:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gb",
        "retmode": "text",
        "tool": NCBI_TOOL_NAME,
    }
    if Entrez.email:
        params["email"] = Entrez.email
    if getattr(Entrez, "api_key", None):
        params["api_key"] = Entrez.api_key

    handle = Entrez.efetch(**params)
    try:
        try:
            record = SeqIO.read(handle, "genbank")
        except ValueError as exc:
            if str(exc) == "No records found in handle":
                raise GenBankRecordUnavailable(accession) from exc
            raise
    finally:
        handle.close()

    sequence = normalize_sequence(str(record.seq))
    validate_hiv1_complete_or_near_complete(record, sequence)
    features = get_gene_locations(record.features)
    return sequence, features or None


def validate_hiv1_complete_or_near_complete(record, sequence: str) -> None:
    sequence_length = len(normalize_sequence(sequence))
    if sequence_length < HIV1_NEAR_COMPLETE_MIN_LENGTH:
        raise GenBankRecordPartial(
            f"{getattr(record, 'id', 'unknown')} length {sequence_length}"
        )


def has_subgenomic_sequence_annotation(text: str) -> bool:
    subgenomic_patterns = (
        r"\bpartial cds\b",
        r"\bfragments?\b",
        r"\bgag[- ]pol\b",
        r"\bpol gene\b",
        r"\benv gene\b",
        r"\bgag gene\b",
        r"\bnef gene\b",
        r"\bpol sequences?\b",
        r"\benv sequences?\b",
        r"\bgag sequences?\b",
        r"\bnef sequences?\b",
        r"\bpol only\b",
        r"\benv only\b",
        r"\bgag only\b",
        r"\bnef only\b",
    )
    return re.search("|".join(f"(?:{pattern})" for pattern in subgenomic_patterns), text) is not None


def partial_genome_annotation_re() -> re.Pattern:
    partial_patterns = (
        r"\bpartial genomes?\b",
        r"\bpartial sequences?\b",
        r"\bpartial cds\b",
        r"\bpartial genomic sequence\b",
        r"\bfragments?\b",
        r"\bgag[- ]pol\b",
        r"\bpol gene\b",
        r"\benv gene\b",
        r"\bgag gene\b",
        r"\bnef gene\b",
        r"\bpol only\b",
        r"\benv only\b",
        r"\bgag only\b",
        r"\bnef only\b",
    )
    return re.compile("|".join(f"(?:{pattern})" for pattern in partial_patterns))


def get_gene_locations(features) -> dict[str, tuple[int, int]]:
    gene_ranges: dict[str, tuple[int, int]] = {}

    for feature in features:
        if feature.type != "gene":
            continue
        gene_name = clean_feature_name(feature.qualifiers.get("gene", ["location"])[0])
        gene_ranges[gene_name] = feature_bounds(feature)

    if gene_ranges:
        return gene_ranges

    for feature in features:
        if feature.type != "CDS":
            continue
        raw_name = (
            feature.qualifiers.get("gene", [None])[0]
            or feature.qualifiers.get("product", [None])[0]
            or feature.qualifiers.get("protein_id", ["CDS_unknown"])[0]
        )
        gene_name = clean_feature_name(raw_name)
        if gene_name not in gene_ranges:
            gene_ranges[gene_name] = feature_bounds(feature)

    return gene_ranges


def clean_feature_name(value) -> str:
    return str(value or "location").replace("/", "_").lower()


def feature_bounds(feature) -> tuple[int, int]:
    return int(feature.location.start) + 1, int(feature.location.end)


def rebuild_sequences_with_locations(
    references_dir: Path,
    sequences_path: Path,
    fetched_features: dict[str, dict[str, tuple[int, int]] | None] | None = None,
    refresh_features: bool = False,
    dry_run: bool = False,
    fetcher: Callable[[str], tuple[str, dict[str, tuple[int, int]] | None]] | None = None,
    sleep_seconds: float = 0.34,
    output_path: Path | None = None,
) -> dict[str, int]:
    fetched_features = fetched_features or {}
    existing_table = read_existing_table(sequences_path)
    columns = list(existing_table.columns) if not existing_table.empty else REFERENCE_TABLE_COLUMNS
    existing_by_accession = index_existing_reference_rows(existing_table)
    duplicate_rows_removed = duplicate_row_count(existing_table)

    rows: list[dict[str, str]] = []
    added_rows = 0
    updated_rows = 0
    updated_feature_rows = 0
    used_accessions: set[str] = set()

    for fasta_path in sorted(references_dir.glob("*.fasta")):
        fasta_info = parse_reference_fasta(fasta_path)
        accession = fasta_info["accession"]
        if accession in used_accessions:
            continue
        used_accessions.add(accession)

        original_row = choose_existing_reference_row(
            existing_by_accession.get(accession, []),
            fasta_info["subtype"],
        )

        if original_row is not None:
            row = original_row.copy()
        else:
            row = default_reference_row(fasta_info)
            added_rows += 1

        row["accession"] = accession
        row["subtype"] = fasta_info["subtype"]
        row["align sequence"] = fasta_info["align sequence"]
        row["sequence"] = fasta_info["sequence"]

        if "group" in columns:
            original_subtype = str(original_row.get("subtype", "")) if original_row else ""
            if not row.get("group") or original_subtype != row["subtype"]:
                row["group"] = infer_group(row["subtype"])

        if needs_feature_update(row.get("features"), refresh_features):
            features = fetched_features.get(accession)
            if features is None and refresh_features and fetcher is not None:
                try:
                    _, features = fetcher(accession)
                except (GenBankRecordUnavailable, GenBankRecordPartial):
                    features = None
                time.sleep(sleep_seconds)
            if refresh_features or accession in fetched_features or original_row is None:
                row["features"] = serialize_features(features)
                if original_row is not None:
                    updated_feature_rows += 1

        row = {column: serialize_cell(row.get(column, "")) for column in columns}
        if original_row is not None:
            original = {column: serialize_cell(original_row.get(column, "")) for column in columns}
            if row != original:
                updated_rows += 1
        rows.append(row)

    if not dry_run:
        new_table = pd.DataFrame(rows, columns=columns)
        validate_final_reference_dataset(new_table, references_dir)
        target_path = output_path or sequences_path
        target_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path is None:
            backup_existing_file(sequences_path)
        new_table.to_csv(target_path, sep="\t", index=False)

    return {
        "added_rows": added_rows,
        "updated_rows": updated_rows,
        "updated_feature_rows": updated_feature_rows,
        "duplicate_rows_removed": duplicate_rows_removed,
        "total_references": len(used_accessions),
    }


def index_existing_reference_rows(
    existing_table: pd.DataFrame,
) -> dict[str, list[dict[str, str]]]:
    by_accession: dict[str, list[dict[str, str]]] = {}
    if existing_table.empty:
        return by_accession

    for row in existing_table.to_dict(orient="records"):
        accession = normalize_accession(serialize_cell(row.get("accession", "")))
        if not accession:
            continue

        if accession not in by_accession:
            by_accession[accession] = []
        by_accession[accession].append(row)

    return by_accession


def choose_existing_reference_row(
    rows: list[dict[str, str]],
    current_subtype: str,
) -> dict[str, str] | None:
    if not rows:
        return None

    return max(
        rows,
        key=lambda row: (
            serialize_cell(row.get("subtype", "")) == current_subtype,
            useful_metadata_score(row),
        ),
    )


def useful_metadata_score(row: dict[str, str]) -> int:
    score = 0
    for key, value in row.items():
        if key in {"subtype", "accession", "align sequence", "sequence"}:
            continue
        text = serialize_cell(value).strip().lower()
        if text and text not in {"na", "none", "nan", "null"}:
            score += 1
    return score


def duplicate_row_count(existing_table: pd.DataFrame) -> int:
    if existing_table.empty or "accession" not in existing_table.columns:
        return 0
    accessions = existing_table["accession"].astype(str).map(normalize_accession)
    return int(accessions.duplicated().sum())


def validate_final_reference_dataset(table: pd.DataFrame, references_dir: Path) -> None:
    assert not table["accession"].duplicated().any()
    fasta_accessions = [
        parse_reference_name(path.stem)["accession"]
        for path in references_dir.glob("*.fasta")
    ]
    unique_accessions = set(table["accession"].astype(str))
    if len(fasta_accessions) != len(table.index) or len(unique_accessions) != len(table.index):
        raise AssertionError(
            "Reference dataset validation failed: FASTA count, TSV rows, and unique accessions differ."
        )
    if set(fasta_accessions) != unique_accessions:
        raise AssertionError(
            "Reference dataset validation failed: FASTA accessions and TSV accessions differ."
        )


def validate_reference_fasta_dir(references_dir: Path) -> None:
    accessions: list[str] = []
    for path in references_dir.glob("*.fasta"):
        fasta_info = parse_reference_fasta(path)
        accessions.append(fasta_info["accession"])
    duplicated = [accession for accession, count in Counter(accessions).items() if count > 1]
    if duplicated:
        raise AssertionError(
            "Reference FASTA validation failed: duplicated accessions "
            + ", ".join(sorted(duplicated))
        )


def commit_staged_reference_dataset(
    staged_references_dir: Path,
    references_dir: Path,
    staged_sequences_path: Path,
    sequences_path: Path,
) -> None:
    table = pd.read_csv(
        staged_sequences_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    validate_reference_fasta_dir(staged_references_dir)
    validate_final_reference_dataset(table, staged_references_dir)
    commit_staged_reference_fastas(staged_references_dir, references_dir)
    backup_existing_file(sequences_path)
    shutil.copy2(staged_sequences_path, sequences_path)


def commit_staged_reference_fastas(
    staged_references_dir: Path,
    references_dir: Path,
) -> None:
    validate_reference_fasta_dir(staged_references_dir)
    backup_existing_dir(references_dir)
    if references_dir.exists():
        shutil.rmtree(references_dir)
    shutil.copytree(staged_references_dir, references_dir)


def read_existing_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=REFERENCE_TABLE_COLUMNS)
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def backup_existing_file(path: Path) -> None:
    if not path.exists():
        return
    backup_path = path.with_suffix(path.suffix + ".bak")
    shutil.copy2(path, backup_path)


def backup_existing_dir(path: Path) -> None:
    if not path.exists():
        return
    backup_path = unique_backup_path(path)
    shutil.copytree(path, backup_path)


def unique_backup_path(path: Path) -> Path:
    base = path.with_name(path.name + ".bak")
    if not base.exists():
        return base
    timestamp = int(time.time())
    candidate = path.with_name(f"{path.name}.bak.{timestamp}")
    suffix = 1
    while candidate.exists():
        suffix += 1
        candidate = path.with_name(f"{path.name}.bak.{timestamp}.{suffix}")
    return candidate


def parse_reference_fasta(path: Path) -> dict[str, str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}.")

    record = records[0]
    parsed_name = parse_reference_name(path.stem)
    sequence = normalize_sequence(str(record.seq))
    aligned = normalize_sequence(str(record.seq), keep_gaps=True)

    return {
        "subtype": parsed_name["subtype"],
        "accession": parsed_name["accession"],
        "align sequence": aligned,
        "sequence": sequence,
    }


def parse_reference_name(stem: str) -> dict[str, str]:
    accession, separator, subtype = stem.partition("-")
    if not separator:
        parts = stem.split(".")
        accession = parts[-1]
        subtype = parts[1] if len(parts) > 1 and parts[0].lower() == "ref" else "Unknown"
    return {
        "accession": accession.split(".")[0],
        "subtype": subtype or "Unknown",
    }


def default_reference_row(fasta_info: dict[str, str]) -> dict[str, str]:
    subtype = fasta_info["subtype"]
    return {
        "subtype": subtype,
        "country": "NA",
        "year": "NA",
        "isolate": "NA",
        "accession": fasta_info["accession"],
        "align sequence": fasta_info["align sequence"],
        "sequence": fasta_info["sequence"],
        "features": "None",
        "group": infer_group(subtype),
    }


def infer_group(subtype: str) -> str:
    normalized = str(subtype).strip().upper()
    if normalized in {"N", "O", "P"}:
        return normalized
    return "M"


def normalize_sequence(sequence: str, keep_gaps: bool = False) -> str:
    sequence = str(sequence).upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if keep_gaps:
        return sequence.replace(".", "-")
    return sequence.replace("-", "").replace(".", "")


def needs_feature_update(features, refresh_features: bool) -> bool:
    if refresh_features:
        return True
    text = str(features or "").strip().lower()
    return text in {"", "none", "nan", "null"}


def serialize_features(features: dict[str, tuple[int, int]] | None) -> str:
    if not features:
        return "None"
    return str(features)


def serialize_cell(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    return str(value)

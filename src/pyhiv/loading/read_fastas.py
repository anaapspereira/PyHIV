import csv
import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


SUPPORTED_FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".ffn")
INPUT_QC_ALLOWED_NT_ALPHABET = set("ACGTRYSWKMBDHVN-")
INPUT_QC_TSV_COLUMNS = [
    "file_name",
    "sequence_id",
    "sequence_length",
    "input_qc_status",
    "invalid_symbols",
    "invalid_symbol_count",
    "invalid_positions",
    "invalid_runs",
    "invalid_run_count",
    "invalid_run_length_class",
]


@dataclass(frozen=True)
class InputSequenceQC:
    file_name: str
    sequence_id: str
    sequence_length: int
    input_qc_status: str
    invalid_symbols: str
    invalid_symbol_count: int
    invalid_positions: str
    invalid_runs: str
    invalid_run_count: int
    invalid_run_length_class: str

    def to_tsv_row(self) -> dict[str, object]:
        return asdict(self)


def discover_fasta_files(input_folder: Path, exclude_dirs: Iterable[Path] = ()) -> List[Path]:
    """
    Recursively find FASTA files with supported extensions under input_folder.

    Parameters
    ----------
    input_folder : Path
        Directory to search, including its subdirectories.
    exclude_dirs : Iterable[Path], optional
        Directories to skip entirely, e.g. a run's own output directory, so
        that previously generated alignment or gene-region FASTA files
        aren't picked up again as new input on a later run.

    Returns
    -------
    List[Path]
        Sorted list of matching FASTA file paths.
    """
    resolved_excludes = [Path(directory).resolve() for directory in exclude_dirs]

    found = []
    for path in input_folder.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in SUPPORTED_FASTA_EXTENSIONS:
            continue
        resolved = path.resolve()
        if any(excluded in resolved.parents for excluded in resolved_excludes):
            continue
        found.append(path)

    return sorted(found)


def read_input_fastas(input_folder: Path, exclude_dirs: Iterable[Path] = ()) -> List[SeqRecord]:
    """
    Reads nucleotide FASTA files (.fasta, .fa, .fna, .ffn) from a specified input
    folder and its subdirectories.

    Parameters
    ----------
    input_folder : Path
        Path to the folder containing the FASTA files.
    exclude_dirs : Iterable[Path], optional
        Directories to skip, e.g. this run's own output directory, so that
        previously generated result FASTAs aren't re-read as new input.

    Returns
    -------
    List[SeqRecord]
        A list of BioPython SeqRecord objects containing sequence IDs and sequences.
        Each record's annotations include "source_file" (the file's path relative
        to input_folder, using "/" separators) and "source_path" (its full path).

    Raises
    ------
    NotADirectoryError
        If the input folder does not exist or is not a directory.
    """
    if not input_folder.is_dir():
        raise NotADirectoryError(f"Input folder {input_folder} is not a directory.")

    fasta_files = discover_fasta_files(input_folder, exclude_dirs=exclude_dirs)

    if not fasta_files:
        logging.warning(f"No FASTA files with supported extensions found in {input_folder}.")

    sequences = []
    for fasta_file in fasta_files:
        try:
            with open(fasta_file, "r") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            if not records:
                logging.warning(f"File {fasta_file} contains no valid sequences.")
            else:
                relative_source = fasta_file.relative_to(input_folder).as_posix()
                for record in records:
                    record.annotations["source_file"] = relative_source
                    record.annotations["source_path"] = str(fasta_file)
                sequences.extend(records)
                logging.info(f"Successfully read {len(records)} sequences from {fasta_file}")
        except Exception as e:
            logging.error(f"Error reading {fasta_file}: {e}")

    return sequences


def validate_input_sequence(
    sequence: str,
    sequence_id: str,
    file_name: str = "",
) -> InputSequenceQC:
    """Validate original nucleotide input without changing pipeline behavior."""
    cleaned = "".join(str(sequence).upper().split())
    invalid_positions: list[int] = []
    invalid_symbols_by_position: list[str] = []
    run_parts: list[str] = []
    run_lengths: list[int] = []
    run_start: int | None = None
    run_symbols: list[str] = []

    def finish_run() -> None:
        nonlocal run_start, run_symbols
        if run_start is None:
            return
        run_end = run_start + len(run_symbols) - 1
        run_text = "".join(run_symbols)
        run_parts.append(f"{run_start}-{run_end}:{run_text}")
        run_lengths.append(len(run_symbols))
        run_start = None
        run_symbols = []

    for index, symbol in enumerate(cleaned, start=1):
        if symbol in INPUT_QC_ALLOWED_NT_ALPHABET:
            finish_run()
            continue
        invalid_positions.append(index)
        invalid_symbols_by_position.append(symbol)
        if run_start is None:
            run_start = index
            run_symbols = [symbol]
        else:
            run_symbols.append(symbol)
    finish_run()

    if not invalid_positions:
        length_class = "none"
    elif any(length >= 3 for length in run_lengths) and any(length < 3 for length in run_lengths):
        length_class = "mixed"
    elif any(length >= 3 for length in run_lengths):
        length_class = "triple_or_longer"
    else:
        length_class = "short_only"

    return InputSequenceQC(
        file_name=file_name,
        sequence_id=sequence_id,
        sequence_length=len(cleaned),
        input_qc_status="WARN" if invalid_positions else "PASS",
        invalid_symbols="".join(sorted(set(invalid_symbols_by_position))),
        invalid_symbol_count=len(invalid_positions),
        invalid_positions=",".join(str(position) for position in invalid_positions),
        invalid_runs=";".join(run_parts),
        invalid_run_count=len(run_parts),
        invalid_run_length_class=length_class,
    )


def input_qc_rows_for_records(records: Iterable[SeqRecord]) -> list[InputSequenceQC]:
    """Return generic input QC rows for already-read FASTA records."""
    return [
        validate_input_sequence(
            sequence=str(record.seq),
            sequence_id=record.id,
            file_name=str(record.annotations.get("source_file", "")),
        )
        for record in records
    ]


def write_input_qc_tsv(rows: Iterable[InputSequenceQC], output_path: Path) -> None:
    """Write one generic input QC row per original sequence."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=INPUT_QC_TSV_COLUMNS, delimiter="	")
        writer.writeheader()
        for row in rows:
            tsv_row = row.to_tsv_row()
            writer.writerow({column: tsv_row[column] for column in INPUT_QC_TSV_COLUMNS})

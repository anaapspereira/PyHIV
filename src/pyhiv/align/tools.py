from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Dict, Tuple

from pyhiv.align.famsa import pyfamsa_align

try:
    from Bio.SeqRecord import SeqRecord
except ImportError:  # pragma: no cover
    raise ImportError(
        "BioPython is required for this module. "
        "Please install it via 'pip install biopython'."
    )


EDLIB_HW = "edlib-HW"
MAFFT = "MAFFT"
PARASAIL_NW = "parasail-NW"
PYFAMSA = "PyFamsa"

# Default backend used when the caller does not explicitly select one.
DEFAULT_ALIGNMENT_TOOL = EDLIB_HW

SUPPORTED_ALIGNMENT_TOOLS = (
    EDLIB_HW,
    MAFFT,
    PARASAIL_NW,
    PYFAMSA,
)

ALIGNMENT_TOOL_ALIASES = {
    "parasail": PARASAIL_NW,
}

ALIGNMENT_TOOL_CHOICES = SUPPORTED_ALIGNMENT_TOOLS + tuple(ALIGNMENT_TOOL_ALIASES)

MAFFT_BIN_ENV = "PYHIV_MAFFT_BIN"

DNA_ALPHABET = "ACGTRYSWKMBDHVN"


def normalize_alignment_tool(alignment_tool: str | None) -> str:
    """Return the canonical alignment tool name."""
    if alignment_tool is None:
        return DEFAULT_ALIGNMENT_TOOL

    for supported_tool in SUPPORTED_ALIGNMENT_TOOLS:
        if alignment_tool.lower() == supported_tool.lower():
            return supported_tool

    alias = ALIGNMENT_TOOL_ALIASES.get(alignment_tool.lower())
    if alias is not None:
        return alias

    supported = ", ".join(ALIGNMENT_TOOL_CHOICES)
    raise ValueError(
        f"Unsupported alignment tool '{alignment_tool}'. "
        f"Supported tools: {supported}."
    )


def align_sequences(
    test_seq: SeqRecord,
    ref_seq: SeqRecord,
    alignment_tool: str | None = DEFAULT_ALIGNMENT_TOOL,
) -> Tuple[str, str]:
    """
    Align a test sequence to a reference sequence using the selected tool.

    Recommended modes:

    - edlib-HW:
        Use when the test sequence is a small region that should be located
        inside a longer reference sequence. Example: HIV-1 protease against
        a full HIV-1 genome.

    Returns:
        Tuple of:
            query_aligned, reference_aligned
    """
    tool = normalize_alignment_tool(alignment_tool)
    aligner = _ALIGNERS[tool]
    return aligner(test_seq, ref_seq)


def validate_alignment_tool_available(alignment_tool: str | None = DEFAULT_ALIGNMENT_TOOL) -> str:
    """
    Validate that the selected alignment backend is usable in this environment.

    Worker-level alignment errors are logged per reference and skipped, so
    missing optional backends must be detected before the pipeline starts.
    """
    tool = normalize_alignment_tool(alignment_tool)

    if tool == EDLIB_HW:
        try:
            import edlib  # noqa: F401
        except ImportError as exc:  # pragma: no cover
            raise RuntimeError(
                "edlib is required for alignment_tool='edlib-HW'. "
                "Install PyHIV dependencies with 'pip install pyhiv-tools' "
                "or install edlib directly with 'pip install edlib'."
            ) from exc

    elif tool == PARASAIL_NW:
        try:
            import parasail  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "parasail is required for alignment_tool='parasail-NW'. "
                "Install it with 'pip install pyhiv-tools[parasail]' "
                "or 'pip install parasail'."
            ) from exc

    elif tool == MAFFT and get_mafft_binary() is None:
        raise RuntimeError(
            "MAFFT support requires a MAFFT executable. "
            f"Install MAFFT and add it to PATH, or set {MAFFT_BIN_ENV}."
        )

    elif tool == PYFAMSA:
        try:
            import pyfamsa  # noqa: F401
        except ImportError as exc:
            raise RuntimeError(
                "PyFamsa is required for alignment_tool='PyFamsa'. "
                "Install it with 'pip install pyhiv-tools[famsa]' "
                "or 'pip install pyfamsa'."
            ) from exc

    return tool


def edlib_hw_align(test_seq: SeqRecord, ref_seq: SeqRecord) -> Tuple[str, str]:
    """
    Align a shorter query sequence inside a longer reference using edlib HW mode.

    This does three things:

    1. Uses edlib HW mode to locate the full query inside the reference.
    2. Re-aligns the query against only that matched reference window using NW.
    3. Projects the query back onto full reference coordinates by padding with
       '-' before and after the matched region.

    The returned alignment has exactly len(ref_seq) columns. This means query
    insertion columns relative to the reference are removed, because they do not
    have reference coordinates.
    """
    try:
        import edlib
    except ImportError as exc:  # pragma: no cover
        raise ImportError("edlib is required for alignment_tool='edlib-HW'.") from exc

    test_sequence = _sequence_to_string(test_seq)
    ref_sequence = _sequence_to_string(ref_seq)

    if not test_sequence:
        raise ValueError("Cannot align an empty test sequence.")

    if not ref_sequence:
        raise ValueError("Cannot align against an empty reference sequence.")

    # Step 1: locate query inside the full reference.
    location_result = edlib.align(
        test_sequence,
        ref_sequence,
        mode="HW",
        task="path",
    )

    if location_result["editDistance"] < 0 or not location_result.get("locations"):
        raise ValueError("Could not align test sequence to reference sequence.")

    start, end = _pick_edlib_location(
        locations=location_result["locations"],
        query_length=len(test_sequence),
    )

    if start < 0 or end < start or end >= len(ref_sequence):
        raise ValueError(
            f"Invalid edlib-HW location {(start, end)} for reference length "
            f"{len(ref_sequence)}."
        )

    ref_window = ref_sequence[start:end + 1]

    # Step 2: clean global alignment only inside the located reference window.
    window_result = edlib.align(
        test_sequence,
        ref_window,
        mode="NW",
        task="path",
    )

    if window_result["editDistance"] < 0:
        raise ValueError("Could not align test sequence to located reference window.")

    window_alignment = edlib.getNiceAlignment(
        window_result,
        test_sequence,
        ref_window,
    )

    query_mid = window_alignment["query_aligned"]
    ref_mid = window_alignment["target_aligned"]

    # Step 3: remove query insertions so the final result stays in reference
    # coordinates and has exactly len(ref_sequence) columns.
    query_mid, ref_mid = _drop_query_insertions(query_mid, ref_mid)

    # _drop_query_insertions always keeps query_mid and ref_mid the same
    # length, and that length always equals the reference window length, so
    # these checks can never fail in practice. Kept as a guard against a
    # future change breaking that invariant.
    expected_window_length = end - start + 1
    if len(query_mid) != expected_window_length:  # pragma: no cover
        raise AssertionError(
            f"Projected query window length {len(query_mid)} does not match "
            f"reference window length {expected_window_length}."
        )

    if len(ref_mid) != expected_window_length:  # pragma: no cover
        raise AssertionError(
            f"Projected reference window length {len(ref_mid)} does not match "
            f"reference window length {expected_window_length}."
        )

    query_aligned = (
        "-" * start
        + query_mid
        + "-" * (len(ref_sequence) - end - 1)
    )

    ref_aligned = (
        ref_sequence[:start]
        + ref_mid
        + ref_sequence[end + 1:]
    )

    # Given the checks above, query_aligned and ref_aligned are always padded
    # to exactly len(ref_sequence), so these can never fail in practice. Kept
    # as a guard against a future change breaking that invariant.
    if len(query_aligned) != len(ref_sequence):  # pragma: no cover
        raise AssertionError(
            f"Projected query alignment length {len(query_aligned)} does not match "
            f"reference length {len(ref_sequence)}."
        )

    if len(ref_aligned) != len(ref_sequence):  # pragma: no cover
        raise AssertionError(
            f"Projected reference alignment length {len(ref_aligned)} does not match "
            f"reference length {len(ref_sequence)}."
        )

    logging.debug(
        "edlib-HW placed query '%s' on reference '%s' at %s-%s with editDistance=%s.",
        getattr(test_seq, "id", "<unknown>"),
        getattr(ref_seq, "id", "<unknown>"),
        start,
        end,
        location_result["editDistance"],
    )

    return query_aligned, ref_aligned


def parasail_nw_align(test_seq: SeqRecord, ref_seq: SeqRecord) -> Tuple[str, str]:
    """Align two sequences using parasail Needleman-Wunsch with traceback."""
    try:
        import parasail
    except ImportError as exc:  # pragma: no cover
        raise ImportError("parasail is required for alignment_tool='parasail-NW'.") from exc

    test_sequence = _sequence_to_string(test_seq)
    ref_sequence = _sequence_to_string(ref_seq)

    matrix = parasail.matrix_create(DNA_ALPHABET, 2, -1)

    result = parasail.nw_trace_scan_16(
        test_sequence,
        ref_sequence,
        5,
        1,
        matrix,
    )

    cigar = result.cigar.decode.decode("ascii")
    return _apply_cigar(test_sequence, ref_sequence, cigar)


def mafft_align(test_seq: SeqRecord, ref_seq: SeqRecord) -> Tuple[str, str]:
    """Align two sequences using the external MAFFT executable."""
    mafft_bin = get_mafft_binary()
    if mafft_bin is None:
        raise RuntimeError(
            "MAFFT support requires a MAFFT executable. "
            f"Set {MAFFT_BIN_ENV} or add 'mafft' to PATH."
        )

    test_sequence = _sequence_to_string(test_seq)
    ref_sequence = _sequence_to_string(ref_seq)

    if not test_sequence:
        raise ValueError("Cannot align an empty test sequence.")

    if not ref_sequence:
        raise ValueError("Cannot align against an empty reference sequence.")

    with tempfile.TemporaryDirectory(prefix="pyhiv-mafft-") as tmp_dir:
        tmp_path = Path(tmp_dir)
        input_fasta = tmp_path / "input.fasta"

        _write_pair_fasta(input_fasta, "query", test_sequence, "reference", ref_sequence)

        command = [
            mafft_bin,
            "--quiet",
            "--thread",
            "1",
            "--auto",
            str(input_fasta),
        ]

        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            output = (result.stderr or result.stdout).strip()
            raise RuntimeError(
                "MAFFT failed with exit code "
                f"{result.returncode}: {output}"
            )

        return _parse_pair_fasta_alignment(result.stdout)


def get_mafft_binary() -> str | None:
    """Resolve MAFFT from an environment override or PATH."""
    env_path = os.environ.get(MAFFT_BIN_ENV)
    if env_path and _is_executable_file(Path(env_path)):
        return env_path

    path_binary = shutil.which("mafft")
    if path_binary:
        return path_binary

    return None


def _is_executable_file(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def _write_pair_fasta(
    path: Path,
    first_id: str,
    first_sequence: str,
    second_id: str,
    second_sequence: str,
) -> None:
    path.write_text(
        f">{first_id}\n{first_sequence}\n"
        f">{second_id}\n{second_sequence}\n"
    )


def _parse_pair_fasta_alignment(text: str) -> Tuple[str, str]:
    records: list[str] = []
    current: list[str] | None = None

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current is not None:
                records.append("".join(current))
            current = []
            continue

        if current is None:
            raise ValueError("Could not parse MAFFT FASTA alignment output.")

        current.append(line.upper())

    if current is not None:
        records.append("".join(current))

    if len(records) != 2 or not records[0] or not records[1]:
        raise ValueError("Could not parse MAFFT FASTA alignment output.")

    if len(records[0]) != len(records[1]):
        raise ValueError("Parsed MAFFT alignment has different sequence lengths.")

    return records[0], records[1]


def _sequence_to_string(seq_record: SeqRecord) -> str:
    """
    Convert a SeqRecord sequence to an uppercase string.

    This intentionally does not remove '-' characters, because doing so would
    silently change coordinates if someone passes a pre-gapped reference.
    Input sequences should normally be ungapped before alignment.
    """
    return str(seq_record.seq).upper()


def _pick_edlib_location(
    locations: list[tuple[int, int]],
    query_length: int,
) -> tuple[int, int]:
    """
    Pick the best edlib location from a list of equally optimal locations.

    Edlib may return multiple equally good placements. Prefer the placement
    whose reference window length is closest to the query length, then the
    earliest start coordinate.
    """
    valid_locations = [
        (start, end)
        for start, end in locations
        if start is not None and end is not None and start >= 0 and end >= start
    ]

    if not valid_locations:
        raise ValueError(f"No valid edlib locations found: {locations!r}")

    return min(
        valid_locations,
        key=lambda loc: (
            abs((loc[1] - loc[0] + 1) - query_length),
            loc[0],
            loc[1],
        ),
    )


def _drop_query_insertions(query_aligned: str, ref_aligned: str) -> Tuple[str, str]:
    """
    Remove columns where the query has an insertion relative to the reference.

    This is required when the desired output must have exactly the same length
    as the ungapped reference.

    Example:

        query: A C T G
        ref:   A - T G

    The C column has no reference coordinate, so it is removed in the
    reference-projected output.
    """
    if len(query_aligned) != len(ref_aligned):
        raise ValueError(
            "Cannot project alignment because query_aligned and ref_aligned "
            "have different lengths."
        )

    projected_query = []
    projected_ref = []
    dropped_insertions = 0

    for query_base, ref_base in zip(query_aligned, ref_aligned):
        if ref_base == "-":
            dropped_insertions += 1
            continue

        projected_query.append(query_base)
        projected_ref.append(ref_base)

    if dropped_insertions:
        logging.debug(
            "Dropped %d query insertion column(s) to keep alignment in reference coordinates.",
            dropped_insertions,
        )

    return "".join(projected_query), "".join(projected_ref)


def _apply_cigar(test_sequence: str, ref_sequence: str, cigar: str) -> Tuple[str, str]:
    """Build aligned query/reference strings from a SAM-style CIGAR."""
    test_aligned = []
    ref_aligned = []
    test_index = 0
    ref_index = 0

    for length, operation in re.findall(r"(\d+)([MID=X])", cigar):
        count = int(length)

        if operation in {"M", "=", "X"}:
            test_aligned.append(test_sequence[test_index:test_index + count])
            ref_aligned.append(ref_sequence[ref_index:ref_index + count])
            test_index += count
            ref_index += count

        elif operation == "I":
            test_aligned.append(test_sequence[test_index:test_index + count])
            ref_aligned.append("-" * count)
            test_index += count

        elif operation == "D":
            test_aligned.append("-" * count)
            ref_aligned.append(ref_sequence[ref_index:ref_index + count])
            ref_index += count

        else:  # pragma: no cover
            logging.warning("Unsupported CIGAR operation '%s' ignored.", operation)

    if test_index < len(test_sequence):
        remaining = len(test_sequence) - test_index
        test_aligned.append(test_sequence[test_index:])
        ref_aligned.append("-" * remaining)

    if ref_index < len(ref_sequence):
        remaining = len(ref_sequence) - ref_index
        test_aligned.append("-" * remaining)
        ref_aligned.append(ref_sequence[ref_index:])

    return "".join(test_aligned), "".join(ref_aligned)


_ALIGNERS: Dict[str, Callable[[SeqRecord, SeqRecord], Tuple[str, str]]] = {
    EDLIB_HW: edlib_hw_align,
    MAFFT: mafft_align,
    PARASAIL_NW: parasail_nw_align,
    PYFAMSA: pyfamsa_align,
}

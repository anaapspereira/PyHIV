"""Offline annotation, independent of alignment and mutation detection.

DRM membership and mutation types (including Accessory) are separate official
HIVDB datasets. Unknown codons are not evidence for absence of resistance.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from types import MappingProxyType
from typing import Mapping

LOCAL_DRM_SOURCE = "HIVDB_HIVFACTS"
LOCAL_DRM_VERSION = "be1c11a5145fea9073fdb71801919d4a34265336"
UNRESOLVED_CODON_STATUSES = frozenset({
    "PARTIAL_CODON", "FRAMESHIFT", "FRAMESHIFT_RESTORED",
    "NOT_COVERED", "PARTIAL_COVERAGE",
})
RESOLVED_AAS = frozenset("ACDEFGHIKLMNPQRSTVWY*-_")
CATALOG_FILES = ("drms_hiv1.json", "mutation-type-pairs_hiv1.json")


def drm_status(is_drm: bool | None) -> str:
    return "UNRESOLVED" if is_drm is None else "DRM" if is_drm else "NON_DRM"


def tri_state_any(values) -> bool | None:
    values = tuple(values)
    if any(value is True for value in values):
        return True
    if any(value is None for value in values):
        return None
    return False


def drm_evaluation_status(components: tuple["DRMComponent", ...]) -> str:
    unresolved = any(component.is_drm is None for component in components)
    resolved = any(component.is_drm is not None for component in components)
    if unresolved and resolved:
        return "PARTIAL"
    if unresolved:
        return "UNRESOLVED"
    return "COMPLETE"


@dataclass(frozen=True)
class DRMComponent:
    aa: str
    is_drm: bool | None
    drm_class: str
    drug_class: str
    mutation_types: str
    is_accessory: bool | None

    @property
    def drm_status(self) -> str:
        return drm_status(self.is_drm)


@dataclass(frozen=True)
class DRMAnnotation:
    is_drm: bool | None
    drm_source: str
    drm_class: str
    drug_class: str
    drm_version: str = LOCAL_DRM_VERSION
    drm_catalog_sha256: str = ""
    mutation_types: str = ""
    is_accessory: bool | None = None
    components: tuple[DRMComponent, ...] = ()

    @property
    def drm_status(self) -> str:
        return drm_status(self.is_drm)

    @property
    def drm_evaluation_status(self) -> str:
        if not self.components:
            return "UNRESOLVED" if self.is_drm is None else "COMPLETE"
        return drm_evaluation_status(self.components)


@dataclass(frozen=True)
class DRMCatalog:
    # DRM entries retain the historical (mutation type, drug class) interface.
    drms: Mapping[tuple[str, int, str], tuple[tuple[str, str], ...]]
    mutation_types: Mapping[tuple[str, int, str], tuple[str, ...]]
    sha256: str


@lru_cache(maxsize=1)
def load_catalog() -> DRMCatalog:
    """Verify both pinned inputs and return an immutable annotation snapshot.

    The bundle hash is SHA-256 of the compact, sorted JSON mapping of filenames
    to verified SHA-256 hashes. It covers membership AND mutation types.
    """
    root = files("pyhiv.mutations").joinpath("data")
    metadata = json.loads(root.joinpath("metadata.json").read_text())
    if metadata["revision"] != LOCAL_DRM_VERSION or metadata["source"] != LOCAL_DRM_SOURCE:
        raise ValueError("DRM catalogue provenance does not match its declared version/source")
    if set(metadata["sha256"]) != set(CATALOG_FILES):
        raise ValueError("DRM catalogue manifest must cover both official datasets")
    decoded, digests = {}, {}
    for name in CATALOG_FILES:
        content = root.joinpath(name).read_bytes()
        digests[name] = hashlib.sha256(content).hexdigest()
        if digests[name] != metadata["sha256"][name]:
            raise ValueError(f"DRM catalogue checksum mismatch: {name}")
        decoded[name] = json.loads(content)
    types: dict[tuple[str, int, str, str], set[str]] = {}
    all_types: dict[tuple[str, int, str], set[str]] = {}
    for row in decoded["mutation-type-pairs_hiv1.json"]:
        for aa in row["aas"]:
            key = (row["gene"], row["position"], aa)
            types.setdefault((*key, row["drugClass"]), set()).add(row["mutationType"])
            all_types.setdefault(key, set()).add(row["mutationType"])
    entries: dict[tuple[str, int, str], set[tuple[str, str]]] = {}
    for drug_class, rows in decoded["drms_hiv1.json"].items():
        for row in rows:
            gene, position, aa = row["gene"], row["position"], row["aa"]
            if not isinstance(position, int) or position < 1 or aa not in RESOLVED_AAS:
                raise ValueError(f"Invalid DRM catalogue entry: {row!r}")
            labels = types.get((gene, position, aa, drug_class), {"Other"})
            entries.setdefault((gene, position, aa), set()).update(
                (label, drug_class) for label in labels
            )
    digest = hashlib.sha256(json.dumps(digests, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return DRMCatalog(
        MappingProxyType({key: tuple(sorted(values)) for key, values in entries.items()}),
        MappingProxyType({key: tuple(sorted(values)) for key, values in all_types.items()}),
        digest,
    )


def load_drm_catalog() -> Mapping[tuple[str, int, str], tuple[tuple[str, str], ...]]:
    """Return verified DRM membership; mutation types are available independently."""
    return load_catalog().drms


def annotate_drm(
    gene: str, position: int, alt_aa: str, *,
    codon_status: str = "", mutation_type: str = "substitution",
) -> DRMAnnotation:
    """Annotate each observed component without inferring unobserved residues.

    is_drm is True / False / None; None means UNRESOLVED, never negative.
    Codon-level uncertainty overrides every component. Otherwise each residue
    is evaluated independently; an unknown component makes the aggregate
    UNRESOLVED only when no resolved component confirms the annotation.
    Insertions are catalogue events ('_'), not mixtures of inserted residues.
    Empty input and unsupported residue symbols are unresolved.
    """
    catalog = load_catalog()
    uncertain = bool(UNRESOLVED_CODON_STATUSES.intersection(codon_status.split("|")))
    uncertain |= mutation_type in {"partial_codon", "frameshift"}
    if mutation_type == "insertion":
        uncertain |= not alt_aa or any(aa not in RESOLVED_AAS for aa in alt_aa)
        residues = ("_",)
    else:
        residues = tuple(sorted(set(alt_aa))) or ("X",)
    components = []
    for residue in residues:
        if uncertain or residue not in RESOLVED_AAS:
            components.append(DRMComponent(residue, None, "", "", "", None))
            continue
        key = (gene, position, residue)
        matches = catalog.drms.get(key, ())
        labels = catalog.mutation_types.get(key, ("Other",))
        components.append(DRMComponent(
            residue, bool(matches),
            ";".join(sorted({label for label, _ in matches})),
            ";".join(sorted({drug for _, drug in matches})),
            ";".join(labels), "Accessory" in labels,
        ))
    def joined(field: str) -> str:
        return ";".join(sorted({value for c in components
                                for value in getattr(c, field).split(";") if value}))
    return DRMAnnotation(
        is_drm=tri_state_any(c.is_drm for c in components),
        drm_source=LOCAL_DRM_SOURCE,
        drm_class=joined("drm_class"), drug_class=joined("drug_class"),
        drm_version=LOCAL_DRM_VERSION, drm_catalog_sha256=catalog.sha256,
        mutation_types=joined("mutation_types"),
        is_accessory=tri_state_any(c.is_accessory for c in components),
        components=tuple(components),
    )

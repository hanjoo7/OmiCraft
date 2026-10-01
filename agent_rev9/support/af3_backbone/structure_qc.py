from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any

AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True)
class AtomRecord:
    atom_name: str
    residue_name: str
    label_chain_id: str
    auth_chain_id: str
    label_seq_id: int | None
    auth_seq_id: int | None
    element: str
    x: float
    y: float
    z: float

    @property
    def coord(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_hydrogen(self) -> bool:
        return self.element.upper() == "H" or self.atom_name.upper().startswith("H")


def _maybe_int(value: str | None) -> int | None:
    if value in {None, "", ".", "?"}:
        return None
    try:
        return int(str(value))
    except ValueError:
        return None


def _maybe_float(value: str | None) -> float | None:
    if value in {None, "", ".", "?"}:
        return None
    try:
        return float(str(value))
    except ValueError:
        return None


def parse_cif_atoms(path: Path) -> list[AtomRecord]:
    lines = path.read_text(encoding="utf-8").splitlines()
    idx = 0
    while idx < len(lines):
        if lines[idx].strip() != "loop_":
            idx += 1
            continue
        idx += 1
        headers: list[str] = []
        while idx < len(lines) and lines[idx].strip().startswith("_"):
            headers.append(lines[idx].strip())
            idx += 1
        if not headers or not all(header.startswith("_atom_site.") for header in headers):
            while idx < len(lines) and lines[idx].strip() not in {"#", "loop_"}:
                idx += 1
            continue
        return _parse_atom_site_rows(headers, lines[idx:])
    return []


def _parse_atom_site_rows(headers: list[str], rows: list[str]) -> list[AtomRecord]:
    header_index = {header: idx for idx, header in enumerate(headers)}

    def get(parts: list[str], header: str) -> str:
        pos = header_index.get(header)
        return parts[pos] if pos is not None and pos < len(parts) else ""

    atoms: list[AtomRecord] = []
    for line in rows:
        stripped = line.strip()
        if not stripped or stripped == "#" or stripped == "loop_" or stripped.startswith("_"):
            break
        parts = shlex.split(stripped, posix=True)
        x = _maybe_float(get(parts, "_atom_site.Cartn_x"))
        y = _maybe_float(get(parts, "_atom_site.Cartn_y"))
        z = _maybe_float(get(parts, "_atom_site.Cartn_z"))
        if x is None or y is None or z is None:
            continue
        atoms.append(
            AtomRecord(
                atom_name=get(parts, "_atom_site.label_atom_id"),
                residue_name=get(parts, "_atom_site.label_comp_id"),
                label_chain_id=get(parts, "_atom_site.label_asym_id"),
                auth_chain_id=get(parts, "_atom_site.auth_asym_id"),
                label_seq_id=_maybe_int(get(parts, "_atom_site.label_seq_id")),
                auth_seq_id=_maybe_int(get(parts, "_atom_site.auth_seq_id")),
                element=get(parts, "_atom_site.type_symbol"),
                x=x,
                y=y,
                z=z,
            )
        )
    return atoms


def _protein_ca_by_residue(
    atoms: list[AtomRecord], protein_chain_id: str
) -> dict[int, tuple[float, float, float]]:
    coords = {}
    for atom in atoms:
        if atom.atom_name != "CA":
            continue
        if atom.label_chain_id != protein_chain_id and atom.auth_chain_id != protein_chain_id:
            continue
        number = atom.label_seq_id or atom.auth_seq_id
        if number is not None:
            coords[number] = atom.coord
    return coords


def _kabsch(points: Any, reference: Any, np: Any) -> tuple[Any, Any, Any]:
    points_centroid = points.mean(axis=0)
    reference_centroid = reference.mean(axis=0)
    centered_points = points - points_centroid
    centered_reference = reference - reference_centroid
    covariance = centered_points.T @ centered_reference
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = vt_matrix.T @ u_matrix.T
    if np.linalg.det(rotation) < 0:
        vt_matrix[-1, :] *= -1
        rotation = vt_matrix.T @ u_matrix.T
    return rotation, points_centroid, reference_centroid

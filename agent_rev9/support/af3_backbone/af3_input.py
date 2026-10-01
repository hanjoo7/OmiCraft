from __future__ import annotations

from typing import Any


def validate_af3_json(payload: dict[str, Any]) -> None:
    if payload.get("dialect") != "alphafold3":
        raise ValueError("AF3 input dialect must be alphafold3")
    if payload.get("version") != 4:
        raise ValueError("AF3 input version must be 4")
    if not payload.get("modelSeeds"):
        raise ValueError("AF3 input must preserve at least one model seed")
    sequences = payload.get("sequences")
    if not isinstance(sequences, list) or len(sequences) != 2:
        raise ValueError("AF3 input requires protein and ligand sequences")
    ligand = sequences[1].get("ligand", {})
    if bool(ligand.get("ccdCodes")) == bool(ligand.get("smiles")):
        raise ValueError("AF3 ligand requires exactly one of ccdCodes or smiles")

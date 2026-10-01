from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class SmilesValidationResult:
    candidate_id: str
    smiles: str
    canonical_smiles: str
    smiles_valid: bool
    rdkit_parse_status: str
    molecular_weight: str
    validation_message: str
    warnings: list[str]

    def as_columns(self) -> dict[str, str]:
        return {
            "candidate_id": self.candidate_id,
            "smiles": self.smiles,
            "canonical_smiles": self.canonical_smiles,
            "smiles_valid": str(self.smiles_valid).lower(),
            "rdkit_parse_status": self.rdkit_parse_status,
            "molecular_weight": self.molecular_weight,
            "validation_message": self.validation_message,
        }


def validate_smiles_rdkit(candidate_id: str, smiles: str) -> SmilesValidationResult:
    original_smiles = smiles or ""
    stripped = original_smiles.strip()
    if not stripped:
        return _invalid(candidate_id, original_smiles, "missing_smiles", "smiles is required")

    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors  # type: ignore
    except Exception:
        return _invalid(
            candidate_id,
            original_smiles,
            "dependency_missing",
            "RDKit is required for SMILES validation",
        )

    mol = Chem.MolFromSmiles(stripped, sanitize=True)
    if mol is None:
        return _invalid(
            candidate_id,
            original_smiles,
            "parse_failed",
            "RDKit could not parse SMILES",
        )

    warnings = _warnings(stripped)
    try:
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return _invalid(
            candidate_id,
            original_smiles,
            "canonicalization_failed",
            "RDKit could not generate canonical SMILES",
        )
    try:
        molecular_weight = f"{Descriptors.MolWt(mol):.6g}"
    except Exception:
        return _invalid(
            candidate_id,
            original_smiles,
            "molecular_weight_failed",
            "RDKit could not calculate molecular weight",
        )
    message = ";".join(warnings)
    return SmilesValidationResult(
        candidate_id=candidate_id,
        smiles=original_smiles,
        canonical_smiles=canonical,
        smiles_valid=True,
        rdkit_parse_status="parsed",
        molecular_weight=molecular_weight,
        validation_message=message,
        warnings=warnings,
    )


def _warnings(smiles: str) -> list[str]:
    warnings: list[str] = []
    if "." in smiles:
        warnings.append("multi_fragment_or_salt_smiles")
    if any(token in smiles for token in ["[Na", "[K", "[Fe", "[Zn", "[Mg", "[Ca"]):
        warnings.append("metal_or_salt_token")
    return warnings


def _invalid(
    candidate_id: str,
    smiles: str,
    status: str,
    message: str,
) -> SmilesValidationResult:
    return SmilesValidationResult(
        candidate_id=candidate_id,
        smiles=smiles,
        canonical_smiles="",
        smiles_valid=False,
        rdkit_parse_status=status,
        molecular_weight="",
        validation_message=message,
        warnings=[message],
    )


def apply_validation_columns(row: dict[str, Any], result: SmilesValidationResult) -> None:
    row["canonical_smiles"] = result.canonical_smiles
    row["smiles_valid"] = str(result.smiles_valid).lower()
    row["rdkit_parse_status"] = result.rdkit_parse_status
    row["molecular_weight"] = result.molecular_weight
    row["validation_message"] = result.validation_message

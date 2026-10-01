from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..pipeline_common.chem_validation import apply_validation_columns, validate_smiles_rdkit
from .io import ordered_fieldnames, read_csv, write_csv
from .status import PredictionStatus

VALIDATION_COLUMNS = [
    "canonical_smiles",
    "smiles_valid",
    "rdkit_parse_status",
    "molecular_weight",
    "validation_message",
]


@dataclass(frozen=True)
class PreparedData:
    prepared: list[dict[str, Any]]
    invalid: list[dict[str, Any]]
    input_fieldnames: list[str]


def validate_smiles(smiles: str) -> list[str]:
    result = validate_smiles_rdkit("", smiles)
    if not result.smiles_valid:
        raise ValueError(f"{result.rdkit_parse_status}: {result.validation_message}")
    return result.warnings


def prepare_candidates(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    candidate_id_column: str = "candidate_id",
    smiles_column: str = "smiles",
    force: bool = False,
) -> PreparedData:
    rows, fieldnames = read_csv(input_path)
    missing = [name for name in [candidate_id_column, smiles_column] if name not in fieldnames]
    if missing:
        raise ValueError(f"Candidate CSV missing required columns: {', '.join(missing)}")

    prepared: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    seen: set[str] = set()
    output = Path(output_dir)

    for index, row in enumerate(rows, start=2):
        candidate_id = (row.get(candidate_id_column) or "").strip()
        smiles = (row.get("model_input_smiles") or row.get(smiles_column) or "").strip()
        handoff_status = (row.get("handoff_status") or "").strip()
        eligible = (row.get("eligible_for_admet") or "").strip().lower()
        if handoff_status and (eligible != "true" or handoff_status != "ready"):
            invalid_row = dict(row)
            invalid_row[candidate_id_column] = candidate_id
            invalid_row[smiles_column] = smiles
            invalid_row["failure_stage"] = "import-handoff"
            invalid_row["error_message"] = ""
            invalid_row["skip_reason"] = row.get("handoff_reason", "")
            invalid_row["prediction_status"] = PredictionStatus.SKIPPED_UPSTREAM.value
            invalid_row.setdefault("canonical_smiles", "")
            invalid_row["smiles_valid"] = "false" if not smiles else ""
            invalid_row["rdkit_parse_status"] = (
                "missing_smiles" if not smiles else "not_evaluated_skipped_upstream"
            )
            invalid_row.setdefault("molecular_weight", "")
            invalid_row.setdefault("validation_message", invalid_row["skip_reason"])
            invalid.append(invalid_row)
            continue
        error: str | None = None
        if not candidate_id:
            error = f"row {index}: candidate_id is required"
        elif candidate_id in seen:
            error = f"duplicate candidate_id: {candidate_id}"
        elif not smiles:
            error = "missing_smiles: smiles is required"
        if error is None:
            validation = validate_smiles_rdkit(candidate_id, smiles)
            if not validation.smiles_valid:
                error = f"{validation.rdkit_parse_status}: {validation.validation_message}"
            else:
                seen.add(candidate_id)
                out_row = dict(row)
                out_row[candidate_id_column] = candidate_id
                out_row[smiles_column] = smiles
                apply_validation_columns(out_row, validation)
                out_row["validation_warnings"] = ";".join(validation.warnings)
                prepared.append(out_row)
                continue
        invalid_row = dict(row)
        invalid_row[candidate_id_column] = candidate_id
        invalid_row[smiles_column] = smiles
        if smiles:
            apply_validation_columns(invalid_row, validate_smiles_rdkit(candidate_id, smiles))
        else:
            invalid_row["canonical_smiles"] = ""
            invalid_row["smiles_valid"] = "false"
            invalid_row["rdkit_parse_status"] = "missing_smiles"
            invalid_row["molecular_weight"] = ""
            invalid_row["validation_message"] = "smiles is required"
        invalid_row["failure_stage"] = "prepare"
        invalid_row["error_message"] = error or "invalid input"
        invalid_row["prediction_status"] = PredictionStatus.INVALID_INPUT.value
        invalid.append(invalid_row)

    prepared_fields = ordered_fieldnames(
        [candidate_id_column, smiles_column, *VALIDATION_COLUMNS, "validation_warnings"],
        prepared,
        fieldnames,
    )
    invalid_fields = ordered_fieldnames(
        [
            candidate_id_column,
            smiles_column,
            *VALIDATION_COLUMNS,
            "failure_stage",
            "error_message",
            "prediction_status",
        ],
        invalid,
        fieldnames,
    )
    write_csv(output / "prepared_candidates.csv", prepared, prepared_fields, force=force)
    write_csv(output / "invalid_candidates.csv", invalid, invalid_fields, force=force)
    return PreparedData(prepared=prepared, invalid=invalid, input_fieldnames=fieldnames)

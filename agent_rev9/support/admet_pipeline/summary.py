from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

from .io import ordered_fieldnames, read_csv, write_csv
from .status import (
    PredictionStatus,
    prediction_status_from_endpoint_counts,
    workflow_status_for_prediction,
)

BASE_SUMMARY_COLUMNS = [
    "candidate_id",
    "smiles",
    "canonical_smiles",
    "pipeline1_run_id",
    "pipeline1_status",
    "structural_call",
    "handoff_status",
    "eligible_for_admet",
    "handoff_reason",
    "selected_structure_path",
    "cross_seed_rmsd",
    "reference_pose_rmsd",
    "pipeline1_warning",
    "smiles_valid",
    "prediction_status",
    "workflow_status",
    "backend",
    "backend_version",
    "model_version",
    "is_mock",
    "available_endpoint_count",
    "missing_endpoint_count",
    "warning_count",
    "error_message",
]
EXCLUDED_EXTRA_FIELDS = {
    "failure_stage",
    "error_message",
    "prediction_status",
    "workflow_status",
}


def summarize_predictions(
    predictions_path: str | Path,
    prepared_path: str | Path,
    invalid_path: str | Path,
    output_dir: str | Path,
    *,
    selected_endpoints: dict[str, dict[str, list[str]]],
    force: bool = False,
) -> list[dict[str, Any]]:
    predictions, prediction_fields = read_csv(predictions_path)
    prepared, prepared_fields = read_csv(prepared_path)
    invalid, invalid_fields = read_csv(invalid_path)
    by_candidate: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in predictions:
        by_candidate[_chemical_key(row)].append(row)

    summaries: list[dict[str, Any]] = []
    selected_names = list(selected_endpoints)
    extra_fields = [
        name
        for name in prepared_fields + invalid_fields + prediction_fields
        if name not in BASE_SUMMARY_COLUMNS
        and name not in selected_names
        and name not in EXCLUDED_EXTRA_FIELDS
        and name not in {"endpoint", "endpoint_name", "value", "endpoint_value"}
    ]

    for candidate in prepared:
        rows = by_candidate.get(_chemical_key(candidate), [])
        failed = [
            row
            for row in rows
            if row["prediction_status"] == PredictionStatus.PREDICTION_FAILED.value
        ]
        successes = [
            row
            for row in rows
            if row["prediction_status"] == PredictionStatus.SUCCESS.value and row["endpoint_name"]
        ]
        status = prediction_status_from_endpoint_counts(
            success_count=len(successes), failed_count=len(failed), total_count=len(rows)
        )
        if not rows:
            error_message = "no prediction rows generated"
        else:
            error_message = "; ".join(
                row["error_message"] for row in failed if row["error_message"]
            )

        summary = _base_summary(candidate, rows)
        summary.update(
            {
                "prediction_status": status,
                "workflow_status": workflow_status_for_prediction(status),
                "backend": rows[0]["backend"] if rows else "",
                "backend_version": rows[0].get("backend_version", "") if rows else "",
                "model_version": rows[0].get("model_version", "") if rows else "",
                "is_mock": rows[0]["is_mock"] if rows else "",
                "available_endpoint_count": len({row["endpoint_name"] for row in successes}),
                "missing_endpoint_count": _missing_count(rows, selected_endpoints),
                "warning_count": _warning_count(candidate),
                "error_message": error_message,
            }
        )
        summary.update(_selected_values(rows, selected_endpoints))
        for field in extra_fields:
            summary[field] = candidate.get(field, "")
        summaries.append(summary)

    for candidate in invalid:
        prediction_status = (
            candidate.get("prediction_status") or PredictionStatus.INVALID_INPUT.value
        )
        summary = _base_summary(candidate, [])
        summary.update(
            {
                "prediction_status": prediction_status,
                "workflow_status": workflow_status_for_prediction(prediction_status),
                "backend": "",
                "backend_version": "",
                "model_version": "",
                "is_mock": "",
                "available_endpoint_count": 0,
                "missing_endpoint_count": len(selected_endpoints),
                "warning_count": 0,
                "error_message": candidate.get("error_message", ""),
            }
        )
        for name in selected_names:
            summary[name] = "not_available"
        for field in extra_fields:
            summary[field] = candidate.get(field, "")
        summaries.append(summary)

    fieldnames = ordered_fieldnames(BASE_SUMMARY_COLUMNS + selected_names, summaries, extra_fields)
    write_csv(Path(output_dir) / "candidate_summary.csv", summaries, fieldnames, force=force)
    return summaries


def _base_summary(candidate: dict[str, str], rows: list[dict[str, str]]) -> dict[str, Any]:
    pipeline1_status = candidate.get("pipeline1_status", "")
    return {
        "candidate_id": candidate.get("candidate_id", ""),
        "variant_id": candidate.get("variant_id", "") or "as_provided",
        "chemical_state_id": candidate.get("chemical_state_id", ""),
        "smiles": candidate.get("smiles", ""),
        "canonical_smiles": candidate.get("canonical_smiles", ""),
        "pipeline1_run_id": candidate.get("pipeline1_run_id", ""),
        "pipeline1_status": pipeline1_status,
        "structural_call": candidate.get("structural_call", "") or pipeline1_status,
        "handoff_status": candidate.get("handoff_status", ""),
        "eligible_for_admet": candidate.get("eligible_for_admet", ""),
        "handoff_reason": candidate.get("handoff_reason", ""),
        "selected_structure_path": candidate.get("selected_structure_path", ""),
        "cross_seed_rmsd": candidate.get("cross_seed_rmsd", ""),
        "reference_pose_rmsd": candidate.get("reference_pose_rmsd", ""),
        "pipeline1_warning": candidate.get("pipeline1_warning", ""),
        "smiles_valid": candidate.get("smiles_valid", ""),
    }


def _chemical_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("candidate_id", ""),
        row.get("variant_id", "") or "as_provided",
        row.get("chemical_state_id", ""),
    )


def _warning_count(candidate: dict[str, str]) -> int:
    warnings = candidate.get("validation_warnings", "")
    return len([item for item in warnings.split(";") if item])


def _selected_values(
    rows: list[dict[str, str]],
    selected_endpoints: dict[str, dict[str, list[str]]],
) -> dict[str, str]:
    endpoint_lookup = {row["endpoint_name"]: row.get("endpoint_value", "") for row in rows}
    values = {}
    for output_name, spec in selected_endpoints.items():
        value = "not_available"
        for alias in spec.get("aliases", []):
            if alias in endpoint_lookup:
                value = endpoint_lookup[alias] or "not_available"
                break
        values[output_name] = value
    return values


def _missing_count(
    rows: list[dict[str, str]], selected_endpoints: dict[str, dict[str, list[str]]]
) -> int:
    values = _selected_values(rows, selected_endpoints)
    return sum(1 for value in values.values() if value in {"", "not_available"})

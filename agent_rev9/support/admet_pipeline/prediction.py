from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .admet_ai_predictor import ADMETAIPredictor
from .io import ordered_fieldnames, read_csv, write_csv
from .mock_predictor import MockADMETPredictor
from .predictor import EndpointPrediction

LONG_COLUMNS = [
    "candidate_id",
    "variant_id",
    "chemical_state_id",
    "smiles",
    "admet_input_smiles",
    "admet_chemical_state_id",
    "admet_chemical_state_source",
    "endpoint_name",
    "endpoint",
    "endpoint_value",
    "value",
    "endpoint_unit",
    "prediction_type",
    "backend",
    "backend_version",
    "model_version",
    "prediction_started_at",
    "prediction_finished_at",
    "prediction_status",
    "error_message",
    "is_mock",
]


def make_predictor(backend: str, *, batch_size: int = 32, device: str = "cpu"):
    if backend == "mock":
        return MockADMETPredictor()
    if backend == "admet_ai":
        return ADMETAIPredictor(batch_size=batch_size, device=device)
    raise ValueError(f"Unsupported ADMET backend: {backend}")


def prediction_to_row(
    prediction: EndpointPrediction,
    predictor: Any,
    *,
    prediction_started_at: str,
    prediction_finished_at: str,
) -> dict[str, Any]:
    return {
        "candidate_id": prediction.candidate_id,
        "smiles": prediction.smiles,
        "endpoint_name": prediction.endpoint_name,
        "endpoint": prediction.endpoint_name,
        "endpoint_value": prediction.endpoint_value,
        "value": prediction.endpoint_value,
        "endpoint_unit": prediction.endpoint_unit or "not_available",
        "prediction_type": prediction.prediction_type,
        "backend": predictor.backend,
        "backend_version": predictor.backend_version,
        "model_version": getattr(predictor, "model_version", "not_available"),
        "prediction_started_at": prediction_started_at,
        "prediction_finished_at": prediction_finished_at,
        "prediction_status": prediction.prediction_status,
        "error_message": prediction.error_message,
        "is_mock": str(bool(predictor.is_mock)).lower(),
    }


def run_predictions(
    prepared_path: str | Path,
    output_dir: str | Path,
    *,
    backend: str,
    batch_size: int,
    device: str,
    force: bool = False,
) -> tuple[list[dict[str, Any]], Any]:
    candidates, input_fields = read_csv(prepared_path)
    canonical = []
    for row in candidates:
        normalized = dict(row)
        normalized["candidate_id"] = row["candidate_id"]
        normalized["smiles"] = row.get("model_input_smiles") or row.get("smiles", "")
        if not normalized["smiles"]:
            continue
        canonical.append(normalized)
    predictor = make_predictor(backend, batch_size=batch_size, device=device)
    prediction_started_at = datetime.now(UTC).isoformat()
    predictions = predictor.predict(canonical)
    prediction_finished_at = datetime.now(UTC).isoformat()
    rows = [
        prediction_to_row(
            prediction,
            predictor,
            prediction_started_at=prediction_started_at,
            prediction_finished_at=prediction_finished_at,
        )
        for prediction in predictions
    ]
    extras = [name for name in input_fields if name not in LONG_COLUMNS]
    metadata = {row["candidate_id"]: row for row in canonical}
    for row in rows:
        source = metadata.get(row["candidate_id"], {})
        row["variant_id"] = source.get("variant_id", "as_provided")
        row["chemical_state_id"] = source.get("chemical_state_id", "")
        row["admet_input_smiles"] = source.get("smiles", "")
        row["admet_chemical_state_id"] = source.get("chemical_state_id", "")
        row["admet_chemical_state_source"] = source.get("chemical_state_source", "")
        for name in extras:
            row[name] = source.get(name, "")
    fieldnames = ordered_fieldnames(LONG_COLUMNS, rows, extras)
    write_csv(Path(output_dir) / "admet_predictions_long.csv", rows, fieldnames, force=force)
    return rows, predictor

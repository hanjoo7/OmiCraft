from __future__ import annotations

import importlib
from importlib import metadata
from typing import Any

from .predictor import EndpointPrediction
from .status import PredictionStatus


class ADMETAIUnavailableError(RuntimeError):
    pass


class ADMETAIPredictor:
    backend = "admet_ai"
    is_mock = False

    def __init__(self, *, batch_size: int = 32, device: str = "cpu") -> None:
        self.batch_size = batch_size
        self.device = device
        try:
            self.module = importlib.import_module("admet_ai")
        except ImportError as exc:
            raise ADMETAIUnavailableError(
                "admet_ai backend was requested, but the admet_ai Python package is not "
                "installed. Install it in this environment and rerun with --backend "
                "admet_ai. The pipeline will not automatically fall back to mock."
            ) from exc
        try:
            self.backend_version = metadata.version("admet_ai")
        except metadata.PackageNotFoundError:
            self.backend_version = "unknown"
        self.model_version = self.backend_version
        self.model = self.module.ADMETModel()

    def predict(self, candidates: list[dict[str, Any]]) -> list[EndpointPrediction]:
        smiles_values = [str(candidate["smiles"]) for candidate in candidates]
        try:
            predictions = self.model.predict(smiles=smiles_values)
        except Exception as exc:
            return [
                EndpointPrediction(
                    candidate_id=str(candidate["candidate_id"]),
                    smiles=str(candidate["smiles"]),
                    endpoint_name="not_available",
                    endpoint_value="",
                    endpoint_unit="not_available",
                    prediction_type="not_available",
                    prediction_status=PredictionStatus.PREDICTION_FAILED.value,
                    error_message=str(exc),
                )
                for candidate in candidates
            ]

        rows: list[EndpointPrediction] = []
        if hasattr(predictions, "iloc"):
            for index, candidate in enumerate(candidates):
                prediction_row = predictions.iloc[index]
                rows.extend(self._to_endpoint_rows(candidate, prediction_row.to_dict()))
        else:
            rows.extend(self._to_endpoint_rows(candidates[0], dict(predictions)))
        return rows

    def _to_endpoint_rows(
        self,
        candidate: dict[str, Any],
        values: dict[str, Any],
    ) -> list[EndpointPrediction]:
        rows: list[EndpointPrediction] = []
        for endpoint_name, value in values.items():
            rows.append(
                EndpointPrediction(
                    candidate_id=str(candidate["candidate_id"]),
                    smiles=str(candidate["smiles"]),
                    endpoint_name=str(endpoint_name),
                    endpoint_value=_format_value(value),
                    endpoint_unit="not_available",
                    prediction_type="admet_ai",
                    prediction_status=PredictionStatus.SUCCESS.value,
                )
            )
        return rows


def _format_value(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except TypeError:
        pass
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)

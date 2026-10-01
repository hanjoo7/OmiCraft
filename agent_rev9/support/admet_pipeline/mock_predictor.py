from __future__ import annotations

import hashlib
from typing import Any

from .predictor import EndpointPrediction
from .status import PredictionStatus


class MockADMETPredictor:
    backend = "mock"
    backend_version = "synthetic-v1"
    model_version = "synthetic-v1"
    is_mock = True

    endpoints = [
        ("Solubility_AqSolDB", "logS", "regression"),
        ("BBB_Martins", "probability", "classification"),
        ("hERG", "probability", "classification"),
        ("AMES", "probability", "classification"),
        ("DILI", "probability", "classification"),
        ("CYP3A4", "probability", "classification"),
    ]

    def predict(self, candidates: list[dict[str, Any]]) -> list[EndpointPrediction]:
        rows: list[EndpointPrediction] = []
        for candidate in candidates:
            candidate_id = str(candidate["candidate_id"])
            smiles = str(candidate["smiles"])
            if candidate_id == "partial_fail":
                rows.append(
                    EndpointPrediction(
                        candidate_id=candidate_id,
                        smiles=smiles,
                        endpoint_name="Solubility_AqSolDB",
                        endpoint_value=self._value(
                            candidate_id, smiles, "Solubility_AqSolDB", "regression"
                        ),
                        endpoint_unit="logS",
                        prediction_type="regression",
                        prediction_status=PredictionStatus.SUCCESS.value,
                    )
                )
                rows.append(
                    EndpointPrediction(
                        candidate_id=candidate_id,
                        smiles=smiles,
                        endpoint_name="hERG",
                        endpoint_value="",
                        endpoint_unit="probability",
                        prediction_type="classification",
                        prediction_status=PredictionStatus.PREDICTION_FAILED.value,
                        error_message="deterministic mock endpoint failure",
                    )
                )
                continue
            if candidate_id == "predict_fail":
                rows.append(
                    EndpointPrediction(
                        candidate_id=candidate_id,
                        smiles=smiles,
                        endpoint_name="not_available",
                        endpoint_value="",
                        endpoint_unit="not_available",
                        prediction_type="not_available",
                        prediction_status=PredictionStatus.PREDICTION_FAILED.value,
                        error_message="deterministic mock failure requested by candidate_id",
                    )
                )
                continue
            for endpoint_name, unit, prediction_type in self.endpoints:
                value = self._value(candidate_id, smiles, endpoint_name, prediction_type)
                rows.append(
                    EndpointPrediction(
                        candidate_id=candidate_id,
                        smiles=smiles,
                        endpoint_name=endpoint_name,
                        endpoint_value=value,
                        endpoint_unit=unit,
                        prediction_type=prediction_type,
                        prediction_status=PredictionStatus.SUCCESS.value,
                    )
                )
        return rows

    @staticmethod
    def _unit_float(candidate_id: str, smiles: str, endpoint_name: str) -> float:
        digest = hashlib.sha256(f"{candidate_id}|{smiles}|{endpoint_name}".encode()).hexdigest()
        return int(digest[:12], 16) / float(0xFFFFFFFFFFFF)

    def _value(
        self, candidate_id: str, smiles: str, endpoint_name: str, prediction_type: str
    ) -> str:
        value = self._unit_float(candidate_id, smiles, endpoint_name)
        if prediction_type == "regression":
            return f"{(-6.0 + value * 5.0):.3f}"
        return f"{value:.3f}"

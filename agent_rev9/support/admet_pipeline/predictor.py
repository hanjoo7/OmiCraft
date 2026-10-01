from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol


@dataclass(frozen=True)
class EndpointPrediction:
    candidate_id: str
    smiles: str
    endpoint_name: str
    endpoint_value: Any
    endpoint_unit: str
    prediction_type: str
    prediction_status: str
    error_message: str = ""


class ADMETPredictor(Protocol):
    backend: str
    backend_version: str
    is_mock: bool

    def predict(self, candidates: list[dict[str, Any]]) -> list[EndpointPrediction]: ...

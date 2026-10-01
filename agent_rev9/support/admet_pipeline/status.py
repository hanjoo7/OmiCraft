from __future__ import annotations

from enum import StrEnum


class PredictionStatus(StrEnum):
    SUCCESS = "success"
    PARTIAL = "partial"
    INVALID_INPUT = "invalid_input"
    PREDICTION_FAILED = "prediction_failed"
    SKIPPED_UPSTREAM = "skipped_upstream"


class WorkflowStatus(StrEnum):
    COMPLETE = "complete"
    PARTIAL = "partial"
    EXCLUDED_BY_STRUCTURE = "excluded_by_structure"
    INVALID_INPUT = "invalid_input"
    FAILED = "failed"


def workflow_status_for_prediction(prediction_status: str) -> str:
    mapping = {
        PredictionStatus.SUCCESS.value: WorkflowStatus.COMPLETE.value,
        PredictionStatus.PARTIAL.value: WorkflowStatus.PARTIAL.value,
        PredictionStatus.SKIPPED_UPSTREAM.value: WorkflowStatus.EXCLUDED_BY_STRUCTURE.value,
        PredictionStatus.INVALID_INPUT.value: WorkflowStatus.INVALID_INPUT.value,
        PredictionStatus.PREDICTION_FAILED.value: WorkflowStatus.FAILED.value,
    }
    return mapping.get(prediction_status, WorkflowStatus.FAILED.value)


def prediction_status_from_endpoint_counts(
    *,
    success_count: int,
    failed_count: int,
    total_count: int,
) -> str:
    if failed_count and success_count:
        return PredictionStatus.PARTIAL.value
    if failed_count and not success_count:
        return PredictionStatus.PREDICTION_FAILED.value
    if total_count == 0:
        return PredictionStatus.PREDICTION_FAILED.value
    return PredictionStatus.SUCCESS.value

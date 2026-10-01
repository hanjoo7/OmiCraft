from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .io import read_csv
from .status import PredictionStatus


def write_report(
    output_dir: str | Path,
    *,
    config: dict[str, Any],
    backend: str,
    backend_version: str,
    model_version: str,
    is_mock: bool,
    started_at: str,
    force: bool = False,
    handoff_manifest: dict[str, Any] | None = None,
    input_path: str | Path | None = None,
) -> dict[str, Any]:
    output = Path(output_dir)
    summary_rows, _ = read_csv(output / "candidate_summary.csv")
    invalid_rows, _ = read_csv(output / "invalid_candidates.csv")
    manifest_path = output / "run_manifest.json"
    report_path = output / "report.md"
    if manifest_path.exists() and not force:
        raise FileExistsError(f"{manifest_path} already exists. Use --force to overwrite.")
    if report_path.exists() and not force:
        raise FileExistsError(f"{report_path} already exists. Use --force to overwrite.")

    ended_at = datetime.now(UTC).isoformat()
    run_id = f"admet_{started_at.replace(':', '').replace('-', '')}"
    input_sha256 = _sha256(Path(input_path)) if input_path else ""
    counts = {
        "total_candidates": len(summary_rows),
        "success_candidates": _count(summary_rows, PredictionStatus.SUCCESS.value),
        "partial_candidates": _count(summary_rows, PredictionStatus.PARTIAL.value),
        "failed_candidates": _count(summary_rows, PredictionStatus.PREDICTION_FAILED.value),
        "invalid_input_candidates": _count(summary_rows, PredictionStatus.INVALID_INPUT.value),
        "skipped_upstream_candidates": _count(
            summary_rows, PredictionStatus.SKIPPED_UPSTREAM.value
        ),
    }
    manifest = {
        "run_id": run_id,
        "created_at": ended_at,
        "project": config.get("project", {}).get("name", "small_molecule_admet"),
        "started_at": started_at,
        "ended_at": ended_at,
        "input_path": str(input_path) if input_path else "",
        "input_sha256": input_sha256,
        "pipeline1_run_id": (handoff_manifest or {}).get("pipeline1_run_id", ""),
        "handoff_schema_version": (handoff_manifest or {}).get("schema_version", ""),
        "backend": backend,
        "backend_version": backend_version,
        "model_version": model_version,
        "is_mock": is_mock,
        "candidate_count": counts["total_candidates"],
        "success_count": counts["success_candidates"],
        "partial_count": counts["partial_candidates"],
        "skipped_count": counts["skipped_upstream_candidates"],
        "failed_count": counts["failed_candidates"],
        "config": config,
        "counts": counts,
        "pipeline1_handoff": handoff_manifest or {},
        "outputs": {
            "prepared_candidates": str(output / "prepared_candidates.csv"),
            "invalid_candidates": str(output / "invalid_candidates.csv"),
            "admet_predictions_long": str(output / "admet_predictions_long.csv"),
            "candidate_summary": str(output / "candidate_summary.csv"),
            "report": str(report_path),
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    report_path.write_text(_render_report(summary_rows, manifest), encoding="utf-8")
    return manifest


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _count(rows: list[dict[str, str]], status: str) -> int:
    return sum(1 for row in rows if row.get("prediction_status") == status)


def _render_report(rows: list[dict[str, str]], manifest: dict[str, Any]) -> str:
    counts = manifest["counts"]
    lines = ["# ADMET Report", ""]
    if manifest["is_mock"]:
        lines.extend(
            [
                "> This report was generated using the mock backend.",
                "> Synthetic values; not scientific ADMET predictions.",
                "",
            ]
        )
    lines.extend(
        [
            f"- Started at: {manifest['started_at']}",
            f"- Ended at: {manifest['ended_at']}",
            f"- Backend: {manifest['backend']} ({manifest['backend_version']})",
            f"- Mock backend: {str(manifest['is_mock']).lower()}",
            f"- Total candidates: {counts['total_candidates']}",
            f"- Success candidates: {counts['success_candidates']}",
            f"- Partial candidates: {counts['partial_candidates']}",
            f"- Failed candidates: {counts['failed_candidates']}",
            f"- Invalid input candidates: {counts['invalid_input_candidates']}",
            f"- Skipped upstream candidates: {counts['skipped_upstream_candidates']}",
        ]
    )
    handoff = manifest.get("pipeline1_handoff") or {}
    if handoff:
        lines.extend(
            [
                f"- Pipeline 1 total candidates: {handoff.get('total_candidate_count', 0)}",
                f"- ADMET eligible candidates: {handoff.get('ready_candidate_count', 0)}",
                f"- Structure-blocked candidates: {handoff.get('blocked_candidate_count', 0)}",
                f"- Pipeline 1 failed candidates: {handoff.get('failed_candidate_count', 0)}",
                f"- Missing SMILES candidates: {handoff.get('missing_smiles_count', 0)}",
            ]
        )
    lines.extend(["", "## Candidate Endpoint Summary", ""])
    endpoint_columns = (
        [
            name
            for name in rows[0].keys()
            if name
            not in {
                "candidate_id",
                "smiles",
                "prediction_status",
                "backend",
                "is_mock",
                "available_endpoint_count",
                "missing_endpoint_count",
                "warning_count",
                "error_message",
                "validation_warnings",
                "canonical_smiles",
                "structural_call",
                "workflow_status",
                "backend_version",
                "model_version",
                "rdkit_parse_status",
                "molecular_weight",
                "validation_message",
                "prediction_started_at",
                "prediction_finished_at",
                "skip_reason",
                "pipeline1_run_id",
                "pipeline1_status",
                "handoff_status",
                "eligible_for_admet",
                "handoff_reason",
                "selected_structure_path",
                "pocket_recovery_fraction",
                "cross_seed_rmsd",
                "reference_pose_rmsd",
                "pipeline1_warning",
            }
        ]
        if rows
        else []
    )
    table_columns = _candidate_table_columns(endpoint_columns)
    lines.extend(_markdown_table(rows, table_columns))
    failures = [
        row
        for row in rows
        if row.get("prediction_status")
        in {
            PredictionStatus.INVALID_INPUT.value,
            PredictionStatus.PREDICTION_FAILED.value,
            PredictionStatus.PARTIAL.value,
            PredictionStatus.SKIPPED_UPSTREAM.value,
        }
    ]
    lines.extend(["", "## Failed Or Partial Candidates", ""])
    if failures:
        lines.extend(
            _markdown_table(
                failures, ["candidate_id", "prediction_status", "handoff_reason", "error_message"]
            )
        )
    else:
        lines.append("No failed or partial candidates.")
    lines.extend(
        [
            "",
            "## Current Limitations",
            "",
            "- Mock values are deterministic synthetic placeholders.",
            "- No docking, selectivity analysis, toxicity adjudication, scoring, or ranking.",
            "- ADMET-AI endpoint names require configured aliases to match returned names.",
            "- Extension point: local ADMET-AI can replace the mock backend via the adapter.",
            "",
        ]
    )
    return "\n".join(lines)


def _candidate_table_columns(endpoint_columns: list[str]) -> list[str]:
    columns = [
        "candidate_id",
        "pipeline1_status",
        "handoff_status",
        "prediction_status",
        *endpoint_columns,
        "handoff_reason",
        "error_message",
    ]
    return columns


def _markdown_table(rows: list[dict[str, str]], columns: list[str]) -> list[str]:
    if not rows:
        return ["No rows."]
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(_clean(row.get(column, "")) for column in columns) + " |")
    return lines


def _clean(value: str) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")

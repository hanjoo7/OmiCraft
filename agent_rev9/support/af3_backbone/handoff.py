from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .reporting import write_csv

SUPPORTED_STATUS = "AF3-supported pose"
UNCERTAIN_STATUS = "structurally uncertain"
REJECTED_STATUS = "rejected by structural QC"
FAILED_STATUS = "execution_or_analysis_failed"
DEFAULT_ELIGIBLE_STATUSES = [SUPPORTED_STATUS]
KNOWN_PIPELINE1_STATUSES = {
    SUPPORTED_STATUS,
    UNCERTAIN_STATUS,
    REJECTED_STATUS,
    FAILED_STATUS,
}
SOURCE_COLUMNS = [
    "candidate_id",
    "parent_id",
    "variant_id",
    "display_name",
    "smiles",
    "ccd_code",
    "role",
    "original_smiles",
    "canonical_smiles",
    "prepared_smiles",
    "model_input_smiles",
    "chemical_state_source",
    "chemical_state_id",
    "formal_charge",
    "protonation_status",
    "protonation_method",
    "target_ph",
    "stereochemistry_status",
    "chemical_state_warning",
]
SUMMARY_COLUMNS = {
    "pipeline1_status": "structural_call",
    "structural_call": "structural_call",
    "pocket_recovery_fraction": "clash_free_prediction_fraction",
    "cross_seed_rmsd": "cross_seed_rmsd",
    "reference_pose_rmsd": "reference_pose_rmsd",
    "pipeline1_warning": "warning",
}
HANDOFF_MIN_COLUMNS = [
    "candidate_id",
    "parent_id",
    "variant_id",
    "smiles",
    "role",
    "pipeline1_status",
    "structural_call",
    "handoff_status",
    "eligible_for_admet",
    "handoff_reason",
    "pipeline1_run_id",
    "selected_structure_path",
    "chemical_state_id",
    "model_input_smiles",
    "formal_charge",
    "protonation_status",
    "protonation_method",
    "target_ph",
    "stereochemistry_status",
    "chemical_state_warning",
]
OPTIONAL_METRIC_COLUMNS = [
    "display_name",
    "ccd_code",
    "pocket_recovery_fraction",
    "cross_seed_rmsd",
    "reference_pose_rmsd",
    "pipeline1_warning",
]


@dataclass(frozen=True)
class HandoffResult:
    rows: list[dict[str, Any]]
    manifest: dict[str, Any]
    handoff_file: Path
    manifest_file: Path


def eligible_statuses_from_config(config: dict[str, Any]) -> list[str]:
    statuses = config.get("handoff", {}).get("eligible_pipeline1_statuses")
    if statuses is None:
        return list(DEFAULT_ELIGIBLE_STATUSES)
    if not isinstance(statuses, list) or not all(isinstance(item, str) for item in statuses):
        raise ValueError("handoff.eligible_pipeline1_statuses must be a list of strings")
    unknown = sorted(set(statuses) - KNOWN_PIPELINE1_STATUSES)
    if unknown:
        raise ValueError(f"Unknown eligible Pipeline 1 statuses: {', '.join(unknown)}")
    return statuses


def export_handoff(
    *,
    candidates_path: str | Path,
    pipeline1_run: str | Path,
    eligible_pipeline1_statuses: list[str],
    force: bool = False,
) -> HandoffResult:
    candidates_path = Path(candidates_path)
    pipeline1_run = Path(pipeline1_run)
    summary_path = pipeline1_run / "candidate_summary.csv"
    run_manifest_path = pipeline1_run / "run_manifest.json"
    handoff_dir = pipeline1_run / "handoff"
    handoff_file = handoff_dir / "admet_candidates.csv"
    handoff_manifest_file = handoff_dir / "handoff_manifest.json"

    if handoff_file.exists() and not force:
        raise FileExistsError(f"{handoff_file} already exists. Use --force to overwrite.")
    if handoff_manifest_file.exists() and not force:
        raise FileExistsError(f"{handoff_manifest_file} already exists. Use --force to overwrite.")

    candidates = _read_csv(candidates_path)
    summaries = _read_csv(summary_path)
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    pipeline1_run_id = _pipeline1_run_id(pipeline1_run, run_manifest)
    _reject_duplicates(candidates, "source candidates")
    _reject_duplicates(summaries, "Pipeline 1 candidate_summary")

    source_by_id = {row["candidate_id"]: row for row in candidates}
    summary_by_id = {row["candidate_id"]: row for row in summaries}
    all_ids = list(source_by_id)
    for candidate_id in summary_by_id:
        if candidate_id not in source_by_id:
            all_ids.append(candidate_id)

    rows = []
    for candidate_id in all_ids:
        source = source_by_id.get(candidate_id, {})
        summary = summary_by_id.get(candidate_id, {})
        row = {name: source.get(name, "") for name in SOURCE_COLUMNS}
        row["candidate_id"] = candidate_id
        for target, source_name in SUMMARY_COLUMNS.items():
            row[target] = summary.get(source_name, "") or "not_available"
        row["selected_structure_path"] = _selected_structure_path(summary)
        row["pipeline1_run_id"] = pipeline1_run_id
        status, eligible, reason = _handoff_state(
            source=source,
            summary=summary,
            eligible_pipeline1_statuses=eligible_pipeline1_statuses,
        )
        row["handoff_status"] = status
        row["eligible_for_admet"] = str(eligible).lower()
        row["handoff_reason"] = reason
        rows.append(row)

    handoff_dir.mkdir(parents=True, exist_ok=True)
    write_csv(handoff_file, rows)
    manifest = _manifest(
        pipeline1_run_id=pipeline1_run_id,
        candidates_path=candidates_path,
        summary_path=summary_path,
        handoff_file=handoff_file,
        eligible_pipeline1_statuses=eligible_pipeline1_statuses,
        rows=rows,
    )
    handoff_manifest_file.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return HandoffResult(rows, manifest, handoff_file, handoff_manifest_file)


def _read_csv(path: Path) -> list[dict[str, str]]:
    import csv

    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _reject_duplicates(rows: list[dict[str, str]], label: str) -> None:
    counts = Counter(row.get("candidate_id", "") for row in rows)
    duplicates = sorted(candidate_id for candidate_id, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate candidate_id in {label}: {', '.join(duplicates)}")


def _pipeline1_run_id(pipeline1_run: Path, run_manifest: dict[str, Any]) -> str:
    return str(
        run_manifest.get("run_id")
        or run_manifest.get("pipeline1_run_id")
        or run_manifest.get("created_at")
        or pipeline1_run.name
    )


def _selected_structure_path(summary: dict[str, str]) -> str:
    for key in ["selected_structure_path", "output_structure_path", "best_structure_path"]:
        if summary.get(key):
            return summary[key]
    return "not_available"


def _handoff_state(
    *,
    source: dict[str, str],
    summary: dict[str, str],
    eligible_pipeline1_statuses: list[str],
) -> tuple[str, bool, str]:
    if not source:
        return "upstream_failed", False, "candidate missing from source candidates"
    if not (source.get("model_input_smiles") or source.get("smiles") or "").strip():
        return "missing_smiles", False, "smiles is required for ADMET"
    if not summary:
        return "upstream_failed", False, "candidate missing from Pipeline 1 summary"
    pipeline1_status = summary.get("structural_call", "")
    if not pipeline1_status:
        return "upstream_failed", False, "Pipeline 1 status missing"
    if pipeline1_status == FAILED_STATUS:
        return "upstream_failed", False, "Pipeline 1 execution or analysis failed"
    if pipeline1_status in eligible_pipeline1_statuses:
        return "ready", True, "eligible Pipeline 1 status"
    if pipeline1_status in {UNCERTAIN_STATUS, REJECTED_STATUS, SUPPORTED_STATUS}:
        return "blocked_by_structure", False, f"Pipeline 1 status not eligible: {pipeline1_status}"
    return "upstream_failed", False, f"unrecognized Pipeline 1 status: {pipeline1_status}"


def _manifest(
    *,
    pipeline1_run_id: str,
    candidates_path: Path,
    summary_path: Path,
    handoff_file: Path,
    eligible_pipeline1_statuses: list[str],
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    counts = Counter(row["handoff_status"] for row in rows)
    return {
        "schema_version": "1.0",
        "source_pipeline": "pipeline1_af3_structure_qc",
        "pipeline1_run_id": pipeline1_run_id,
        "source_candidates_path": str(candidates_path),
        "source_candidate_summary": str(summary_path),
        "source_candidate_summary_path": str(summary_path),
        "handoff_csv": str(handoff_file),
        "handoff_file": str(handoff_file),
        "eligible_pipeline1_statuses": eligible_pipeline1_statuses,
        "candidate_count": len(rows),
        "total_candidate_count": len(rows),
        "eligible_count": counts.get("ready", 0),
        "ready_candidate_count": counts.get("ready", 0),
        "blocked_candidate_count": counts.get("blocked_by_structure", 0),
        "failed_candidate_count": counts.get("upstream_failed", 0),
        "missing_smiles_count": counts.get("missing_smiles", 0),
        "skipped_count": len(rows) - counts.get("ready", 0),
        "handoff_csv_sha256": _sha256(handoff_file),
        "created_at": datetime.now(UTC).isoformat(),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


HANDOFF_FIELDNAMES = HANDOFF_MIN_COLUMNS + OPTIONAL_METRIC_COLUMNS

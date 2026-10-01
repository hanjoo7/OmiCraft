from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .io import read_csv

REQUIRED_HANDOFF_COLUMNS = {
    "candidate_id",
    "smiles",
    "pipeline1_status",
    "handoff_status",
    "eligible_for_admet",
    "pipeline1_run_id",
}
REQUIRED_MANIFEST_FIELDS = {
    "schema_version",
    "pipeline1_run_id",
    "created_at",
    "source_candidate_summary",
    "handoff_csv",
    "candidate_count",
    "eligible_count",
    "skipped_count",
    "handoff_csv_sha256",
}
VALID_ELIGIBLE_VALUES = {"true", "false"}
READY_STATUS = "ready"
UPSTREAM_COLUMNS = [
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
]


@dataclass(frozen=True)
class ImportedHandoff:
    input_path: Path
    manifest_path: Path
    manifest: dict[str, Any]
    rows: list[dict[str, str]]

    @property
    def pipeline1_run_id(self) -> str:
        return str(self.manifest["pipeline1_run_id"])


def resolve_handoff_input(pipeline1_run: str | Path) -> ImportedHandoff:
    pipeline1_run = Path(pipeline1_run)
    input_path = pipeline1_run / "handoff" / "admet_candidates.csv"
    manifest_path = pipeline1_run / "handoff" / "handoff_manifest.json"
    return load_handoff(input_path=input_path, manifest_path=manifest_path)


def load_handoff(input_path: str | Path, manifest_path: str | Path) -> ImportedHandoff:
    input_path = Path(input_path)
    manifest_path = Path(manifest_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Pipeline 1 handoff CSV not found: {input_path}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Pipeline 1 handoff manifest not found: {manifest_path}")
    rows, fieldnames = read_csv(input_path)
    missing = REQUIRED_HANDOFF_COLUMNS - set(fieldnames)
    if missing:
        raise ValueError(f"Handoff CSV missing required columns: {', '.join(sorted(missing))}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _validate_manifest_fields(manifest)
    _validate_handoff(rows, manifest, input_path, manifest_path)
    return ImportedHandoff(
        input_path=input_path, manifest_path=manifest_path, manifest=manifest, rows=rows
    )


def _validate_manifest_fields(manifest: dict[str, Any]) -> None:
    missing = REQUIRED_MANIFEST_FIELDS - set(manifest)
    if missing:
        raise ValueError(f"Handoff manifest missing required fields: {', '.join(sorted(missing))}")
    if manifest.get("schema_version") != "1.0":
        raise ValueError(f"Unsupported handoff schema_version: {manifest.get('schema_version')}")


def _validate_handoff(
    rows: list[dict[str, str]], manifest: dict[str, Any], input_path: Path, manifest_path: Path
) -> None:
    run_id = str(manifest.get("pipeline1_run_id", ""))
    if not run_id:
        raise ValueError("Handoff manifest missing pipeline1_run_id")
    manifest_csv = _resolve_manifest_path(str(manifest["handoff_csv"]), manifest_path)
    if manifest_csv.exists() and _sha256(manifest_csv) != _sha256(input_path):
        raise ValueError("Handoff manifest handoff_csv points to different CSV content")
    if str(manifest["handoff_csv_sha256"]) != _sha256(input_path):
        raise ValueError("Handoff CSV checksum does not match manifest")
    counts = Counter(row.get("candidate_id", "") for row in rows)
    duplicates = sorted(candidate_id for candidate_id, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate candidate_id in handoff CSV: {', '.join(duplicates)}")
    candidate_count = int(manifest["candidate_count"])
    if candidate_count != len(rows):
        raise ValueError(
            "Handoff manifest candidate_count does not match CSV row count: "
            f"{candidate_count} != {len(rows)}"
        )
    eligible_count = int(manifest["eligible_count"])
    actual_eligible = sum(1 for row in rows if row.get("eligible_for_admet", "").lower() == "true")
    if eligible_count != actual_eligible:
        raise ValueError(
            "Handoff manifest eligible_count does not match CSV eligible rows: "
            f"{eligible_count} != {actual_eligible}"
        )
    skipped_count = int(manifest["skipped_count"])
    if skipped_count != len(rows) - actual_eligible:
        raise ValueError(
            "Handoff manifest skipped_count does not match CSV skipped rows: "
            f"{skipped_count} != {len(rows) - actual_eligible}"
        )
    for index, row in enumerate(rows, start=2):
        if row.get("pipeline1_run_id") != run_id:
            raise ValueError(f"row {index}: pipeline1_run_id does not match manifest")
        eligible = (row.get("eligible_for_admet") or "").lower()
        if eligible not in VALID_ELIGIBLE_VALUES:
            raise ValueError(f"row {index}: eligible_for_admet must be true or false")
        if row.get("handoff_status") == READY_STATUS and eligible != "true":
            raise ValueError(f"row {index}: ready handoff row must be eligible_for_admet=true")
        if eligible == "true" and row.get("handoff_status") != READY_STATUS:
            raise ValueError(f"row {index}: eligible handoff row must have handoff_status=ready")
        if row.get("handoff_status") == READY_STATUS and not row.get("smiles", "").strip():
            raise ValueError(f"row {index}: ready handoff row is missing SMILES")


def _resolve_manifest_path(path_text: str, manifest_path: Path) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return (manifest_path.parent / path).resolve()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_ready_for_admet(row: dict[str, str]) -> bool:
    return (
        row.get("eligible_for_admet", "").lower() == "true"
        and row.get("handoff_status") == READY_STATUS
    )

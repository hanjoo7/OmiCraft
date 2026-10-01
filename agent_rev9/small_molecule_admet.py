"""ADMET execution with a checksum-validated structure handoff."""

import csv
import json
import subprocess
from pathlib import Path

from .model_runtime import gpu_environment
from .screening_model_tools import ADMETTool
from .small_molecule_io import (
    finite,
    now,
    relative_path,
    sha256,
    unavailable,
    write_csv,
    write_json,
)
from .support.admet_pipeline.handoff import (
    is_ready_for_admet,
    load_handoff,
)
from .support.af3_backbone.handoff import export_handoff


def eligible_for_property_prediction(candidate):
    workflow = candidate.get("design_workflow") or {}
    return (
        candidate["route"] == "small_molecule"
        and candidate["rdkit_validation"]["smiles_valid"]
        and candidate.get("qualification", {}).get("safety_verdict") != "REJECT"
        and not workflow.get("safety_veto")
        and (candidate["validation"]["eligible_for_admet"]
             or candidate.get("admet_execution_policy") == "valid_input")
    )


def export_candidate_handoff(candidate, directory, run_id):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    source = write_csv(
        directory / "input_candidates.csv",
        [
            {
                k: candidate.get(k)
                for k in (
                    "candidate_id",
                    "smiles",
                    "canonical_smiles",
                    "prepared_smiles",
                    "display_name",
                )
            }
        ],
    )
    eligible = eligible_for_property_prediction(candidate)
    structurally_supported = candidate["validation"]["eligible_for_admet"]
    policy = candidate.get("admet_execution_policy", "structural_pass")
    write_csv(
        directory / "candidate_summary.csv",
        [
            {
                "candidate_id": candidate["candidate_id"],
                "structural_call": "AF3-supported pose" if structurally_supported else "structurally uncertain",
            }
        ],
    )
    write_json(directory / "run_manifest.json", {"run_id": run_id, "route": candidate["route"]})
    result = export_handoff(
        candidates_path=source,
        pipeline1_run=directory,
        eligible_pipeline1_statuses=(["AF3-supported pose", "structurally uncertain"]
                                     if eligible and policy == "valid_input" else ["AF3-supported pose"]),
    )
    rows = [{**row, "route": candidate["route"]} for row in result.rows]
    write_csv(result.handoff_file, rows)
    manifest = {
        **result.manifest,
        "route": candidate["route"],
        "admet_execution_policy": policy,
        "structural_eligible_for_admet": structurally_supported,
        "handoff_csv_sha256": sha256(result.handoff_file),
    }
    for field in (
        "source_candidates_path",
        "source_candidate_summary",
        "source_candidate_summary_path",
        "handoff_csv",
        "handoff_file",
    ):
        if manifest.get(field):
            manifest[field] = relative_path(manifest[field], result.manifest_file)
    write_json(result.manifest_file, manifest)

    imported = load_handoff(result.handoff_file, result.manifest_file)
    if (
        imported.rows[0]["route"] != "small_molecule"
        or is_ready_for_admet(imported.rows[0]) != eligible
    ):
        raise ValueError("ADMET handoff route/eligibility mismatch")
    return result


class SmallMoleculeADMETTool:
    def __init__(self, options=None, *, gpu_device=0):
        self.existing = ADMETTool(options, gpu_device=gpu_device)
        self.gpu_device = gpu_device

    def run(self, candidate, output_dir, *, run_id, is_mock=False):
        if candidate["route"] != "small_molecule":
            return unavailable("skipped_route", is_mock=is_mock)
        if not eligible_for_property_prediction(candidate):
            return unavailable("skipped_structural_screening", is_mock=is_mock)
        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        handoff = export_candidate_handoff(candidate, directory / "structure_handoff", run_id)
        if is_mock:
            from .support.admet_pipeline.mock_predictor import (
                MockADMETPredictor,
            )

            predictor = MockADMETPredictor()
            records = predictor.predict(
                [
                    {
                        "candidate_id": candidate["candidate_id"],
                        "smiles": candidate["smiles"],
                    }
                ]
            )
            rows = [{**record.__dict__, "is_mock": True} for record in records]
            return {
                "status": "success",
                "backend": "mock",
                "backend_version": predictor.backend_version,
                "model_version": predictor.model_version,
                "is_mock": True,
                "endpoints": rows,
                "metrics": {
                    row["endpoint_name"]: finite(float(row["endpoint_value"]))
                    for row in rows
                    if row["endpoint_value"] not in (None, "")
                },
                "available_endpoint_count": len(rows),
                "missing_endpoint_count": 0,
                "handoff_manifest": str(handoff.manifest_file),
                "warning": ["SYNTHETIC_ADMET_NOT_A_PREDICTION"],
                "started_at": now(),
                "finished_at": now(),
                "error_message": "",
            }
        configured, reason = self.existing.is_configured()
        if not configured:
            return unavailable(
                "dependency_missing",
                reason,
                handoff_manifest=str(handoff.manifest_file),
            )
        c = self.existing.config
        predictions = directory / "predictions"
        command = [
            c.python_path,
            c.script_path,
            "--config",
            c.config_path,
            "--input-mode",
            "pipeline1_handoff",
            "--pipeline1-run",
            str((directory / "structure_handoff").resolve()),
            "--output",
            str(predictions.resolve()),
            "--backend",
            "admet_ai",
            "--stage",
            "all",
        ]
        started = now()
        try:
            with (directory / "admet.log").open("w") as handle:
                result = subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=c.timeout_seconds,
                    env=gpu_environment(self.gpu_device, use_gpu=c.use_gpu),
                )
            if result.returncode:
                raise ValueError(f"ADMET exited {result.returncode}")
            parsed = self.parse_output(predictions, candidate["candidate_id"])
        except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as exc:
            parsed = unavailable("prediction_failed", str(exc))
        return {
            **parsed,
            "command": command,
            "started_at": started,
            "finished_at": now(),
            "handoff_manifest": str(handoff.manifest_file),
            "is_mock": False,
        }

    @staticmethod
    def parse_output(directory, candidate_id):
        directory = Path(directory)
        manifest = json.loads((directory / "run_manifest.json").read_text())
        if manifest.get("is_mock") is not False or manifest.get("backend") != "admet_ai":
            raise ValueError("ADMET backend provenance mismatch")
        with (directory / "candidate_summary.csv").open() as handle:
            summaries = [r for r in csv.DictReader(handle) if r["candidate_id"] == candidate_id]
        if len(summaries) != 1:
            raise ValueError("ADMET candidate summary missing or duplicated")
        summary = summaries[0]
        with (directory / "admet_predictions_long.csv").open() as handle:
            rows = [r for r in csv.DictReader(handle) if r["candidate_id"] == candidate_id]
        metrics = {}
        warnings = []
        for row in rows:
            if row.get("is_mock", "").lower() != "false":
                raise ValueError("Mock ADMET endpoint in real run")
            try:
                value = finite(float(row.get("endpoint_value", "")))
            except (TypeError, ValueError):
                value = None
            row["endpoint_value"] = value
            if row.get("prediction_status") == "success" and value is not None:
                if row["endpoint_name"] in metrics:
                    raise ValueError("Duplicate ADMET endpoint")
                metrics[row["endpoint_name"]] = value
        available = len(metrics)
        if (
            summary.get("available_endpoint_count")
            and int(summary["available_endpoint_count"]) != available
        ):
            warnings.append("ADMET_ENDPOINT_COUNT_MISMATCH")
        missing = (
            int(summary["missing_endpoint_count"])
            if summary.get("missing_endpoint_count")
            else None
        )
        status = summary.get("prediction_status", "prediction_failed")
        if not metrics:
            status = "prediction_failed"
        elif missing is None or missing or status == "partial" or warnings:
            status = "partial"
        elif status not in {"success", "partial"}:
            status = "prediction_failed"
        return {
            "status": status,
            "backend": "admet_ai",
            "backend_version": manifest.get("backend_version"),
            "model_version": manifest.get("model_version"),
            "is_mock": False,
            "metrics": metrics,
            "endpoints": rows,
            "available_endpoint_count": available,
            "missing_endpoint_count": missing,
            "existing_summary": summary,
            "warning": warnings,
            "error_message": summary.get("error_message", ""),
        }

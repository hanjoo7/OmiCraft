from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from ..pipeline_common.preparation import validate_preparation
from .config import load_config
from .handoff import ImportedHandoff, resolve_handoff_input
from .input import prepare_candidates
from .prediction import make_predictor, run_predictions
from .reporting import write_report
from .summary import summarize_predictions

STAGES = {"prepare", "predict", "summarize", "report", "all"}
INPUT_MODES = {"standalone", "shared_preparation", "pipeline1_handoff"}


def run_pipeline(
    *,
    config_path: str | Path,
    input_path: str | Path | None,
    output_dir: str | Path,
    input_mode: str = "standalone",
    pipeline1_run: str | Path | None = None,
    preparation_manifest: str | Path | None = None,
    backend_override: str | None,
    stage: str,
    force: bool = False,
) -> dict[str, Any]:
    if stage not in STAGES:
        raise ValueError(f"Unsupported stage: {stage}")
    if input_mode not in INPUT_MODES:
        raise ValueError(f"Unsupported ADMET input mode: {input_mode}")
    imported_handoff: ImportedHandoff | None = None
    if input_mode == "shared_preparation":
        if input_path is None or preparation_manifest is None:
            raise ValueError("shared_preparation requires input_path and preparation_manifest")
        validate_preparation(input_path, preparation_manifest)
    elif input_mode == "pipeline1_handoff":
        if input_path is not None:
            raise ValueError("--input and --pipeline1-run cannot be used together")
        if pipeline1_run is None:
            raise ValueError("--pipeline1-run is required when --input-mode pipeline1_handoff")
        imported_handoff = resolve_handoff_input(pipeline1_run)
        input_path = imported_handoff.input_path
    elif input_path is None:
        raise ValueError("--input is required when --input-mode standalone")
    config = load_config(config_path)
    backend = backend_override or config.backend
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    started_at = datetime.now(UTC).isoformat()
    stages = ["prepare", "predict", "summarize", "report"] if stage == "all" else [stage]
    predictor = None

    if "prepare" in stages and input_mode != "shared_preparation":
        prepared = prepare_candidates(
            input_path,
            output,
            candidate_id_column=config.candidate_id_column,
            smiles_column=config.smiles_column,
            force=force,
        )
    else:
        prepared = None
        if input_mode == "shared_preparation":
            target = output / "prepared_candidates.csv"
            target.write_bytes(Path(input_path).read_bytes())
            invalid_target = output / "invalid_candidates.csv"
            if not invalid_target.exists() or force:
                invalid_target.write_text(
                    "candidate_id,variant_id,chemical_state_id\n", encoding="utf-8"
                )

    if "predict" in stages:
        predictions, predictor = run_predictions(
            output / "prepared_candidates.csv",
            output,
            backend=backend,
            batch_size=config.batch_size,
            device=config.device,
            force=force,
        )
    else:
        predictions = []
        if "report" in stages:
            predictor = make_predictor(backend, batch_size=config.batch_size, device=config.device)

    if "summarize" in stages:
        summaries = summarize_predictions(
            output / "admet_predictions_long.csv",
            output / "prepared_candidates.csv",
            output / "invalid_candidates.csv",
            output,
            selected_endpoints=config.selected_endpoints,
            force=force,
        )
    else:
        summaries = []

    if "report" in stages:
        if predictor is None:
            predictor = make_predictor(backend, batch_size=config.batch_size, device=config.device)
        manifest = write_report(
            output,
            config=config.raw,
            backend=predictor.backend,
            backend_version=predictor.backend_version,
            model_version=getattr(predictor, "model_version", "not_available"),
            is_mock=bool(predictor.is_mock),
            started_at=started_at,
            force=force,
            handoff_manifest=imported_handoff.manifest if imported_handoff else None,
            input_path=input_path,
        )
    else:
        manifest = {}

    return {
        "stage": stage,
        "backend": backend,
        "output_dir": str(output.resolve()),
        "input_mode": input_mode,
        "pipeline1_run": str(pipeline1_run) if pipeline1_run else None,
        "prepared_candidates": len(prepared.prepared) if prepared else None,
        "invalid_candidates": len(prepared.invalid) if prepared else None,
        "prediction_rows": len(predictions),
        "summary_rows": len(summaries),
        "manifest": manifest,
    }

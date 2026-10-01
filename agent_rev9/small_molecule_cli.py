"""Command-line interface for small-molecule screening."""

import argparse
import importlib.util
import json
from dataclasses import replace
from pathlib import Path

from .af3_ligand_tool import AF3LigandTool
from .docking_backend import resolve_executable, tool_readiness
from .input_paths import local_paths
from .ligand_preparation import validate_candidate
from .routes import SMALL_MOLECULE, validate_route
from .screening_tools import ScreeningOptions
from .small_molecule_admet import SmallMoleculeADMETTool
from .small_molecule_config import SmallMoleculeConfig
from .small_molecule_io import redact, write_json
from .small_molecule_workflow import run_small_molecule_workflow


def load_job(input_path, config_path):
    payload = local_paths(json.loads(input_path.read_text()), input_path.resolve().parent)
    if validate_route(payload.get("route")) != SMALL_MOLECULE:
        raise ValueError("Use the existing binder CLI for protein_binder_rfd3")
    candidates = payload["candidates"]
    for candidate in candidates:
        if validate_route(candidate.get("route", payload["route"])) != SMALL_MOLECULE:
            raise ValueError("Mixed routes are not accepted by this ligand CLI")
    if not candidates:
        raise ValueError("At least one candidate is required")
    ids = [c.get("candidate_id") for c in candidates]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate candidate_id")
    raw = local_paths(json.loads(config_path.read_text()), config_path.resolve().parent)
    unknown = set(raw) - {"small_molecule", "screening_options"}
    if unknown:
        raise ValueError("Unknown CLI configuration keys: " + ", ".join(sorted(unknown)))
    config = SmallMoleculeConfig(**raw.get("small_molecule", {}))
    if config.structural_backend != "af3":
        raise ValueError(
            "This CLI implements AF3; legacy Boltz remains available through Screening"
        )
    options = ScreeningOptions.from_state({"screening_options": raw.get("screening_options", {})})
    return payload, config, options


def preflight(payload, config, options):
    """Read-only input checks and executable version checks; no model inference."""
    context = payload["context"]
    prepared = context.get("prepared_receptor_path")
    available, reason = SmallMoleculeADMETTool(options.admet_options).existing.is_configured()
    return {
        "route": SMALL_MOLECULE,
        "inference_started": False,
        "is_mock": config.docking.backend == "mock",
        "docking": "explicit_mock"
        if config.docking.backend == "mock"
        else "available"
        if resolve_executable(config.docking.backend, configured=config.docking.executable)["path"]
        else "dependency_missing",
        "docking_readiness": ({"tool_status": "explicit_mock", "is_mock": True}
            if config.docking.backend == "mock" else tool_readiness(
                config.docking.backend, configured=config.docking.executable,
                environment=config.docking.environment)),
        "prepared_receptor": "available"
        if prepared and Path(prepared).is_file()
        else "dependency_missing",
        "af3": AF3LigandTool(config.af3).dependency_status(),
        "admet": "configured_presence_only" if available else "dependency_missing",
        "admet_message": reason,
        "optional_tools": {
            name: "available" if importlib.util.find_spec(name) else "dependency_missing"
            for name in ("rdkit", "posebusters", "prolif")
        },
        "candidates": [
            validate_candidate(c, seed=config.docking.seed) for c in payload["candidates"]
        ],
        "threshold_status": config.thresholds.threshold_status,
        "warning": [
            "PRESENCE_CHECK_ONLY_NOT_BACKEND_EXECUTION",
            "NULL_THRESHOLDS_DO_NOT_ESTABLISH_STRUCTURAL_PASS",
        ],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--execute",
        action="store_true",
        help="Run configured real backends, with no downloads",
    )
    mode.add_argument(
        "--mock",
        action="store_true",
        help="Explicit synthetic backend; never real success",
    )
    mode.add_argument(
        "--preflight", action="store_true", help="Read-only dependency and RDKit checks"
    )
    parser.add_argument("--gpu-device", type=int)
    parser.add_argument("--prepared-receptor", type=Path)
    parser.add_argument("--docking-executable", "--gnina-executable", dest="docking_executable")
    args = parser.parse_args(argv)
    record = None
    try:
        payload, config, options = load_job(args.input, args.config)
        if args.gpu_device is not None:
            if args.gpu_device < 0:
                raise ValueError("GPU device must be nonnegative")
            options = replace(options, gpu_device=args.gpu_device)
        if args.prepared_receptor:
            payload["context"]["prepared_receptor_path"] = str(args.prepared_receptor.resolve())
        if args.docking_executable:
            config.docking.executable = args.docking_executable
        if args.mock:
            config.docking.backend = "mock"
        if args.execute and config.docking.backend == "mock":
            raise ValueError("Use --mock explicitly for synthetic execution")
        if args.preflight:
            print(
                json.dumps(
                    redact(preflight(payload, config, options)),
                    indent=2,
                    ensure_ascii=False,
                )
            )
            return 0
        if args.output is None:
            raise ValueError("--output is required except with --preflight")
        if args.output.exists() and any(args.output.iterdir()):
            raise ValueError(
                "Output directory must be fresh; existing artifacts are never overwritten"
            )
        from .run_records import RunRecord
        record = RunRecord(args.output, "SMALL_MOLECULE")
        result = run_small_molecule_workflow(
            payload["candidates"],
            payload["context"],
            config,
            args.output.resolve(),
            options=options,
            dry_run=not (args.execute or args.mock),
            qualification=payload.get("qualification"),
        )
        from .critic_agent import review_modality
        from .modality_dispatch import result as modality_result
        from .validation_agent import validate_modality
        complete = bool(result["candidates"]) and all(
            c.get("docking", {}).get("status") == "success"
            and (not config.run_af3 or c.get("af3_ligand", {}).get("status") == "success")
            for c in result["candidates"])
        stored = modality_result("SMALL_MOLECULE", "NOT_RUN" if not (args.execute or args.mock)
            else "COMPLETED" if complete else "PARTIAL", "small_molecule", summary=result,
            artifacts=[{"path": str(p.resolve()), "label": p.name} for p in args.output.rglob("*")
                       if p.is_file() and p.name not in {"events.jsonl", "run_summary.json"}],
            is_mock=result["is_mock"])
        stored["validation"] = validate_modality(stored)
        stored["validation_decision"] = stored["validation"]["decision"]
        stored["critic"] = review_modality(stored)
        record.finish(stored)
        write_json(args.output / "input_manifest.json", redact(payload))
        print(
            json.dumps(
                {
                    "candidate_summary_path": result["candidate_summary_path"],
                    "is_mock": result["is_mock"],
                    "candidates": [
                        {
                            "candidate_id": c["candidate_id"],
                            "workflow_status": c["workflow_status"],
                            "critic_verdict": c["critic"]["verdict"],
                        }
                        for c in result["candidates"]
                    ],
                },
                indent=2,
            )
        )
        return (
            0
            if all(c["workflow_status"] in {"complete", "partial"} for c in result["candidates"])
            or not args.execute
            else 2
        )
    except (OSError, ValueError, TypeError, KeyError) as exc:
        if record is not None:
            from .modality_dispatch import result as modality_result
            record.finish(modality_result("SMALL_MOLECULE", "FAILED", "execution", blockers=[str(exc)]))
        parser.exit(2, f"Input/configuration error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())

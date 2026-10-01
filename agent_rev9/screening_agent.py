"""Agent ⑦ Screening: one orchestrator, two deterministic tool cascades.

Binder: RFD3 → ProteinMPNN → AF3 → interface QC → Validation → Critic.
Small molecule: RDKit → Vina/GNINA → pose QC → AF3 → consensus → ADMET.
An explicit boltz2_legacy option retains the previous small-molecule cascade.
No LLM verdicts, fabricated model scores, or additional agents.
"""

import json
import logging
import math
import re
import tempfile
from dataclasses import asdict, replace
from pathlib import Path

from .configuration import get_config
from .input_paths import local_paths
from .protein_mpnn_tool import ProteinMPNNTool
from .screening_gates import PoseScreeningTool, verdict_result
from .screening_model_tools import ADMETTool, Boltz2Tool, GNINATool, skipped
from .screening_tools import (
    BINDER_MODALITY,
    PIPELINE,
    SMALL_MOLECULE_PIPELINE,
    AF3Tool,
    RFD3Tool,
    ScreeningOptions,
    StructuralScreeningTool,
    execution_status,
    utc_now,
)
from .state import OmiCraftState

logger = logging.getLogger(__name__)


# ── 자원 프로파일 정의 ───────────────────────────────────
RESOURCE_PROFILES = {
    "A": {
        "name": "80GB GPU × 1",
        "vram_gb": 80,
        "gpu_budget_h": 6.0,
        "can_boltz2": True,
        "can_boltzina": True,
        "can_protenix": True,
        "can_rfdiffusion": True,
        "can_alphafold": True,
        "batch_size": "large",
    },
    "B": {
        "name": "48GB GPU × 2",
        "vram_gb": 48,
        "gpu_budget_h": 6.0,
        "can_boltz2": True,
        "can_boltzina": True,
        "can_protenix": False,
        "can_rfdiffusion": True,
        "can_alphafold": True,
        "batch_size": "reduced",
    },
    "C": {
        "name": "GPU 없음 (CPU only)",
        "vram_gb": 0,
        "gpu_budget_h": 0,
        "can_boltz2": False,
        "can_boltzina": False,
        "can_protenix": False,
        "can_rfdiffusion": False,
        "can_alphafold": False,
        "batch_size": "cpu_only",
    },
}


def _design_handoff(state: dict) -> tuple[list[dict], list[str]]:
    """Read the Design agent's existing file schema; select current targets only."""
    path = Path(
        state.get("design_results_path")
        or (Path(get_config().data.agent_results) / "design_output.json")
    )
    errors = []
    designs = []
    if path.is_file():
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            designs = local_paths(data["designs"], path.resolve().parent)
            if not isinstance(designs, list) or not all(isinstance(d, dict) for d in designs):
                raise ValueError("designs must be a list of objects")
        except (ValueError, KeyError, TypeError, OSError) as exc:
            errors.append(f"Design handoff invalid: {exc}")
            designs = []
    elif state.get("design_results_path"):
        errors.append(f"Design handoff not found: {path}")
    advance = state.get("advance_targets", [])
    if advance:
        by_gene = {d.get("gene_name"): d for d in designs}
        designs = [
            {
                **target,
                **by_gene.get(target.get("gene_name"), {}),
                "gene_name": target.get("gene_name"),
                "modality": target.get("modality")
                or by_gene.get(target.get("gene_name"), {}).get("modality"),
            }
            for target in advance
        ]
    return designs, errors


def _empty_screening(state, *, status="NOT_EXECUTED", errors=()):
    root = Path(get_config().data.agent_results) / "screening"
    root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="screening_", dir=root))
    dry_run = state.get("dry_run", True)
    return {
        "status": status,
        "modality": "NONE",
        "target": "",
        "execution_mode": "mock" if dry_run else "real",
        "dry_run": dry_run,
        "is_mock": dry_run,
        "candidates": [],
        "targets": [],
        "pipeline": [],
        "summary": {"pass": 0, "review": 0, "fail": 0},
        "mock_summary": {"pass": 0, "review": 0, "fail": 0},
        "accepted_candidate_ids": [],
        "profiles": {},
        "binder_profiles": {},
        "profiles_execution_mode": "not_executed",
        "generated_designs": 0,
        "generated_rfd3": 0,
        "mpnn_completed": 0,
        "af3_evaluated": 0,
        "provenance": {
            "run_id": state.get("run_id") or directory.name,
            "output_dir": str(directory),
            "started_at": utc_now(),
            "ended_at": utc_now(),
        },
        "configuration_errors": list(errors),
    }


def screening_node(state: OmiCraftState) -> dict:
    """Single Agent: read Design, route, allocate, call tools, aggregate, persist."""
    from .routes import normalize_route_state, validate_handoff_routes

    state = normalize_route_state(state)
    if state.get("route") in {"adc", "degrader"}:
        from .modality_dispatch import execute
        cfg = {**state.get("modality_config", {}), "output_dir": state["output_dir"],
               "dry_run": state.get("dry_run", True)}
        return {"modality_result": execute({"modality": state["route"], "input": state.get("modality_input", {})}, cfg)}
    from .protein_design_route import protein_route_selected

    designs, handoff_errors = _design_handoff(state)
    if not state.get("advance_targets") and not state.get("design_results_path"):
        path = Path(get_config().data.agent_results) / "final_dossier.json"
        if path.is_file():
            try:
                restored = json.loads(path.read_text()).get("advance_targets", [])
                if restored:
                    state = {**state, "advance_targets": restored}
                    designs, handoff_errors = _design_handoff(state)
            except (ValueError, OSError) as exc:
                handoff_errors.append(str(exc))
    if not designs and protein_route_selected(state):
        designs.append({"gene_name": "configured_target", "modality": BINDER_MODALITY})
    try:
        options = ScreeningOptions.from_state(state)
    except (ValueError, TypeError) as exc:
        output = _empty_screening(state, status="NOT_CONFIGURED", errors=[str(exc)])
        return _persist_screening(
            state,
            output,
            _screening_evidence(output, Path(output["provenance"]["output_dir"]).name),
            handoff_errors + [str(exc)],
            ScreeningOptions(),
        )
    validate_handoff_routes(designs, state.get("route"))
    if any(d.get("modality") in {"ADC", "DEGRADER"} for d in designs):
        raise ValueError("Unsupported modality in legacy handoff; use explicit ADC/Degrader executor")

    # Validate file-based candidates too, before either route can launch a model.
    for design in designs:
        if design.get("modality") == "SMALL_MOLECULE":
            _, ligand_inputs = _small_molecule_inputs(design)
            validate_handoff_routes([{**design, "candidates": ligand_inputs}], state.get("route"))
    selected, skipped_targets = [], []
    for design in designs:
        modality = design.get("modality")
        if modality not in (BINDER_MODALITY, "SMALL_MOLECULE"):
            skipped_targets.append(
                {
                    "gene_name": design.get("gene_name"),
                    "modality": modality,
                    "status": "SKIPPED_MODALITY",
                }
            )
        elif len(selected) >= options.max_targets:
            skipped_targets.append(
                {
                    "gene_name": design.get("gene_name"),
                    "modality": modality,
                    "status": "SKIPPED_RESOURCE",
                }
            )
        else:
            selected.append(design)
    branches, errors = {}, list(handoff_errors)
    binders = [d for d in selected if d.get("modality") == BINDER_MODALITY]
    molecules = [d for d in selected if d.get("modality") == "SMALL_MOLECULE"]
    if binders:
        binder_output, branch_errors, _ = _screen_binders(state, binders, handoff_errors, {})
        branches["protein_binder"] = binder_output
        errors.extend(branch_errors)
    if molecules:
        molecule_options = options
        ceiling = state.get("budget", {}).get("max_candidates")
        if ceiling is not None and binders:
            allocated = sum(
                t.get("allocated_designs", 0) for t in branches["protein_binder"]["targets"]
            )
            molecule_options = replace(
                options,
                max_gnina_candidates=min(options.max_gnina_candidates, max(0, ceiling - allocated)),
            )
        molecule_output, branch_errors = _screen_small_molecules(
            state, molecules, handoff_errors, molecule_options
        )
        branches["small_molecule"] = molecule_output
        errors.extend(branch_errors)
    if not branches:
        output = _empty_screening(
            state,
            status="NOT_CONFIGURED"
            if handoff_errors
            else "SKIPPED_RESOURCE"
            if selected == [] and any(t["status"] == "SKIPPED_RESOURCE" for t in skipped_targets)
            else "NOT_EXECUTED",
        )
    else:
        output = dict(next(iter(branches.values())))
        output["branch_results"] = branches
        output["modality"] = "MIXED" if len(branches) > 1 else output["modality"]
        output["candidates"] = [c for branch in branches.values() for c in branch["candidates"]]
        output["targets"] = [t for branch in branches.values() for t in branch["targets"]]
        output["target"] = ", ".join(t["gene_name"] for t in output["targets"])
        output["summary"] = {
            v: sum(c["verdict"] == v.upper() and not c["is_mock"] for c in output["candidates"])
            for v in ("pass", "review", "fail")
        }
        output["mock_summary"] = {
            v: sum(c["verdict"] == v.upper() and c["is_mock"] for c in output["candidates"])
            for v in ("pass", "review", "fail")
        }
        output["accepted_candidate_ids"] = [
            c["candidate_id"]
            for c in output["candidates"]
            if c["verdict"] == "PASS" and not c["is_mock"]
        ]
        output["is_mock"] = bool(
            output["dry_run"] or any(c["is_mock"] for c in output["candidates"])
        )
        output["execution_mode"] = (
            "mock"
            if output["dry_run"]
            or (output["candidates"] and all(c["is_mock"] for c in output["candidates"]))
            else "real"
        )
        states = {b["status"] for b in branches.values()}
        if len(states) > 1:
            output["status"] = (
                "PARTIAL_EXECUTION"
                if any(b["status"] in ("COMPLETED", "PARTIAL_EXECUTION") for b in branches.values())
                else "FAILED_EXECUTION"
                if "FAILED_EXECUTION" in states
                else "NOT_CONFIGURED"
                if "NOT_CONFIGURED" in states
                else "NOT_EXECUTED"
            )
        output["pipelines"] = {name: branch["pipeline"] for name, branch in branches.items()}
        if "small_molecule" in branches:
            for key in (
                "validation_results",
                "candidate_rankings",
                "candidate_summary_paths",
            ):
                if key in branches["small_molecule"]:
                    output[key] = branches["small_molecule"][key]
            for key in (
                "gnina_attempted",
                "gnina_completed",
                "boltz2_attempted",
                "boltz2_completed",
                "admet_attempted",
                "admet_completed",
                "final_small_molecule_hits",
            ):
                output[key] = branches["small_molecule"][key]
    for branch_name, branch in branches.items():
        route = "protein_binder_rfd3" if branch_name == "protein_binder" else "small_molecule"
        branch["route"] = route
        for candidate in branch["candidates"]:
            candidate["route"] = route
        if route == "protein_binder_rfd3":
            from .small_molecule_io import write_json

            write_json(
                Path(branch["provenance"]["output_dir"]) / "run_manifest.json",
                {
                    "route": route,
                    "run_id": branch["provenance"]["run_id"],
                    "is_mock": branch["is_mock"],
                    "candidate_count": len(branch["candidates"]),
                    "pipeline": branch["pipeline"],
                    "admet_status": "skipped_route",
                },
            )
    if len(branches) == 1:
        output["route"] = next(iter(branches.values()))["route"]
    elif len(branches) > 1:
        output.pop("route", None)
        output["routes"] = [branch["route"] for branch in branches.values()]
    output["skipped_targets"] = skipped_targets
    output.update(
        pass_count=output["summary"]["pass"],
        review_count=output["summary"]["review"],
        fail_count=output["summary"]["fail"],
    )
    output["provenance"] = {**output["provenance"], "ended_at": utc_now()}
    output["configuration"] = asdict(options)
    cards = _screening_evidence(output, Path(output["provenance"]["output_dir"]).name)
    return _persist_screening(state, output, cards, list(dict.fromkeys(errors)), options)


def _small_molecule_inputs(design):
    """Use explicit Design values; never use Design's expected library/pass counts."""
    values = dict(design["target"]) if isinstance(design.get("target"), dict) else {}
    for step in design.get("design_steps", []):
        if step.get("step") in ("pocket_validation", "virtual_screening", "structure_prediction"):
            values.update(step.get("params", {}))
    values.update({k: v for k, v in design.items() if k not in ("target", "design_steps")})
    candidates = values.get("candidates", values.get("ligands", []))
    if values.get("candidates_path"):
        import csv

        with Path(values["candidates_path"]).open(newline="") as handle:
            candidates = list(csv.DictReader(handle))
    if not isinstance(candidates, list) or not all(isinstance(c, dict) for c in candidates):
        raise ValueError("Design candidates must be a list of ligand objects")
    context = {
        k: values.get(k)
        for k in (
            "receptor_path",
            "raw_receptor_path",
            "prepared_receptor_path",
            "target_sequence",
            "target_chain",
            "receptor_prepared",
            "receptor_preparation",
            "target_msa_path",
            "docking_box",
            "expected_pocket_residues",
            "hotspot_residues",
        )
    }
    context["docking_box"] = context["docking_box"] or values.get("binding_site")
    if design.get("workflow_schema") == "therapeutic-workflow-v1":
        context["design_workflow"] = {
            key: design.get(key)
            for key in (
                "workflow_schema",
                "entry_criteria",
                "safety_veto",
                "is_mock",
                "final_status",
                "evidence",
            )
        }
    return context, candidates


def _stage_call(tool, **kwargs):
    """Keep execution errors local to a candidate; tools own computation."""
    try:
        return tool.run(**kwargs)
    except Exception as exc:
        return skipped(str(exc), status="failed", is_mock=bool(kwargs.get("dry_run")))


def _screen_small_molecules(state, designs, handoff_errors, options):
    config = get_config().small_molecule
    if config.structural_backend == "boltz2_legacy":
        return _screen_small_molecules_legacy(state, designs, handoff_errors, options)
    from .small_molecule_workflow import run_small_molecule_workflow

    output = _empty_screening(state)
    output.update(
        modality="SMALL_MOLECULE",
        route="small_molecule",
        pipeline=[
            "RDKit",
            "preparation",
            "docking",
            "pose_qc",
            "AlphaFold3 ligand",
            "structural_consensus",
            "Validation",
            "ADMET",
            "deterministic_ranking",
            "Critic",
        ],
    )
    errors = list(handoff_errors)
    counters = {
        f"{stage}_{kind}": 0
        for stage in ("docking", "af3", "admet")
        for kind in ("tool_calls", "attempted", "completed")
    }
    summaries, validations, rankings = [], [], []
    directory = Path(output["provenance"]["output_dir"])
    remaining_docking, remaining_af3, remaining_admet = (
        options.max_gnina_candidates,
        options.max_af3_candidates,
        options.max_admet_candidates,
    )
    for ti, design in enumerate(designs):
        record = {
            "gene_name": design.get("gene_name", "unknown"),
            "modality": "SMALL_MOLECULE",
            "route": "small_molecule",
            "status": "NOT_EXECUTED",
            "failure_reasons": [],
        }
        output["targets"].append(record)
        try:
            if handoff_errors:
                raise ValueError("; ".join(handoff_errors))
            context, candidates = _small_molecule_inputs(design)
            if not candidates:
                raise ValueError("No small-molecule candidates supplied")
            context["gene_name"] = record["gene_name"]
            target_options = replace(
                options,
                max_gnina_candidates=remaining_docking,
                max_af3_candidates=remaining_af3,
                max_admet_candidates=remaining_admet,
            )
            result = run_small_molecule_workflow(
                candidates,
                context,
                config,
                directory / f"target_{ti:03d}",
                options=target_options,
                dry_run=state.get("dry_run", True),
                qualification=design,
                run_id=output["provenance"]["run_id"],
            )
            output["candidates"].extend(result["candidates"])
            summaries.append(result["candidate_summary_path"])
            validations.extend(result["validation_results"])
            rankings.extend(result["candidate_rankings"])
            for key, value in result["counts"].items():
                counters[key] += value
            remaining_docking -= result["counts"]["docking_tool_calls"]
            remaining_af3 -= result["counts"]["af3_tool_calls"]
            remaining_admet -= result["counts"]["admet_tool_calls"]
            record["status"] = "NOT_EXECUTED" if result["is_mock"] else "EXECUTED"
        except (OSError, ValueError, TypeError, KeyError) as exc:
            record.update(status="NOT_CONFIGURED", failure_reasons=[str(exc)])
            errors.append(f"[Screening:{record['gene_name']}] {exc}")
    real = [c for c in output["candidates"] if not c["is_mock"]]
    output.update(counters)
    output.update(
        summary={
            v: sum(c["verdict"] == v.upper() for c in real) for v in ("pass", "review", "fail")
        },
        mock_summary={
            v: sum(c["verdict"] == v.upper() for c in output["candidates"] if c["is_mock"])
            for v in ("pass", "review", "fail")
        },
        validation_results=validations,
        candidate_rankings=rankings,
        candidate_summary_paths=summaries,
        accepted_candidate_ids=[c["candidate_id"] for c in real if c["verdict"] == "PASS"],
        final_small_molecule_hits=sum(c["verdict"] == "PASS" for c in real),
        gnina_attempted=counters["docking_attempted"],
        gnina_completed=counters["docking_completed"],
        boltz2_attempted=0,
        boltz2_completed=0,
    )
    output["status"] = (
        "NOT_CONFIGURED"
        if errors and not output["candidates"]
        else "NOT_EXECUTED"
        if not real
        else "COMPLETED"
        if all(c["workflow_status"] == "complete" for c in real)
        else "PARTIAL_EXECUTION"
    )
    return output, errors


def _screen_small_molecules_legacy(state, designs, handoff_errors, options):
    from .protein_design_route import protein_route_options

    dry_run = protein_route_options(state)["dry_run"]
    output = _empty_screening({**state, "dry_run": dry_run})
    output.update(modality="SMALL_MOLECULE", pipeline=SMALL_MOLECULE_PIPELINE)
    directory = Path(output["provenance"]["output_dir"])
    gnina_options = dict(options.gnina_options)
    if options.resource_profile == "C":
        gnina_options["use_gpu"] = False
    gnina = GNINATool(gnina_options, gpu_device=options.gpu_device)
    boltz = Boltz2Tool(options.boltz2_options, gpu_device=options.gpu_device)
    admet_options = dict(options.admet_options)
    if options.resource_profile == "C":
        admet_options["use_gpu"] = False
    admet = ADMETTool(admet_options, gpu_device=options.gpu_device)
    gates = PoseScreeningTool(options.cascade_thresholds, options.thresholds)
    remaining = {stage: options.stage_limit(stage) for stage in ("gnina", "boltz2", "admet")}
    counts = {
        f"{stage}_{suffix}": 0
        for stage in remaining
        for suffix in ("tool_calls", "attempted", "completed")
    }
    errors, eligible = list(handoff_errors), []
    resamples_remaining = options.max_boltz2_resamples
    for ti, design in enumerate(designs):
        record = {
            "gene_name": design.get("gene_name", "unknown"),
            "modality": "SMALL_MOLECULE",
            "status": "NOT_EXECUTED",
            "failure_reasons": [],
        }
        output["targets"].append(record)
        try:
            context, ligands = _small_molecule_inputs(design)
            if handoff_errors or not ligands:
                record.update(
                    status="NOT_CONFIGURED",
                    failure_reasons=handoff_errors or ["MISSING_LIGAND_CANDIDATES"],
                )
                continue
            workflow = context.get("design_workflow")
            candidate_mock = bool(dry_run or (workflow and workflow.get("is_mock")))
            for ci, ligand in enumerate(ligands):
                candidate = {
                    "candidate_id": f"sm_target_{ti:03d}_candidate_{ci:04d}",
                    "input_candidate_id": ligand.get("candidate_id"),
                    "gene_name": record["gene_name"],
                    "modality": "SMALL_MOLECULE",
                    "smiles": ligand.get("smiles", ""),
                    "ligand_path": ligand.get("ligand_path"),
                    "is_mock": candidate_mock,
                    "verdict": "REVIEW",
                    "final_verdict": "REVIEW",
                    "status": "NOT_EXECUTED",
                    "failure_reasons": [],
                }
                output["candidates"].append(candidate)
                candidate_context = {**context, "smiles": candidate["smiles"]}
                if workflow:
                    candidate["design_workflow"] = workflow
                if workflow and (
                    workflow["safety_veto"] or not workflow["entry_criteria"]["eligible"]
                ):
                    dock = skipped("DESIGN_ENTRY_CRITERIA_NOT_MET", is_mock=candidate_mock)
                elif not remaining["gnina"]:
                    dock = skipped("GNINA_BUDGET_EXHAUSTED", is_mock=candidate_mock)
                    candidate["status"] = "SKIPPED_RESOURCE"
                else:
                    remaining["gnina"] -= 1
                    counts["gnina_tool_calls"] += 1
                    dock = _stage_call(
                        gnina,
                        candidate=candidate,
                        context=candidate_context,
                        output_dir=str(directory / candidate["candidate_id"] / "gnina"),
                        dry_run=candidate_mock,
                    )
                    counts["gnina_attempted"] += int(
                        bool(dock.get("provenance", {}).get("inference_started"))
                        and not dock.get("is_mock")
                    )
                candidate["gnina"] = dock
                candidate["is_mock"] = bool(candidate_mock or dock.get("is_mock"))
                counts["gnina_completed"] += int(
                    dock["status"] == "success" and not candidate["is_mock"]
                )
                try:
                    gate = gates.run(dock, candidate_context)
                except Exception as exc:
                    gate = verdict_result("REVIEW", ["POSE_SCREENING_ERROR: " + str(exc)])
                candidate["pose_screening"] = gate
                candidate.update(
                    verdict=gate["verdict"],
                    final_verdict=gate["verdict"],
                    failure_reasons=gate["failure_reasons"],
                )
                if candidate["status"] != "SKIPPED_RESOURCE":
                    candidate["status"] = {
                        "failed": "FAILED_GNINA",
                        "dependency_not_found": "NOT_CONFIGURED",
                        "invalid_input": "NOT_CONFIGURED",
                        "dry_run": "NOT_EXECUTED",
                    }.get(dock["status"], gate["verdict"])
                candidate["boltz2"] = skipped("GNINA_GATE_NOT_PASSED", is_mock=candidate["is_mock"])
                candidate["admet"] = skipped("UPSTREAM_NOT_PASSED", is_mock=candidate["is_mock"])
                candidate["cross_model_validation"] = verdict_result(
                    gate["verdict"], gate["failure_reasons"], status="NOT_EVALUATED"
                )
                if dock["status"] == "success" and gate["verdict"] == "PASS":
                    eligible.append((candidate, candidate_context))
            record["status"] = (
                "EXECUTED"
                if any(
                    c["gene_name"] == record["gene_name"] and c["gnina"]["status"] == "success"
                    for c in output["candidates"]
                )
                else "NOT_EXECUTED"
            )
        except Exception as exc:
            record.update(status="NOT_CONFIGURED", failure_reasons=[str(exc)])
            errors.append(f"[Screening:{record['gene_name']}] {exc}")
    eligible.sort(
        key=lambda item: item[0]["gnina"].get("scores", {}).get("CNNscore", -math.inf), reverse=True
    )
    for candidate, context in eligible:
        mock = candidate["is_mock"]
        if not remaining["boltz2"]:
            candidate["boltz2"] = skipped("BOLTZ2_BUDGET_EXHAUSTED", is_mock=mock)
            candidate.update(
                status="SKIPPED_RESOURCE",
                verdict="REVIEW",
                final_verdict="REVIEW",
                failure_reasons=["BOLTZ2_BUDGET_EXHAUSTED"],
            )
            continue
        remaining["boltz2"] -= 1
        counts["boltz2_tool_calls"] += 1
        result = _stage_call(
            boltz,
            target_sequence=context["target_sequence"],
            smiles=candidate["smiles"],
            msa_path=context.get("target_msa_path"),
            output_dir=str(directory / candidate["candidate_id"] / "boltz2"),
            dry_run=mock,
        )
        candidate["boltz2"] = result
        counts["boltz2_attempted"] += int(
            bool(result.get("provenance", {}).get("inference_started"))
            and not result.get("is_mock")
        )
        candidate["is_mock"] = bool(mock or result.get("is_mock"))
        counts["boltz2_completed"] += int(
            result["status"] == "success" and not candidate["is_mock"]
        )
        try:
            cross = gates.boltz_consistency(candidate["pose_screening"], result, context)
        except Exception as exc:
            cross = verdict_result("REVIEW", ["BOLTZ2_CONSISTENCY_ERROR: " + str(exc)])

        # Optional fallback uses the existing Boltz adapter, a distinct seed and
        # the SAME total call budget. Conflicting samples remain REVIEW.
        disagreement = cross["verdict"] == "REVIEW" and any(
            "DISAGREEMENT" in reason for reason in cross["failure_reasons"]
        )
        if disagreement and options.max_boltz2_resamples:
            candidate["consensus_sampling"] = {
                "status": "SKIPPED_RESOURCE",
                "initial_validation": cross,
                "samples": [],
            }
            if remaining["boltz2"] and resamples_remaining:
                remaining["boltz2"] -= 1
                resamples_remaining -= 1
                counts["boltz2_tool_calls"] += 1
                sampling_options = {**options.boltz2_options, "seed": boltz.config.seed + 1}
                additional = Boltz2Tool(sampling_options, gpu_device=options.gpu_device)
                sampled = _stage_call(
                    additional,
                    target_sequence=context["target_sequence"],
                    smiles=candidate["smiles"],
                    msa_path=context.get("target_msa_path"),
                    output_dir=str(directory / candidate["candidate_id"] / "boltz2_resample"),
                    dry_run=candidate["is_mock"],
                )
                candidate["is_mock"] = bool(candidate["is_mock"] or sampled.get("is_mock"))
                counts["boltz2_attempted"] += int(
                    bool(sampled.get("provenance", {}).get("inference_started"))
                    and not sampled.get("is_mock")
                )
                counts["boltz2_completed"] += int(
                    sampled.get("status") == "success" and not candidate["is_mock"]
                )
                try:
                    additional_gate = gates.boltz_consistency(
                        candidate["pose_screening"], sampled, context
                    )
                except Exception as exc:
                    additional_gate = verdict_result("REVIEW", ["RESAMPLING_ERROR: " + str(exc)])
                candidate["consensus_sampling"].update(
                    status=sampled["status"], samples=[sampled], validation=additional_gate
                )
                cross = {
                    **cross,
                    "verdict": "REVIEW",
                    "failure_reasons": list(
                        dict.fromkeys(
                            cross["failure_reasons"]
                            + additional_gate["failure_reasons"]
                            + ["ADDITIONAL_SAMPLING_REQUIRES_REVIEW"]
                        )
                    ),
                }
        candidate["cross_model_validation"] = cross
        candidate.update(
            verdict=cross["verdict"],
            final_verdict=cross["verdict"],
            failure_reasons=cross["failure_reasons"],
        )
        candidate["status"] = {
            "failed": "FAILED_BOLTZ2",
            "dependency_not_found": "NOT_CONFIGURED",
            "dry_run": "NOT_EXECUTED",
        }.get(result["status"], cross["verdict"])
        if result["status"] != "success" or cross["verdict"] == "FAIL":
            continue
        if not remaining["admet"]:
            candidate["admet"] = skipped("ADMET_BUDGET_EXHAUSTED", is_mock=candidate["is_mock"])
            candidate.update(
                status="SKIPPED_RESOURCE",
                verdict="REVIEW",
                final_verdict="REVIEW",
                failure_reasons=cross["failure_reasons"] + ["ADMET_BUDGET_EXHAUSTED"],
            )
            continue
        remaining["admet"] -= 1
        counts["admet_tool_calls"] += 1
        predicted = _stage_call(
            admet,
            candidate=candidate,
            output_dir=str(directory / candidate["candidate_id"] / "admet"),
            dry_run=candidate["is_mock"],
        )
        candidate["admet"] = predicted
        counts["admet_attempted"] += int(
            bool(predicted.get("provenance", {}).get("inference_started"))
            and not predicted.get("is_mock")
        )
        candidate["is_mock"] = bool(candidate["is_mock"] or predicted.get("is_mock"))
        counts["admet_completed"] += int(
            predicted["status"] == "success" and not candidate["is_mock"]
        )
        admet_gate = gates.admet(predicted)
        candidate["admet_screening"] = admet_gate
        final = (
            "FAIL"
            if admet_gate["verdict"] == "FAIL"
            else "PASS"
            if cross["verdict"] == admet_gate["verdict"] == "PASS"
            else "REVIEW"
        )
        candidate.update(
            verdict=final,
            final_verdict=final,
            failure_reasons=cross["failure_reasons"] + admet_gate["failure_reasons"],
        )
        candidate["status"] = {
            "failed": "FAILED_ADMET",
            "dependency_not_found": "NOT_CONFIGURED",
            "dry_run": "NOT_EXECUTED",
        }.get(predicted["status"], final)
    for candidate in output["candidates"]:
        workflow = candidate.get("design_workflow")
        if not workflow:
            continue
        selectivity = next((e for e in workflow["evidence"] if e["name"] == "selectivity"), {})
        if workflow["safety_veto"]:
            candidate.update(status="FAIL", verdict="FAIL", final_verdict="FAIL")
            candidate["failure_reasons"].append("SAFETY_VETO")
        elif selectivity.get("type") == "CONTRADICTORY":
            candidate.update(status="FAIL", verdict="FAIL", final_verdict="FAIL")
            candidate["failure_reasons"].append("SELECTIVITY_CONTRADICTORY")
        elif candidate["verdict"] == "PASS" and selectivity.get("type") != "SUPPORTIVE":
            candidate.update(status="REVIEW", verdict="REVIEW", final_verdict="REVIEW")
            candidate["failure_reasons"].append("SELECTIVITY_NOT_CONFIRMED")
    output.update(counts)
    candidates = output["candidates"]
    real = [c for c in candidates if not c["is_mock"]]
    output["summary"] = {
        v: sum(c["verdict"] == v.upper() for c in real) for v in ("pass", "review", "fail")
    }
    output["mock_summary"] = {
        v: sum(c["verdict"] == v.upper() for c in candidates if c["is_mock"])
        for v in ("pass", "review", "fail")
    }
    output["accepted_candidate_ids"] = [c["candidate_id"] for c in real if c["verdict"] == "PASS"]
    output["final_small_molecule_hits"] = len(output["accepted_candidate_ids"])
    if dry_run or (candidates and not real):
        output["status"] = "NOT_EXECUTED"
    elif counts["gnina_completed"]:
        output["status"] = (
            "COMPLETED"
            if all(c["status"] in ("PASS", "REVIEW", "FAIL") for c in real)
            else "PARTIAL_EXECUTION"
        )
    elif any(c["status"].startswith("FAILED_") for c in real):
        output["status"] = "FAILED_EXECUTION"
    elif any(t["status"] == "NOT_CONFIGURED" for t in output["targets"]) or any(
        c["status"] == "NOT_CONFIGURED" for c in real
    ):
        output["status"] = "NOT_CONFIGURED"
    elif candidates and all(c["status"] == "SKIPPED_RESOURCE" for c in candidates):
        output["status"] = "SKIPPED_RESOURCE"
    output["provenance"]["ended_at"] = utc_now()
    return output, errors


def _design_inputs(design: dict, state: dict) -> tuple[dict, list[str]]:
    """Translate data fields only; never execute Design's legacy tool names."""
    config = get_config()
    values = config.rfd3.model_dump()
    recognized = set(values)
    aliases = {
        "hotspots": "hotspot_residues",
        "hotspot": "hotspot_residues",
        "binder_length": "design_length",
    }
    conditions = []
    sources = [design.get("target", {})]
    sources += [
        step.get("params", {})
        for step in design.get("design_steps", [])
        if step.get("step") in ("rfdiffusion_backbone", "rfd3_generation")
    ]
    sources += [design, design.get("rfd3_input", {}), state.get("protein_design_input", {})]
    for source in sources:
        if not isinstance(source, dict):
            continue
        for key, value in source.items():
            key = aliases.get(key, key)
            if key not in recognized:
                continue
            if key in ("hotspot_residues", "motif_residues", "binding_site"):
                if isinstance(value, str):
                    tokens = [token.strip() for token in value.split(",") if token.strip()]
                    if not all(re.fullmatch(r"[A-Za-z]+[0-9]+", token) for token in tokens):
                        conditions.append(key)
                        continue
                    value = tokens
                if not isinstance(value, (list, tuple)) or not all(
                    isinstance(v, str) for v in value
                ):
                    conditions.append(key)
                    continue
                value = tuple(value)
                if value and key in conditions:
                    conditions = [c for c in conditions if c != key]
            values[key] = value
    return values, [key for key in conditions if not values.get(key)]


def _import_previous_candidates(
    previous: dict, target_id: str, binder_chain: str, gene: str, expected_length: int
) -> list[dict]:
    """Reuse RFD3 artifacts; old AF3 predictions cannot bypass sequence design."""
    from .protein_design_route import extract_protein_sequence

    candidates = []
    for index, old in enumerate(previous.get("candidates", [])):
        if old.get("gene_name") and old["gene_name"] != gene:
            continue
        source = previous.get("rfd3", {})
        mock = bool(
            previous.get("is_mock")
            or old.get("is_mock")
            or source.get("status") == "dry_run"
            or source.get("metadata", {}).get("is_mock")
        )
        path = old.get("rfd3_structure_path") or old.get("structure_path", "")
        chain = (
            old.get("binder_chain") or old.get("metadata", {}).get("binder_chain") or binder_chain
        )
        try:
            sequence = extract_protein_sequence(path, chain)
        except Exception:
            sequence = ""
        valid = source.get("status") in ("success", "dry_run") and len(sequence) == expected_length
        candidate = {
            "candidate_id": old.get("candidate_id") or f"{target_id}_candidate_{index:04d}",
            "rfd3_candidate_id": old.get("rfd3_candidate_id")
            or old.get("candidate_id")
            or f"{target_id}_candidate_{index:04d}",
            "rfd3_sequence": sequence,
            "rfd3_structure_path": path,
            "binder_chain": chain,
            "generated": source.get("status") == "success"
            and valid
            and not source.get("metadata", {}).get("is_mock", False),
            "is_mock": mock,
            "status": "RFD3_COMPLETED" if valid else "FAILED_RFD3",
            "failure_reasons": []
            if valid
            else ["RFD3 artifact or successful execution record missing"],
            "metadata": {"generated_by": "RFdiffusion3", "source": "protein_design_result"},
            "rfd3": old.get(
                "rfd3",
                {"structure_path": path, "status": "RFD3_COMPLETED" if valid else "FAILED_RFD3"},
            ),
        }
        if old.get("protein_mpnn"):
            candidate["protein_mpnn"] = old["protein_mpnn"]
            candidate["af3_validation"] = old.get("af3_validation", {"status": "not_run"})
        candidates.append(candidate)
    return candidates


def _mpnn_plan(design: dict):
    """Reuse Design's model/temperature/omit-AA fields and total sequence count."""
    from .configuration import ProteinMPNNConfig

    configured = get_config().protein_mpnn
    values = configured.model_dump()
    cap = None
    for step in design.get("design_steps", []):
        if step.get("step") != "proteinmpnn_sequence":
            continue
        params = step.get("params", {})
        for key, dest in (
            ("model_version", "model_name"),
            ("sampling_temp", "sampling_temp"),
            ("omit_AAs", "omit_AAs"),
        ):
            if key in params and dest not in configured.model_fields_set:
                values[dest] = params[key]
        if "num_seqs" in params:
            cap = params["num_seqs"]
            if type(cap) is not int or cap < 0:
                raise ValueError(
                    "Design proteinmpnn_sequence.num_seqs must be a nonnegative integer"
                )
    return ProteinMPNNConfig(**values), cap


def _design_sequences(batch, *, tool, output_dir, remaining, target_cap, dry_run):
    """Candidate-level orchestration only; ProteinMPNN owns sequence parsing/design."""
    results = []
    counts = {"mpnn_requested": 0, "mpnn_tool_calls": 0, "mpnn_attempted": 0}
    target_remaining = remaining if target_cap is None else min(target_cap, remaining)
    for backbone in batch:
        backbone["rfd3_candidate_id"] = (
            backbone.get("rfd3_candidate_id") or backbone["candidate_id"]
        )
        backbone["is_mock"] = bool(dry_run or backbone["is_mock"])
        backbone["generated"] = bool(backbone.get("generated") and not dry_run)
        cached = backbone.get("protein_mpnn", {})
        if backbone["status"] != "FAILED_RFD3" and tool.reusable(
            cached, len(backbone["rfd3_sequence"])
        ):
            backbone.update(
                sequence=cached["sequence"],
                mpnn_status="MPNN_COMPLETED",
                is_mock=bool(backbone["is_mock"] or cached.get("is_mock")),
                mpnn_reused=True,
            )
            results.append(backbone)
            continue

        # An AF3 result for the former RFD3 sequence is invalid after MPNN redesign.
        backbone.pop("af3_validation", None)
        allocation = min(tool.config.num_sequences_per_backbone, target_remaining)
        if backbone["status"] == "FAILED_RFD3":
            output = {"status": "not_run", "error_message": "RFD3_FAILED"}
        elif allocation == 0:
            output = {"status": "not_run", "error_message": "MPNN_BUDGET_EXHAUSTED"}
        else:
            remaining -= allocation
            target_remaining -= allocation
            counts["mpnn_requested"] += allocation
            counts["mpnn_tool_calls"] += 1
            counts["mpnn_attempted"] += int(not backbone["is_mock"])
            try:
                output = tool.run(
                    backbone=backbone,
                    output_dir=str(output_dir / backbone["candidate_id"] / "protein_mpnn"),
                    num_sequences=allocation,
                    dry_run=backbone["is_mock"],
                )
                if output.get("provenance", {}).get("execution_mode") == "REUSED_RESULT":
                    counts["mpnn_tool_calls"] -= 1
                    counts["mpnn_attempted"] -= int(not backbone["is_mock"])
            except Exception as exc:
                output = {"status": "failed", "error_message": str(exc)}
        sequences = (
            output.get("sequences", [])[:allocation]
            if allocation and output.get("status") in ("success", "dry_run")
            else []
        )
        if not sequences:
            phase = (
                "NOT_CONFIGURED"
                if output.get("status") in ("dependency_not_found", "invalid_input")
                else "SKIPPED_RESOURCE"
                if output.get("error_message") == "MPNN_BUDGET_EXHAUSTED"
                else "NOT_EXECUTED"
                if output.get("status") in ("not_run", "dry_run")
                else "FAILED_MPNN"
            )
            sequences = [
                {
                    "candidate_id": backbone["candidate_id"],
                    "sequence": "",
                    "score": None,
                    "sequence_path": "",
                    "status": phase,
                    "failure_reasons": [output.get("error_message") or "MPNN_RETURNED_NO_SEQUENCE"],
                }
            ]
        for sequence in sequences:
            sequence = {
                **sequence,
                "provenance": output.get("provenance", {}),
                "execution_status": output.get(
                    "execution_status", execution_status(output.get("status", "not_run"))
                ),
            }
            mock = bool(backbone["is_mock"] or output.get("is_mock") or sequence.get("is_mock"))
            sequence["is_mock"] = mock
            if sequence["status"] == "MPNN_COMPLETED" and not tool.reusable(
                sequence, len(backbone["rfd3_sequence"])
            ):
                sequence.update(
                    status="FAILED_MPNN",
                    failure_reasons=["MPNN sampled sequence artifact is missing or invalid"],
                )
            candidate = {
                **backbone,
                "candidate_id": sequence["candidate_id"],
                "sequence": sequence["sequence"],
                "protein_mpnn": sequence,
                "mpnn_status": sequence["status"],
                "is_mock": mock,
            }
            results.append(candidate)
    return results, remaining, counts


def _screen_binders(
    state: dict, designs: list[dict], handoff_errors: list[str], profiles: dict
) -> dict:
    from .configuration import RFD3Config
    from .protein_design_route import AlphaFold3BinderAdapter, protein_route_options
    from .rfdiffusion3_stage import RFdiffusion3Input, RFdiffusion3Runner

    config = get_config()
    progress = state.get('_progress_callback', lambda *args: None)
    route = protein_route_options(state)
    dry_run = route["dry_run"]
    started = utc_now()
    errors = list(handoff_errors)
    candidates, targets = [], []
    options_error = None
    try:
        options = ScreeningOptions.from_state(state)
    except (TypeError, ValueError) as exc:
        options = ScreeningOptions(max_rfd3_candidates=0, max_af3_candidates=0)
        options_error = f"Invalid Screening configuration: {exc}"
        errors.append(options_error)
    remaining_rfd3, remaining_mpnn, remaining_af3 = options.limits()
    protenix_counts = {"protenix_tool_calls": 0, "protenix_attempted": 0, "protenix_completed": 0}
    root = Path(
        state.get("protein_design_input", {}).get("output_dir")
        or config.rfd3.output_dir
        or (Path(config.data.agent_results) / "binder_screening")
    )
    root.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix="screening_", dir=root))
    run_id = state.get("run_id") or run_dir.name
    af3_tool = AF3Tool(
        AlphaFold3BinderAdapter(
            **config.af3_binder.model_dump(exclude={"seeds"}), gpu_device=options.gpu_device
        )
    )
    structural_tool = StructuralScreeningTool(options.thresholds)
    requested = requested_rfd3 = attempted_af3 = evaluated_af3 = af3_tool_calls = 0
    mpnn_counts = {"mpnn_requested": 0, "mpnn_tool_calls": 0, "mpnn_attempted": 0}
    first_rfd3 = {"status": "not_run"}
    previous = state.get("protein_design_result")
    for target_index, design in enumerate(designs):
        target_id = f"target_{target_index:03d}"
        gene = design.get("gene_name", "unknown")
        record = {
            "target_id": target_id,
            "gene_name": gene,
            "status": "NOT_EXECUTED",
            "failure_reasons": [],
        }
        targets.append(record)
        try:
            raw, unresolved = _design_inputs(design, state)
            try:
                settings = RFD3Config(**raw)
                mpnn_settings, mpnn_target_cap = _mpnn_plan(design)
            except (TypeError, ValueError) as exc:
                record.update(status="NOT_CONFIGURED", failure_reasons=[str(exc)])
                continue
            mpnn_tool = ProteinMPNNTool(mpnn_settings, gpu_device=options.gpu_device)
            record["requested_designs"] = settings.num_designs
            requested += settings.num_designs
            if not settings.backbone_paths:
                requested_rfd3 += settings.num_designs
            context = {
                "target_structure": settings.target_structure,
                "target_chain": settings.target_chain,
                "target_sequence": settings.target_sequence,
                "ligand": settings.ligand,
                "hotspot_residues": tuple(
                    dict.fromkeys(settings.hotspot_residues + settings.binding_site)
                ),
                "unresolved_design_conditions": unresolved,
            }
            record["design_conditions"] = context
            record["input_provenance"] = {
                "source": state.get("design_results_path")
                or str(Path(config.data.agent_results) / "design_output.json"),
                "rfd3_input": settings.model_dump(),
            }
            if options_error or handoff_errors:
                record.update(
                    status="NOT_CONFIGURED",
                    failure_reasons=[options_error] if options_error else handoff_errors,
                )
                continue
            if target_index >= options.max_targets:
                record.update(status="SKIPPED_RESOURCE", failure_reasons=["MAX_TARGETS_REACHED"])
                continue
            if settings.backbone_paths:
                from .binder_redesign import import_backbones
                limit = min(len(settings.backbone_paths), settings.num_designs, options.max_rfd3_candidates)
                if not limit:
                    record.update(status="SKIPPED_RESOURCE", failure_reasons=["BACKBONE_BUDGET_EXHAUSTED"])
                    continue
                progress("backbone_import", "RUNNING", gene)
                batch = import_backbones(settings, run_dir / target_id / "backbone_import", target_id, limit)
                progress("backbone_import", "COMPLETED", gene)
                record.update(status="EXECUTED", design_protocol="sequence_redesign",
                              rfd3={"status":"not_run", "reason":"supplied_backbone"})
            elif (
                previous
                and previous.get("rfd3", {}).get("status") != "not_run"
                and target_index == 0
            ):
                batch = _import_previous_candidates(
                    previous, target_id, settings.binder_chain, gene, settings.design_length
                )
                record["status"] = execution_status(
                    previous.get("rfd3", {}).get("status", "not_run")
                )
                record["reused_existing_output"] = True
                record["rfd3"] = previous.get("rfd3", {"status": "not_run"})
            else:
                allowed = min(settings.num_designs, remaining_rfd3)
                record["allocated_designs"] = allowed
                if not allowed:
                    record.update(
                        status="SKIPPED_RESOURCE", failure_reasons=["RFD3_BUDGET_EXHAUSTED"]
                    )
                    continue
                if not settings.target_sequence or set(settings.target_sequence) - set(
                    "ACDEFGHIKLMNPQRSTVWY"
                ):
                    record.update(
                        status="NOT_CONFIGURED",
                        failure_reasons=["Valid target_sequence is required from Design/config"],
                    )
                    continue
                remaining_rfd3 -= allowed
                request = RFdiffusion3Input(
                    **{
                        key: getattr(settings, key)
                        for key in RFdiffusion3Input.__dataclass_fields__
                        if key not in {"design_type", "output_dir", "num_designs"}
                    },
                    design_type="protein_binder",
                    output_dir=str(run_dir / target_id / "rfd3"),
                    num_designs=allowed,
                )
                progress("rfd3", "RUNNING", gene)
                output = RFD3Tool(
                    RFdiffusion3Runner(
                        settings.executable,
                        timeout_seconds=settings.timeout_seconds,
                        gpu_device=options.gpu_device,
                        environment=settings.environment,
                    )
                ).run(request, dry_run=dry_run)
                progress("rfd3", "COMPLETED" if output["rfd3"]["status"] == "success" else "FAILED", gene)
                batch = output["candidates"]
                for candidate in batch:
                    candidate["candidate_id"] = target_id + "_" + candidate["candidate_id"]
                record.update(
                    status=output["status"],
                    rfd3=output["rfd3"],
                    provenance=output.get("provenance", {}),
                )
                if output["rfd3"].get("error_message"):
                    record["failure_reasons"].append(output["rfd3"]["error_message"])
            if target_index == 0:
                first_rfd3 = record.get("rfd3", first_rfd3)
            progress("protein_mpnn", "RUNNING", gene)
            batch, remaining_mpnn, target_mpnn_counts = _design_sequences(
                batch,
                tool=mpnn_tool,
                output_dir=run_dir / target_id,
                remaining=remaining_mpnn,
                target_cap=mpnn_target_cap,
                dry_run=dry_run,
            )
            progress("protein_mpnn", "COMPLETED" if any(c.get("mpnn_status") == "MPNN_COMPLETED" for c in batch) else "FAILED", gene)
            for key, value in target_mpnn_counts.items():
                mpnn_counts[key] += value
            record.update(target_mpnn_counts)
            for candidate in batch:
                candidate.update(
                    gene_name=gene,
                    target_id=target_id,
                    modality=BINDER_MODALITY,
                    is_mock=bool(dry_run or candidate["is_mock"]),
                )
                candidates.append(candidate)
                inherited_af3 = candidate.get("af3_validation")
                was_failed = candidate["status"] == "FAILED_RFD3"
                mpnn_failed = candidate["mpnn_status"] == "FAILED_MPNN"
                can_fold = candidate["mpnn_status"] == "MPNN_COMPLETED" or (
                    candidate["is_mock"]
                    and bool(candidate["sequence"])
                    and candidate["mpnn_status"] == "NOT_EXECUTED"
                )
                if inherited_af3 is not None and not af3_tool.input_matches(
                    inherited_af3,
                    settings.target_sequence,
                    candidate["sequence"],
                    expected_payload=af3_tool.adapter.input_payload(
                        target_sequence=settings.target_sequence,
                        binder_sequence=candidate["sequence"],
                        seeds=config.af3_binder.seeds,
                    ),
                ):
                    inherited_af3 = None
                if was_failed or not can_fold:
                    reason = "; ".join(
                        candidate["failure_reasons"]
                        if was_failed
                        else candidate["protein_mpnn"].get("failure_reasons", [])
                    )
                    af3 = {"status": "not_run", "error_message": reason, "metrics": {}}
                elif inherited_af3 is not None:
                    af3 = inherited_af3
                    candidate["af3_status"] = "REUSED_OUTPUT"
                elif remaining_af3 == 0:
                    af3 = {
                        "status": "not_run",
                        "error_message": "AF3_BUDGET_EXHAUSTED",
                        "metrics": {},
                    }
                    candidate["af3_status"] = "SKIPPED_RESOURCE"
                else:
                    remaining_af3 -= 1
                    af3_tool_calls += 1
                    attempted_af3 += int(not dry_run and not candidate["is_mock"])
                    try:
                        progress("af3", "RUNNING", candidate["candidate_id"])
                        af3 = af3_tool.run(
                            target_sequence=settings.target_sequence,
                            binder_sequence=candidate["sequence"],
                            output_dir=str(run_dir / target_id / candidate["candidate_id"] / "af3"),
                            seeds=config.af3_binder.seeds,
                            dry_run=dry_run or candidate["is_mock"],
                        )
                        if af3.get("runtime_metadata", {}).get("execution_mode") == "REUSED_RESULT":
                            candidate["af3_status"] = "REUSED_OUTPUT"
                            af3_tool_calls -= 1
                            attempted_af3 -= int(not dry_run and not candidate["is_mock"])
                    except Exception as exc:
                        af3 = {"status": "failed", "error_message": str(exc), "metrics": {}}
                if af3.get("status") == "success":
                    af3_metrics = af3.get("metrics", {})
                    model_paths = af3_metrics.get("model_paths") or [af3_metrics.get("model_path")]
                    if not any(path and Path(path).is_file() for path in model_paths):
                        af3 = {
                            **af3,
                            "status": "failed",
                            "error_message": "AF3 reported success but model artifact is missing",
                        }
                candidate["is_mock"] = bool(
                    candidate["is_mock"]
                    or af3.get("status") == "dry_run"
                    or af3.get("metrics", {}).get("is_mock")
                )
                candidate["af3_validation"] = af3
                candidate["af3_output_path"] = af3.get("metrics", {}).get("model_path") or af3.get(
                    "output_dir", ""
                )
                if af3.get("status") == "success" and not candidate["is_mock"]:
                    evaluated_af3 += 1
                try:
                    screening_context = {**context}
                    if was_failed or mpnn_failed:
                        screening_context["execution_failure"] = (
                            af3.get("error_message") or candidate["mpnn_status"]
                        )
                    screened = structural_tool.run(
                        af3, screening_context, is_mock=candidate["is_mock"]
                    )
                except Exception as exc:
                    screened = {
                        "verdict": "REVIEW",
                        "failure_reasons": [f"STRUCTURAL_SCREENING_ERROR: {exc}"],
                        "af3_metrics": af3.get("metrics", {}),
                        "structural_metrics": {"status": "FAILED_EXECUTION"},
                        "hotspot_preserved": None,
                        "is_mock": candidate["is_mock"],
                    }
                candidate.update(screened)
                phase = candidate.get("af3_status")
                if phase != "SKIPPED_RESOURCE":
                    phase = {
                        "success": "AF3_COMPLETED",
                        "failed": "FAILED_AF3",
                        "dependency_not_found": "NOT_CONFIGURED",
                    }.get(af3.get("status"), "NOT_EXECUTED")
                progress("af3", "COMPLETED" if phase == "AF3_COMPLETED" else "NOT_RUN" if not can_fold else "FAILED", candidate["candidate_id"])
                candidate["af3_status"] = phase
                candidate["status"] = (
                    "FAILED_RFD3"
                    if was_failed
                    else "FAILED_MPNN"
                    if mpnn_failed
                    else candidate["mpnn_status"]
                    if not can_fold
                    else candidate["verdict"]
                    if phase == "AF3_COMPLETED"
                    else phase
                )
                candidate["af3"] = af3
                candidate["primary_structural_screening"] = {
                    **screened,
                    "interface_metrics": screened["structural_metrics"],
                }
                secondary = {"status": "not_used", "is_mock": False,
                             "cross_validation_backend": "protenix-v2",
                             "cross_validation_required": False}
                candidate["protenix_v2"] = secondary
                candidate.update(cross_validation_backend="protenix-v2",
                                 cross_validation_status="not_used",
                                 cross_validation_required=False)
                primary_call = screened["verdict"]
                if candidate["is_mock"] and primary_call == "PASS":
                    primary_call = "REVIEW"
                cross = verdict_result(primary_call, screened["failure_reasons"])
                cross.update(status="not_used", required=False,
                             basis="AF3 and interface QC only")
                candidate["cross_model_validation"] = cross
                candidate.update(
                    verdict=cross["verdict"],
                    final_verdict=cross["verdict"],
                    failure_reasons=cross["failure_reasons"],
                )
                if phase == "AF3_COMPLETED":
                    candidate["status"] = candidate["verdict"]

                # Compatibility with the existing Critic's deterministic protein gate.
                candidate["structure_path"] = candidate["rfd3_structure_path"]
                candidate["binder_sequence"] = candidate["sequence"]
                candidate["structural_validation"] = {
                    "status": candidate["verdict"].lower(),
                    "passed": candidate["verdict"] == "PASS" and not candidate["is_mock"],
                    "reason": "; ".join(candidate["failure_reasons"]),
                    "is_mock": candidate["is_mock"],
                    "metrics": candidate["structural_metrics"],
                }
                candidate["admet_status"] = "skipped_not_applicable"
        except Exception as exc:
            record.update(status="FAILED_EXECUTION", failure_reasons=[str(exc)])
            errors.append(f"[Screening:{gene}] {exc}")
    real = [c for c in candidates if not c["is_mock"]]
    mock = [c for c in candidates if c["is_mock"]]
    summary = {v.lower(): sum(c["verdict"] == v for c in real) for v in ("PASS", "REVIEW", "FAIL")}
    accepted = [c["candidate_id"] for c in real if c["verdict"] == "PASS"]
    statuses = {t["status"] for t in targets}
    if not candidates and "NOT_CONFIGURED" in statuses:
        status = "NOT_CONFIGURED"
    elif not candidates and "FAILED_EXECUTION" in statuses:
        status = "FAILED_EXECUTION"
    elif not candidates and "SKIPPED_RESOURCE" in statuses:
        status = "SKIPPED_RESOURCE"
    elif (
        any(c["status"] in ("FAILED_AF3", "FAILED_MPNN", "FAILED_RFD3") for c in candidates)
        and not evaluated_af3
    ):
        status = "FAILED_EXECUTION"
    elif dry_run or (mock and not real):
        status = "NOT_EXECUTED"
    elif evaluated_af3:
        status = (
            "COMPLETED"
            if all(
                c["af3_status"] == "AF3_COMPLETED"

                for c in real
            )
            and statuses <= {"EXECUTED"}
            else "PARTIAL_EXECUTION"
        )
    elif "FAILED_EXECUTION" in statuses or any(
        c["status"] in ("FAILED_AF3", "FAILED_MPNN", "FAILED_RFD3") for c in real
    ):
        status = "FAILED_EXECUTION"
    elif "NOT_CONFIGURED" in statuses or any(c["status"] == "NOT_CONFIGURED" for c in real):
        status = "NOT_CONFIGURED"
    elif "SKIPPED_RESOURCE" in statuses or any(c["status"] == "SKIPPED_RESOURCE" for c in real):
        status = "SKIPPED_RESOURCE"
    else:
        status = "NOT_EXECUTED"
    legacy_status = (
        "success"
        if accepted
        else "dry_run"
        if dry_run or mock
        else "failed"
        if status == "FAILED_EXECUTION"
        else "not_run"
    )
    first = candidates[0] if candidates else {}
    legacy = {
        "rfd3_candidate_id": first.get("candidate_id"),
        "af3_validation": first.get("af3_validation", {"status": "not_run"}),
        "structural_validation": first.get(
            "structural_validation", {"status": "not_run", "passed": None}
        ),
        "status": legacy_status,
        "design_mode": "protein_binder",
        "generated_by": "RFdiffusion3",
        "is_mock": bool(dry_run or mock),
        "rfd3": first_rfd3,
        "candidates": candidates,
        "candidate_count": len(candidates),
        "validated_candidate_count": len(accepted),
        "admet_status": "skipped_not_applicable",
        "pipeline_version": "agent_rev9-screening-tools-v3",
    }
    output = {
        "modality": BINDER_MODALITY,
        "design_mode": "protein_binder",
        "generated_by": "RFdiffusion3",
        "target": ", ".join(t["gene_name"] for t in targets),
        "status": status,
        "execution_mode": "mock" if dry_run or (mock and not real) else "real",
        "dry_run": dry_run,
        "is_mock": bool(dry_run or mock),
        "pipeline": PIPELINE,
        "requested_designs": requested,
        "generated_designs": len(
            {c["rfd3_structure_path"] for c in candidates if c.get("generated")}
        ),
        **mpnn_counts,
        **protenix_counts,
        "mpnn_completed": sum(
            c.get("mpnn_status") == "MPNN_COMPLETED"
            and not dry_run
            and not c.get("protein_mpnn", {}).get("is_mock", True)
            for c in candidates
        ),
        "af3_attempted": attempted_af3,
        "af3_tool_calls": af3_tool_calls,
        "af3_evaluated": evaluated_af3,
        "candidates": candidates,
        "summary": summary,
        "mock_summary": {
            v.lower(): sum(c["verdict"] == v for c in mock) for v in ("PASS", "REVIEW", "FAIL")
        },
        "targets": targets,
        "accepted_candidate_ids": accepted,
        "final_binder_hits": len(accepted),
        "admet_status": "skipped_not_applicable",
        "profiles": profiles,
        "binder_profiles": {},
        "profiles_execution_mode": "not_executed",
        "protein_design_result": legacy,
        "configuration": asdict(options),
        "provenance": {
            "run_id": run_id,
            "started_at": started,
            "ended_at": utc_now(),
            "output_dir": str(run_dir),
        },
    }
    if any(t.get("design_protocol") == "sequence_redesign" for t in targets):
        output["design_protocol"] = legacy["design_protocol"] = "sequence_redesign"
        output["generated_by"] = legacy["generated_by"] = "ProteinMPNN"
        output["pipeline"] = ["Provided backbone", *PIPELINE[1:]]
        output["provided_backbones"] = len({c["rfd3_structure_path"] for c in candidates if c.get("metadata", {}).get("generated_by") == "USER_BACKBONE"})
    if config.binder_analysis.enabled:
        try:
            from .binder_analysis import analyze_binders
            progress("binder_analysis", "RUNNING", ", ".join(t["gene_name"] for t in targets))
            output["binder_analysis"] = analyze_binders(output, run_dir / "analysis", config.binder_analysis)
            progress("binder_analysis", "COMPLETED", output["target"])
        except Exception as exc:
            output["binder_analysis"] = {"status":"FAILED", "warnings":[str(exc)]}
            progress("binder_analysis", "FAILED", output["target"])
    output.update(
        requested_rfd3=requested_rfd3,
        generated_rfd3=output["generated_designs"],
        pass_count=summary["pass"],
        review_count=summary["review"],
        fail_count=summary["fail"],
    )
    return output, errors, options


def _screening_evidence(output: dict, evidence_prefix: str) -> list[dict]:
    """Keep the existing EvidenceStore API; store per-tool evidence and failures."""
    cards, metric_cards = [], []
    summary_payload = {
        "execution_status": output["status"],
        "execution_mode": output["execution_mode"],
        "summary": output["summary"],
        "mock_summary": output["mock_summary"],
        "targets": output["targets"],
        "skipped_targets": output.get("skipped_targets", []),
    }
    cards.append(
        {
            "id": f"EV-{evidence_prefix}-SUMMARY",
            "source": "Screening",
            "type": "UNKNOWN",
            "content": json.dumps(_json_safe(summary_payload), ensure_ascii=False),
            "gene": "",
            "confidence": 0.0,
        }
    )
    for c in output["candidates"]:
        is_mock = bool(c["is_mock"])
        ev_type = (
            "UNKNOWN"
            if is_mock or c["verdict"] == "REVIEW"
            else "SUPPORTIVE"
            if c["verdict"] == "PASS"
            else "CONTRADICTORY"
        )
        if c["modality"] == BINDER_MODALITY:
            stages = [
                ("RFdiffusion3", c.get("rfd3", {})),
                ("ProteinMPNN", c.get("protein_mpnn", {})),
                ("AlphaFold3", c.get("af3_validation", {})),
                ("Protenix-v2", c.get("protenix_v2", {})),
            ]
        else:
            stages = (
                [
                    ("Docking", c.get("docking", {})),
                    ("AF3 ligand", c.get("af3_ligand", {})),
                    ("ADMET", c.get("admet", {})),
                ]
                if "af3_ligand" in c
                else [
                    ("GNINA", c.get("gnina", {})),
                    ("Boltz-2", c.get("boltz2", {})),
                    ("ADMET", c.get("admet", {})),
                ]
            )
            stages += [
                ("Boltz-2 additional sampling", sample)
                for sample in c.get("consensus_sampling", {}).get("samples", [])
            ]

        # Critic's existing text context truncates each content field. Start with
        # the verdict/tool statuses; full structured details follow in the card.
        compact = (
            c["verdict"]
            + " "
            + ",".join(f"{tool}:{result.get('status', 'not_run')}" for tool, result in stages)
        )
        detail = {
            "candidate": c["candidate_id"],
            "modality": c["modality"],
            "final_verdict": c["verdict"],
            "candidate_status": c["status"],
            "failure_reasons": c["failure_reasons"],
            "cross_model_validation": c.get("cross_model_validation"),
            "admet_screening": c.get("admet_screening"),
            "is_mock": is_mock,
            "run_id": output["provenance"]["run_id"],
        }
        if c["modality"] == "SMALL_MOLECULE":
            detail.update(
                design_workflow=c.get("design_workflow"),
                consensus_sampling=c.get("consensus_sampling"),
                validation=c.get("validation"),
                candidate_ranking=c.get("ranking"),
                candidate_critic=c.get("critic"),
            )
        cards.append(
            {
                "id": f"EV-{evidence_prefix}-{c['candidate_id']}-FINAL",
                "source": "Screening",
                "type": ev_type,
                "content": json.dumps(
                    _json_safe({"summary": compact, **detail}), ensure_ascii=False
                ),
                "gene": c["gene_name"],
                "confidence": 0.0 if ev_type == "UNKNOWN" else 1.0,
            }
        )
        gate_names = (
            ("AF3 structural gate", "primary_structural_screening"),
            ("Pose gate", "pose_screening"),
            ("Cross-model gate", "cross_model_validation"),
            ("ADMET gate", "admet_screening"),
        )
        for label, field in gate_names:
            gate = c.get(field)
            if gate and gate.get("status") == "not_used":
                continue
            if gate is None:
                continue
            gate_type = (
                "UNKNOWN"
                if is_mock or gate["verdict"] == "REVIEW"
                else "SUPPORTIVE"
                if gate["verdict"] == "PASS"
                else "CONTRADICTORY"
            )
            payload = {
                "verdict": gate["verdict"],
                "reasons": gate["failure_reasons"],
                "candidate": c["candidate_id"],
                "is_mock": is_mock,
                "details": gate,
            }
            cards.append(
                {
                    "id": f"EV-{evidence_prefix}-{c['candidate_id']}-{field}",
                    "source": label,
                    "type": gate_type,
                    "content": json.dumps(_json_safe(payload), ensure_ascii=False),
                    "gene": c["gene_name"],
                    "confidence": 0.0 if gate_type == "UNKNOWN" else 1.0,
                }
            )
        for tool, result in stages:
            metrics = result.get("metrics") or result.get("scores") or {}
            if tool == "ProteinMPNN":
                metrics = {"mpnn_status": c.get("mpnn_status"), "mpnn_score": result.get("score")}
            keys = [
                key
                for key in metrics
                if key not in ("samples", "model_paths", "plddt_summary")
                and not key.endswith("_path")
            ]
            if "plddt_summary" in metrics:
                keys.append("plddt_summary")
            for key, value in [("status", result.get("status", "not_run"))] + [
                (key, metrics[key]) for key in keys
            ]:
                payload = {
                    "candidate": c["candidate_id"],
                    "tool": tool,
                    "status": result.get("status", "not_run"),
                    "metric": key,
                    "value": value,
                    "final_verdict": c["verdict"],
                    "failure_reasons": c["failure_reasons"],
                    "error_message": result.get("error_message"),
                    "provenance": result.get("provenance", {}),
                    "output_path": result.get("sequence_path")
                    or result.get("pose_path")
                    or result.get("output_path")
                    or metrics.get("model_path")
                    or result.get("output_dir"),
                    "execution_mode": "mock" if is_mock or result.get("is_mock") else "real",
                    "is_mock": is_mock,
                    "candidate_status": c["status"],
                    "run_id": output["provenance"]["run_id"],
                }
                metric_type = (
                    "UNKNOWN"
                    if key == "status"
                    or value is None
                    or result.get("status") not in ("success", "MPNN_COMPLETED")
                    else ev_type
                )
                metric_cards.append(
                    {
                        "id": f"EV-{evidence_prefix}-{c['candidate_id']}-{tool}-{key}",
                        "source": tool,
                        "type": metric_type,
                        "content": json.dumps(_json_safe(payload), ensure_ascii=False),
                        "gene": c["gene_name"],
                        "confidence": 0.0 if metric_type == "UNKNOWN" else 1.0,
                    }
                )

    # Candidate decisions precede verbose raw metrics in the legacy Critic window.
    return cards + metric_cards


def _json_safe(value):
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _persist_screening(state, output, cards, errors, options):
    from .evidence_store import EvidenceStore

    config = get_config()
    directory = Path(config.data.agent_results)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "screening_output.json"
    output = _json_safe(output)
    db_path = str(
        Path(
            options.evidence_db_path or (Path(config.data.base_dir) / "evidence_store.db")
        ).resolve()
    )
    store = None
    try:
        store = EvidenceStore(db_path)
        store.add_many(cards, run_id=output["provenance"]["run_id"], agent="screening")
        output["evidence_store"] = {"status": "SAVED", "path": db_path}
    except Exception as exc:
        errors.append(f"Evidence persistence failed: {exc}")
        output["evidence_store"] = {"status": "FAILED_EXECUTION", "error": str(exc)}
    finally:
        if store is not None:
            store.close()
    text = json.dumps(output, ensure_ascii=False, indent=2, allow_nan=False)
    (Path(output["provenance"]["output_dir"]) / "screening_output.json").write_text(
        text, encoding="utf-8"
    )
    path.write_text(text, encoding="utf-8")
    update = {
        "current_agent": "screening",
        "screening_result": output,
        "screening_results_path": str(path),
        "evidence_cards": cards,
        "errors": errors,
        "messages": [
            f"[Screening] 스크리닝 대상 표적: {output.get('target', '')}",
            f"[Screening] {output['status']}: mode={output['execution_mode']}, PASS/REVIEW/FAIL={output['summary']}",
        ],
    }
    molecule_branch = output.get("branch_results", {}).get("small_molecule", output)
    update.update(
        validation_results=molecule_branch.get("validation_results", []),
        candidate_rankings=molecule_branch.get("candidate_rankings", []),
        candidate_summary_path=(molecule_branch.get("candidate_summary_paths") or [""])[0],
    )
    if output.get("protein_design_result"):
        update.update(
            use_rfd3=True,
            design_mode="protein_binder",
            protein_design_result=output["protein_design_result"],
        )
    return update


def main(argv=None) -> int:
    """Run only Screening. Default is dry-run; the demo cannot run real models."""
    import argparse

    from .configuration import OmiCraftConfig, set_config

    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--state", type=Path, help="JSON with the existing OmiCraftState fields")
    source.add_argument("--demo", action="store_true", help="Synthetic target; dry-run only")
    parser.add_argument(
        "--config", type=Path, help="OmiCraftConfig JSON; existing RFD3/AF3 settings"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--execute", action="store_true", help="Explicitly allow configured model execution"
    )
    args = parser.parse_args(argv)
    if args.demo and args.execute:
        parser.error("--demo always uses dry-run; --execute is not allowed")
    config = (
        OmiCraftConfig(
            **local_paths(json.loads(args.config.read_text()), args.config.resolve().parent)
        )
        if args.config
        else OmiCraftConfig()
    )
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    config.data.base_dir = str(directory)
    if not config.rfd3.output_dir:
        config.rfd3.output_dir = str(directory / "binder_screening")
    previous = get_config()
    try:
        set_config(config)
        if args.demo:
            from .rfdiffusion3_stage import _mock_cif

            target = directory / "mock_target.cif"
            target.write_text(_mock_cif(10, "T"))
            design = directory / "mock_design_output.json"
            design.write_text(
                json.dumps(
                    {
                        "designs": [
                            {
                                "gene_name": "MOCK_TARGET",
                                "modality": BINDER_MODALITY,
                                "target_structure": str(target),
                                "target_chain": "T",
                                "target_sequence": "A" * 10,
                                "hotspots": ["T3"],
                                "binder_length": 8,
                                "num_designs": 2,
                            }
                        ]
                    }
                )
            )
            state = {
                "design_results_path": str(design),
                "advance_targets": [{"gene_name": "MOCK_TARGET", "modality": BINDER_MODALITY}],
                "use_rfd3": False,
            }
        else:
            state = local_paths(json.loads(args.state.read_text()), args.state.resolve().parent)
        state["dry_run"] = not args.execute
        from .modality_dispatch import result as modality_result
        from .orchestrator import _run_screening_cli
        from .run_records import RunRecord

        record = RunRecord(directory, "DE_NOVO_BINDER")
        update = _run_screening_cli(state)
        raw = update.get("screening_result", {})
        status = {"COMPLETED": "COMPLETED", "NOT_EXECUTED": "NOT_RUN",
                  "PARTIAL_EXECUTION": "PARTIAL", "FAILED_EXECUTION": "FAILED"}.get(raw.get("status"), "BLOCKED")
        stored = modality_result("DE_NOVO_BINDER", status, "binder", summary=raw,
            artifacts=[{"path": update["dossier_path"], "label": "Screening dossier"}],
            is_mock=state["dry_run"])
        from .critic_agent import review_modality
        from .validation_agent import validate_modality
        stored["validation"] = validate_modality(stored)
        stored["validation_decision"] = stored["validation"]["decision"]
        stored["critic"] = review_modality(stored)
        record.finish(stored)
        result = update.get("screening_result", {})
        print(
            json.dumps(
                {
                    key: result.get(key)
                    for key in (
                        "status",
                        "execution_mode",
                        "requested_designs",
                        "generated_designs",
                        "mpnn_completed",
                        "af3_evaluated",
                        "summary",
                    )
                },
                indent=2,
            )
        )
        print("Result:", update.get("screening_results_path", ""))
        return (
            0
            if result.get("status") == "COMPLETED"
            or (
                not args.execute
                and result.get("status") == "NOT_EXECUTED"
                and result.get("candidates")
            )
            else 2
        )
    finally:
        set_config(previous)


if __name__ == "__main__":
    raise SystemExit(main())

"""Ligand preparation, docking, AF3 comparison and ADMET workflow."""

import hashlib
from dataclasses import asdict
from pathlib import Path

from .af3_ligand_tool import AF3LigandTool
from .candidate_ranking import rank_candidates
from .docking_backend import DockingRequest, docking_backend, run_docking
from .ligand_preparation import prepare_ligand, prepare_receptor, validate_candidate
from .pose_qc import PoseQCTool
from .routes import SMALL_MOLECULE, validate_route
from .small_molecule_admet import SmallMoleculeADMETTool, eligible_for_property_prediction
from .small_molecule_critic import audit_candidate
from .small_molecule_io import config_hash, now, unavailable, write_json
from .small_molecule_mock import SyntheticAF3Tool
from .small_molecule_reporting import workflow_status, write_artifacts
from .structural_consensus import StructuralConsensusTool
from .validation_agent import validate_small_molecule


def run_small_molecule_workflow(
    candidates,
    context,
    config,
    directory,
    *,
    options,
    dry_run=True,
    run_id=None,
    qualification=None,
    tools=None,
    progress=None,
):
    """dry_run plans native jobs; backend=mock explicitly executes synthetic fixtures."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    run_id = run_id or directory.name
    def stage(name, status, candidate_id=None):
        from .small_molecule_io import now
        entry = dict(stage=name, execution_status=status, candidate_id=candidate_id, timestamp=now())
        write_json(directory / "stage_status.json", entry)
        if progress:
            progress(name, status, candidate_id)

    def tool_status(result):
        return {"success": "COMPLETED", "partial": "PARTIAL", "not_run": "NOT_RUN"}.get(result.get("status"), "FAILED")
    started = now()
    tools = tools or {}
    seen = set()
    for c in candidates:
        if validate_route(c.get("route", SMALL_MOLECULE)) != SMALL_MOLECULE:
            raise ValueError("Binder candidate cannot enter small-molecule workflow")
        candidate_id = c.get("candidate_id")
        if candidate_id and candidate_id in seen:
            raise ValueError("Duplicate candidate_id: " + candidate_id)
        if candidate_id:
            seen.add(candidate_id)
    config = config.model_copy(deep=True)
    if options.resource_profile == "C":
        config.docking.use_gpu = False

    def effective_config():
        return {
            "small_molecule": config.model_dump(),
            "screening_options": asdict(options),
        }

    before = config_hash(effective_config())
    mock = config.docking.backend == "mock"
    backend = tools.get("docking") or docking_backend(config.docking, gpu_device=options.gpu_device)
    af3 = tools.get("af3") or (SyntheticAF3Tool if mock else AF3LigandTool)(
        config.af3, gpu_device=options.gpu_device
    )
    qc_tool = tools.get("pose_qc") or PoseQCTool(config.pose_qc)
    admet_options = dict(options.admet_options)
    if options.resource_profile == "C":
        admet_options["use_gpu"] = False
    admet = tools.get("admet") or SmallMoleculeADMETTool(
        admet_options, gpu_device=options.gpu_device
    )
    consensus_tool = StructuralConsensusTool(config.thresholds, config.pose_qc)
    context = dict(context)
    context["raw_receptor_path"] = context.get("raw_receptor_path") or context.get("receptor_path")
    stage("receptor_preparation", "RUNNING")
    receptor = prepare_receptor(context, directory / "receptor_preparation", is_mock=mock)
    stage("receptor_preparation", tool_status(receptor))
    context["prepared_receptor_path"] = receptor.get("prepared_receptor_path")
    source_path = write_json(directory / "input_candidates.json", candidates)
    remaining = {
        "docking": options.max_gnina_candidates,
        "af3": options.limits()[2],
        "admet": options.max_admet_candidates,
    }
    counts = {
        f"{stage}_{kind}": 0
        for stage in remaining
        for kind in ("tool_calls", "attempted", "completed")
    }
    output = []
    for index, source in enumerate(candidates):
        cid = source.get("candidate_id")
        safe = f"candidate_{index:04d}_" + hashlib.sha256(str(cid).encode()).hexdigest()[:8]
        work = directory / safe
        work.mkdir()
        stage("ligand_preparation", "RUNNING", cid)
        chemical = validate_candidate(source, seed=config.docking.seed)
        c = {
            **source,
            "route": SMALL_MOLECULE,
            "gene_name": context.get("gene_name", "unknown"),
            "modality": "SMALL_MOLECULE",
            "is_mock": bool(
                mock
                or dry_run
                or source.get("is_mock")
                or (qualification or {}).get("is_mock")
                or (context.get("design_workflow") or {}).get("is_mock")
            ),
            "admet_execution_policy": config.admet_policy,
            "canonical_smiles": chemical["canonical_smiles"],
            "rdkit_validation": chemical,
            "prepared_smiles": chemical["prepared_smiles"],
            "binding_site": context.get("docking_box"),
            "design_workflow": context.get("design_workflow"),
            "qualification": qualification or {},
            "pose_qc": [],
            "config_sha256": before,
            "input_candidate_id": cid,
            "warning": [],
        }
        write_json(work / "rdkit_validation.json", chemical)
        ligand = prepare_ligand(
            source, chemical, work / "ligand_preparation", seed=config.docking.seed
        )
        c["preparation"] = {"receptor": receptor, "ligand": ligand}
        stage("ligand_preparation", tool_status(ligand), cid)
        stage("docking", "RUNNING", cid)
        c["prepared_smiles"] = chemical["prepared_smiles"] = ligand.get("prepared_smiles")
        blocked = (
            not config.enabled
            or not chemical["docking_eligible"]
            or (qualification or {}).get("safety_verdict") == "REJECT"
        )
        upstream_mock = bool(
            source.get("is_mock")
            or (qualification or {}).get("is_mock")
            or (context.get("design_workflow") or {}).get("is_mock")
        )
        if upstream_mock and not mock and not dry_run:
            blocked = True
        workflow = context.get("design_workflow")
        if workflow and (
            workflow.get("safety_veto") or not workflow.get("entry_criteria", {}).get("eligible")
        ):
            blocked = True
        dock = unavailable(
            "invalid_input" if not chemical["smiles_valid"] else "not_run",
            poses=[],
            is_mock=mock,
            backend=backend.name,
        )
        if blocked:
            dock["error_message"] = (
                "MOCK_INPUT_REQUIRES_EXPLICIT_MOCK_BACKEND"
                if upstream_mock and not mock
                else "INPUT_OR_QUALIFICATION_NOT_ELIGIBLE"
            )
        elif receptor["status"] != "success" or ligand["status"] != "success":
            dock = unavailable(
                receptor["status"] if receptor["status"] != "success" else ligand["status"],
                "Preparation is incomplete",
                poses=[],
                backend=backend.name,
                is_mock=mock,
            )
        elif config.run_docking and remaining["docking"]:
            remaining["docking"] -= 1
            counts["docking_tool_calls"] += 1
            if dry_run and not mock:
                dock = unavailable(
                    "not_run",
                    "Dry-run: no docking inference",
                    poses=[],
                    backend=backend.name,
                )
            else:
                try:
                    box = context["docking_box"]
                    request = DockingRequest(
                        str(cid),
                        receptor["prepared_receptor_path"],
                        ligand["prepared_ligand_path"],
                        tuple(box["center"]),
                        tuple(box["size"]),
                        c["canonical_smiles"],
                        str(work / "docking"),
                        source.get("prepared_ligand_pdbqt_path"),
                        context.get("reference_ligand_path"),
                    )
                    dock = run_docking(backend, request)
                    counts["docking_attempted"] += int(not mock and bool(dock.get("command")))
                    counts["docking_completed"] += int(not mock and dock["status"] == "success")
                except (KeyError, TypeError, ValueError) as exc:
                    dock = unavailable(
                        "invalid_input",
                        str(exc),
                        poses=[],
                        backend=backend.name,
                        is_mock=mock,
                    )
        else:
            dock["error_message"] = "DOCKING_DISABLED_OR_BUDGET_EXHAUSTED"
        dock.setdefault("raw_pose_count", len(dock.get("poses", [])))
        c["docking"] = dock
        stage("docking", tool_status(dock), cid)
        stage("pose_qc", "RUNNING", cid)
        if config.run_pose_qc:
            c["pose_qc"] = [
                qc_tool.run(p, context, c["canonical_smiles"]) for p in dock.get("poses", [])
            ]
        good = any(q["pose_qc_status"] in {"pass", "pass_with_warning"} for q in c["pose_qc"])
        c["pose_qc_status"] = (
            "pass_with_warning"
            if good and any(q["warning"] for q in c["pose_qc"])
            else "pass"
            if good
            else "fail"
            if c["pose_qc"]
            else "not_run"
        )
        stage("pose_qc", "COMPLETED" if c["pose_qc"] else "NOT_RUN", cid)
        stage("af3", "RUNNING", cid)
        predicted = unavailable("not_run", "Pose QC has not passed", samples=[], is_mock=mock)
        if (
            config.run_af3
            and not blocked
            and chemical["af3_eligible"]
            and remaining["af3"]
            and (good or dry_run and not mock)
        ):
            remaining["af3"] -= 1
            counts["af3_tool_calls"] += 1
            try:
                predicted = af3.run(
                    candidate=c,
                    context=context,
                    output_dir=str(work / "af3"),
                    dry_run=dry_run and not mock,
                    is_mock=mock,
                )
            except (OSError, ValueError, RuntimeError, KeyError) as exc:
                predicted = unavailable("backend_failed", str(exc), samples=[], is_mock=mock)
            counts["af3_attempted"] += int(not mock and bool(predicted.get("command")))
            counts["af3_completed"] += int(not mock and predicted["status"] == "success")
        c["af3_ligand"] = predicted
        stage("af3", tool_status(predicted), cid)
        c["is_mock"] = c["is_mock"] or any(
            stage.get("is_mock", False) for stage in (dock, predicted)
        )
        stage("structural_comparison", "RUNNING", cid)
        if config.run_consensus:
            try:
                consensus = consensus_tool.run(c, dock, c["pose_qc"], predicted, context)
            except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
                consensus = {
                    "consensus_status": "INDETERMINATE",
                    "reason_codes": ["CONSENSUS_CALCULATION_FAILED"],
                    "warning": [str(exc)],
                    "comparisons": [],
                    "missing_metrics": ["consensus_geometry"],
                }
        else:
            consensus = {
                "consensus_status": "INDETERMINATE",
                "reason_codes": ["CONSENSUS_DISABLED"],
                "warning": [],
                "comparisons": [],
            }
        c["consensus"] = consensus
        stage("structural_comparison", "COMPLETED" if consensus.get("comparisons") else "NOT_RUN", cid)
        c["validation"] = validate_small_molecule(c, config.validation)
        stage("admet", "RUNNING", cid)
        property_eligible = not blocked and eligible_for_property_prediction(c)
        c["admet"] = unavailable("skipped_structural_screening", is_mock=mock)
        if (
            config.run_admet
            and property_eligible
            and remaining["admet"]
            and (mock or not dry_run)
        ):
            remaining["admet"] -= 1
            counts["admet_tool_calls"] += 1
            try:
                c["admet"] = admet.run(c, str(work / "admet"), run_id=run_id, is_mock=mock)
            except (OSError, ValueError, RuntimeError, KeyError) as exc:
                c["admet"] = unavailable("prediction_failed", str(exc), is_mock=mock)
            counts["admet_attempted"] += int(not mock and bool(c["admet"].get("command")))
            counts["admet_completed"] += int(
                not mock and c["admet"]["status"] in {"success", "partial"}
            )
        elif property_eligible:
            c["admet"] = unavailable("not_run", "ADMET_DISABLED_OR_BUDGET_EXHAUSTED", is_mock=mock)
        c["admet"]["execution_policy"] = config.admet_policy
        stage("admet", tool_status(c["admet"]), cid)
        output.append(c)
    stage("validation_and_critic", "RUNNING")
    rankings = rank_candidates(output)
    for c in output:
        c["config_unchanged"] = before == config_hash(effective_config())
        c["critic"] = audit_candidate(c, config_unchanged=c["config_unchanged"])
        c["workflow_status"] = workflow_status(c)
        c["verdict"] = c["final_verdict"] = (
            "FAIL"
            if c["critic"]["verdict"] == "NO_GO"
            else "PASS"
            if c["critic"]["verdict"] in {"GO", "GO_WITH_WARNING"} and not c["is_mock"]
            else "REVIEW"
        )
        c["status"] = c["verdict"]
        c["failure_reasons"] = list(
            dict.fromkeys(c["validation"]["reason_codes"] + c["critic"]["reason_codes"])
        )
        c["cross_model_validation"] = {
            "verdict": c["verdict"],
            "failure_reasons": c["failure_reasons"],
            "consensus_status": c["consensus"]["consensus_status"],
        }
    stage("validation_and_critic", "COMPLETED")
    stage("report", "RUNNING")
    paths = [
        str(source_path),
        context.get("raw_receptor_path"),
        context.get("target_msa_path") or config.af3.target_msa_path,
        context.get("prepared_receptor_path"),
    ] + [c.get("raw_ligand_path") or c.get("ligand_path") for c in candidates]
    paths.extend(
        c.get(key) for c in candidates for key in ("ccd_path", "prepared_ligand_pdbqt_path")
    )
    summary = write_artifacts(
        directory,
        output,
        run_id=run_id,
        config=effective_config(),
        config_sha256=before,
        input_paths=paths,
        started_at=started,
    )
    stage("report", "COMPLETED")
    return {
        "candidates": output,
        "candidate_rankings": rankings,
        "candidate_summary_path": str(summary),
        "validation_results": [c["validation"] for c in output],
        "counts": counts,
        "run_id": run_id,
        "started_at": started,
        "finished_at": now(),
        "is_mock": bool(mock or dry_run or any(c["is_mock"] for c in output)),
        "config_sha256": before,
    }

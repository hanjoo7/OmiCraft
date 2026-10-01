from collections import Counter
from pathlib import Path

from .small_molecule_io import now, redact, relative_path, sha256, write_csv, write_json
from .support.admet_pipeline.summary import BASE_SUMMARY_COLUMNS

SUMMARY_COLUMNS = """candidate_id display_name route smiles canonical_smiles prepared_smiles smiles_valid rdkit_parse_status molecular_weight formal_charge fragment_count docking_backend docking_backend_version docking_status docking_score cnn_score cnn_affinity best_docking_pose_path pose_qc_status pocket_occupancy severe_clash_count chirality_valid af3_status af3_successful_seed_count af3_ligand_plddt_mean af3_ligand_chain_iptm af3_protein_ligand_pae af3_cross_seed_pose_rmsd af3_cross_seed_contact_reproducibility docking_af3_pose_rmsd contact_residue_jaccard interaction_fingerprint_similarity consensus_status screening_call eligible_for_admet admet_prediction_status backend backend_version model_version is_mock available_endpoint_count missing_endpoint_count ranking_status rank critic_verdict workflow_status warning error_message""".split()


def workflow_status(c):
    if not c["rdkit_validation"]["smiles_valid"]:
        return "invalid_input"
    if c["docking"]["status"] in {"backend_failed", "failed"}:
        return "docking_failed"
    if c["pose_qc_status"] == "fail":
        return "structurally_unsupported"
    if c["af3_ligand"]["status"] in {"backend_failed", "failed"}:
        return "af3_failed"
    if c["admet"]["status"] in {"prediction_failed", "backend_failed", "failed"}:
        return "admet_failed"
    if c["admet"]["status"] == "partial":
        return "partial"
    if c["critic"]["verdict"] in {"GO", "GO_WITH_WARNING"} and c["admet"]["status"] == "success":
        return "complete"
    return "needs_review"


def candidate_summary(c):
    d, a, consensus, admet = c["docking"], c["af3_ligand"], c["consensus"], c["admet"]
    pose = next(
        (p for p in d.get("poses", []) if p["pose_rank"] == consensus.get("docking_pose_rank")),
        {},
    )
    qc = next((q for q in c["pose_qc"] if q["pose_rank"] == pose.get("pose_rank")), {})
    sample = next(
        (
            s
            for s in a.get("samples", [])
            if (s["seed"], s["sample_id"])
            == (consensus.get("af3_seed"), consensus.get("af3_sample_id"))
        ),
        {},
    )
    row = {
        **admet.get("existing_summary", {}),
        **{k: c.get(k) for k in ("candidate_id", "display_name", "route", "smiles", "is_mock")},
        **{
            k: c["rdkit_validation"].get(k)
            for k in (
                "canonical_smiles",
                "prepared_smiles",
                "smiles_valid",
                "rdkit_parse_status",
                "molecular_weight",
                "formal_charge",
                "fragment_count",
            )
        },
        "docking_backend": d.get("backend"),
        "docking_backend_version": d.get("backend_version"),
        "docking_status": d["status"],
        **{k: pose.get(k) for k in ("docking_score", "cnn_score", "cnn_affinity")},
        "best_docking_pose_path": pose.get("pose_path"),
        "pose_qc_status": c["pose_qc_status"],
        **{k: qc.get(k) for k in ("pocket_occupancy", "severe_clash_count", "chirality_valid")},
        "af3_status": a["status"],
        "af3_successful_seed_count": a.get("successful_seed_count"),
        "af3_ligand_plddt_mean": sample.get("ligand_atom_plddt_mean"),
        "af3_ligand_chain_iptm": sample.get("ligand_chain_iptm"),
        "af3_protein_ligand_pae": sample.get("protein_ligand_pae_min"),
        **{
            k: consensus.get(k)
            for k in (
                "af3_cross_seed_pose_rmsd",
                "af3_cross_seed_contact_reproducibility",
                "docking_af3_pose_rmsd",
                "contact_residue_jaccard",
                "interaction_fingerprint_similarity",
                "consensus_status",
            )
        },
        **{k: c["validation"][k] for k in ("screening_call", "eligible_for_admet")},
        "admet_prediction_status": admet["status"],
        **{
            k: admet.get(k)
            for k in (
                "backend",
                "backend_version",
                "model_version",
                "available_endpoint_count",
                "missing_endpoint_count",
            )
        },
        "ranking_status": c["ranking"]["ranking_status"],
        "rank": c["ranking"]["rank"],
        "critic_verdict": c["critic"]["verdict"],
        "workflow_status": workflow_status(c),
        "warning": list(dict.fromkeys(c["validation"]["warnings"] + c["critic"]["warnings"])),
        "error_message": "; ".join(
            str(stage.get("error_message", ""))
            for stage in (d, a, admet)
            if stage.get("error_message")
        ),
    }
    return row


def write_artifacts(
    directory,
    candidates,
    *,
    run_id,
    config,
    config_sha256,
    input_paths=(),
    started_at=None,
):
    directory = Path(directory)
    rows = [candidate_summary(c) for c in candidates]
    path = write_csv(
        directory / "candidate_summary.csv",
        rows,
        list(dict.fromkeys(BASE_SUMMARY_COLUMNS + SUMMARY_COLUMNS)),
    )
    write_json(directory / "candidates.json", redact(candidates))
    inputs = [
        {"path": relative_path(p, directory / "run_manifest.json"), "sha256": sha256(p)}
        for p in input_paths
        if p and Path(p).is_file()
    ]
    manifest = {
        "schema_version": "small-molecule-v1",
        "run_id": run_id,
        "created_at": now(),
        "route": "small_molecule",
        "config_sha256": config_sha256,
        "configuration": redact(config),
        "inputs": inputs,
        "is_mock": any(c["is_mock"] for c in candidates) if candidates else False,
        "started_at": started_at,
        "finished_at": now(),
        "candidate_count": len(candidates),
        "warning": [],
        "summary_path": relative_path(path, directory / "run_manifest.json"),
        "summary_sha256": sha256(path),
        "success_count": sum(r["workflow_status"] == "complete" for r in rows),
        "partial_count": sum(r["workflow_status"] == "partial" for r in rows),
        "failed_count": sum(
            r["workflow_status"]
            in {
                "invalid_input",
                "docking_failed",
                "af3_failed",
                "structurally_unsupported",
                "admet_failed",
                "failed",
            }
            for r in rows
        ),
        "skipped_count": sum(r["workflow_status"] == "needs_review" for r in rows),
    }
    if len({c["is_mock"] for c in candidates}) > 1:
        manifest["warning"].append("MOCK_REAL_ROWS_MIXED")
    write_json(directory / "run_manifest.json", manifest)
    for name, key in [
        ("docking", "docking"),
        ("af3", "af3_ligand"),
        ("consensus", "consensus"),
        ("admet", "admet"),
    ]:
        records = [
            {
                "candidate_id": c["candidate_id"],
                "route": c["route"],
                "is_mock": c["is_mock"],
                **redact(c[key]),
            }
            for c in candidates
        ]
        statuses = Counter(r.get("status", r.get("consensus_status", "not_run")) for r in records)
        write_json(
            directory / f"{name}_manifest.json",
            {
                **manifest,
                "stage": name,
                "success_count": sum(statuses[s] for s in ("success", "CONSENSUS_SUPPORTED")),
                "partial_count": sum(
                    statuses[s]
                    for s in (
                        "partial",
                        "POCKET_SUPPORTED_POSE_UNCERTAIN",
                        "AF3_SUPPORTED_ONLY",
                        "DOCKING_SUPPORTED_ONLY",
                        "DISCORDANT",
                        "INDETERMINATE",
                    )
                ),
                "failed_count": sum(
                    statuses[s]
                    for s in (
                        "failed",
                        "backend_failed",
                        "prediction_failed",
                        "invalid_input",
                        "STRUCTURALLY_UNSUPPORTED",
                    )
                ),
                "skipped_count": sum(
                    n
                    for s, n in statuses.items()
                    if s
                    in {
                        "not_run",
                        "dependency_missing",
                        "weights_missing",
                        "database_missing",
                        "skipped_structural_screening",
                        "skipped_route",
                    }
                ),
                "dependency_status": dict(statuses),
                "backend": sorted({str(r.get("backend")) for r in records if r.get("backend")}),
                "records": records,
            },
        )
    write_dossier(directory, candidates, run_id=run_id)
    return path


def write_dossier(directory, candidates, *, run_id):
    """Write the complete candidate audit and ranking results."""
    return write_json(
        Path(directory) / "final_dossier.json",
        {
            "schema_version": "small-molecule-dossier-v1",
            "run_id": run_id,
            "route": "small_molecule",
            "interpretation": "Structural evidence and ADMET predictions only. Rank is evidence ordering, not affinity, efficacy or clinical success.",
            "candidate_summary_path": "candidate_summary.csv",
            "candidates": [
                {
                    "candidate_id": c["candidate_id"],
                    "route": c["route"],
                    "is_mock": c["is_mock"],
                    "validation": c["validation"],
                    "ranking": c["ranking"],
                    "critic": c["critic"],
                    "workflow_status": c["workflow_status"],
                }
                for c in candidates
            ],
        },
    )


def ligand_dossier_node(state, *, scenario):
    from .configuration import get_config

    directory = Path(get_config().data.agent_results)
    path = directory / ("final_dossier.json" if scenario == 1 else "final_dossier_s2.json")
    screening = state.get("screening_result", {})
    write_json(
        path,
        redact(
            {
                "research_question": state.get("research_question"),
                "scenario": scenario,
                "advance_targets": state.get("advance_targets", []),
                "screening_result": screening,
                "candidate_rankings": screening.get("candidate_rankings", []),
                "candidate_summary_paths": screening.get("candidate_summary_paths", []),
                "critic_verdict": state.get("critic_verdict"),
                "critic_reason": state.get("critic_reason"),
                "dossier_text": "Docking and AF3 are independent structural predictions. Confidence and evidence ranks do not estimate binding affinity or clinical success.",
            }
        ),
    )
    return {
        "current_agent": "dossier",
        "messages": [f"[Dossier] Structural evidence dossier: {path}"],
    }

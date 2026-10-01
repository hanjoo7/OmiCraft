EVIDENCE_ORDER = {
    "CONSENSUS_SUPPORTED": 0,
    "POCKET_SUPPORTED_POSE_UNCERTAIN": 1,
    "AF3_SUPPORTED_ONLY": 2,
    "DOCKING_SUPPORTED_ONLY": 3,
    "DISCORDANT": 4,
    "INDETERMINATE": 5,
    "STRUCTURALLY_UNSUPPORTED": 6,
}


def rank_candidates(candidates):
    def warning_count(c):
        return len(set(c["validation"]["warnings"] + c["admet"].get("warning", [])))

    def key(c):
        valid = c["rdkit_validation"]["smiles_valid"]
        failed = c["admet"].get("status") in {
            "prediction_failed",
            "backend_failed",
            "failed",
        }
        complete = c["admet"].get("available_endpoint_count") or 0
        execution_failed = max(
            {"success": 0, "partial": 1}.get(c[stage].get("status"), 2)
            for stage in ("docking", "af3_ligand")
        )
        missing = c["admet"].get("missing_endpoint_count")
        return (
            not valid,
            execution_failed,
            failed,
            EVIDENCE_ORDER[c["consensus"]["consensus_status"]],
            {
                "pass": 0,
                "pass_with_warning": 1,
                "not_run": 2,
                "dependency_missing": 3,
                "fail": 4,
            }.get(c.get("pose_qc_status"), 2),
            {"consistent": 0, "not_assessed": 1, "inconsistent": 2}.get(
                c["consensus"].get("af3_consistency_status"), 1
            ),
            missing is None,
            missing if missing is not None else 0,
            -complete,
            warning_count(c),
            str(c["candidate_id"]),
        )

    rows = []
    for rank, c in enumerate(sorted(candidates, key=key), 1):
        admet_status = c["admet"]["status"]
        status = (
            "invalid_input"
            if not c["rdkit_validation"]["smiles_valid"]
            else "admet_failed"
            if admet_status in {"prediction_failed", "backend_failed", "failed"}
            else "mock_only"
            if c["is_mock"]
            else "partial_evidence"
            if admet_status != "success"
            else "ranked"
        )
        row = {
            "candidate_id": c["candidate_id"],
            "route": c["route"],
            "structural_evidence_tier": c["consensus"]["consensus_status"],
            "consensus_status": c["consensus"]["consensus_status"],
            "pose_qc_status": c.get("pose_qc_status"),
            "af3_consistency_status": c["consensus"].get("af3_consistency_status"),
            "admet_status": admet_status,
            "ranking_status": status,
            "unresolved_warning_count": warning_count(c),
            "rank": rank,
            "ranking_meaning": "evidence_order_only",
        }
        c["ranking"] = row
        rows.append(row)
    return rows

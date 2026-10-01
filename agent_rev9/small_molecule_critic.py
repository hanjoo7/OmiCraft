def audit_candidate(candidate, *, config_unchanged=True):
    missing, contradictions, warnings = [], [], []
    v = candidate["validation"]
    a = candidate["af3_ligand"]
    d = candidate["docking"]
    if not config_unchanged or candidate.get("config_unchanged") is False:
        contradictions.append("CONFIG_CHANGED_DURING_RUN")
    stages = [d, a, candidate["admet"]]
    flags = {
        stage.get("is_mock", False)
        for stage in stages
        if stage.get("status") in {"success", "partial"}
    }
    row_flags = {
        r.get("is_mock", candidate.get("is_mock", False))
        for r in a.get("samples", [])
        if r.get("af3_status") == "success"
    }
    flags.update(row_flags)
    if len(flags) > 1:
        contradictions.append("MOCK_REAL_MIXTURE")
    if candidate["route"] != "small_molecule" and candidate["admet"].get("status") not in {
        "skipped_route",
        "not_run",
    }:
        contradictions.append("ADMET_APPLIED_TO_BINDER")
    if a.get("observed_sample_count", 0) < a.get("expected_sample_count", 0):
        missing.append("AF3_SAMPLES_INCOMPLETE")
    if len(a.get("samples", [])) < a.get("observed_sample_count", 0):
        contradictions.append("AF3_SAMPLE_SELECTION")
    if len(d.get("poses", [])) < d.get("raw_pose_count", len(d.get("poses", []))):
        contradictions.append("DOCKING_POSES_HIDDEN")
    method = candidate["consensus"].get("metric_method", {})
    if any(
        row.get("docking_af3_pose_rmsd") is not None
        for row in candidate["consensus"].get("comparisons", [])
    ):
        if method.get("symmetry_aware") is not True:
            contradictions.append("NAIVE_ATOM_INDEX_RMSD")
        if method.get("receptor_aligned") is not True:
            contradictions.append("RECEPTOR_NOT_ALIGNED")
    prep = candidate.get("preparation", {})
    receptor = prep.get("receptor", {})
    ligand = prep.get("ligand", {})
    if receptor.get("raw_receptor_path") and receptor.get("raw_receptor_path") == receptor.get(
        "prepared_receptor_path"
    ):
        contradictions.append("RAW_PREPARED_RECEPTOR_MIXED")
    if ligand.get("stereochemistry_preservation") is False:
        contradictions.append("STEREOCHEMISTRY_CHANGED")
    if candidate.get("qualification", {}).get("safety_verdict") == "REJECT":
        contradictions.append("QUALIFICATION_SAFETY_VETO")
    for tag in (
        "af3_as_affinity",
        "combined_drug_score",
        "mixed_docking_score_scales",
        "binder_thresholds_applied",
        "missing_values_filled_zero",
    ):
        if candidate.get("interpretation", {}).get(tag):
            contradictions.append(tag.upper())
    missing += v.get("missing_evidence", [])
    if not candidate.get("binding_site"):
        missing.append("BINDING_SITE_MISSING")
    if candidate["is_mock"]:
        warnings.append("MOCK_EVIDENCE_ONLY")
    if candidate["admet"]["status"] in {
        "prediction_failed",
        "backend_failed",
        "failed",
    }:
        warnings.append("ADMET_FAILED")
        if candidate["ranking"]["ranking_status"] == "ranked":
            contradictions.append("FAILED_ADMET_RANKED_NORMAL")
    workflow = candidate.get("design_workflow") or {}
    if workflow.get("safety_veto"):
        contradictions.append("QUALIFICATION_SAFETY_VETO")
    if (
        not any(
            e.get("name") == "selectivity" and e.get("type") == "SUPPORTIVE"
            for e in workflow.get("evidence", [])
        )
        and workflow
    ):
        missing.append("SELECTIVITY_EVIDENCE_MISSING")
    if (
        candidate.get("qualification", {}).get("safety_verdict")
        not in {"PASS", "ADVANCE", "REJECT"}
        and not candidate["is_mock"]
    ):
        missing.append("QUALIFICATION_SAFETY_NOT_ASSESSED")
    warnings += v["warnings"] + candidate["admet"].get("warning", [])
    if contradictions or v["screening_call"] == "FAIL":
        verdict = "NO_GO"
    elif missing:
        verdict = "INSUFFICIENT_EVIDENCE"
    elif (
        candidate["is_mock"]
        or v["screening_call"] in {"REVISE", "INDETERMINATE"}
        or candidate["admet"]["status"] not in {"success", "partial"}
    ):
        verdict = "NEEDS_REVIEW"
    else:
        verdict = (
            "GO_WITH_WARNING" if warnings or candidate["admet"]["status"] == "partial" else "GO"
        )
    return {
        "candidate_id": candidate["candidate_id"],
        "verdict": verdict,
        "reason_codes": list(
            dict.fromkeys(
                contradictions
                + missing
                + (["ADMET_PARTIAL"] if candidate["admet"]["status"] == "partial" else [])
            )
        ),
        "contradictions": contradictions,
        "missing_evidence": missing,
        "warnings": list(dict.fromkeys(warnings)),
        "rerun_stage": "af3" if any("AF3" in reason for reason in missing) else None,
        "human_review_required": verdict in {"NO_GO", "NEEDS_REVIEW", "INSUFFICIENT_EVIDENCE"},
    }

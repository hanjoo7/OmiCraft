from itertools import combinations

import numpy as np

from .ligand_geometry import (
    aligned_ligand,
    jaccard,
    ligand_from_cif,
    pose_geometry,
    receptor_alignment,
    symmetry_rmsd,
)
from .ligand_preparation import read_ligand

CONSENSUS_STATUSES = (
    "CONSENSUS_SUPPORTED",
    "POCKET_SUPPORTED_POSE_UNCERTAIN",
    "AF3_SUPPORTED_ONLY",
    "DOCKING_SUPPORTED_ONLY",
    "DISCORDANT",
    "STRUCTURALLY_UNSUPPORTED",
    "INDETERMINATE",
)


def classify_consensus(metrics, thresholds):
    """Null thresholds calculate/report metrics but never create a hard gate."""
    reasons, missing, checks = [], [], {}
    rules = [
        ("minimum_successful_af3_seeds", "af3_successful_seed_count", "min"),
        ("maximum_cross_seed_pose_rmsd", "af3_cross_seed_pose_rmsd", "max"),
        ("maximum_docking_af3_pose_rmsd", "docking_af3_pose_rmsd", "max"),
        ("minimum_contact_residue_jaccard", "contact_residue_jaccard", "min"),
        ("minimum_pocket_reproducibility", "af3_pocket_reproducibility", "min"),
        ("minimum_ligand_plddt", "ligand_atom_plddt_mean", "min"),
        ("minimum_ligand_chain_iptm", "ligand_chain_iptm", "min"),
        ("maximum_protein_ligand_pae", "protein_ligand_pae_min", "max"),
    ]
    for setting, metric, direction in rules:
        cutoff = getattr(thresholds, setting)
        if cutoff is None:
            continue
        value = metrics.get(metric)
        if value is None:
            missing.append(metric)
            checks[setting] = None
        else:
            checks[setting] = value >= cutoff if direction == "min" else value <= cutoff
    failures = [name for name, value in checks.items() if value is False]
    docking, af3 = metrics.get("docking_supported"), metrics.get("af3_supported")
    if not docking and not af3:
        status = (
            "STRUCTURALLY_UNSUPPORTED"
            if metrics.get("both_structurally_failed")
            else "INDETERMINATE"
        )
    elif docking and not af3:
        status = "DOCKING_SUPPORTED_ONLY"
    elif af3 and not docking:
        status = "AF3_SUPPORTED_ONLY"
    elif missing or metrics.get("mapping_status") != "success":
        status = "INDETERMINATE"
    elif any(
        name in failures
        for name in ("maximum_docking_af3_pose_rmsd", "minimum_contact_residue_jaccard")
    ):
        status = "DISCORDANT"
    elif failures:
        status = "POCKET_SUPPORTED_POSE_UNCERTAIN"
    elif not checks:
        status = "INDETERMINATE"
        reasons.append("THRESHOLDS_NOT_CONFIGURED")
    elif (
        thresholds.maximum_docking_af3_pose_rmsd is None
        or thresholds.maximum_cross_seed_pose_rmsd is None
        or thresholds.minimum_successful_af3_seeds is None
    ):
        status = "POCKET_SUPPORTED_POSE_UNCERTAIN"
        reasons.append("POSE_OR_SEED_THRESHOLDS_NOT_CONFIGURED")
    else:
        status = "CONSENSUS_SUPPORTED"
    if thresholds.threshold_status == "provisional":
        reasons.append("PROVISIONAL_THRESHOLDS")
    if failures:
        reasons += ["THRESHOLD_NOT_MET:" + name for name in failures]
    if missing:
        reasons += ["METRIC_MISSING:" + name for name in missing]
    return {
        "consensus_status": status,
        "reason_codes": reasons,
        "threshold_checks": checks,
        "missing_metrics": missing,
    }


class StructuralConsensusTool:
    def __init__(self, thresholds, pose_config):
        self.thresholds, self.pose_config = thresholds, pose_config

    def run(self, candidate, docking, qcs, af3, context):
        samples = af3.get("samples", [])
        rows, aligned = [], {}
        valid_poses = {
            p["pose_rank"]: p
            for p, q in zip(docking.get("poses", []), qcs)
            if q["pose_qc_status"] in {"pass", "pass_with_warning"}
        }
        qc_by_rank = {q["pose_rank"]: q for q in qcs}
        for sample in samples:
            if not sample.get("model_path"):
                continue
            try:
                alignment = receptor_alignment(
                    context["raw_receptor_path"],
                    sample["model_path"],
                    context["target_chain"],
                    sample["protein_chain_id"],
                )
                raw_mol = ligand_from_cif(
                    sample["model_path"],
                    sample["ligand_chain_id"],
                    candidate["canonical_smiles"],
                )
                from .ligand_geometry import chirality_valid

                sample["chirality_valid"] = chirality_valid(raw_mol, candidate["canonical_smiles"])
                mol = aligned_ligand(raw_mol, alignment)
                rotation, center, target_center, mapping, ref, moving, rmsd = alignment
                transformed = [
                    {
                        **r,
                        "atoms": [
                            (np.asarray(x) - center) @ rotation.T + target_center
                            for x in r["atoms"]
                        ],
                    }
                    for r in moving
                ]
                labels = [
                    ref[mapping[i]]["label"] if i in mapping else None for i in range(len(moving))
                ]
                geometry = pose_geometry(
                    mol,
                    transformed,
                    context["docking_box"],
                    contact_distance=self.pose_config.contact_distance_angstrom,
                    clash_distance=self.pose_config.severe_clash_distance_angstrom,
                    labels=labels,
                )
                sample.update(geometry, receptor_alignment_rmsd=rmsd, geometry_status="success")
                if self.pose_config.run_prolif:
                    from .pose_qc import prolif_fingerprint

                    fingerprint = prolif_fingerprint(mol, context["prepared_receptor_path"])
                    sample["warning"].extend(fingerprint.pop("warning"))
                    sample.update(fingerprint)
                aligned[(sample["seed"], sample["sample_id"])] = mol
            except (OSError, ValueError, RuntimeError, KeyError) as exc:
                sample.update(geometry_status="mapping_failed", geometry_error=str(exc))
        usable = [
            s
            for s in samples
            if s.get("geometry_status") == "success" and s["af3_status"] == "success"
        ]
        cross_rmsd, cross_contact = [], []
        cross_rows = []
        for left, right in combinations(usable, 2):
            if left["seed"] == right["seed"]:
                continue
            record = {
                "left": [left["seed"], left["sample_id"]],
                "right": [right["seed"], right["sample_id"]],
                "rmsd_angstrom": None,
                "mapping_status": "success",
            }
            try:
                value = symmetry_rmsd(
                    aligned[tuple(record["left"])], aligned[tuple(record["right"])]
                )
                cross_rmsd.append(value)
                record["rmsd_angstrom"] = value
            except ValueError as exc:
                record.update(mapping_status="mapping_failed", error_message=str(exc))
            overlap = jaccard(left.get("contact_residues"), right.get("contact_residues"))
            if overlap is not None:
                cross_contact.append(overlap)
            record["contact_jaccard"] = overlap
            cross_rows.append(record)
        successful_seeds = {s["seed"] for s in usable}
        shared = {
            "af3_successful_seed_count": len(successful_seeds),
            "af3_cross_seed_pose_rmsd": max(cross_rmsd)
            if cross_rmsd and all(r["mapping_status"] == "success" for r in cross_rows)
            else None,
            "af3_cross_seed_contact_reproducibility": min(cross_contact) if cross_contact else None,
            "af3_pocket_reproducibility": sum(s.get("pocket_occupancy") == 1 for s in usable)
            / len(usable)
            if usable
            else None,
        }
        for pose in docking.get("poses", []) or [{"pose_rank": None}]:
            qc = qc_by_rank.get(pose["pose_rank"], {})
            for sample in samples or [{}]:
                row = {
                    "candidate_id": candidate["candidate_id"],
                    "docking_backend": docking.get("backend"),
                    "docking_pose_rank": pose["pose_rank"],
                    "af3_seed": sample.get("seed"),
                    "af3_sample_id": sample.get("sample_id"),
                    "docking_af3_pose_rmsd": None,
                    "contact_residue_jaccard": None,
                    "interaction_fingerprint_similarity": None,
                    "hotspot_contact_overlap": None,
                    "binding_site_center_distance": None,
                    "mapping_status": "not_run",
                    "warning": [],
                    "error_message": "",
                    **shared,
                    "docking_supported": pose["pose_rank"] in valid_poses,
                    "af3_supported": sample.get("geometry_status") == "success"
                    and sample.get("pocket_occupancy") == 1
                    and sample.get("has_clash") is False
                    and sample.get("chirality_valid") is True
                    and sample.get("severe_clash_count") == 0,
                    "ligand_atom_plddt_mean": _worst(samples, "ligand_atom_plddt_mean", min),
                    "ligand_chain_iptm": _worst(samples, "ligand_chain_iptm", min),
                    "protein_ligand_pae_min": _worst(samples, "protein_ligand_pae_min", max),
                }
                row["both_structurally_failed"] = (
                    qc.get("pose_qc_status") == "fail"
                    and sample.get("geometry_status") == "success"
                    and not row["af3_supported"]
                )
                key = sample.get("seed"), sample.get("sample_id")
                if pose.get("pose_path") and key in aligned:
                    try:
                        row["docking_af3_pose_rmsd"] = symmetry_rmsd(
                            read_ligand(pose["pose_path"]), aligned[key]
                        )
                        row["mapping_status"] = "success"
                        row["contact_residue_jaccard"] = jaccard(
                            qc.get("contact_residues"), sample.get("contact_residues")
                        )
                        row["binding_site_center_distance"] = sample.get(
                            "binding_site_center_distance"
                        )
                        row["interaction_fingerprint_similarity"] = jaccard(
                            qc.get("interaction_fingerprint"),
                            sample.get("interaction_fingerprint"),
                        )
                        hotspots = set(context.get("hotspot_residues") or [])
                        row["hotspot_contact_overlap"] = (
                            len(
                                hotspots
                                & set(qc.get("contact_residues") or [])
                                & set(sample.get("contact_residues") or [])
                            )
                            / len(hotspots)
                            if hotspots
                            else None
                        )
                    except (OSError, ValueError, RuntimeError) as exc:
                        row.update(mapping_status="mapping_failed", error_message=str(exc))
                if row["interaction_fingerprint_similarity"] is None:
                    row["warning"].append("INTERACTION_FINGERPRINT_UNAVAILABLE")
                row.update(classify_consensus(row, self.thresholds))
                rows.append(row)
        order = {name: i for i, name in enumerate(CONSENSUS_STATUSES)}
        best = min(
            rows,
            key=lambda r: (
                order[r["consensus_status"]],
                r["docking_af3_pose_rmsd"]
                if r["docking_af3_pose_rmsd"] is not None
                else float("inf"),
                r["docking_pose_rank"] or 0,
            ),
        )
        seed_consistency = "not_assessed"
        cutoff = self.thresholds.maximum_cross_seed_pose_rmsd
        if cutoff is not None and shared["af3_cross_seed_pose_rmsd"] is not None:
            seed_consistency = (
                "consistent" if shared["af3_cross_seed_pose_rmsd"] <= cutoff else "inconsistent"
            )
        return {
            **best,
            "comparisons": rows,
            "cross_seed_comparisons": cross_rows,
            "af3_consistency_status": seed_consistency,
            "metric_method": {
                "receptor_aligned": True,
                "symmetry_aware": True,
                "rmsd_units": "angstrom",
                "hydrogens": "excluded",
                "rmsd_method": "RDKit CalcRMS, graph automorphisms, no ligand refit",
                "cross_seed_summary": "maximum RMSD and minimum contact Jaccard across all different-seed sample pairs",
            },
            "thresholds": self.thresholds.model_dump(),
        }


def _worst(samples, field, aggregate):
    values = [s.get(field) for s in samples]
    return aggregate(values) if values and all(v is not None for v in values) else None

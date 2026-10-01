"""Binder tools and screening policies."""

from __future__ import annotations

import json
import math
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from .protein_design_route import AlphaFold3BinderAdapter, extract_protein_sequence
from .rfdiffusion3_stage import RFdiffusion3Input, RFdiffusion3Runner
from .screening_gates import CascadeThresholds, cif_atoms, contacts
from .structure_io import protein_residues

PIPELINE = [
    "RFdiffusion3",
    "ProteinMPNN",
    "AlphaFold3",
    "structural_screening",
    "Validation",
    "Critic",
]
SMALL_MOLECULE_PIPELINE = ["GNINA", "pose_screening", "Boltz-2", "cross_model_validation", "ADMET"]
BINDER_MODALITY = "DE_NOVO_BINDER"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def package_version(name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


@dataclass(frozen=True)
class StructuralThresholds:
    """Project heuristics, not biological success criteria. All gates live here."""

    iptm_fail_below: float = 0.4
    iptm_pass_min: float = 0.6
    ptm_fail_below: float = 0.3
    ptm_pass_min: float = 0.5
    pair_iptm_fail_below: float = 0.4
    pair_iptm_pass_min: float = 0.6
    pair_pae_pass_max: float = 10.0
    pair_pae_fail_above: float = 20.0
    binder_plddt_fail_below: float = 50.0
    binder_plddt_pass_min: float = 70.0
    contact_distance_angstrom: float = 5.0
    severe_clash_distance_angstrom: float = 1.5
    min_interface_contacts: int = 1
    require_hotspot_contact: bool = True

    def __post_init__(self):
        for name, value in asdict(self).items():
            if name == "require_hotspot_contact":
                if type(value) is not bool:
                    raise ValueError(f"{name} must be boolean")
            elif type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid threshold: {name}")
        for low, high in (
            (self.iptm_fail_below, self.iptm_pass_min),
            (self.ptm_fail_below, self.ptm_pass_min),
            (self.pair_iptm_fail_below, self.pair_iptm_pass_min),
        ):
            if not 0 <= low <= high <= 1:
                raise ValueError("TM thresholds must satisfy 0 <= fail <= pass <= 1")
        if not self.pair_pae_pass_max <= self.pair_pae_fail_above:
            raise ValueError("PAE pass threshold must not exceed fail threshold")
        if not 0 <= self.binder_plddt_fail_below <= self.binder_plddt_pass_min <= 100:
            raise ValueError("pLDDT thresholds must use the 0..100 scale")
        if not 0 < self.severe_clash_distance_angstrom < self.contact_distance_angstrom:
            raise ValueError("Clash distance must be positive and below contact distance")
        if type(self.min_interface_contacts) is not int or self.min_interface_contacts < 1:
            raise ValueError("min_interface_contacts must be a positive integer")


@dataclass(frozen=True)
class ScreeningOptions:
    max_rfd3_candidates: int = 8
    max_mpnn_candidates: int = 8
    max_af3_candidates: int = 8
    max_protenix_candidates: int = 8
    max_gnina_candidates: int = 100
    max_boltz2_candidates: int = 8
    max_boltz2_resamples: int = 0
    max_admet_candidates: int = 8
    max_targets: int = 1
    resource_profile: str | None = None
    gpu_device: int = 0
    evidence_db_path: str = ""
    thresholds: StructuralThresholds = field(default_factory=StructuralThresholds)
    cascade_thresholds: CascadeThresholds = field(default_factory=CascadeThresholds)
    protenix_options: dict = field(default_factory=dict)
    gnina_options: dict = field(default_factory=dict)
    boltz2_options: dict = field(default_factory=dict)
    admet_options: dict = field(default_factory=dict)

    def __post_init__(self):
        for name in (
            "max_rfd3_candidates",
            "max_mpnn_candidates",
            "max_af3_candidates",
            "max_protenix_candidates",
            "max_gnina_candidates",
            "max_boltz2_candidates",
            "max_admet_candidates",
            "max_boltz2_resamples",
            "max_targets",
            "gpu_device",
        ):
            if type(getattr(self, name)) is not int or getattr(self, name) < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.resource_profile not in (None, "A", "B", "C"):
            raise ValueError("resource_profile must be A, B, C or null")
        from .screening_model_tools import ToolSettings

        for name in ("protenix_options", "gnina_options", "boltz2_options", "admet_options"):
            if not isinstance(getattr(self, name), dict):
                raise ValueError(f"{name} must be an object")
            ToolSettings(**getattr(self, name))

    @classmethod
    def from_state(cls, state: dict) -> ScreeningOptions:

        values = {
            **state.get("budget", {}).get("screening", {}),
            **state.get("screening_options", {}),
        }

        ceiling = state.get("budget", {}).get("max_candidates")
        if ceiling is not None:
            if type(ceiling) is not int or ceiling < 0:
                raise ValueError("budget.max_candidates must be a nonnegative integer")
            for key in (
                "max_rfd3_candidates",
                "max_mpnn_candidates",
                "max_af3_candidates",
                "max_protenix_candidates",
                "max_gnina_candidates",
                "max_boltz2_candidates",
                "max_admet_candidates",
            ):
                values[key] = min(
                    values.get(key, 100 if key == "max_gnina_candidates" else 8), ceiling
                )
        thresholds = StructuralThresholds(**values.pop("thresholds", {}))
        cascade_thresholds = CascadeThresholds(**values.pop("cascade_thresholds", {}))
        return cls(**values, thresholds=thresholds, cascade_thresholds=cascade_thresholds)

    def limits(self) -> tuple[int, int, int]:
        caps = (self.max_rfd3_candidates, self.max_mpnn_candidates, self.max_af3_candidates)
        profile_cap = {None: max(caps), "A": 8, "B": 4, "C": 0}[self.resource_profile]
        return tuple(min(value, profile_cap) for value in caps)

    def stage_limit(self, stage: str) -> int:
        cap = getattr(self, f"max_{stage}_candidates")
        if stage == "protenix" and self.resource_profile in ("B", "C"):
            return 0
        if stage == "boltz2":
            return min(cap, {None: cap, "A": 8, "B": 4, "C": 0}[self.resource_profile])
        return cap


def execution_status(status: str) -> str:
    return {
        "success": "EXECUTED",
        "dry_run": "NOT_EXECUTED",
        "not_run": "NOT_EXECUTED",
        "dependency_not_found": "NOT_CONFIGURED",
        "invalid_input": "NOT_CONFIGURED",
        "failed": "FAILED_EXECUTION",
    }.get(status, "NOT_EXECUTED")


class RFD3Tool:
    """Reuse Foundry; validate backbones before ProteinMPNN sequence design."""

    def __init__(self, runner: RFdiffusion3Runner):
        self.runner = runner

    def run(self, request: RFdiffusion3Input, *, dry_run: bool) -> dict:
        started = utc_now()
        output = self.runner.run(request, dry_run=dry_run)
        is_mock = dry_run or output.status == "dry_run" or bool(output.metadata.get("is_mock"))
        candidates = []
        if output.status in ("success", "dry_run"):
            for index, path in enumerate(output.generated_structures):
                candidate = {
                    "candidate_id": f"candidate_{index:04d}",
                    "sequence": "",
                    "rfd3_structure_path": path,
                    "is_mock": is_mock,
                    "generated": Path(path).is_file() and not is_mock,
                    "status": "RFD3_COMPLETED",
                    "failure_reasons": [],
                    "target_chain": request.target_chain,
                    "binder_chain": request.binder_chain,
                    "metadata": {
                        "generated_by": "RFdiffusion3",
                        "binder_chain": request.binder_chain,
                    },
                }
                try:
                    sequence = extract_protein_sequence(path, request.binder_chain)
                    if len(sequence) != request.design_length:
                        raise ValueError(
                            "Binder sequence length differs from design_length; check binder_chain"
                        )
                    candidate["sequence"] = sequence
                except Exception as exc:
                    candidate.update(
                        status="FAILED_RFD3", generated=False, failure_reasons=[str(exc)]
                    )
                candidates.append(candidate)
        for candidate in candidates:
            candidate["rfd3_sequence"] = candidate.pop("sequence")
            candidate["rfd3"] = {
                "structure_path": candidate["rfd3_structure_path"],
                "status": candidate["status"],
                "target_chain": request.target_chain,
                "binder_chain": request.binder_chain,
                "metadata": output.metadata,
            }
        return {
            "status": execution_status(output.status),
            "rfd3": output.as_dict(),
            "candidates": candidates,
            "is_mock": is_mock,
            "provenance": {
                "started_at": started,
                "ended_at": utc_now(),
                "tool_version": package_version("rc-foundry"),
                "config": asdict(request),
                "seed": request.seed,
                "timeout_seconds": self.runner.timeout_seconds,
                "command": output.command_or_config.get("command"),
                "gpu_device": self.runner.gpu_device,
                "CUDA_VISIBLE_DEVICES": output.metadata.get("CUDA_VISIBLE_DEVICES"),
                "parent_CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
            },
        }


class AF3Tool:
    """Predict a target–binder complex from the two sequences."""

    def __init__(self, adapter: AlphaFold3BinderAdapter):
        self.adapter = adapter

    @staticmethod
    def input_matches(
        output: dict,
        target_sequence: str,
        binder_sequence: str,
        *,
        expected_payload: dict | None = None,
    ) -> bool:
        if output.get("status") != "success":
            return False
        try:
            payload = json.loads(Path(output["input_json"]).read_text())
            if expected_payload is not None and any(
                payload.get(key) != expected_payload.get(key)
                for key in ("sequences", "modelSeeds", "dialect", "version")
            ):
                return False
            proteins = {item["protein"]["id"]: item["protein"] for item in payload["sequences"]}
            return (
                proteins["A"]["sequence"] == target_sequence
                and proteins["B"]["sequence"] == binder_sequence
                and not any(protein.get("templates") for protein in proteins.values())
            )
        except (OSError, KeyError, TypeError, ValueError):
            return False

    def run(
        self,
        *,
        target_sequence: str,
        binder_sequence: str,
        output_dir: str,
        seeds: tuple[int, ...],
        dry_run: bool,
    ) -> dict:
        started = utc_now()
        output = self.adapter.run(
            target_sequence=target_sequence,
            binder_sequence=binder_sequence,
            output_dir=output_dir,
            seeds=seeds,
            dry_run=dry_run,
        )
        return {
            **output.as_dict(),
            "execution_status": execution_status(output.status),
            "provenance": {
                "started_at": started,
                "ended_at": utc_now(),
                "tool_version": package_version("alphafold3"),
                "seeds": list(seeds),
                "command": list(output.command),
                "input_json": output.input_json,
                "gpu_device": self.adapter.gpu_device,
                "timeout_seconds": self.adapter.timeout_seconds,
                "CUDA_VISIBLE_DEVICES": output.runtime_metadata.get("CUDA_VISIBLE_DEVICES"),
                "parent_CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
            },
        }


class StructuralScreeningTool:
    """Deterministic confidence, geometry and explicitly mapped hotspot gates."""

    def __init__(self, thresholds: StructuralThresholds | None = None):
        self.thresholds = thresholds or StructuralThresholds()

    def run(self, af3: dict, context: dict, *, is_mock: bool = False) -> dict:
        metrics = af3.get("metrics", {})
        failures, review = [], []
        geometry = {
            "status": "NOT_EXECUTED",
            "hotspot_status": "NOT_EXECUTED",
            "hotspot_preserved": None,
        }
        if context.get("execution_failure"):
            return self._result("FAIL", [context["execution_failure"]], metrics, geometry, is_mock)
        if af3.get("status") not in ("success", "dry_run"):
            reason = af3.get("error_message") or "AF3_NOT_EXECUTED"
            return self._result(
                "FAIL" if af3.get("status") == "failed" else "REVIEW",
                [reason],
                metrics,
                geometry,
                is_mock,
            )
        if af3.get("status") == "dry_run":
            return self._result("REVIEW", ["AF3_NOT_EXECUTED"], metrics, geometry, True)
        t = self.thresholds
        if metrics.get("has_clash") is True:
            failures.append("AF3_HAS_CLASH")
        elif metrics.get("has_clash") is not False:
            review.append("AF3_CLASH_STATUS_MISSING")
        self._score(
            metrics.get("iptm"), t.iptm_fail_below, t.iptm_pass_min, 1, "IPTM", failures, review
        )
        self._score(
            metrics.get("ptm"), t.ptm_fail_below, t.ptm_pass_min, 1, "PTM", failures, review
        )
        pair_iptm = self._pair_values(metrics, "chain_pair_iptm")
        pair_pae = self._pair_values(metrics, "chain_pair_pae_min")
        self._score(
            min(pair_iptm) if pair_iptm else None,
            t.pair_iptm_fail_below,
            t.pair_iptm_pass_min,
            1,
            "PAIR_IPTM",
            failures,
            review,
        )
        if not pair_pae or any(v < 0 for v in pair_pae):
            review.append("PAIR_PAE_MISSING_OR_INVALID")
        elif max(pair_pae) > t.pair_pae_fail_above:
            failures.append("PAIR_PAE_ABOVE_HARD_LIMIT")
        elif max(pair_pae) > t.pair_pae_pass_max:
            review.append("PAIR_PAE_BORDERLINE")
        plddt = (metrics.get("plddt_summary") or {}).get("by_chain_mean", {}).get("B")
        self._score(
            plddt,
            t.binder_plddt_fail_below,
            t.binder_plddt_pass_min,
            100,
            "BINDER_PLDDT",
            failures,
            review,
        )
        geometry = self._geometry(metrics.get("model_path"), context)
        if geometry["status"] != "COMPLETED":
            review.append(geometry["status"])
        else:
            if geometry["min_interchain_distance_angstrom"] < t.severe_clash_distance_angstrom:
                failures.append("SEVERE_INTERCHAIN_CLASH")
            if geometry["interface_contact_count"] < t.min_interface_contacts:
                failures.append("INSUFFICIENT_INTERFACE_CONTACTS")
        if t.require_hotspot_contact and context.get("hotspot_residues"):
            if geometry["hotspot_preserved"] is False:
                failures.append("HOTSPOT_CONTACT_LOST")
            elif geometry["hotspot_preserved"] is None:
                review.append(geometry["hotspot_status"])
        if context.get("unresolved_design_conditions"):
            review.append("DESIGN_CONDITIONS_UNRESOLVED")
        if context.get("ligand"):
            review.append("LIGAND_CONTEXT_NOT_VALIDATED_BY_PROTEIN_ONLY_AF3")
        if is_mock or metrics.get("is_mock"):
            review.append("MOCK_RESULT_NOT_BIOLOGICAL_EVIDENCE")
        verdict = "FAIL" if failures else "REVIEW" if review else "PASS"
        return self._result(verdict, failures + review, metrics, geometry, is_mock)

    @staticmethod
    def _score(value, hard_min, pass_min, maximum, label, failures, review):
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= maximum:
            review.append(label + "_MISSING_OR_INVALID")
        elif value < hard_min:
            failures.append(label + "_BELOW_HARD_LIMIT")
        elif value < pass_min:
            review.append(label + "_BORDERLINE")

    @staticmethod
    def _pair_values(metrics, key):
        chains = metrics.get("chain_ids") or []
        matrix = metrics.get(key)
        if "A" not in chains or "B" not in chains or not isinstance(matrix, list):
            return []
        a, b = chains.index("A"), chains.index("B")
        try:
            values = [matrix[a][b], matrix[b][a]]
            return (
                values if all(type(v) in (int, float) and math.isfinite(v) for v in values) else []
            )
        except (IndexError, TypeError):
            return []

    @staticmethod
    def _result(verdict, reasons, af3_metrics, geometry, is_mock):
        return {
            "verdict": verdict,
            "failure_reasons": reasons,
            "af3_metrics": af3_metrics,
            "structural_metrics": geometry,
            "hotspot_preserved": geometry["hotspot_preserved"],
            "is_mock": is_mock,
        }

    def _geometry(self, model_path: str | None, context: dict) -> dict:
        result = {
            "status": "MODEL_MISSING",
            "hotspot_status": "NOT_REQUESTED",
            "hotspot_preserved": None,
        }
        if not model_path or not Path(model_path).is_file():
            return result
        try:
            atoms = cif_atoms(model_path)
            target = [a for a in atoms if a.label_chain_id == "A" and not a.is_hydrogen]
            binder = [a for a in atoms if a.label_chain_id == "B" and not a.is_hydrogen]
            if not target or not binder:
                return {
                    **result,
                    "status": "AF3_CHAINS_MISSING",
                    "hotspot_status": "AF3_CHAINS_MISSING",
                }
            geometry = contacts(
                [a.coord for a in target],
                [a.coord for a in binder],
                [a.label_seq_id for a in target],
                self.thresholds.contact_distance_angstrom,
            )
            by_residue = geometry.pop("residue_min_distance_angstrom")
            result.update(
                status="COMPLETED",
                **geometry,
                interface_target_residue_count=len(geometry["interface_target_residue_ids"]),
            )
            hotspots = context.get("hotspot_residues", ())
            if hotspots:
                try:
                    residues = [
                        (n, aa)
                        for c, n, aa in protein_residues(context["target_structure"])
                        if c == context["target_chain"]
                    ]
                except Exception as exc:
                    result.update(
                        hotspot_status="HOTSPOT_MAPPING_UNAVAILABLE", hotspot_error=str(exc)
                    )
                    return result
                if "".join(aa for _, aa in residues) != context.get("target_sequence"):
                    result["hotspot_status"] = "HOTSPOT_MAPPING_UNAVAILABLE_SEQUENCE_MISMATCH"
                    return result
                mapping = {context["target_chain"] + n: i + 1 for i, (n, _) in enumerate(residues)}
                if any(h not in mapping or mapping[h] not in by_residue for h in hotspots):
                    result["hotspot_status"] = "HOTSPOT_MAPPING_UNAVAILABLE_RESIDUE_MISSING"
                    return result
                values = {h: by_residue[mapping[h]] for h in hotspots}
                result.update(
                    hotspot_status="MAPPED_BY_EXACT_TARGET_SEQUENCE",
                    hotspot_distances_angstrom=values,
                    hotspot_preserved=all(
                        d <= self.thresholds.contact_distance_angstrom for d in values.values()
                    ),
                )
        except Exception as exc:
            result.update(
                status="GEOMETRY_UNAVAILABLE",
                hotspot_status="HOTSPOT_MAPPING_UNAVAILABLE",
                error=str(exc),
            )
        return result

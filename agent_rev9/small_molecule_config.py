from __future__ import annotations
"""Configuration and screening thresholds for small molecules."""

from typing import Literal

from pydantic import AliasChoices, BaseModel, ConfigDict, Field, model_validator


class StrictConfig(BaseModel):
    model_config = ConfigDict(extra="forbid")


class ConsensusThresholds(StrictConfig):
    threshold_status: Literal["provisional", "validated"] = "provisional"
    minimum_successful_af3_seeds: int | None = Field(default=None, ge=1)
    maximum_cross_seed_pose_rmsd: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    maximum_docking_af3_pose_rmsd: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    minimum_contact_residue_jaccard: float | None = Field(default=None, ge=0, le=1)
    minimum_pocket_reproducibility: float | None = Field(default=None, ge=0, le=1)
    minimum_ligand_plddt: float | None = Field(default=None, ge=0, le=100)
    minimum_ligand_chain_iptm: float | None = Field(default=None, ge=0, le=1)
    maximum_protein_ligand_pae: float | None = Field(default=None, ge=0, allow_inf_nan=False)


class PoseQCConfig(StrictConfig):
    contact_distance_angstrom: float = Field(default=5.0, gt=0, allow_inf_nan=False)
    severe_clash_distance_angstrom: float = Field(default=1.5, gt=0, allow_inf_nan=False)
    run_posebusters: bool = True
    run_prolif: bool = True


class ValidationPolicy(StrictConfig):
    discordant_call: Literal["INDETERMINATE", "REVISE", "FAIL"] = "INDETERMINATE"
    eligible_consensus_statuses: list[str] = Field(
        default_factory=lambda: [
            "CONSENSUS_SUPPORTED",
            "POCKET_SUPPORTED_POSE_UNCERTAIN",
        ]
    )
    eligible_screening_calls: list[str] = Field(
        default_factory=lambda: ["PASS", "PASS_WITH_WARNING"]
    )
    require_all_af3_seeds: bool = True

    @model_validator(mode="after")
    def check_statuses(self):
        valid = {
            "CONSENSUS_SUPPORTED",
            "POCKET_SUPPORTED_POSE_UNCERTAIN",
            "AF3_SUPPORTED_ONLY",
            "DOCKING_SUPPORTED_ONLY",
            "DISCORDANT",
            "STRUCTURALLY_UNSUPPORTED",
            "INDETERMINATE",
        }
        if set(self.eligible_consensus_statuses) - valid:
            raise ValueError("Unknown ADMET consensus eligibility status")
        if set(self.eligible_screening_calls) - {
            "PASS",
            "PASS_WITH_WARNING",
            "REVISE",
            "FAIL",
            "INDETERMINATE",
        }:
            raise ValueError("Unknown screening call")
        return self


class DockingConfig(StrictConfig):
    backend: Literal["vina", "gnina", "mock"] = "gnina"
    executable: str = ""
    version: str | None = None
    seed: int = Field(default=11, ge=0)
    exhaustiveness: int = Field(default=8, ge=1)
    num_poses: int = Field(default=9, ge=1, validation_alias=AliasChoices("num_poses", "num_modes"))
    cnn_scoring: Literal["rescore", "refinement", "all", "none"] = "rescore"
    custom_cnn_model: str | None = None
    environment: dict[str, str] = Field(default_factory=dict)
    device: int | None = Field(default=None, ge=0)
    cpu: int = Field(default=4, ge=1)
    use_gpu: bool = True
    timeout_seconds: float = Field(default=3600, gt=0, allow_inf_nan=False)


class LigandAF3Config(StrictConfig):
    python_path: str = ""
    script_path: str = ""
    model_dir: str = ""
    database_dir: str = ""
    hmmer_bin_dir: str = ""
    target_msa_path: str = ""
    target_msa_cache_dir: str = ""
    seeds: tuple[int, ...] = (11, 22, 33)
    num_samples: int = Field(default=5, ge=1)
    input_version: Literal[4] = 4
    protein_chain_id: str = "A"
    ligand_chain_id: str = "B"
    use_templates: bool = False
    timeout_seconds: float = Field(default=21600, gt=0, allow_inf_nan=False)
    data_pipeline_timeout_seconds: float = Field(default=5400, gt=0, allow_inf_nan=False)
    inference_timeout_seconds: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    environment: dict[str, str] = Field(default_factory=dict)

    @model_validator(mode="after")
    def valid_seeds_chains(self):
        if (
            not self.seeds
            or len(set(self.seeds)) != len(self.seeds)
            or any(type(x) is not int or x < 0 for x in self.seeds)
        ):
            raise ValueError("AF3 seeds must be unique nonnegative integers")
        if self.protein_chain_id == self.ligand_chain_id or not all(
            c.isalpha() and c.isupper() for c in (self.protein_chain_id, self.ligand_chain_id)
        ):
            raise ValueError(
                "AF3 protein and ligand chain IDs must be distinct uppercase identifiers"
            )
        return self


class SmallMoleculeConfig(StrictConfig):
    enabled: bool = True
    structural_backend: Literal["af3", "boltz2_legacy"] = "af3"
    run_docking: bool = True
    run_pose_qc: bool = True
    run_af3: bool = True
    run_consensus: bool = True
    run_admet: bool = True
    admet_policy: Literal["structural_pass", "valid_input"] = "structural_pass"
    docking: DockingConfig = Field(default_factory=DockingConfig)
    af3: LigandAF3Config = Field(default_factory=LigandAF3Config)
    pose_qc: PoseQCConfig = Field(default_factory=PoseQCConfig)
    thresholds: ConsensusThresholds = Field(default_factory=ConsensusThresholds)
    validation: ValidationPolicy = Field(default_factory=ValidationPolicy)

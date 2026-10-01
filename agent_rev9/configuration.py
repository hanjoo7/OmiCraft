from __future__ import annotations
"""
OmiCraft rev_1.0 — Pydantic 설정 중앙화
Robin 패턴 적용: 모든 파라미터를 Pydantic BaseModel로 관리한다.
"""

import os
import sys
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from .small_molecule_config import SmallMoleculeConfig
from .runtime_paths import load_json


class ConfigModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class DataPaths(ConfigModel):
    """데이터 경로 설정"""

    base_dir: str = Field(default_factory=lambda: str(Path.cwd() / "runs"))

    counts_rds: str | None = None
    metadata_rds: str | None = None
    annotation_rds: str | None = None
    clinical_rds: str | None = None
    gene_sets_rds: str | None = None
    sample_manifest: dict | str | None = None
    chronos_path: str | None = None
    demeter2_path: str | None = None
    model_meta_path: str | None = None
    hpa_path: str | None = None
    expr_rds: str | None = None
    meta_rds: str | None = None

    @property
    def agent_results(self) -> str:
        return os.path.join(self.base_dir, "agent_results")


class RFD3Config(ConfigModel):
    enabled: bool = False
    dry_run: bool = True
    executable: str = "rfd3"
    target_structure: str = ""
    backbone_paths: tuple[str, ...] = ()
    target_chain: str = "A"
    binder_chain: str = "A"  # Binder-first contig; verify for the installed RFD3 version.
    target_sequence: str = ""
    output_dir: str = ""
    ligand: str | None = None
    hotspot_residues: tuple[str, ...] = ()
    motif_residues: tuple[str, ...] = ()
    binding_site: tuple[str, ...] = ()
    design_length: int = Field(default=30, ge=5)
    num_designs: int = Field(default=8, ge=1)
    seed: int | None = Field(default=11, ge=0)
    checkpoint_path: str | None = None
    is_non_loopy: bool | None = None
    step_scale: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    gamma_0: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    environment: dict[str, str] = Field(
        default_factory=lambda: {"DEBUG": "false", "WANDB_MODE": "disabled"}
    )
    timeout_seconds: float = Field(default=3600, gt=0, allow_inf_nan=False)


class ProteinMPNNConfig(ConfigModel):
    """Official ProteinMPNN checkout; installation paths must be supplied."""

    python_path: str = sys.executable
    script_path: str = ""
    model_weights_dir: str = ""
    model_name: str = "v_48_020"
    sampling_temp: float = Field(default=0.1, gt=0, allow_inf_nan=False)
    omit_AAs: str = "CX"
    num_sequences_per_backbone: int = Field(default=1, ge=1, strict=True)
    seed: int = Field(default=11, ge=1, strict=True)  # Upstream seed=0 is random.
    timeout_seconds: float = Field(default=3600, gt=0, allow_inf_nan=False)


class AF3BinderConfig(ConfigModel):
    python_path: str = sys.executable
    script_path: str = ""
    model_dir: str = ""
    database_dir: str = ""
    seeds: tuple[int, ...] = (11,)
    num_samples: int = Field(default=5, ge=1)
    target_msa_path: str = ""
    target_msa_cache_dir: str = ""
    binder_msa_mode: Literal["search", "single_sequence"] = "search"
    use_templates: bool = True
    hmmer_bin_dir: str = ""
    environment: dict[str, str] = Field(default_factory=dict)
    timeout_seconds: float = Field(default=21600, gt=0, allow_inf_nan=False)


class ProtenixConfig(ConfigModel):
    enabled: Literal[False] = False
    required: Literal[False] = False


class CompetitionDemoConfig(ConfigModel):
    source_run: str = ""
    selection_mode: Literal["explicit_user", "showcase_all"] = "explicit_user"
    linker_payload_template: str | None = None
    degrader_component_config: str | None = None
    degrader_reference_config: str = str(Path(__file__).parent / "configs/degrader_sjf1528_reference.json")


class UpstreamConfig(ConfigModel):
    rscript: str | None = None
    max_candidates: int = Field(default=100, ge=1, le=500)
    report_reference_run: str = ""
    cpu_threads: int = Field(default=4, ge=1, le=32)
    r_timeout: int = Field(default=14400, ge=60)
    cache_dir: str | None = None
    depmap_release: str = "unknown"
    census_version: str | None = None
    census_provenance_path: str | None = None


class BinderAnalysisConfig(ConfigModel):
    enabled: bool = True
    min_length: int = Field(default=1, ge=1)
    max_length: int = Field(default=10000, ge=1)
    plddt_min: float = Field(default=80, ge=0, le=100, allow_inf_nan=False)
    ipae_max: float = Field(default=10, ge=0, allow_inf_nan=False)
    rmsd_max: float = Field(default=50, ge=0, allow_inf_nan=False)
    rank_by: Literal["i_pae", "binder_rmsd", "binder_plddt", "iptm", "ptm"] = "i_pae"


class OmiCraftConfig(ConfigModel):
    """OmiCraft 전체 설정"""

    execution_profile: Literal["standard", "competition_demo"] = "standard"
    competition_demo: CompetitionDemoConfig = Field(default_factory=CompetitionDemoConfig)
    protenix: ProtenixConfig = Field(default_factory=ProtenixConfig)
    small_molecule: SmallMoleculeConfig = Field(default_factory=SmallMoleculeConfig)
    rfd3: RFD3Config = Field(default_factory=RFD3Config)
    protein_mpnn: ProteinMPNNConfig = Field(default_factory=ProteinMPNNConfig)
    af3_binder: AF3BinderConfig = Field(default_factory=AF3BinderConfig)
    binder_analysis: BinderAnalysisConfig = Field(default_factory=BinderAnalysisConfig)
    data: DataPaths = Field(default_factory=DataPaths)
    upstream: UpstreamConfig = Field(default_factory=UpstreamConfig)
    agent_logging: bool = False


_default_config: OmiCraftConfig | None = None
_scoped_config = ContextVar("omicraft_config", default=None)


def get_config() -> OmiCraftConfig:
    global _default_config
    if _scoped_config.get() is not None:
        return _scoped_config.get()
    if _default_config is None:
        config_path = os.environ.get("OMICRAFT_CONFIG")
        _default_config = (OmiCraftConfig.model_validate(load_json(config_path))
                           if config_path else OmiCraftConfig())
    return _default_config


def set_config(config: OmiCraftConfig) -> None:
    global _default_config
    _default_config = config


@contextmanager
def configuration_scope(config):
    token = _scoped_config.set(config)
    try:
        yield config
    finally:
        _scoped_config.reset(token)

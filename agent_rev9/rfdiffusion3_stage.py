"""RFdiffusion3 input preparation and execution."""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from .result_reuse import cached_tool
from typing import Any, Literal

from .model_runtime import gpu_environment
from .structure_io import protein_residues

RFD3Status = Literal[
    "not_run", "dry_run", "success", "dependency_not_found", "failed", "invalid_input"
]


@dataclass(frozen=True)
class RFdiffusion3Input:
    target_structure: str
    target_chain: str
    design_type: str
    output_dir: str
    binder_chain: str = "A"
    ligand: str | None = None
    hotspot_residues: tuple[str, ...] = ()
    motif_residues: tuple[str, ...] = ()
    binding_site: tuple[str, ...] = ()
    design_length: int = 30
    num_designs: int = 8
    seed: int | None = None
    checkpoint_path: str | None = None
    is_non_loopy: bool | None = None
    step_scale: float | None = None
    gamma_0: float | None = None


@dataclass(frozen=True)
class RFdiffusion3Output:
    status: RFD3Status
    generated_structures: tuple[str, ...] = ()
    metadata: dict[str, Any] = field(default_factory=dict)
    command_or_config: dict[str, Any] = field(default_factory=dict)
    error_message: str = ""

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


class RFdiffusion3Runner:
    """Prepare, execute, and parse RFD3 through its public Foundry CLI."""

    def __init__(
        self,
        executable: str = "rfd3",
        *,
        timeout_seconds: float | None = None,
        gpu_device: int = 0,
        environment: dict[str, str] | None = None,
    ) -> None:
        self.executable = executable
        self.timeout_seconds = timeout_seconds
        self.gpu_device = gpu_device
        self.environment = {"DEBUG": "false", "WANDB_MODE": "disabled", **(environment or {})}

    def available(self) -> bool:
        return shutil.which(self.executable) is not None

    def prepare_input(self, request: RFdiffusion3Input) -> tuple[Path, dict[str, Any]]:
        self._validate(request)
        output_dir = Path(request.output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        residues = [
            n for c, n, _ in protein_residues(request.target_structure) if c == request.target_chain
        ]
        if not residues or any(not n.isdigit() for n in residues):
            raise ValueError("Target chain needs standard residues with positive integer numbering")
        target_contig = ",".join(f"{request.target_chain}{n}" for n in residues)
        allowed = {f"{request.target_chain}{n}" for n in residues}
        for token in request.hotspot_residues + request.motif_residues + request.binding_site:
            if token not in allowed:
                raise ValueError(f"Selection is not a target residue: {token}")
        specification: dict[str, Any] = {
            "input": str(Path(request.target_structure).resolve()),
            "contig": f"{request.design_length},/0,{target_contig}",
        }
        if request.is_non_loopy is not None:
            specification["is_non_loopy"] = request.is_non_loopy
        if request.ligand:
            specification["ligand"] = request.ligand
        selected = tuple(dict.fromkeys(request.hotspot_residues + request.binding_site))
        if selected:
            specification["select_hotspots"] = ",".join(selected)
            specification["infer_ori_strategy"] = "hotspots"

        specification["select_fixed_atoms"] = True
        payload = {"protein_binder": specification}
        input_path = output_dir / "rfd3_input.json"
        input_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return input_path, payload

    def build_command(self, request: RFdiffusion3Input, input_path: Path) -> list[str]:
        command = [
            self.executable,
            "design",
            f"out_dir={Path(request.output_dir).resolve()}",
            f"inputs={input_path.resolve()}",
            "n_batches=1",
            f"diffusion_batch_size={request.num_designs}",
            "prevalidate_inputs=True",
        ]
        if request.seed is not None:
            command.append(f"seed={request.seed}")
        if request.checkpoint_path:
            command.append(f"ckpt_path={Path(request.checkpoint_path).resolve()}")
        for name in ("step_scale", "gamma_0"):
            value = getattr(request, name)
            if value is not None:
                command.append(f"inference_sampler.{name}={value}")
        return command

    @cached_tool('rfd3', ['rfdiffusion3_stage.py', 'structure_io.py'])
    def run(self, request: RFdiffusion3Input, *, dry_run: bool = False) -> RFdiffusion3Output:
        try:
            input_path, payload = self.prepare_input(request)
        except (OSError, ValueError, KeyError, ImportError) as exc:
            return RFdiffusion3Output(status="invalid_input", error_message=str(exc))
        command = self.build_command(request, input_path)
        common = {
            "input_path": str(input_path),
            "official_cli": True,
            "is_mock": dry_run,
            "gpu_device": self.gpu_device,
            "timeout_seconds": self.timeout_seconds,
        }
        command_info = {"command": command, "input_specification": payload}
        if dry_run:
            mock = Path(request.output_dir).resolve() / "mock_rfd3_candidate_000.cif"
            mock.write_text(
                _mock_cif(request.design_length, request.binder_chain), encoding="utf-8"
            )
            return RFdiffusion3Output(
                status="dry_run",
                generated_structures=(str(mock),),
                metadata={**common, "candidate_count": 1},
                command_or_config=command_info,
            )
        if not self.available():
            return RFdiffusion3Output(
                status="dependency_not_found",
                metadata=common,
                command_or_config=command_info,
                error_message=(
                    "RFdiffusion3 CLI 'rfd3' is unavailable. Install rc-foundry[rfd3] "
                    "and a checkpoint with 'foundry install rfd3'."
                ),
            )
        log_path = Path(request.output_dir).resolve() / "rfd3.log"
        existing = set(self.parse_output(request.output_dir))
        if existing:
            return RFdiffusion3Output(
                status="invalid_input",
                metadata=common,
                command_or_config=command_info,
                error_message="Output directory contains structures; use a fresh directory",
            )
        try:
            env = gpu_environment(self.gpu_device, environment=self.environment)
            common["CUDA_VISIBLE_DEVICES"] = env["CUDA_VISIBLE_DEVICES"]
            with log_path.open("w", encoding="utf-8") as log:
                completed = subprocess.run(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=self.timeout_seconds,
                    env=env,
                )
        except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
            return RFdiffusion3Output(
                status="failed",
                metadata=common,
                command_or_config=command_info,
                error_message=str(exc),
            )
        structures = tuple(
            p
            for p in self.parse_output(request.output_dir)
            if Path(p).resolve() != Path(request.target_structure).resolve()
        )
        if completed.returncode != 0 or not structures:
            return RFdiffusion3Output(
                status="failed",
                metadata={**common, "return_code": completed.returncode, "log_path": str(log_path)},
                command_or_config=command_info,
                error_message="RFD3 returned no generated structure; inspect rfd3.log",
            )
        return RFdiffusion3Output(
            status="success",
            generated_structures=structures,
            metadata={
                **common,
                "return_code": 0,
                "log_path": str(log_path),
                "candidate_count": len(structures),
            },
            command_or_config=command_info,
        )

    @staticmethod
    def parse_output(output_dir: str) -> tuple[str, ...]:
        root = Path(output_dir).resolve()
        paths = sorted(
            p
            for p in root.rglob("*")
            if p.is_file() and p.name.endswith((".cif", ".cif.gz", ".pdb", ".pdb.gz"))
        )
        return tuple(
            str(path)
            for path in paths
            if not path.name.startswith("mock_")
            and "af3_validation" not in path.relative_to(root).parts
            and not any("trajector" in part for part in path.relative_to(root).parts)
        )

    @staticmethod
    def _validate(request: RFdiffusion3Input) -> None:
        if request.is_non_loopy is not None and type(request.is_non_loopy) is not bool:
            raise ValueError("is_non_loopy must be a boolean")
        for name in ("step_scale", "gamma_0"):
            value = getattr(request, name)
            if value is not None and (
                type(value) not in (float, int)
                or not math.isfinite(value)
                or value < 0
                or (name == "step_scale" and value == 0)
            ):
                raise ValueError(f"Invalid {name}")
        target = Path(request.target_structure)
        if not target.is_file():
            raise ValueError(f"Target structure not found: {target}")
        if not re.fullmatch(r"[A-Za-z]+", request.target_chain):
            raise ValueError("target_chain is required")
        if request.design_type != "protein_binder":
            raise ValueError("RFD3 design_type must be protein_binder")
        if request.design_length < 5 or request.num_designs < 1:
            raise ValueError("design_length >= 5 and num_designs >= 1 are required")


def _mock_cif(length: int = 30, chain: str = "A") -> str:
    header = """data_mock_rfd3_candidate
_audit_conform.dict_name mock_dry_run_only
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
"""
    return (
        header
        + "".join(
            f"ATOM {i} C CA ALA {chain} {i} {i * 3.8:.1f} 0.0 0.0\n" for i in range(1, length + 1)
        )
        + "#\n"
    )

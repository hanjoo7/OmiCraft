from __future__ import annotations
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from .model_runtime import gpu_environment
from .screening_model_tools import GNINATool
from .small_molecule_config import DockingConfig
from .small_molecule_io import now, unavailable, write_json


def resolve_executable(backend="gnina", *, cli=None, configured=None, environ=None):
    env = os.environ if environ is None else environ
    if backend not in {"gnina", "vina"}:
        raise ValueError("Unsupported docking backend")
    path = shutil.which(backend, path=env.get("PATH", ""))
    choices = ([("cli", cli), ("config", configured),
                ("environment", env.get("GNINA_EXECUTABLE")), ("PATH", path)]
               if backend == "gnina" else
               [("environment", env.get("OMICRAFT_VINA")), ("PATH", path),
                ("config", cli or configured)])
    source, requested = next(((k, v) for k, v in choices if v), ("none", None))
    resolved = shutil.which(str(requested), path=env.get("PATH", "")) if requested else None
    return {"backend": backend, "path": resolved, "source": source,
            "tool_status": "available" if resolved else "not_available",
            "blocker": None if resolved else backend + "_executable_missing"}


def tool_readiness(backend="gnina", *, cli=None, configured=None, environment=None):
    result = resolve_executable(backend, cli=cli, configured=configured)
    result.update(version=None, is_mock=False)
    if not result["path"]:
        return result
    try:
        check = subprocess.run([result["path"], "--version"], capture_output=True, text=True,
                               timeout=30, check=False, env={**os.environ, **(environment or {})})
        result.update(exit_code=check.returncode, version=check.stdout.strip() or check.stderr.strip())
        result["tool_status"] = "ready" if check.returncode == 0 else "not_available"
        result["blocker"] = None if check.returncode == 0 else backend + "_version_failed"
    except (OSError, subprocess.TimeoutExpired) as exc:
        result.update(tool_status="not_available", blocker=backend + "_version_failed", error=str(exc))
    return result


@dataclass(frozen=True)
class DockingRequest:
    candidate_id: str
    prepared_receptor_path: str
    prepared_ligand_path: str
    binding_site_center: tuple
    binding_site_size: tuple
    smiles: str
    output_dir: str
    prepared_ligand_pdbqt_path: str | None = None
    autobox_ligand_path: str | None = None


@dataclass(frozen=True)
class PreparedDockingJob:
    request: DockingRequest
    command: tuple[str, ...]
    output_path: str


class DockingBackend(Protocol):
    name: str
    version: str | None
    is_mock: bool

    def prepare(self, request: DockingRequest) -> PreparedDockingJob: ...
    def run(self, job: PreparedDockingJob) -> dict: ...
    def parse(self, output_dir: Path, smiles: str) -> list[dict]: ...


class NativeDockingBackend:
    def __init__(self, config: DockingConfig, *, gpu_device=0):
        self.config = config.model_copy()
        self.executable_resolution = resolve_executable(config.backend, configured=config.executable)
        self.config.executable = self.executable_resolution["path"]
        self.name = config.backend
        self.version = config.version
        self.gpu_device = config.device if config.device is not None else gpu_device
        self.is_mock = False
        if self.name not in {"vina", "gnina"}:
            raise ValueError("Native docking backend must be vina or gnina")

    def prepare(self, request):
        if not (self.name == "gnina" and request.autobox_ligand_path):
            GNINATool.validate_box(
                {"center": request.binding_site_center, "size": request.binding_site_size}
            )
        directory = Path(request.output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        if not self.config.executable or not shutil.which(self.config.executable):
            raise FileNotFoundError(f"{self.name}_executable_missing")
        ready = tool_readiness(self.name, configured=self.config.executable, environment=self.config.environment)
        if ready["tool_status"] != "ready":
            raise ValueError(ready["blocker"])
        self.version = ready["version"]
        ligand = request.prepared_ligand_path
        if not Path(request.prepared_receptor_path).is_file() or not Path(ligand).is_file():
            raise ValueError("Prepared receptor or ligand missing")
        if self.name == "vina":
            if Path(request.prepared_receptor_path).suffix.lower() != ".pdbqt":
                raise ValueError("Vina requires an explicitly prepared receptor PDBQT")
            if (
                not request.prepared_ligand_pdbqt_path
                or not Path(request.prepared_ligand_pdbqt_path).is_file()
            ):
                raise FileNotFoundError(
                    "Vina requires ligand PDBQT prepared with Meeko; no automatic format or charge guessing"
                )
            ligand = request.prepared_ligand_pdbqt_path
        output = directory / ("poses.sdf" if self.name == "gnina" else "poses.pdbqt")
        if output.exists():
            raise ValueError("Docking output already exists; use a fresh directory")
        command = [
            self.config.executable,
            "--receptor",
            request.prepared_receptor_path,
            "--ligand",
            ligand,
            "--out",
            str(output),
            "--seed",
            str(self.config.seed),
            "--cpu",
            str(self.config.cpu),
            "--exhaustiveness",
            str(self.config.exhaustiveness),
            "--num_modes",
            str(self.config.num_poses),
        ]
        for key, values in ([] if self.name == "gnina" and request.autobox_ligand_path else [
            ("center", request.binding_site_center),
            ("size", request.binding_site_size),
        ]):
            for axis, value in zip("xyz", values):
                command.extend([f"--{key}_{axis}", str(value)])
        if self.name == "gnina":
            if request.autobox_ligand_path:
                if not Path(request.autobox_ligand_path).is_file():
                    raise FileNotFoundError("reference_ligand_missing")
                command += ["--autobox_ligand", request.autobox_ligand_path]
            command += ["--cnn_scoring", self.config.cnn_scoring] + (
                ["--device", "0"] if self.config.use_gpu else ["--no_gpu"]
            )
            if self.config.custom_cnn_model:
                command += ["--cnn_model", self.config.custom_cnn_model]
        return PreparedDockingJob(request, tuple(command), str(output))

    def run(self, job):
        result = {
            "status": "backend_failed",
            "backend": self.name,
            "backend_version": self.version,
            "is_mock": False,
            "poses": [],
            "command": list(job.command),
            "seed": self.config.seed,
            "binary_path": self.config.executable,
            "started_at": now(),
            "warning": [],
            "error_message": "",
        }
        directory = Path(job.request.output_dir)
        try:
            with (directory / "docking.log").open("w") as log:
                completed = subprocess.run(
                    job.command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=self.config.timeout_seconds,
                    env=gpu_environment(
                        self.gpu_device,
                        use_gpu=self.name == "gnina" and self.config.use_gpu,
                        environment=self.config.environment,
                    ),
                )
            result["exit_code"] = completed.returncode
            poses = self.parse(directory, job.request.smiles)
            result["poses"] = poses
            if completed.returncode != 0:
                raise ValueError(f"Docking exited {completed.returncode}")
            if not any(p["docking_status"] == "success" for p in poses):
                raise ValueError("No valid docking poses")
            result["status"] = "success"
        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            result["error_message"] = str(exc)
        result["finished_at"] = now()
        for row in result["poses"]:
            row.update(
                candidate_id=job.request.candidate_id,
                backend=self.name,
                backend_version=self.version,
                is_mock=False,
                started_at=result["started_at"],
                finished_at=result["finished_at"],
            )
        write_json(directory / "docking_result.json", result)
        return result

    def parse(self, output_dir, smiles):
        from rdkit import Chem

        directory = Path(output_dir)
        if self.name == "vina":
            return self._parse_vina(directory, smiles)
        path = directory / "poses.sdf"
        if not path.is_file():
            raise ValueError("Docking SDF missing")
        poses = []
        for index, block in enumerate(_sdf_records(path.read_text())):
            supplier = Chem.SDMolSupplier()
            supplier.SetData(block, removeHs=False)
            try:
                mol = supplier[0]
            except (IndexError, ValueError):
                mol = None
            row = self._pose(index + 1)
            if mol is None:
                row["error_message"] = "Malformed SDF pose"
                poses.append(row)
                continue
            for output, source in [
                ("docking_score", "minimizedAffinity"),
                ("cnn_score", "CNNscore"),
                ("cnn_affinity", "CNNaffinity"),
            ]:
                if mol.HasProp(source):
                    try:
                        from .small_molecule_io import finite

                        row[output] = finite(float(mol.GetProp(source)))
                    except ValueError:
                        row["warning"].append("MALFORMED_" + source)
            destination = directory / f"pose_{index + 1:03d}.sdf"
            with Chem.SDWriter(str(destination)) as writer:
                writer.write(mol)
            row.update(pose_path=str(destination.resolve()), docking_status="success",
                       raw_score_fields={key: mol.GetProp(key) for key in mol.GetPropNames()})
            if row["docking_score"] is None:
                row.update(docking_status="failed", error_message="docking_score_missing")
            poses.append(row)
        return poses

    def _parse_vina(self, directory, smiles):
        from rdkit import Chem

        from .ligand_geometry import molecule_from_atoms

        text = (directory / "poses.pdbqt").read_text()
        blocks = re.findall(r"MODEL\s+\d+\s*\n(.*?)ENDMDL", text, re.S) or (
            [text] if "ATOM" in text else []
        )
        poses = []
        atom_types = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H", "HS": "H"}
        for index, block in enumerate(blocks):
            row = self._pose(index + 1)
            raw_path = directory / f"pose_{index + 1:03d}.pdbqt"
            raw_path.write_text(block)
            row["raw_pose_path"] = str(raw_path.resolve())
            try:
                lines = [line for line in block.splitlines() if line.startswith(("ATOM", "HETATM"))]
                mol = molecule_from_atoms(
                    [
                        atom_types.get(line.split()[-1], line.split()[-1].capitalize())
                        for line in lines
                    ],
                    [
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                        for line in lines
                    ],
                    smiles,
                )
                score = re.search(r"REMARK VINA RESULT:\s*([-+\d.eE]+)", block)
                if score:
                    from .small_molecule_io import finite

                    row["docking_score"] = finite(float(score.group(1)))
                path = directory / f"pose_{index + 1:03d}.sdf"
                with Chem.SDWriter(str(path)) as writer:
                    writer.write(mol)
                row.update(pose_path=str(path.resolve()), docking_status="success")
                if row["docking_score"] is None:
                    row.update(docking_status="failed", error_message="docking_score_missing")
            except (ValueError, RuntimeError) as exc:
                row["error_message"] = str(exc)
            poses.append(row)
        return poses

    @staticmethod
    def _pose(rank):
        return {
            "pose_rank": rank,
            "docking_score": None,
            "cnn_score": None,
            "cnn_affinity": None,
            "pose_path": None,
            "docking_status": "backend_failed",
            "warning": [],
            "error_message": "",
        }


class MockDockingBackend:
    name = "mock"
    version = "synthetic-v1"
    is_mock = True

    def prepare(self, request):
        return PreparedDockingJob(request, (), str(Path(request.output_dir) / "mock_pose.sdf"))

    def run(self, job):
        import numpy as np
        from rdkit import Chem

        from .ligand_preparation import read_ligand

        directory = Path(job.request.output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        mol = read_ligand(job.request.prepared_ligand_path)
        conf = mol.GetConformer()
        xyz = conf.GetPositions()
        xyz += np.asarray(job.request.binding_site_center) - Chem.RemoveHs(
            Chem.Mol(mol)
        ).GetConformer().GetPositions().mean(0)
        for i, point in enumerate(xyz):
            conf.SetAtomPosition(i, tuple(point))
        with Chem.SDWriter(job.output_path) as writer:
            writer.write(mol)
        row = NativeDockingBackend._pose(1)
        row.update(
            candidate_id=job.request.candidate_id,
            backend="mock",
            backend_version=self.version,
            is_mock=True,
            pose_path=job.output_path,
            docking_status="success",
            started_at=now(),
            finished_at=now(),
        )
        return {
            "status": "success",
            "backend": "mock",
            "backend_version": self.version,
            "is_mock": True,
            "poses": [row],
            "started_at": row["started_at"],
            "finished_at": row["finished_at"],
            "warning": ["SYNTHETIC_COORDINATES_NOT_A_PREDICTION"],
            "error_message": "",
        }

    def parse(self, output_dir, smiles):
        raise ValueError("Mock output is generated explicitly by run")


def docking_backend(config, *, gpu_device=0):
    return (
        MockDockingBackend()
        if config.backend == "mock"
        else NativeDockingBackend(config, gpu_device=gpu_device)
    )


def run_docking(backend, request):
    try:
        return backend.run(backend.prepare(request))
    except FileNotFoundError as exc:
        return unavailable(
            "dependency_missing",
            str(exc),
            tool_status="not_available", blocker=str(exc),
            backend=backend.name,
            backend_version=backend.version,
            is_mock=backend.is_mock,
            poses=[],
        )
    except (OSError, ValueError, RuntimeError) as exc:
        return unavailable(
            "backend_failed",
            str(exc),
            backend=backend.name,
            backend_version=backend.version,
            is_mock=backend.is_mock,
            poses=[],
        )


def _sdf_records(text):

    lines = []
    for line in text.splitlines(keepends=True):
        lines.append(line)
        if line.strip() == "$$$$":
            yield "".join(lines)
            lines = []
    if any(line.strip() for line in lines):
        yield "".join(lines) + "\n$$$$\n"

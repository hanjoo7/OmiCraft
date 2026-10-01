"""Subprocess adapters for structure, docking and ADMET tools."""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from .model_runtime import gpu_environment
from .screening_gates import finite, ligand_pose
from .screening_tools import utc_now


from .resource_usage import measured_tool


@dataclass(frozen=True)
class ToolSettings:
    executable: str = ""
    python_path: str = ""
    script_path: str = ""
    config_path: str = ""
    checkpoint_path: str = ""
    cache_dir: str = ""
    required_files: tuple[str, ...] = ()
    timeout_seconds: float = 3600
    seed: int = 11
    num_samples: int = 1
    use_msa: bool = True
    allow_msa_server: bool = False
    use_gpu: bool = True
    cpu: int = 4
    exhaustiveness: int = 8
    num_modes: int = 9

    def __post_init__(self):
        for name in ("seed", "num_samples", "cpu", "exhaustiveness", "num_modes"):
            value = getattr(self, name)
            if type(value) is not int or value < (0 if name == "seed" else 1):
                raise ValueError(f"Invalid {name}")
        if not finite(self.timeout_seconds) or self.timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive and finite")
        for name in ("use_msa", "allow_msa_server", "use_gpu"):
            if type(getattr(self, name)) is not bool:
                raise ValueError(f"{name} must be boolean")
        if not isinstance(self.required_files, (tuple, list)) or any(
            not isinstance(p, str) for p in self.required_files
        ):
            raise ValueError("required_files must be explicit file paths")


def skipped(reason, *, status="not_run", is_mock=False):
    return {
        "status": status,
        "error_message": reason,
        "metrics": {},
        "is_mock": is_mock,
        "execution_mode": "mock" if is_mock else "real",
        "provenance": {"inference_started": False},
    }


def valid_sequence(sequence):
    if not isinstance(sequence, str) or not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("A nonempty standard amino-acid sequence is required")


class ModelTool:
    name = "model"

    def __init__(self, options=None, *, gpu_device=0):
        self.config = ToolSettings(**(options or {}))
        self.gpu_device = gpu_device

    def is_configured(self):
        executable = self.config.executable
        if not executable or not (
            shutil.which(executable)
            or (Path(executable).is_file() and os.access(executable, os.X_OK))
        ):
            return False, f"{self.name} executable is not configured"
        for path in self.config.required_files:
            if not Path(path).is_file():
                return False, f"{self.name} required asset missing: {path}"
        return True, ""

    def _execute(self, command, directory, *, dry_run, input_path=None):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        provenance = {
            "started_at": utc_now(),
            "command": command,
            "config": asdict(self.config),
            "input_path": str(input_path) if input_path else None,
            "inference_started": False,
            "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
        }
        result = {
            "status": "dry_run" if dry_run else "dependency_not_found",
            "metrics": {},
            "is_mock": dry_run,
            "execution_mode": "mock" if dry_run else "real",
            "output_dir": str(directory),
            "provenance": provenance,
        }
        try:
            if not dry_run and (directory / "tool_execution.json").exists():
                raise ValueError(
                    "Tool output directory was already used; select a fresh output directory"
                )
            configured, reason = self.is_configured()
            if dry_run:
                result["error_message"] = "NOT_EXECUTED: dry-run"
            elif not configured:
                result["error_message"] = reason
            else:
                env = gpu_environment(self.gpu_device, use_gpu=self.config.use_gpu)
                provenance["subprocess_CUDA_VISIBLE_DEVICES"] = env["CUDA_VISIBLE_DEVICES"]
                provenance["inference_started"] = True
                log_path = directory / (self.name + ".log")
                result["log_path"] = str(log_path)
                with log_path.open("w") as log:
                    completed = subprocess.run(
                        command,
                        env=env,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        timeout=self.config.timeout_seconds,
                        check=False,
                    )
                provenance["return_code"] = completed.returncode
                result["status"] = "success" if completed.returncode == 0 else "failed"
                if completed.returncode:
                    result["error_message"] = f"{self.name} failed; inspect {log_path.name}"
        except (OSError, subprocess.TimeoutExpired, ValueError) as exc:
            result.update(status="failed", error_message=str(exc))
        provenance["ended_at"] = utc_now()
        (directory / "tool_execution.json").write_text(json.dumps(result, indent=2))
        return result

    @staticmethod
    def _parsed(result, parser):
        if result["status"] == "success":
            try:
                result.update(parser())
            except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
                result.update(status="failed", error_message=str(exc))
        if result.get("output_dir"):
            (Path(result["output_dir"]) / "tool_result.json").write_text(
                json.dumps(result, indent=2)
            )
        return result


class ProtenixTool(ModelTool):
    name = "protenix_v2"

    def is_configured(self):
        configured, reason = super().is_configured()
        if configured and not self.config.required_files:
            return (
                False,
                "Protenix-v2 checkpoint/assets must be explicitly listed in required_files",
            )
        return configured, reason

    @measured_tool
    def run(self, *, target_sequence, binder_sequence, output_dir, dry_run, msa_paths=None):
        valid_sequence(target_sequence)
        valid_sequence(binder_sequence)
        directory = Path(output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        proteins = [
            {"proteinChain": {"sequence": sequence, "count": 1, "id": [chain]}}
            for chain, sequence in (("A", target_sequence), ("B", binder_sequence))
        ]
        msa_paths = msa_paths or {}
        for entry in proteins:
            protein = entry["proteinChain"]
            path = msa_paths.get(protein["id"][0])
            if path:
                if not Path(path).is_file():
                    return skipped(
                        "Protenix precomputed MSA missing",
                        status="dependency_not_found",
                        is_mock=dry_run,
                    )
                protein["unpairedMsaPath"] = str(Path(path).resolve())
        input_path = directory / "input.json"
        input_path.write_text(
            json.dumps([{"name": "binder_validation", "sequences": proteins}], indent=2)
        )
        c = self.config
        if (
            c.use_msa
            and not c.allow_msa_server
            and not all(p["proteinChain"].get("unpairedMsaPath") for p in proteins)
            and not dry_run
        ):
            return skipped(
                "Protenix requires precomputed MSAs or explicit allow_msa_server",
                status="dependency_not_found",
            )
        command = (
            [
                c.executable,
                "pred",
                "-i",
                str(input_path),
                "-o",
                str(directory / "predictions"),
                "-n",
                "protenix-v2",
                "-s",
                str(c.seed),
                "-e",
                str(c.num_samples),
                "--use_msa",
                str(c.use_msa).lower(),
                "--use_template",
                "false",
            ]
            if c.executable
            else []
        )
        result = self._execute(command, directory, dry_run=dry_run, input_path=input_path)
        result["input_json"] = str(input_path)
        return self._parsed(
            result, lambda: {"metrics": self.parse_output(directory / "predictions")}
        )

    @staticmethod
    def parse_output(directory):
        samples = []
        for path in Path(directory).rglob("*_summary_confidence_sample_*.json"):
            model = path.with_name(
                path.name.replace("_summary_confidence_sample_", "_sample_").replace(
                    ".json", ".cif"
                )
            )
            if not model.is_file():
                continue
            raw = json.loads(path.read_text())
            if not isinstance(raw, dict):
                raise ValueError("Protenix summary must be an object")
            metrics = {
                k: raw[k]
                for k in (
                    "ptm",
                    "iptm",
                    "plddt",
                    "gpde",
                    "ranking_score",
                    "has_clash",
                    "chain_ptm",
                    "chain_iptm",
                    "chain_pair_iptm",
                    "chain_plddt",
                    "chain_pair_plddt",
                )
                if k in raw
            }
            metrics.update(model_path=str(model), summary_path=str(path))
            samples.append(metrics)
        if not samples:
            raise ValueError("Protenix returned no matching model and confidence artifacts")
        best = max(
            samples,
            key=lambda x: x["ranking_score"] if finite(x.get("ranking_score")) else -math.inf,
        )
        return {**best, "samples": samples, "sample_count": len(samples)}


class GNINATool(ModelTool):
    name = "gnina"

    @staticmethod
    def validate_box(box):
        if not isinstance(box, dict) or set(box) != {"center", "size"}:
            raise ValueError("MISSING_BINDING_SITE: explicit docking_box center and size required")
        for key in ("center", "size"):
            if (
                not isinstance(box[key], (list, tuple))
                or len(box[key]) != 3
                or any(not finite(x) for x in box[key])
            ):
                raise ValueError("Invalid docking box")
        if any(x <= 0 for x in box["size"]):
            raise ValueError("Docking box size must be positive")

    @staticmethod
    def prepare_ligand(candidate, directory, seed):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        from .support.pipeline_common.chem_validation import validate_smiles_rdkit

        validation = validate_smiles_rdkit(candidate["candidate_id"], candidate["smiles"])
        if not validation.smiles_valid:
            raise ValueError(validation.validation_message)
        if "." in validation.canonical_smiles:
            raise ValueError("Explicit single-component chemical state required for docking")
        if candidate.get("ligand_path"):
            mol, _ = ligand_pose(candidate["ligand_path"], candidate["smiles"])
            source = "provided_3d_ligand"
        else:
            mol = Chem.AddHs(Chem.MolFromSmiles(validation.canonical_smiles))
            if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
                raise ValueError("RDKit conformer generation failed")
            if not AllChem.UFFHasAllMoleculeParams(mol):
                raise ValueError("RDKit UFF parameters unavailable")
            if AllChem.UFFOptimizeMolecule(mol, maxIters=1000) != 0:
                raise ValueError("RDKit conformer optimization did not converge")
            source = "rdkit_conformer"
        AllChem.ComputeGasteigerCharges(mol)
        if any(not math.isfinite(float(a.GetProp("_GasteigerCharge"))) for a in mol.GetAtoms()):
            raise ValueError("Non-finite ligand partial charges")
        path = Path(directory) / "ligand.sdf"
        with Chem.SDWriter(str(path)) as writer:
            writer.write(mol)
        return path, {
            "source": source,
            "canonical_smiles": validation.canonical_smiles,
            "formal_charge": Chem.GetFormalCharge(mol),
            "protonation": "as_supplied; not enumerated",
            "partial_charge_method": "RDKit Gasteiger",
            "validation": validation.as_columns(),
        }

    @measured_tool
    def run(self, *, candidate, context, output_dir, dry_run):
        directory = Path(output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        try:
            self.validate_box(context.get("docking_box"))
            receptor = Path(context["receptor_path"])
            if not receptor.is_file() or receptor.suffix.lower() not in {".pdb", ".pdbqt"}:
                return skipped(
                    "Prepared receptor PDB/PDBQT missing",
                    status="dependency_not_found",
                    is_mock=dry_run,
                )
            ligand, preparation = self.prepare_ligand(candidate, directory, self.config.seed)
        except (ValueError, OSError, KeyError, ImportError) as exc:
            return skipped(str(exc), status="invalid_input", is_mock=dry_run)
        c, box = self.config, context["docking_box"]
        command = [
            c.executable,
            "--receptor",
            str(receptor.resolve()),
            "--ligand",
            str(ligand),
            "--out",
            str(directory / "poses.sdf"),
            "--seed",
            str(c.seed),
            "--cpu",
            str(c.cpu),
            "--exhaustiveness",
            str(c.exhaustiveness),
            "--num_modes",
            str(c.num_modes),
            "--cnn_scoring",
            "rescore",
        ]
        for key in ("center", "size"):
            for axis, value in zip("xyz", box[key]):
                command.extend([f"--{key}_{axis}", str(value)])
        command += ["--device", "0"] if c.use_gpu else ["--no_gpu"]
        result = self._execute(
            command if c.executable else [], directory, dry_run=dry_run, input_path=ligand
        )
        result["preparation"] = preparation
        return self._parsed(result, lambda: self.parse_output(directory / "poses.sdf"))

    @staticmethod
    def parse_output(path):
        from rdkit import Chem

        molecules = [
            mol for mol in Chem.SDMolSupplier(str(path), removeHs=False) if mol is not None
        ]
        if not molecules:
            raise ValueError("GNINA returned no valid SDF pose")
        records = []
        for index, mol in enumerate(molecules):
            scores = {}
            for name in ("minimizedAffinity", "CNNscore", "CNNaffinity", "CNNaffinity_variance"):
                if mol.HasProp(name):
                    try:
                        number = float(mol.GetProp(name))
                        if math.isfinite(number):
                            scores[name] = number
                    except ValueError:
                        pass
            records.append({"rank": index + 1, "scores": scores})
        index = max(
            range(len(records)), key=lambda i: records[i]["scores"].get("CNNscore", -math.inf)
        )
        selected = Path(path).with_name("selected_pose.sdf")
        with Chem.SDWriter(str(selected)) as writer:
            writer.write(molecules[index])
        return {
            **records[index],
            "pose_path": str(selected),
            "poses_path": str(path),
            "poses": records,
        }


class Boltz2Tool(ModelTool):
    name = "boltz2"

    def is_configured(self):
        configured, reason = super().is_configured()
        c = self.config
        if configured and (
            not c.checkpoint_path
            or not Path(c.checkpoint_path).is_file()
            or not c.cache_dir
            or not Path(c.cache_dir).is_dir()
        ):
            return False, "Boltz-2 checkpoint and populated cache directory must be configured"
        return configured, reason

    @measured_tool
    def run(self, *, target_sequence, smiles, output_dir, dry_run, msa_path=None):
        valid_sequence(target_sequence)
        directory = Path(output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        protein = {"id": "A", "sequence": target_sequence}
        if msa_path:
            if not Path(msa_path).is_file():
                return skipped(
                    "Boltz precomputed MSA missing", status="dependency_not_found", is_mock=dry_run
                )
            protein["msa"] = str(Path(msa_path).resolve())
        elif not self.config.allow_msa_server and not dry_run:
            return skipped(
                "Boltz requires precomputed MSA or explicit allow_msa_server",
                status="dependency_not_found",
            )
        import yaml

        path = directory / "input.yaml"
        path.write_text(
            yaml.safe_dump(
                {
                    "version": 1,
                    "sequences": [{"protein": protein}, {"ligand": {"id": "B", "smiles": smiles}}],
                    "properties": [{"affinity": {"binder": "B"}}],
                },
                sort_keys=False,
            )
        )
        c = self.config
        command = [
            c.executable,
            "predict",
            str(path),
            "--out_dir",
            str(directory / "predictions"),
            "--seed",
            str(c.seed),
            "--diffusion_samples",
            str(c.num_samples),
            "--devices",
            "1",
            "--accelerator",
            "gpu" if c.use_gpu else "cpu",
            "--output_format",
            "mmcif",
        ]
        if c.checkpoint_path:
            command += ["--checkpoint", str(Path(c.checkpoint_path).resolve())]
        if c.cache_dir:
            command += ["--cache", str(Path(c.cache_dir).resolve())]
        if c.allow_msa_server:
            command += ["--use_msa_server"]
        result = self._execute(
            command if c.executable else [], directory, dry_run=dry_run, input_path=path
        )
        result["input_yaml"] = str(path)
        return self._parsed(
            result, lambda: {"metrics": self.parse_output(directory / "predictions")}
        )

    @staticmethod
    def parse_output(directory):
        samples = []
        for path in Path(directory).rglob("confidence_*_model_*.json"):
            stem = path.stem.removeprefix("confidence_")
            model = path.with_name(stem + ".cif")
            if not model.is_file():
                continue
            raw = json.loads(path.read_text())
            metrics = {
                k: raw[k]
                for k in (
                    "confidence_score",
                    "ptm",
                    "iptm",
                    "ligand_iptm",
                    "protein_iptm",
                    "complex_plddt",
                    "complex_iplddt",
                    "chains_ptm",
                    "pair_chains_iptm",
                )
                if k in raw
            }
            affinity_path = path.with_name("affinity_" + stem.rsplit("_model_", 1)[0] + ".json")
            if affinity_path.is_file():
                affinity = json.loads(affinity_path.read_text())
                metrics.update({k: v for k, v in affinity.items() if k.startswith("affinity_")})
                metrics["affinity_path"] = str(affinity_path)
            metrics.update(model_path=str(model), summary_path=str(path))
            samples.append(metrics)
        if not samples:
            raise ValueError("Boltz returned no matching model and confidence artifacts")
        best = max(
            samples,
            key=lambda s: s["confidence_score"] if finite(s.get("confidence_score")) else -math.inf,
        )
        return {
            **best,
            "samples": samples,
            "sample_count": len(samples),
            "affinity_units": "log10(IC50 in micromolar); model prediction",
        }


class ADMETTool(ModelTool):
    name = "admet"

    def __init__(self, options=None, *, gpu_device=0):
        super().__init__({"use_gpu": False, **(options or {})}, gpu_device=gpu_device)

    def is_configured(self):
        c = self.config
        if (
            not c.python_path
            or not Path(c.python_path).is_file()
            or not c.script_path
            or not Path(c.script_path).is_file()
            or not c.config_path
            or not Path(c.config_path).is_file()
        ):
            return False, "Existing ADMET pipeline Python, script and config must be configured"
        root = Path(c.python_path).absolute().parent.parent
        package_found = any(root.glob("lib/python*/site-packages/admet_ai"))
        if Path(c.python_path).absolute() == Path(sys.executable).absolute():
            import importlib.util

            package_found = package_found or importlib.util.find_spec("admet_ai") is not None
        if not package_found:
            return False, "admet_ai package is not installed in configured Python environment"
        return True, ""

    @measured_tool
    def run(self, *, candidate, output_dir, dry_run):
        directory = Path(output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "candidate.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["candidate_id", "smiles"])
            writer.writeheader()
            writer.writerow({k: candidate[k] for k in writer.fieldnames})
        c = self.config
        command = (
            [
                c.python_path,
                c.script_path,
                "--config",
                c.config_path,
                "--input",
                str(path),
                "--output",
                str(directory / "predictions"),
                "--backend",
                "admet_ai",
                "--stage",
                "all",
            ]
            if c.python_path and c.script_path
            else []
        )
        result = self._execute(command, directory, dry_run=dry_run, input_path=path)
        return self._parsed(
            result, lambda: self.parse_output(directory / "predictions", candidate["candidate_id"])
        )

    @staticmethod
    def parse_output(directory, candidate_id):
        directory = Path(directory)
        manifest = json.loads((directory / "run_manifest.json").read_text())
        if manifest.get("is_mock") is not False or manifest.get("backend") != "admet_ai":
            raise ValueError("ADMET result is not a real admet_ai execution")
        with (directory / "admet_predictions_long.csv").open() as handle:
            rows = [
                row for row in csv.DictReader(handle) if row.get("candidate_id") == candidate_id
            ]
        metrics = {}
        for row in rows:
            if (
                row.get("is_mock", "").lower() != "false"
                or row.get("prediction_status") != "success"
            ):
                continue
            try:
                value = float(row["endpoint_value"])
                if math.isfinite(value):
                    metrics[row["endpoint_name"]] = value
            except (ValueError, KeyError):
                pass
        if not metrics:
            raise ValueError("ADMET has no successful numeric endpoints for this candidate")
        return {
            "metrics": metrics,
            "endpoints": rows,
            "manifest": manifest,
            "output_path": str(directory / "admet_predictions_long.csv"),
        }

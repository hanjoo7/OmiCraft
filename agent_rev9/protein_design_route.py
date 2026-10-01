"""AF3 binder inputs, output parsing and route state."""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from .result_reuse import cached_tool
from typing import Any, Literal

from .model_runtime import gpu_environment
from .af3_msa_cache import TargetMSACache
from .rfdiffusion3_stage import RFdiffusion3Output


@dataclass(frozen=True)
class AlphaFold3ValidationOutput:
    status: Literal["not_run", "dry_run", "success", "dependency_not_found", "failed"]
    input_json: str = ""
    output_dir: str = ""
    command: tuple[str, ...] = ()
    metrics: dict[str, Any] = field(default_factory=dict)
    error_message: str = ""
    runtime_metadata: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


class AlphaFold3BinderAdapter:
    """AlphaFold 3 binder inference adapter."""

    def __init__(
        self,
        *,
        python_path: str,
        script_path: str,
        model_dir: str,
        database_dir: str,
        gpu_device: int = 0,
        timeout_seconds: float | None = None,
        num_samples: int = 5,
        target_msa_path: str = "",
        target_msa_cache_dir: str = "",
        binder_msa_mode: str = "search",
        use_templates: bool = True,
        hmmer_bin_dir: str = "",
        environment: dict[str, str] | None = None,
    ) -> None:
        self.python_path = python_path
        self.script_path = script_path
        self.model_dir = model_dir
        self.database_dir = database_dir
        self.gpu_device = gpu_device
        self.timeout_seconds = timeout_seconds
        if type(num_samples) is not int or num_samples < 1:
            raise ValueError("num_samples must be a positive integer")
        self.num_samples = num_samples
        self.target_msa_path = target_msa_path
        self.target_msa_cache_dir = target_msa_cache_dir
        self.binder_msa_mode = binder_msa_mode
        self.use_templates = use_templates
        self.hmmer_bin_dir = hmmer_bin_dir
        self.environment = dict(environment or {})

    def input_payload(
        self, *, target_sequence: str, binder_sequence: str, seeds: tuple[int, ...]
    ) -> dict:
        """Build the AF3 input used for execution and cache checks."""
        alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        if (
            not target_sequence
            or not binder_sequence
            or set(target_sequence + binder_sequence) - alphabet
        ):
            raise ValueError("AF3 requires nonempty standard amino-acid sequences")
        if not seeds or any(type(seed) is not int or seed < 0 for seed in seeds):
            raise ValueError("AF3 requires nonnegative integer seeds")
        payload = {
            "name": "rfd3_binder_validation",
            "modelSeeds": list(seeds),
            "sequences": [
                {"protein": {"id": "A", "sequence": target_sequence}},
                {"protein": {"id": "B", "sequence": binder_sequence}},
            ],
            "dialect": "alphafold3",
            "version": 4,
        }
        proteins = [item["protein"] for item in payload["sequences"]]
        if self.target_msa_path:
            cache = json.loads(Path(self.target_msa_path).read_text(encoding="utf-8"))
            if not isinstance(cache, dict) or cache.get("sequence") != target_sequence:
                raise ValueError("Cached target MSA sequence mismatch")
            for key in ("unpairedMsa", "pairedMsa"):
                if not isinstance(cache.get(key), str):
                    raise ValueError(f"Cached target MSA requires string {key}")
                proteins[0][key] = cache[key]
        elif self.target_msa_cache_dir:
            cache = self.msa_cache().load(target_sequence)
            if cache:
                proteins[0].update({key: cache[key] for key in ('unpairedMsa', 'pairedMsa')})
        if self.binder_msa_mode not in {"search", "single_sequence"}:
            raise ValueError("Invalid binder_msa_mode")
        if self.binder_msa_mode == "single_sequence":
            proteins[1].update(unpairedMsa="", pairedMsa="")
        if type(self.use_templates) is not bool:
            raise ValueError("use_templates must be a boolean")
        if not self.use_templates:
            for protein in proteins:
                protein["templates"] = []
        return payload

    def prepare_input(
        self, *, target_sequence: str, binder_sequence: str, output_dir: str, seeds: tuple[int, ...]
    ) -> Path:
        payload = self.input_payload(
            target_sequence=target_sequence, binder_sequence=binder_sequence, seeds=seeds
        )
        path = Path(output_dir).resolve() / "af3_binder_input.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return path

    def msa_cache(self):
        return TargetMSACache(self.target_msa_cache_dir, self.script_path, self.database_dir)

    @cached_tool('af3_binder', ['protein_design_route.py', 'af3_msa_cache.py'])
    def run(self, *, target_sequence, binder_sequence, output_dir, seeds=(11,), dry_run=False):
        cache = self.msa_cache()
        with cache.execution(target_sequence, output_dir, Path(output_dir) / 'af3_predictions', 'A',
                             explicit=bool(self.target_msa_path), dry_run=dry_run) as audit:
            result = self._run(target_sequence=target_sequence, binder_sequence=binder_sequence,
                               output_dir=output_dir, seeds=seeds, dry_run=dry_run)
        result.runtime_metadata['target_msa_cache'] = audit
        return result

    def _run(
        self,
        *,
        target_sequence: str,
        binder_sequence: str,
        output_dir: str,
        seeds: tuple[int, ...] = (11,),
        dry_run: bool = False,
    ) -> AlphaFold3ValidationOutput:
        try:
            input_path = self.prepare_input(
                target_sequence=target_sequence,
                binder_sequence=binder_sequence,
                output_dir=output_dir,
                seeds=seeds,
            )
        except (OSError, ValueError) as exc:
            return AlphaFold3ValidationOutput(status="failed", error_message=str(exc))
        prediction_dir = Path(output_dir).resolve() / "af3_predictions"
        command = (
            self.python_path,
            self.script_path,
            f"--json_path={input_path}",
            f"--model_dir={Path(self.model_dir).resolve()}",
            f"--db_dir={Path(self.database_dir).resolve()}",
            f"--output_dir={prediction_dir}",
            "--gpu_device=0",
            f"--num_diffusion_samples={self.num_samples}",
        )
        if dry_run:
            return AlphaFold3ValidationOutput(
                status="dry_run",
                input_json=str(input_path),
                output_dir=str(prediction_dir),
                command=command,
                metrics={"is_mock": True},
            )
        if not Path(self.python_path).is_file() or not Path(self.script_path).is_file():
            return AlphaFold3ValidationOutput(
                status="dependency_not_found",
                input_json=str(input_path),
                command=command,
                error_message="AlphaFold3 Python executable or run_alphafold.py not found",
            )
        if (
            not Path(self.model_dir).is_dir()
            or not self.model_dir
            or not self.database_dir
            or not Path(self.database_dir).is_dir()
        ):
            return AlphaFold3ValidationOutput(
                status="dependency_not_found",
                input_json=str(input_path),
                command=command,
                error_message="AF3 model_dir and database_dir must be existing directories",
            )
        if prediction_dir.exists() and any(prediction_dir.iterdir()):
            return AlphaFold3ValidationOutput(
                status="failed",
                input_json=str(input_path),
                command=command,
                error_message="AF3 prediction directory must be empty",
            )
        prediction_dir.mkdir(parents=True, exist_ok=True)
        log_path = Path(output_dir).resolve() / "af3_binder.log"
        try:
            env = gpu_environment(self.gpu_device, environment=self.environment)
            if self.hmmer_bin_dir:
                env["PATH"] = (
                    str(Path(self.hmmer_bin_dir).resolve()) + os.pathsep + env.get("PATH", "")
                )
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
            return AlphaFold3ValidationOutput(
                status="failed", input_json=str(input_path), command=command, error_message=str(exc)
            )
        metrics = self.parse_output(prediction_dir)
        status = "success" if completed.returncode == 0 and metrics.get("model_paths") else "failed"
        return AlphaFold3ValidationOutput(
            status=status,
            input_json=str(input_path),
            output_dir=str(prediction_dir),
            command=command,
            metrics=metrics,
            runtime_metadata={
                "gpu_device": self.gpu_device,
                "CUDA_VISIBLE_DEVICES": env["CUDA_VISIBLE_DEVICES"],
            },
            error_message=""
            if status == "success"
            else "AF3 returned no model; inspect af3_binder.log",
        )

    @staticmethod
    def parse_output(output_dir: Path) -> dict[str, Any]:

        models = sorted(
            str(path)
            for path in output_dir.rglob("*model.cif*")
            if path.is_file() and path.name.endswith(("model.cif", "model.cif.zst"))
        )
        samples = []
        summaries = sorted(p for p in output_dir.rglob("*summary_confidences.json*")
                           if p.name.endswith((".json", ".json.gz", ".json.zst")))
        sample_dirs = {p.parent.parent for p in summaries
                       if re.fullmatch(r"seed-\d+_sample-\d+", p.parent.name)}
        seen = set()
        for path in summaries:
            # The top-level ranked model repeats one of the seed/sample results.
            if path.parent in sample_dirs:
                continue
            prefix = path.name.split("summary_confidences.json", 1)[0]
            key = (path.parent, prefix)
            if key in seen:
                continue
            model = next((path.with_name(prefix + "model.cif" + suffix)
                          for suffix in ("", ".zst")
                          if path.with_name(prefix + "model.cif" + suffix).is_file()), None)
            if model is None:
                continue
            try:
                data = _read_af3_json(path)
            except (OSError, ValueError):
                continue
            seen.add(key)
            full_path = next((path.with_name(prefix + "confidences.json" + suffix)
                              for suffix in ("", ".gz", ".zst")
                              if path.with_name(prefix + "confidences.json" + suffix).is_file()),
                             path.with_name(prefix + "confidences.json"))
            full = {}
            warnings = []
            try:
                full = _read_af3_json(full_path)
            except (OSError, ValueError, ImportError):
                warnings.append("full_confidences_unavailable")
            sample = {
                "model_path": str(model),
                "summary_path": str(path),
                "confidences_path": str(full_path) if full_path.is_file() else "",
                "parse_warnings": warnings,
            }
            identity = re.fullmatch(r"seed-(\d+)_sample-(\d+)", path.parent.name)
            sample.update(seed=int(identity[1]) if identity else None,
                          sample_id=int(identity[2]) if identity else None)
            for key in ("ranking_score", "ptm", "iptm", "fraction_disordered"):
                sample[key] = _finite_number(data.get(key))
            # This checkout serializes has_clash as float 0.0/1.0.
            clash = data.get("has_clash")
            sample["has_clash"] = (
                bool(clash) if type(clash) in (bool, int, float) and clash in (0, 1) else None
            )
            for key in ("chain_pair_iptm", "chain_pair_pae_min", "chain_ptm", "chain_iptm"):
                sample[key] = _numeric_array(data.get(key))
            chains = data.get("chain_ids") or full.get("token_chain_ids", [])
            sample["chain_ids"] = (
                list(dict.fromkeys(chains))
                if isinstance(chains, list) and all(isinstance(c, str) for c in chains)
                else []
            )
            # Derived summary of the official full-confidence atom_plddts field.
            scores = full.get("atom_plddts", [])
            atom_chains = full.get("atom_chain_ids", [])
            values = [_finite_number(v) for v in scores] if isinstance(scores, list) else []
            valid = [v for v in values if v is not None and 0 <= v <= 100]
            means = {}
            if isinstance(atom_chains, list) and len(atom_chains) == len(values):
                for chain in sample["chain_ids"]:
                    chain_values = [
                        v
                        for c, v in zip(atom_chains, values)
                        if c == chain and v is not None and 0 <= v <= 100
                    ]
                    means[chain] = sum(chain_values) / len(chain_values) if chain_values else None
            sample["plddt_summary"] = {
                "source": "atom_plddts",
                "mean": sum(valid) / len(valid) if valid else None,
                "by_chain_mean": means,
                "atom_count": len(valid),
            }
            samples.append(sample)
        best = max(
            samples,
            key=lambda s: s["ranking_score"] if s["ranking_score"] is not None else -math.inf,
            default={},
        )
        return {
            "model_paths": models,
            "sample_count": len(samples),
            "samples": samples,
            **{
                key: best.get(key)
                for key in (
                    "model_path",
                    "summary_path",
                    "confidences_path",
                    "ranking_score",
                    "ptm",
                    "iptm",
                    "has_clash",
                    "fraction_disordered",
                    "chain_ids",
                    "chain_pair_iptm",
                    "chain_pair_pae_min",
                    "chain_ptm",
                    "chain_iptm",
                    "plddt_summary",
                    "parse_warnings",
                )
            },
        }


def _read_af3_json(path: Path) -> dict:
    path = Path(path)
    if path.suffix == ".zst":
        import zstandard

        with (
            path.open("rb") as handle,
            zstandard.ZstdDecompressor().stream_reader(handle) as reader,
        ):
            data = json.loads(reader.read())
    elif path.suffix == ".gz":
        import gzip

        with gzip.open(path, "rt", encoding="utf-8") as reader:
            data = json.load(reader)
    else:
        data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("AF3 confidence file must contain an object")
    return data


def _finite_number(value):
    return float(value) if type(value) in (int, float) and math.isfinite(value) else None


def _numeric_array(value):
    if not isinstance(value, list):
        return None
    return [_numeric_array(v) if isinstance(v, list) else _finite_number(v) for v in value]


def extract_protein_sequence(structure_path: str, chain: str = "A") -> str:
    from .structure_io import protein_residues

    return "".join(aa for c, _, aa in protein_residues(structure_path) if c == chain)


def _route_result(rfd3: RFdiffusion3Output, candidates: list[dict]) -> dict[str, Any]:
    first = candidates[0] if candidates else {}
    passed = sum(
        c["structural_validation"].get("passed") is True and not c["is_mock"] for c in candidates
    )
    is_mock = (
        rfd3.status == "dry_run"
        or rfd3.metadata.get("is_mock", False)
        or any(c["is_mock"] for c in candidates)
    )
    if not candidates or rfd3.status not in {"success", "dry_run"}:
        status = rfd3.status if rfd3.status not in {"success", "dry_run"} else "failed"
    elif is_mock and all(c["af3_validation"]["status"] == "dry_run" for c in candidates):
        status = "dry_run"
    else:
        status = "success" if passed else "failed"
    return {
        "status": status,
        "design_mode": "protein_binder",
        "generated_by": "RFdiffusion3",
        "is_mock": is_mock,
        "rfd3": rfd3.as_dict(),
        "candidates": candidates,
        "candidate_count": len(candidates),
        "validated_candidate_count": passed,
        "rfd3_candidate_id": first.get("rfd3_candidate_id"),
        "af3_validation": first.get("af3_validation", {"status": "not_run"}),
        "structural_validation": first.get(
            "structural_validation", {"status": "not_run", "passed": None}
        ),
        "admet_status": "skipped_not_applicable",
        "pipeline_version": "agent_rev9+rfd3-optional-v1",
    }


def protein_route_options(state: dict) -> dict:
    """Explicit state selection wins over the default, which keeps RFD3 disabled."""
    from .configuration import get_config

    config = get_config()
    enabled = state.get("use_rfd3", config.rfd3.enabled)
    mode = state.get("design_mode", "protein_binder" if enabled else "small_molecule")
    if type(enabled) is not bool or type(state.get("dry_run", config.rfd3.dry_run)) is not bool:
        raise ValueError("use_rfd3 and dry_run must be booleans")
    if mode not in {"small_molecule", "protein_binder"}:
        raise ValueError(f"Unsupported design_mode: {mode}")
    if enabled and mode != "protein_binder":
        raise ValueError("RFD3 requires design_mode=protein_binder")
    return {
        "design_mode": mode,
        "use_rfd3": enabled,
        "dry_run": state.get("dry_run", config.rfd3.dry_run),
    }


def protein_route_selected(state: dict) -> bool:
    return protein_route_options(state)["use_rfd3"]


def protein_design_node(state: dict) -> dict:
    """Prepare binder state for the budgeted Screening workflow."""
    options = protein_route_options(state)
    if not options["use_rfd3"]:
        raise ValueError("protein_design requires explicit use_rfd3=True")
    result = state.get("protein_design_result") or _route_result(
        RFdiffusion3Output(status="not_run"), []
    )
    return {
        **options,
        "protein_design_result": result,
        "current_agent": "protein_design",
        "messages": [
            "[RFD3] Handoff prepared; Screening will run RFD3 → ProteinMPNN → AF3 within configured budgets"
        ],
        "errors": [],
        "evidence_cards": [],
    }

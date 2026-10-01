import os
import re
import subprocess
from pathlib import Path

from .model_runtime import gpu_environment
from .af3_msa_cache import TargetMSACache
from .protein_design_route import AlphaFold3BinderAdapter, _read_af3_json
from .small_molecule_io import finite, now, unavailable, write_json
from .support.af3_backbone.af3_input import validate_af3_json


def confidence_chain_ids(data, details=None):
    """Map chain-level arrays using stable unique chain order, not token offsets."""
    raw = data.get('chain_ids') or (details or {}).get('token_chain_ids', [])
    if not isinstance(raw, list) or any(not isinstance(v, str) or not v for v in raw):
        raise ValueError('INVALID_CONFIDENCE_CHAIN_IDS')
    ids = list(dict.fromkeys(raw))
    n = len(ids)
    for key in ('chain_iptm', 'chain_ptm'):
        values = data.get(key)
        if values is not None and (not isinstance(values, list) or len(values) != n):
            raise ValueError('CONFIDENCE_CHAIN_DIMENSION_MISMATCH:' + key)
    for key in ('chain_pair_iptm', 'chain_pair_pae_min'):
        values = data.get(key)
        if values is not None and (not isinstance(values, list) or len(values) != n
                or any(not isinstance(row, list) or len(row) != n for row in values)):
            raise ValueError('CONFIDENCE_CHAIN_DIMENSION_MISMATCH:' + key)
    observed = (details or {}).get('token_chain_ids')
    if observed and set(observed) != set(ids):
        raise ValueError('CONFIDENCE_CHAIN_SET_MISMATCH')
    return ids


def ccd_smiles(path, code):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(path))
    if data.get("_chem_comp.id", [""])[0].upper() != code.upper():
        raise ValueError("CCD code differs from local component")
    kinds = data.get("_pdbx_chem_comp_descriptor.type", [])
    values = data.get("_pdbx_chem_comp_descriptor.descriptor", [])
    options = [
        value
        for kind, value in zip(kinds, values)
        if kind.upper() in {"SMILES_CANONICAL", "SMILES"}
    ]
    if not options:
        raise ValueError("CCD has no SMILES descriptor for identity validation")
    return options


from .resource_usage import measured_tool

class AF3LigandTool:
    def __init__(self, config, *, gpu_device=0):
        self.config = config
        self.gpu_device = gpu_device

    def prepare_input(self, candidate, context, directory):
        from rdkit import Chem

        config = self.config
        common = AlphaFold3BinderAdapter(
            python_path=config.python_path,
            script_path=config.script_path,
            model_dir=config.model_dir,
            database_dir=config.database_dir,
            target_msa_path=context.get("target_msa_path") or config.target_msa_path,
            target_msa_cache_dir=config.target_msa_cache_dir,
            use_templates=config.use_templates,
        )
        payload = common.input_payload(
            target_sequence=context["target_sequence"],
            binder_sequence="A",
            seeds=config.seeds,
        )
        protein = payload["sequences"][0]["protein"]
        protein["id"] = config.protein_chain_id
        smiles = candidate["canonical_smiles"]
        if not smiles or Chem.MolFromSmiles(smiles) is None:
            raise ValueError("Validated ligand SMILES required")
        ligand = {"id": config.ligand_chain_id}
        if candidate.get("ccd_code"):
            if not candidate.get("ccd_path"):
                raise ValueError("A local CCD file is required to verify ligand identity")
            identity = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
            identities = [
                Chem.MolToSmiles(m, isomericSmiles=True)
                for text in ccd_smiles(candidate["ccd_path"], candidate["ccd_code"])
                if (m := Chem.MolFromSmiles(text)) is not None
            ]
            if identity not in identities:
                raise ValueError("CCD_SMILES_IDENTITY_MISMATCH")
            ligand["ccdCodes"] = [candidate["ccd_code"]]
        else:
            ligand["smiles"] = smiles
        payload.update(
            name="ligand_" + re.sub(r"[^A-Za-z0-9_]", "_", candidate["candidate_id"]),
            sequences=[{"protein": protein}, {"ligand": ligand}],
            version=config.input_version,
        )
        validate_af3_json(payload)

        return write_json(Path(directory) / "af3_ligand_input.json", payload)

    def dependency_status(self):
        c = self.config
        if (
            not c.python_path
            or not Path(c.python_path).is_file()
            or not c.script_path
            or not Path(c.script_path).is_file()
        ):
            return "dependency_missing"
        if not c.model_dir or not any(Path(c.model_dir).glob("af3.bin*")):
            return "weights_missing"
        if (
            not c.database_dir
            or not Path(c.database_dir).is_dir()
            or not any(Path(c.database_dir).iterdir())
        ):
            return "database_missing"
        return "available"

    @measured_tool
    def run(self, *, candidate, context, output_dir, dry_run=False, is_mock=False):
        c = self.config
        cache = TargetMSACache(c.target_msa_cache_dir, c.script_path, c.database_dir)
        with cache.execution(context.get('target_sequence', ''), output_dir,
                             Path(output_dir) / 'predictions', c.protein_chain_id,
                             explicit=bool(context.get('target_msa_path') or c.target_msa_path),
                             dry_run=dry_run or is_mock) as audit:
            result = self._run(candidate=candidate, context=context, output_dir=output_dir,
                               dry_run=dry_run, is_mock=is_mock)
        result['target_msa_cache'] = audit
        return result

    def _run(self, *, candidate, context, output_dir, dry_run=False, is_mock=False):
        directory = Path(output_dir)
        directory.mkdir(parents=True, exist_ok=True)
        c = self.config
        try:
            input_path = self.prepare_input(candidate, context, directory)
        except (OSError, ValueError, KeyError) as exc:
            return unavailable("invalid_input", str(exc), samples=[], is_mock=is_mock)
        if dry_run:
            return unavailable(
                "not_run",
                "AF3 inference was not requested",
                samples=[],
                input_json=str(input_path),
                is_mock=is_mock,
            )
        if is_mock:
            return unavailable(
                "not_run",
                "Use the explicit synthetic AF3 fixture backend",
                samples=[],
                input_json=str(input_path),
                is_mock=True,
            )
        dep = self.dependency_status()
        if dep != "available":
            return unavailable(
                dep,
                "AF3 local runtime is not ready; no downloads attempted",
                samples=[],
                input_json=str(input_path),
            )
        predictions = directory / "predictions"
        if predictions.exists() and any(predictions.iterdir()):
            return unavailable(
                "backend_failed", "AF3 prediction directory must be fresh", samples=[]
            )
        command = [
            c.python_path,
            c.script_path,
            f"--json_path={input_path.resolve()}",
            f"--model_dir={Path(c.model_dir).resolve()}",
            f"--db_dir={Path(c.database_dir).resolve()}",
            f"--output_dir={predictions.resolve()}",
            "--gpu_device=0",
            f"--num_diffusion_samples={c.num_samples}",
        ]
        started = now()
        try:
            env = gpu_environment(self.gpu_device, environment=c.environment)
            if c.hmmer_bin_dir:
                env["PATH"] = (
                    str(Path(c.hmmer_bin_dir).resolve()) + os.pathsep + env.get("PATH", "")
                )
            with (directory / "af3.log").open("w") as handle:
                result = subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=c.timeout_seconds,
                    env=env,
                )
            parsed = self.parse_output(predictions, candidate["candidate_id"])
            if result.returncode:
                parsed.update(
                    status="backend_failed",
                    error_message=f"AF3 exited {result.returncode}",
                )
            if parsed.get("is_mock"):
                parsed.update(
                    status="backend_failed",
                    error_message="Synthetic AF3 outputs in a real backend run",
                )
        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            parsed = unavailable("backend_failed", str(exc), samples=[])
        from .small_molecule_io import redact

        return {
            **parsed,
            "input_json": str(input_path),
            "command": redact(command),
            "started_at": started,
            "finished_at": now(),
            "is_mock": parsed.get("is_mock", False),
            "backend": "alphafold3",
            "seeds": list(c.seeds),
            "expected_sample_count": len(c.seeds) * c.num_samples,
            "gpu_device": self.gpu_device,
        }

    def parse_output(self, directory, candidate_id):
        root = Path(directory)
        c = self.config
        locations = {}
        for path in root.rglob("*"):
            match = re.fullmatch(r"seed-(\d+)_sample-(\d+)", path.name)
            if path.is_dir() and match:
                key = tuple(map(int, match.groups()))
                if key in locations:
                    raise ValueError("Duplicate seed/sample outputs; use a fresh run directory")
                locations[key] = path
        expected = {(seed, sample) for seed in c.seeds for sample in range(c.num_samples)}
        samples = [
            self._sample(locations.get(key), candidate_id, *key)
            for key in sorted(expected | locations.keys())
        ]
        successful = {row["seed"] for row in samples if row["af3_status"] == "success"}
        complete = (
            bool(samples)
            and all(row["af3_status"] == "success" for row in samples)
            and set(locations) == expected
        )
        return {
            "status": "success" if complete else "partial" if successful else "backend_failed",
            "samples": samples,
            "successful_seed_count": len(successful),
            "successful_seeds": sorted(successful),
            "expected_seeds": list(c.seeds),
            "expected_sample_count": len(expected),
            "observed_sample_count": len(locations),
            "warning": [] if complete else ["INCOMPLETE_OR_UNEXPECTED_AF3_SAMPLES"],
            "error_message": "",
            "is_mock": any(row.get("is_mock", False) for row in samples),
            "backend": "alphafold3",
        }

    def _sample(self, directory, candidate_id, seed, sample):
        c = self.config
        row = {
            key: None
            for key in (
                "model_path",
                "ranking_score",
                "ptm",
                "iptm",
                "ligand_chain_iptm",
                "protein_ligand_chain_pair_iptm",
                "protein_ligand_pae_min",
                "ligand_atom_plddt_mean",
                "ligand_atom_plddt_min",
                "pocket_occupancy",
                "contact_probability_summary",
                "has_clash",
                "chirality_valid",
            )
        }
        row.update(
            candidate_id=candidate_id,
            seed=seed,
            sample_id=sample,
            ligand_chain_id=c.ligand_chain_id,
            protein_chain_id=c.protein_chain_id,
            af3_status="not_run",
            warning=[],
            error_message="",
        )
        if directory is None:
            row["warning"] = ["MISSING_AF3_SAMPLE"]
            return row

        def find(suffix):
            paths = [p for p in directory.iterdir() if p.name.endswith((suffix, suffix + ".zst"))]
            if suffix == "confidences.json":
                paths = [p for p in paths if "summary_confidences" not in p.name]
            return paths[0] if len(paths) == 1 else None

        model, summary, full = (
            find("model.cif"),
            find("summary_confidences.json"),
            find("confidences.json"),
        )
        row.update(
            model_path=str(model) if model else None,
            summary_path=str(summary) if summary else None,
            confidences_path=str(full) if full else None,
        )
        if model is None:
            row["warning"].append("MISSING_MODEL")
        for path, label in [(summary, "SUMMARY"), (full, "FULL_CONFIDENCES")]:
            if path is None:
                row["warning"].append("MISSING_" + label)
        try:
            data = _read_af3_json(summary) if summary else {}
            details = _read_af3_json(full) if full else {}
            row["is_mock"] = bool(data.get("is_mock") or details.get("is_mock"))
            for key in ("ranking_score", "ptm", "iptm"):
                row[key] = finite(data.get(key))
            for key in ("ptm", "iptm"):
                if row[key] is not None and not 0 <= row[key] <= 1:
                    row[key] = None
                    row["warning"].append("INVALID_" + key.upper())
            clash = data.get("has_clash")
            row["has_clash"] = (
                bool(clash) if type(clash) in (bool, int, float) and clash in (0, 1) else None
            )
            ids = confidence_chain_ids(data, details)
            if c.ligand_chain_id in ids and c.protein_chain_id in ids:
                li, pi = ids.index(c.ligand_chain_id), ids.index(c.protein_chain_id)
                for name, value in [
                    ("ligand_chain_iptm", lambda: data["chain_iptm"][li]),
                    (
                        "protein_ligand_chain_pair_iptm",
                        lambda: data["chain_pair_iptm"][pi][li],
                    ),
                    (
                        "protein_ligand_pae_min",
                        lambda: data["chain_pair_pae_min"][pi][li],
                    ),
                ]:
                    try:
                        row[name] = finite(value())
                    except (KeyError, IndexError, TypeError):
                        row["warning"].append("MISSING_" + name.upper())
            else:
                row["warning"].append("CHAIN_IDS_MISSING_OR_MISMATCH")
            for key in (
                "ligand_chain_iptm",
                "protein_ligand_chain_pair_iptm",
                "protein_ligand_pae_min",
            ):
                value = row[key]
                if value is not None and (
                    value < 0 or key != "protein_ligand_pae_min" and value > 1
                ):
                    row[key] = None
                    row["warning"].append("INVALID_" + key.upper())
            chains, values = (
                details.get("atom_chain_ids", []),
                details.get("atom_plddts", []),
            )
            scores = (
                [finite(v) for chain, v in zip(chains, values) if chain == c.ligand_chain_id]
                if len(chains) == len(values)
                else []
            )
            from .screening_gates import cif_atoms

            ligand_count = (
                sum(a.label_chain_id == c.ligand_chain_id for a in cif_atoms(str(model)))
                if model
                else None
            )
            if (
                scores
                and len(scores) == ligand_count
                and all(v is not None and 0 <= v <= 100 for v in scores)
            ):
                row["ligand_atom_plddt_mean"] = sum(scores) / len(scores)
                row["ligand_atom_plddt_min"] = min(scores)
            else:
                row["warning"].append("LIGAND_PLDDT_MISSING_OR_INVALID")
            tokens = details.get("token_chain_ids", [])
            probs = details.get("contact_probs", [])
            pi = [i for i, ch in enumerate(tokens) if ch == c.protein_chain_id]
            li = [i for i, ch in enumerate(tokens) if ch == c.ligand_chain_id]
            if (
                pi
                and li
                and len(probs) == len(tokens)
                and all(len(x) == len(tokens) for x in probs)
            ):
                numbers = [finite(probs[i][j]) for i in pi for j in li]
                if all(v is not None and 0 <= v <= 1 for v in numbers):
                    row["contact_probability_summary"] = {
                        "mean": sum(numbers) / len(numbers),
                        "max": max(numbers),
                        "pair_count": len(numbers),
                    }
            if row["contact_probability_summary"] is None:
                row["warning"].append("CONTACT_PROBABILITIES_UNAVAILABLE")
            row["af3_status"] = "success" if model and summary and full else "partial"
        except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
            row.update(af3_status="backend_failed", error_message=str(exc))
        return row

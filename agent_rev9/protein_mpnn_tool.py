"""ProteinMPNN sequence design for binder backbones."""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import math
import os
import re
import string
import subprocess
from pathlib import Path
from .result_reuse import cached_tool

from .configuration import ProteinMPNNConfig
from .screening_tools import execution_status, utc_now


class ProteinMPNNTool:
    def __init__(self, config: ProteinMPNNConfig, *, gpu_device: int = 0):
        self.config = config
        self.gpu_device = gpu_device

    def is_configured(self) -> tuple[bool, str]:
        c = self.config
        if not c.script_path or not Path(c.script_path).is_file():
            return False, "Configure protein_mpnn.script_path with the official protein_mpnn_run.py"
        if not Path(c.python_path).is_file() or not os.access(c.python_path, os.X_OK):
            return False, "ProteinMPNN Python executable is unavailable"
        if not re.fullmatch(r"[A-Za-z0-9_]+", c.model_name):
            return False, "Invalid ProteinMPNN model_name"
        if (
            not c.model_weights_dir
            or not (Path(c.model_weights_dir) / (c.model_name + ".pt")).is_file()
        ):
            return False, "ProteinMPNN model_weights_dir/model_name.pt is unavailable"
        if set(c.omit_AAs) - set("ACDEFGHIKLMNPQRSTVWYX") or not set("ACDEFGHIKLMNPQRSTVWY") - set(
            c.omit_AAs
        ):
            return False, "Invalid ProteinMPNN omit_AAs"
        return True, ""

    @staticmethod
    def prepare_backbone(
        structure_path: str, binder_chain: str, output_dir: Path
    ) -> tuple[Path, dict]:
        """Use Biopython parsers/writer; export protein coordinates with explicit chain mapping."""
        from Bio.Data.PDBData import protein_letters_3to1
        from Bio.PDB import PDBIO, Chain, MMCIFParser, Model, PDBParser, Structure

        source = Path(structure_path)
        opener = gzip.open if source.suffix == ".gz" else open
        cif = source.name.lower().endswith((".cif", ".cif.gz"))
        parser = (
            MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
            if cif
            else PDBParser(QUIET=True)
        )
        with opener(source, "rt") as stream:
            model = next(parser.get_structure("source", stream).get_models())
        protein = [
            (chain, [r for r in chain if r.resname in protein_letters_3to1 and r.id[0] == " "])
            for chain in model
        ]
        protein = [(chain, residues) for chain, residues in protein if residues]
        if binder_chain not in {chain.id for chain, _ in protein}:
            raise ValueError(f"Binder chain {binder_chain!r} missing from RFD3 backbone")
        alphabet = string.ascii_uppercase + string.ascii_lowercase + string.digits
        preserved = {
            chain.id for chain, _ in protein if len(chain.id) == 1 and chain.id in alphabet
        }
        available = iter(letter for letter in alphabet if letter not in preserved)
        mapping = {
            chain.id: chain.id if chain.id in preserved else next(available) for chain, _ in protein
        }
        exported = Structure.Structure("backbone")
        exported_model = Model.Model(0)
        exported.add(exported_model)
        residue_map = {}
        for chain, residues in protein:
            if len(residues) > 9999:
                raise ValueError("ProteinMPNN PDB input exceeds 9999 residues per chain")
            output_chain = Chain.Chain(mapping[chain.id])
            exported_model.add(output_chain)
            residue_map[chain.id] = []
            for index, residue in enumerate(residues, 1):
                if not {"N", "CA", "C", "O"}.issubset(residue.child_dict):
                    raise ValueError(f"Incomplete N/CA/C/O backbone: {chain.id}:{residue.id}")
                for atom in residue:
                    if not all(
                        math.isfinite(float(v)) and -999.999 <= v <= 9999.999 for v in atom.coord
                    ):
                        raise ValueError("Invalid or out-of-range PDB coordinates")
                cloned = residue.copy()
                cloned.detach_parent()
                cloned.id = (" ", index, " ")
                output_chain.add(cloned)
                residue_map[chain.id].append(
                    {
                        "source_residue": str(residue.id[1]) + residue.id[2].strip(),
                        "pdb_residue": index,
                    }
                )
        path = output_dir / "backbone.pdb"
        writer = PDBIO()
        writer.set_structure(exported)
        writer.save(str(path))
        return path, {
            "chain_mapping": mapping,
            "binder_pdb_chain": mapping[binder_chain],
            "residue_mapping": residue_map,
            "nonprotein_atoms": "excluded",
            "source_structure": str(source.resolve()),
            "model_index": 0,
        }

    def build_command(
        self, pdb_path: Path, binder_chain: str, output_dir: Path, num_sequences: int
    ) -> list[str]:
        c = self.config
        return [
            str(Path(c.python_path).absolute()),
            str(Path(c.script_path).resolve()),
            "--pdb_path",
            str(pdb_path),
            "--pdb_path_chains",
            binder_chain,
            "--out_folder",
            str(output_dir),
            "--path_to_model_weights",
            str(Path(c.model_weights_dir).resolve()),
            "--model_name",
            c.model_name,
            "--num_seq_per_target",
            str(num_sequences),
            "--batch_size",
            "1",
            "--sampling_temp",
            str(c.sampling_temp),
            "--seed",
            str(c.seed),
            "--omit_AAs",
            c.omit_AAs,
        ]

    @staticmethod
    def parse_output(
        path: Path,
        *,
        candidate_id: str,
        expected_length: int,
        binder_chain: str,
        max_sequences: int,
    ) -> list[dict]:
        from Bio import SeqIO

        records = list(SeqIO.parse(str(path), "fasta"))
        if not records:
            raise ValueError("ProteinMPNN returned empty FASTA")
        designed = re.search(r"designed_chains=(\[[^\]]*\])", records[0].description)
        if not designed or ast.literal_eval(designed.group(1)) != [binder_chain]:
            raise ValueError(
                "ProteinMPNN FASTA must identify only the requested designed binder chain"
            )
        results = []
        seen = set()
        for record in records[1:]:
            sample = re.search(r"(?:^|,\s*)sample=(\d+)(?:,|$)", record.description)
            if not sample:
                continue
            sample_id = int(sample.group(1))
            sequence = str(record.seq)
            match = re.search(r"(?:^|,\s*)score=([^,\s]+)", record.description)
            try:
                score = float(match.group(1)) if match else None
            except ValueError:
                score = None
            reasons = []
            if sample_id in seen:
                reasons.append("DUPLICATE_MPNN_SAMPLE_ID")
            seen.add(sample_id)
            if len(sequence) != expected_length or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
                reasons.append("INVALID_MPNN_BINDER_SEQUENCE")
            if score is None or not math.isfinite(score) or score < 0:
                score = None
                reasons.append("MPNN_SCORE_MISSING_OR_INVALID")
            results.append(
                {
                    "candidate_id": f"{candidate_id}_seq_{len(results) + 1:04d}",
                    "sequence": sequence,
                    "score": score,
                    "sequence_path": str(path),
                    "sample_id": sample_id,
                    "status": "FAILED_MPNN" if reasons else "MPNN_COMPLETED",
                    "is_mock": False,
                    "failure_reasons": reasons,
                    "header": record.description,
                    "model_header": records[0].description,
                }
            )
            if len(results) >= max_sequences:
                break
        if not results:
            raise ValueError("ProteinMPNN returned no sampled sequences")
        return results

    @staticmethod
    def reusable(result: dict, expected_length: int) -> bool:
        """Require a completed sampled FASTA artifact, not just a cached sequence string."""
        try:
            from Bio import SeqIO

            sequence = result.get("sequence", "")
            score = result.get("score")
            return (
                result.get("status") == "MPNN_COMPLETED"
                and type(score) in (int, float)
                and math.isfinite(score)
                and score >= 0
                and len(sequence) == expected_length
                and not set(sequence) - set("ACDEFGHIKLMNPQRSTVWY")
                and any(
                    str(record.seq) == sequence
                    and re.search(r"(?:^|,\s*)sample=\d+", record.description)
                    for record in SeqIO.parse(result["sequence_path"], "fasta")
                )
            )
        except (OSError, KeyError, TypeError, ValueError):
            return False

    @cached_tool('protein_mpnn', ['protein_mpnn_tool.py', 'structure_io.py'])
    def run(self, *, backbone: dict, output_dir: str, num_sequences: int, dry_run: bool) -> dict:
        c = self.config
        directory = Path(output_dir).resolve()
        directory.mkdir(parents=True, exist_ok=True)
        mock = bool(dry_run or backbone.get("is_mock"))
        provenance = {
            "started_at": utc_now(),
            "config": c.model_dump(),
            "seed": c.seed,
            "source_structure": backbone["rfd3_structure_path"],
            "command": [],
            "output_dir": str(directory),
            "tool_version": None,
            "gpu_device": self.gpu_device,
            "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "inference_started": False,
        }

        def finish(status, sequences=None, error=""):
            provenance["ended_at"] = utc_now()
            result = {
                "status": status,
                "execution_status": execution_status(status),
                "sequences": sequences or [],
                "is_mock": mock,
                "error_message": error,
                "provenance": provenance,
            }
            (directory / "protein_mpnn_run.json").write_text(
                json.dumps(result, indent=2), encoding="utf-8"
            )
            return result

        if type(num_sequences) is not int or num_sequences < 1:
            return finish("invalid_input", error="num_sequences must be a positive integer")
        if mock:
            sequence = backbone.get("rfd3_sequence", "")
            path = directory / "mock_sequence.fasta"
            path.write_text(f">MOCK_NOT_PROTEINMPNN\n{sequence}\n")
            return finish(
                "dry_run",
                [
                    {
                        "candidate_id": backbone["candidate_id"] + "_seq_0001",
                        "sequence": sequence,
                        "score": None,
                        "sequence_path": str(path),
                        "status": "NOT_EXECUTED",
                        "is_mock": True,
                        "failure_reasons": ["MPNN_NOT_EXECUTED"],
                    }
                ],
            )
        configured, reason = self.is_configured()
        if not configured:
            return finish("dependency_not_found", error=reason)
        try:
            prediction_dir = directory / "mpnn_predictions"
            if prediction_dir.exists() and any(prediction_dir.iterdir()):
                return finish("failed", error="ProteinMPNN prediction directory must be empty")
            pdb_path, mapping = self.prepare_backbone(
                backbone["rfd3_structure_path"], backbone["binder_chain"], directory
            )
            command = self.build_command(
                pdb_path, mapping["binder_pdb_chain"], prediction_dir, num_sequences
            )
            provenance.update(
                command=command,
                input_pdb=str(pdb_path),
                mapping=mapping,
                script_sha256=hashlib.sha256(Path(c.script_path).read_bytes()).hexdigest(),
            )
            from .screening_model_tools import gpu_environment

            try:
                env = gpu_environment(self.gpu_device)
            except ValueError:
                return finish("dependency_not_found", error="ProteinMPNN GPU index is not visible")
            provenance["subprocess_CUDA_VISIBLE_DEVICES"] = env["CUDA_VISIBLE_DEVICES"]
            provenance["inference_started"] = True
            with (directory / "protein_mpnn.log").open("w") as log:
                completed = subprocess.run(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                    timeout=c.timeout_seconds,
                    env=env,
                )
            provenance["return_code"] = completed.returncode
            if completed.returncode:
                return finish("failed", error="ProteinMPNN failed; inspect protein_mpnn.log")
            results = self.parse_output(
                prediction_dir / "seqs" / (pdb_path.stem + ".fa"),
                candidate_id=backbone["candidate_id"],
                expected_length=len(backbone["rfd3_sequence"]),
                binder_chain=mapping["binder_pdb_chain"],
                max_sequences=num_sequences,
            )
            native_header = results[0].get("model_header", "")
            version = re.search(r"(?:^|,\s*)git_hash=([^,\s]+)", native_header)
            provenance["tool_version"] = (
                version.group(1) if version and version.group(1) != "unknown" else None
            )
            for result in results:
                if set(result["sequence"]) & set(c.omit_AAs):
                    result.update(
                        status="FAILED_MPNN",
                        failure_reasons=result["failure_reasons"]
                        + ["MPNN_OMITTED_AMINO_ACID_GENERATED"],
                    )
            return finish("success", results)
        except (
            OSError,
            ValueError,
            KeyError,
            StopIteration,
            ImportError,
            subprocess.TimeoutExpired,
        ) as exc:
            return finish("failed", error=f"{type(exc).__name__}: {exc}")

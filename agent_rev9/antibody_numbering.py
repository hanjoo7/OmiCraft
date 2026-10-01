from __future__ import annotations
"""IMGT antibody numbering through an isolated ANARCII runtime."""

import importlib.metadata
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Protocol


@dataclass
class NumberingResult:
    sequence_id: str
    backend: str
    backend_version: str | None = None
    scheme: str = "imgt"
    chain_type: str | None = None
    score: float | None = None
    query_start: int | None = None
    query_end: int | None = None
    is_mock: bool = False
    numbering_status: str = "failed"
    numbering: list = field(default_factory=list)
    error_message: str | None = None

    def to_dict(self):
        return asdict(self)


class AntibodyNumberingBackend(Protocol):
    def number(self, sequence_id: str, sequence: str, scheme: str = "imgt") -> NumberingResult: ...


def normalize(sequence_id, sequence, raw, backend, version, scheme="imgt"):
    result = NumberingResult(sequence_id, backend, version, scheme.lower())
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        result.error_message = "invalid_sequence"
        return result
    result.error_message = raw.get("error")
    if result.error_message or not raw.get("numbering"):
        result.error_message = str(result.error_message or "numbering_missing")
        return result
    if str(raw.get("scheme", "")).lower() != scheme.lower():
        result.error_message = "numbering_scheme_mismatch"
        return result
    try:
        result.chain_type = raw["chain_type"]
        if result.chain_type not in {"H", "K", "L"}:
            raise ValueError("not_an_antibody_chain")
        result.score = float(raw["score"]) if raw.get("score") is not None else None
        result.query_start, result.query_end = int(raw["query_start"]), int(raw["query_end"])
        index = result.query_start
        seen = set()
        for (position, insertion), aa in raw["numbering"]:
            if aa == "-":
                continue
            key = (int(position), str(insertion).strip())
            if key in seen or index >= len(sequence) or index < 0 or sequence[index] != aa:
                raise ValueError("numbering_sequence_or_duplicate_mismatch")
            seen.add(key)
            result.numbering.append(
                {
                    "sequence_index": index + 1,
                    "imgt_position": key[0],
                    "imgt_insertion_code": key[1],
                    "amino_acid": aa,
                }
            )
            index += 1
        if not result.numbering:
            raise ValueError("numbering_empty")
        result.numbering_status = "success"
    except (ValueError, KeyError, TypeError) as exc:
        result.error_message = str(exc)
        result.numbering = []
    return result


class AnarciiBackend:
    def __init__(self, python_path=None, timeout=300):
        local = Path(__file__).resolve().parents[2] / "tools/anarcii/.venv/bin/python"
        self.python_path = str(
            python_path
            or os.environ.get("ANARCII_PYTHON")
            or (local if local.is_file() else sys.executable)
        )
        self.timeout = timeout
        self.last_raw = {}

    def availability(self):
        try:
            result = subprocess.run(
                [
                    self.python_path,
                    "-c",
                    'import anarcii,importlib.metadata;print(importlib.metadata.version("anarcii"))',
                ],
                text=True,
                capture_output=True,
                timeout=30,
                env={**os.environ, "CUDA_VISIBLE_DEVICES": ""},
                check=False,
            )
            return {
                "status": "ready" if result.returncode == 0 else "not_installed",
                "version": result.stdout.strip() if result.returncode == 0 else None,
                "python_path": self.python_path,
                "error": result.stderr if result.returncode else None,
            }
        except (OSError, subprocess.TimeoutExpired) as exc:
            return {"status": "not_installed", "version": None, "error": str(exc)}

    def number_many(self, sequences, scheme="imgt"):
        availability = self.availability()
        if availability["status"] != "ready":
            return {key: NumberingResult(key, "ANARCII", scheme=scheme,
                        numbering_status="backend_unavailable", error_message=availability.get("error"))
                    for key in sequences}
        valid = {
            k: v for k, v in sequences.items() if v and not set(v) - set("ACDEFGHIKLMNPQRSTVWY")
        }
        results = {
            k: NumberingResult(k, "ANARCII", scheme=scheme, error_message="invalid_sequence")
            for k in sequences
            if k not in valid
        }
        if not valid:
            return results
        try:
            with tempfile.TemporaryDirectory(prefix="anarcii_") as temporary:
                source, target = Path(temporary) / "input.json", Path(temporary) / "output.json"
                source.write_text(json.dumps({"sequences": valid, "scheme": scheme}))
                process = subprocess.run(
                    [self.python_path, str(Path(__file__).resolve()), str(source), str(target)],
                    capture_output=True,
                    text=True,
                    timeout=self.timeout,
                    check=False,
                    env={**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "4"},
                )
                if process.returncode or not target.is_file():
                    raise RuntimeError(
                        process.stderr.strip() or f"ANARCII exited {process.returncode}"
                    )
                self.last_raw = json.loads(target.read_text())
                for key, sequence in valid.items():
                    raw = self.last_raw["results"].get(key, {"error": "sequence_result_missing"})
                    results[key] = normalize(
                        key, sequence, raw, "ANARCII", self.last_raw["version"], scheme
                    )
        except (OSError, RuntimeError, ValueError, subprocess.TimeoutExpired) as exc:
            for key in valid:
                results[key] = NumberingResult(
                    key, "ANARCII", scheme=scheme, error_message=str(exc)
                )
        return results

    def number(self, sequence_id, sequence, scheme="imgt"):
        return self.number_many({sequence_id: sequence}, scheme)[sequence_id]


class LegacyAnarciBackend:
    def number(self, sequence_id, sequence, scheme="imgt"):
        try:
            from anarci import anarci

            version = importlib.metadata.version("anarci")
            domains, details, _ = anarci([(sequence_id, sequence)], scheme=scheme)
            if not domains[0] or len(domains[0]) != 1:
                raise ValueError("numbering_failed_or_multidomain")
            numbering, start, end = domains[0][0]
            raw = {
                "numbering": numbering,
                "chain_type": details[0][0]["chain_type"],
                "score": details[0][0].get("bitscore"),
                "query_start": start,
                "query_end": end,
                "error": None,
                "scheme": scheme,
            }
            return normalize(sequence_id, sequence, raw, "ANARCI", version, scheme)
        except ImportError as exc:
            return NumberingResult(sequence_id, "ANARCI", scheme=scheme,
                numbering_status="backend_unavailable", error_message=str(exc))
        except (OSError, ValueError, IndexError, KeyError, TypeError) as exc:
            return NumberingResult(sequence_id, "ANARCI", scheme=scheme, error_message=str(exc))


def cdr_region(position):
    for start, end, name in ((27, 38, "CDR1"), (56, 65, "CDR2"), (105, 117, "CDR3")):
        if start <= position <= end:
            return name
    return "framework"


def _worker(source, target):
    from anarcii import Anarcii

    request = json.loads(Path(source).read_text())
    model = Anarcii(cpu=True, ncpu=4)
    results = model.number(request["sequences"])
    if request["scheme"].lower() != "imgt":
        results = model.to_scheme(request["scheme"])
    Path(target).write_text(
        json.dumps(
            {"version": importlib.metadata.version("anarcii"), "results": results},
            default=lambda v: v.item() if hasattr(v, "item") else str(v),
        )
    )


if __name__ == "__main__":
    _worker(*sys.argv[1:])

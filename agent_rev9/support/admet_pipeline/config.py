from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class ADMETConfig:
    root: Path
    raw: dict[str, Any]

    @property
    def candidate_id_column(self) -> str:
        return self.raw["input"].get("candidate_id_column", "candidate_id")

    @property
    def smiles_column(self) -> str:
        return self.raw["input"].get("smiles_column", "smiles")

    @property
    def backend(self) -> str:
        return self.raw["predictor"].get("backend", "mock")

    @property
    def batch_size(self) -> int:
        return int(self.raw["predictor"].get("batch_size", 32))

    @property
    def device(self) -> str:
        return self.raw["predictor"].get("device", "cpu")

    @property
    def selected_endpoints(self) -> dict[str, dict[str, list[str]]]:
        return self.raw.get("report", {}).get("selected_endpoints", {})


REQUIRED_SECTIONS = ("project", "input", "predictor", "report")


def load_config(path: str | Path) -> ADMETConfig:
    config_path = Path(path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    missing = [section for section in REQUIRED_SECTIONS if section not in raw]
    if missing:
        raise ValueError(f"Missing required config sections: {', '.join(missing)}")
    return ADMETConfig(root=config_path.parent.parent, raw=raw)

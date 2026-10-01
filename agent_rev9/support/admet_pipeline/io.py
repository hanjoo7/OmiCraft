from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def read_csv(path: str | Path) -> tuple[list[dict[str, str]], list[str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def write_csv(
    path: str | Path, rows: list[dict[str, Any]], fieldnames: list[str], *, force: bool
) -> None:
    output_path = Path(path)
    if output_path.exists() and not force:
        raise FileExistsError(f"{output_path} already exists. Use --force to overwrite.")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def ordered_fieldnames(
    required: list[str], rows: list[dict[str, Any]], input_order: list[str] | None = None
) -> list[str]:
    names: list[str] = []
    for name in required + list(input_order or []):
        if name not in names:
            names.append(name)
    for row in rows:
        for name in row:
            if name not in names:
                names.append(name)
    return names

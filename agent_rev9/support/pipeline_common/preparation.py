from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_preparation(prepared_path: str | Path, manifest_path: str | Path) -> dict[str, Any]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    if manifest.get("schema_version") != "1.0":
        raise ValueError("Unsupported preparation manifest schema_version")
    if sha256_file(Path(prepared_path)) != manifest.get("output_sha256"):
        raise ValueError("Prepared candidate checksum mismatch")
    return manifest

import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def config_hash(config):
    return hashlib.sha256(
        json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def resolve_path(value, manifest):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (Path(manifest).parent / path).resolve()


def relative_path(path, manifest):
    return os.path.relpath(Path(path).resolve(), Path(manifest).parent.resolve())


def redact(value):
    if isinstance(value, dict):
        return {
            key: (
                "<redacted>"
                if any(
                    part in key.lower()
                    for part in (
                        "api_key",
                        "api_token",
                        "access_token",
                        "authorization",
                        "credential",
                        "password",
                        "secret",
                        "model_dir",
                        "checkpoint",
                        "weights",
                        "environment",
                    )
                )
                else redact(item)
            )
            for key, item in value.items()
        }
    if isinstance(value, (tuple, list)):
        return [redact(item) for item in value]
    if isinstance(value, str) and any(
        flag in value for flag in ("--model_dir=", "ckpt_path=", "--checkpoint=")
    ):
        return value.split("=", 1)[0] + "=<redacted>"
    return value


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    import tempfile
    serialized = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix='.' + path.name, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(serialized)
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path


def write_csv(path, rows, columns=()):
    fields = list(dict.fromkeys([*columns, *(key for row in rows for key in row)]))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, ensure_ascii=False)
                    if isinstance(value, (dict, list, tuple))
                    else value
                    for key, value in row.items()
                }
            )
    return path


def finite(value):
    return float(value) if type(value) in (int, float) and math.isfinite(value) else None


def unavailable(status="not_run", message="", **extra):
    return {
        "status": status,
        "error_message": message,
        "warning": [],
        "is_mock": False,
        **extra,
    }

#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

if not __package__:
    import importlib.util
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("_omicraft_launch", root / "launch.py")
    launcher = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(launcher)
    launcher.load_package()
    __package__ = "agent_rev9.scripts"

from ..support.admet_pipeline.pipeline import INPUT_MODES, STAGES, run_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the minimal ADMET pipeline.")
    parser.add_argument("--config", default="configs/admet.yaml")
    parser.add_argument("--input", default=None)
    parser.add_argument("--input-mode", choices=sorted(INPUT_MODES), default="standalone")
    parser.add_argument("--pipeline1-run", default=None)
    parser.add_argument("--preparation-manifest", default=None)
    parser.add_argument("--output", required=True)
    parser.add_argument("--backend", default=None)
    parser.add_argument("--stage", choices=sorted(STAGES), default="all")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.input_mode == "pipeline1_handoff" and args.input:
        parser.error("--input and --pipeline1-run cannot be used together")
    if args.input_mode == "pipeline1_handoff" and not args.pipeline1_run:
        parser.error("--pipeline1-run is required with --input-mode pipeline1_handoff")
    if args.input_mode == "standalone" and not args.input:
        parser.error("--input is required with --input-mode standalone")
    if args.input_mode == "shared_preparation" and (
        not args.input or not args.preparation_manifest
    ):
        parser.error("shared_preparation requires --input and --preparation-manifest")

    result = run_pipeline(
        config_path=args.config,
        input_path=args.input,
        output_dir=args.output,
        input_mode=args.input_mode,
        pipeline1_run=args.pipeline1_run,
        preparation_manifest=args.preparation_manifest,
        backend_override=args.backend,
        stage=args.stage,
        force=args.force,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

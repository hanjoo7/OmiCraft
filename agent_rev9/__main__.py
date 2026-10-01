import argparse
import sys


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        description="OmiCraft rev9 — Binder / Small-molecule / Qualification pipeline"
    )
    parser.add_argument(
        "route",
        choices=("binder", "small-molecule", "adc", "degrader", "qualify", "full"),
        help=(
            "binder: 단백질 binder 설계 (WF 09 protein binder route); "
            "small-molecule: 저분자 설계 (WF 09 small molecule route); "
            "qualify: WF 01-08 qualification pipeline (표적 선정까지); "
            "full: WF 01-09 전체 파이프라인 (gate 통과 시 설계 포함)"
        ),
    )
    if not arguments or arguments[0] in {"-h", "--help"}:
        parser.print_help()
        return 0

    route = parser.parse_args(arguments[:1]).route

    if route == "binder":
        from .screening_agent import main as run
        return run(arguments[1:])

    if route == "small-molecule":
        from .small_molecule_cli import main as run
        return run(arguments[1:])

    if route in {"adc", "degrader"}:
        import json
        from pathlib import Path

        from .input_paths import local_paths
        from .modality_dispatch import execute
        command = argparse.ArgumentParser(prog=f"omicraft {route}")
        command.add_argument("--input", type=Path, required=True)
        command.add_argument("--config", type=Path)
        command.add_argument("--output", type=Path, required=True)
        command.add_argument("--execute", action="store_true")
        args = command.parse_args(arguments[1:])
        data = local_paths(json.loads(args.input.read_text()), args.input.resolve().parent)
        config = local_paths(json.loads(args.config.read_text()), args.config.resolve().parent) if args.config else {}
        config.update(output_dir=str(args.output), dry_run=not args.execute)
        print(json.dumps(execute({"modality": route, "input": data}, config), indent=2))
        return 0

    # qualify / full — 새 orchestrator.py 기반 파이프라인
    qual_parser = argparse.ArgumentParser(prog=f"omicraft {route}")
    qual_parser.add_argument(
        "--question", "-q",
        default="삼중음성 유방암(TNBC)에서 개발 가능한 치료 표적과 모달리티를 제안하라",
        help="연구 질문",
    )
    qual_parser.add_argument(
        "--disease", default="breast_cancer", help="질환 (기본: breast_cancer)"
    )
    qual_parser.add_argument(
        "--subtype", default="TNBC", help="아형 (기본: TNBC)"
    )
    qual_parser.add_argument(
        "--dry-run", action="store_true", help="API 호출 없이 구조만 실행"
    )
    qual_parser.add_argument(
        "--auto-select", action="store_true", default=False,
        help="테스트용 자동 선택; 실제 사용자 선택을 대신하지 않음"
    )
    qual_parser.add_argument(
        "--verbose", "-v", action="store_true", default=True
    )
    args = qual_parser.parse_args(arguments[1:])

    state = {
        "research_question": args.question,
        "disease": args.disease,
        "subtype": args.subtype,
        "dry_run": args.dry_run,
        "budget": {"auto_select": args.auto_select},
    }

    import logging
    if args.verbose:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s %(name)s — %(message)s",
        )

    if route == "qualify":
        from .orchestrator import run_qualification_pipeline
        result = run_qualification_pipeline(state)
        _print_summary(result, verbose=args.verbose)
    else:  # full
        from .orchestrator import run_full_pipeline
        result = run_full_pipeline(state)
        _print_summary(result, verbose=args.verbose)

    return 0


def _print_summary(state: dict, verbose: bool = True) -> None:
    print("\n" + "=" * 60)
    print("OmiCraft rev9 — 실행 요약")
    print("=" * 60)
    print(f"Gate status   : {state.get('gate_status', 'N/A')}")
    print(f"Gate reason   : {state.get('gate_reason', '')}")

    advance = state.get("advance_targets") or []
    print(f"ADVANCE 표적  : {len(advance)}개")
    for t in advance[:5]:
        tier = t.get("Tier") or t.get("tier") or t.get("best_tier") or "?"
        b = t.get("bio_tier") or t.get("bio_confidence") or "?"
        print(
            f"  {t.get('gene_name',''):12s}  {t.get('modality',''):18s}"
            f"  Tier={tier}  B={b}"
        )

    output_tables = state.get("output_tables") or {}
    if output_tables:
        print("\n출력 테이블:")
        for name, path in output_tables.items():
            print(f"  {name}: {path}")

    if verbose:
        errors = state.get("errors") or []
        if errors:
            print(f"\n오류 {len(errors)}건:")
            for e in errors[:5]:
                print(f"  {e}")


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations
"""
R Script Tool Wrapper
R 스크립트를 subprocess로 호출하고 결과 TSV/JSON을 읽어 반환한다.

제안서: "역할은 계획·해석·재계획을, 수치 계산은 검증된 deterministic module이 담당"
→ Agent(LLM)는 계획/해석, R은 정확한 통계 계산

사용법:
    from r_tool import run_r_analysis
    result = run_r_analysis(counts_path, metadata_path, output_dir)
"""

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from typing import Optional

# R 스크립트 경로
R_SCRIPT = os.path.join(
    os.path.dirname(__file__),
    "20260911_tcga_brca_tnbc_analytic_01_07_hr.R"
)


CONDA_EXE = None
for p in ["/c/ProgramData/anaconda3/Scripts/conda.exe",
          os.path.expanduser("~/anaconda3/Scripts/conda.exe"),
          os.path.expanduser("~/miniconda3/Scripts/conda.exe")]:
    if os.path.exists(p):
        CONDA_EXE = p
        break

R_CONDA_ENV = "r_env"  # R 전용 conda 환경


def find_rscript() -> Optional[str]:
    """Rscript 실행 방법을 반환한다.
    Returns: "conda_run" | Rscript 경로 | None
    """
    from .configuration import get_config
    configured = os.environ.get("OMICRAFT_RSCRIPT") or get_config().upstream.rscript
    if configured:
        return configured if Path(configured).is_file() and os.access(configured, os.X_OK) else None

    # 1. conda run 방식 (Windows에서 가장 안정적)
    if CONDA_EXE:
        try:
            result = subprocess.run(
                [CONDA_EXE, "run", "-n", R_CONDA_ENV, "Rscript", "--version"],
                capture_output=True, text=True, timeout=15
            )
            if result.returncode == 0:
                return "conda_run"
        except Exception:
            pass

    # 2. PATH에서 찾기
    rscript = shutil.which("Rscript")
    if rscript:
        return rscript

    # 3. Windows 기본 설치 경로
    for root in ["C:/Program Files/R", "C:/Program Files (x86)/R"]:
        if os.path.exists(root):
            for ver in sorted(os.listdir(root), reverse=True):
                candidate = os.path.join(root, ver, "bin", "Rscript.exe")
                if os.path.exists(candidate):
                    return candidate

    return None


def run_r_analysis(
    counts_rds: str = None,
    metadata_rds: str = None,
    annotation_rds: str = None,
    clinical_rds: str = None,
    gene_sets_rds: str = None,
    output_dir: str = None,
    r_script: str = None,
    timeout: int | None = None,
    progress_callback=None,
    reuse_results=True,
) -> dict:
    """R 스크립트를 실행하고 결과를 반환한다.

    Args:
        counts_rds: raw count matrix RDS 경로
        metadata_rds: sample metadata RDS 경로
        annotation_rds: gene annotation RDS 경로
        clinical_rds: clinical data RDS 경로
        gene_sets_rds: MSigDB gene sets RDS 경로
        output_dir: 결과 출력 디렉토리
        r_script: R 스크립트 경로 (기본: 동봉된 스크립트)
        timeout: 실행 제한 시간 (초)

    Returns:
        {
            "success": bool,
            "output_dir": str,
            "de_summary": dict,      # DE 요약 (DEG 수, 방향 등)
            "gsea_summary": dict,     # GSEA 요약
            "survival_summary": dict, # Survival 요약
            "report_path": str,       # HTML 리포트 경로
            "error": str,             # 실패 시 에러 메시지
        }
    """
    from .configuration import get_config
    config = get_config().upstream
    timeout = timeout or config.r_timeout
    r_script = os.path.abspath(r_script or R_SCRIPT)
    if not os.path.isfile(r_script):
        return {"success": False, "execution_status": "BLOCKED_INPUT",
                "error": f"R script not found: {r_script}"}
    cache = Path(r_script).parent.parent.parent / "omics_pipeline/output/omicraft_core/20260825_tcga_brca_tnbc_vs_non_tnbc"
    inputs = {
        "HR_COUNTS_RDS": (counts_rds, "03_comparison/20260825_counts_selected_raw_hr.rds"),
        "HR_METADATA_RDS": (metadata_rds, "03_comparison/20260825_metadata_selected_hr.rds"),
        "HR_ANNOTATION_RDS": (annotation_rds, "03_comparison/20260825_gene_annotation_selected_hr.rds"),
        "HR_CLINICAL_RDS": (clinical_rds, "01_download/20260825_BRCA_receptor_clinical_raw_hr.rds"),
        "HR_GENE_SETS_RDS": (gene_sets_rds, "08_gsea/20260825_MSigDB_gene_sets_snapshot_hr.rds"),
    }
    env = os.environ.copy()
    missing = []
    for name, (provided, default) in inputs.items():
        path = Path(provided or env.get(name) or cache / default).expanduser().resolve()
        env[name] = str(path)
        if not path.is_file():
            missing.append(name)
    rscript_path = find_rscript()
    blockers = (["Missing RDS inputs: " + ", ".join(missing)] if missing else [])
    if not rscript_path:
        blockers.append("Rscript not found in server environment; configure OMICRAFT_RSCRIPT or PATH")
    if blockers:
        return {"success": False,
                "execution_status": "BLOCKED_INPUT" if missing else "BLOCKED_ENVIRONMENT",
                "error": "; ".join(blockers), "blockers": blockers,
                "missing_inputs": missing, "output_dir": output_dir}
    if output_dir:
        env["HR_OUTPUT_DIR"] = str(Path(output_dir).expanduser().resolve())

    for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        env[name] = str(config.cpu_threads)
    env["HR_CPU_THREADS"] = str(config.cpu_threads)
    output = Path(env["HR_OUTPUT_DIR"]) if output_dir else None
    if output:
        output.mkdir(parents=True, exist_ok=True)
    fingerprint = hashlib.sha256()
    for path in [r_script, *[env[name] for name in inputs]]:
        fingerprint.update(Path(path).read_bytes())
    parameters = {k: v for k, v in env.items() if k.startswith("HR_") and k != "HR_OUTPUT_DIR"}
    fingerprint.update(json.dumps([parameters, rscript_path], sort_keys=True).encode())
    cache_key = fingerprint.hexdigest()
    cached = Path(config.cache_dir) / cache_key if config.cache_dir else None
    from .resource_usage import record_usage
    if reuse_results and output and cached and (cached / "cache_manifest.json").is_file():
        try:
            manifest = json.loads((cached / "cache_manifest.json").read_text())
            files = manifest.get("files", {})
            valid = bool(files) and all(
                (cached / name).resolve().is_relative_to(cached.resolve())
                and (cached / name).is_file()
                and hashlib.sha256((cached / name).read_bytes()).hexdigest() == digest
                for name, digest in files.items())
            status = json.loads((cached / "report_status.json").read_text())
            valid = valid and status.get("status") == "SUCCEEDED" and "03_de/DE_all_genes.tsv" in files
            if valid:
                shutil.copytree(cached, output, dirs_exist_ok=True)
                parsed = _parse_r_results(str(output))
                parsed.update(execution_mode="VERIFIED_CACHE_REUSE", cache_key=cache_key,
                              source_output_dir=manifest["source_output_dir"])
                record_usage('cache_hit', tool='R analysis', source=manifest['source_output_dir'])
                if progress_callback:
                    progress_callback("R:동일 입력·설정의 완료 결과 재사용")
                return parsed
        except (OSError, ValueError, KeyError, TypeError):
            pass  # Corrupt or unfinished cache entries require a fresh analysis.

    log_dir = Path(tempfile.mkdtemp(prefix=output.name + "-logs-", dir=output.parent)) if output and progress_callback else None
    started = time.time()
    stop_progress = threading.Event()
    def monitor():
        previous = None
        previous_message = None
        while not stop_progress.wait(3):
            path = output / "progress.log" if output else None
            if path and path.is_file():
                stages = path.read_text().strip().splitlines()
                stage = stages[-1] if stages else None
                if stage and stage != previous:
                    previous = stage
                    if progress_callback:
                        progress_callback("R:" + stage)
            if log_dir and (log_dir / "stderr.log").is_file():
                messages = (log_dir / "stderr.log").read_text(errors="replace").strip().splitlines()
                message = messages[-1] if messages else None
                if message and message != previous_message:
                    previous_message = message
                    progress_callback("R:" + message[:240])
    watcher = threading.Thread(target=monitor, daemon=True)
    watcher.start()
    # R 실행 (conda run 또는 직접 호출)
    try:
        if rscript_path == "conda_run":
            cmd = [CONDA_EXE, "run", "-n", R_CONDA_ENV, "Rscript", r_script]
        else:
            cmd = [rscript_path, r_script]

        record_usage('tool_call', tool='R analysis')
        options = {"env": env, "text": True, "timeout": timeout, "cwd": os.path.dirname(r_script)}
        if log_dir:
            with (log_dir / "stdout.log").open("w") as stdout, (log_dir / "stderr.log").open("w") as stderr:
                result = subprocess.run(cmd, stdout=stdout, stderr=stderr, **options)
            result.stdout = (log_dir / "stdout.log").read_text(errors="replace")
            result.stderr = (log_dir / "stderr.log").read_text(errors="replace")
        else:
            result = subprocess.run(cmd, capture_output=True, **options)

        if output:
            (output / "r_stdout.log").write_text(result.stdout or "", encoding="utf-8")
            (output / "r_stderr.log").write_text(result.stderr or "", encoding="utf-8")
        if result.returncode != 0:
            return {
                "success": False,
                "error": result.stderr[-500:] if result.stderr else "R script failed",
                "stdout": result.stdout[-200:] if result.stdout else "",
            }

    except subprocess.TimeoutExpired:
        return {"success": False, "error": f"R script timeout ({timeout}s)"}
    except Exception as e:
        return {"success": False, "error": str(e)}

    finally:
        stop_progress.set()
        watcher.join(timeout=5)
        if log_dir:
            for name in ("stdout", "stderr"):
                source = log_dir / (name + ".log")
                if source.is_file():
                    shutil.copy2(source, output / ("r_" + name + ".log"))
            shutil.rmtree(log_dir)

    # 결과 파싱
    if output_dir is None:
        # R 스크립트가 자동 생성한 output_dir 찾기
        # 최신 디렉토리를 찾는다
        import glob
        parent = os.path.dirname(r_script)
        candidates = sorted(glob.glob(os.path.join(parent, "output", "*_tcga_brca_*")), reverse=True)
        if candidates:
            output_dir = candidates[0]

    if not output_dir or not os.path.exists(output_dir):
        return {"success": False, "error": "Output directory not found", "output_dir": output_dir}

    parsed = _parse_r_results(output_dir)
    parsed.update(execution_mode="FRESH_ANALYSIS", cache_key=cache_key, elapsed_seconds=round(time.time()-started, 2))
    if not parsed.get("output_files", {}).get("de_all"):
        return {**parsed, "success": False, "error": "DE_all_genes.tsv missing"}
    try:
        completed = json.loads((output / "report_status.json").read_text()).get("status") == "SUCCEEDED" if output else False
    except (OSError, ValueError):
        completed = False
    if cached and output and completed:
        cached.mkdir(parents=True, exist_ok=True)
        files = {}
        for source in output.rglob("*"):
            if not source.is_file() or source.suffix == ".rds":
                continue
            relative = source.relative_to(output)
            target = cached / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            files[str(relative)] = hashlib.sha256(target.read_bytes()).hexdigest()
        (cached / "cache_manifest.json").write_text(json.dumps({"files": files, "source_output_dir": str(output),
            "cache_key": cache_key, "created_at": time.time()}, indent=2))
    return parsed


def _parse_r_results(output_dir: str) -> dict:
    """R 스크립트 출력 디렉토리에서 결과를 읽는다."""
    result = {
        "success": True,
        "output_dir": output_dir,
    }

    # DE summary
    de_summary_path = os.path.join(output_dir, "03_de", "de_summary.json")
    if os.path.exists(de_summary_path):
        with open(de_summary_path, "r", encoding="utf-8") as f:
            result["de_summary"] = json.load(f)

    summary = result.get('de_summary', {})
    directions = summary.get('direction_counts')
    if isinstance(directions, dict):
        summary.update(up_in_tnbc=directions.get('Up_in_TNBC', 0),
                       down_in_tnbc=directions.get('Down_in_TNBC', 0))
        summary['n_significant'] = summary['up_in_tnbc'] + summary['down_in_tnbc']

    # GSEA summary
    gsea_summary_path = os.path.join(output_dir, "05_gsea", "gsea_summary.json")
    if os.path.exists(gsea_summary_path):
        with open(gsea_summary_path, "r", encoding="utf-8") as f:
            result["gsea_summary"] = json.load(f)

    # Survival summary
    survival_summary_path = os.path.join(output_dir, "06_survival", "survival_summary.json")
    if os.path.exists(survival_summary_path):
        with open(survival_summary_path, "r", encoding="utf-8") as f:
            result["survival_summary"] = json.load(f)

    # Input manifest
    manifest_path = os.path.join(output_dir, "input_manifest.json")
    if os.path.exists(manifest_path):
        with open(manifest_path, "r", encoding="utf-8") as f:
            result["input_manifest"] = json.load(f)

    # Analysis parameters
    params_path = os.path.join(output_dir, "analysis_parameters.json")
    if os.path.exists(params_path):
        with open(params_path, "r", encoding="utf-8") as f:
            result["analysis_parameters"] = json.load(f)

    # Batch selection
    batch_notes_path = os.path.join(output_dir, "02_batch", "batch_selection_notes.txt")
    if os.path.exists(batch_notes_path):
        with open(batch_notes_path, "r") as f:
            result["batch_selection"] = f.read()

    # Report
    report_path = os.path.join(output_dir, "index.html")
    if os.path.exists(report_path):
        result["report_path"] = report_path

    # Report status
    status_path = os.path.join(output_dir, "report_status.json")
    if os.path.exists(status_path):
        with open(status_path, "r") as f:
            result["report_status"] = json.load(f)

    # Key output files
    result["output_files"] = {
        "deg_table": os.path.join(output_dir, "03_de", "DEG.tsv"),
        "de_all": os.path.join(output_dir, "03_de", "DE_all_genes.tsv"),
        "integrated": os.path.join(output_dir, "07_merge", "DE_survival_all_genes.tsv"),
        "volcano": os.path.join(output_dir, "04_de_plots", "volcano.png"),
        "heatmap": os.path.join(output_dir, "04_de_plots", "DEG_heatmap.png"),
        "pca_before": os.path.join(output_dir, "01_qc", "before_batch_pca_group.png"),
        "pca_after": os.path.join(output_dir, "02_batch", "after_batch_pca_group.png"),
        "hallmark_gsea": os.path.join(output_dir, "05_gsea", "Hallmark_NES_top.png"),
        "gobp_gsea": os.path.join(output_dir, "05_gsea", "GOBP_NES_top.png"),
    }

    # 존재하는 파일만 유지
    result["output_files"] = {k: v for k, v in result["output_files"].items() if os.path.exists(v)}

    return result


if __name__ == "__main__":
    rscript = find_rscript()
    if rscript:
        print(f"Rscript found: {rscript}")
    else:
        print("Rscript NOT found — R not installed")
        print("Install via: conda install -c conda-forge r-base r-deseq2 r-sva")
    print(f"R script: {R_SCRIPT}")
    print(f"R script exists: {os.path.exists(R_SCRIPT)}")

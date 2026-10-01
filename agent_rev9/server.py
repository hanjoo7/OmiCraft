from __future__ import annotations
"""
OmiCraft FastAPI + WebSocket Server
브라우저에서 실시간으로 에이전트 실행 상태를 확인할 수 있다.

엔드포인트:
  GET  /             → 대시보드 HTML
  POST /api/run      → 파이프라인 실행 (비동기)
  WS   /ws           → 실시간 실행 로그 스트리밍
  GET  /api/dossier  → 최종 Dossier 조회
"""

import asyncio
import json
import os
import re
import sys
import threading
import time
import uuid as _uuid
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

if sys.stdout.encoding != "utf-8":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

from fastapi import FastAPI, HTTPException, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from .graph import build_graph
from .execution_control import start_worker, cancel_event
from .candidate_display import state_hidden_genes, visible_gene

# ── 입력 검증 ────────────────────────────────────────────────────────

_BIOMEDICAL_TERMS = {
    # English — cancers
    "cancer", "tumor", "tumour", "carcinoma", "sarcoma", "lymphoma", "leukemia",
    "melanoma", "glioma", "adenocarcinoma", "tnbc", "triple negative",
    "breast", "lung", "colon", "colorectal", "ovarian", "pancreatic", "prostate",
    "hepatocellular", "gastric", "cervical", "endometrial", "renal", "bladder",
    # English — biology / drug discovery
    "gene", "protein", "pathway", "mutation", "biomarker", "expression",
    "target", "therapy", "therapeutic", "treatment", "drug", "inhibitor",
    "antibody", "kinase", "receptor", "oncogene", "suppressor", "ligand",
    "erbb2", "her2", "egfr", "kras", "braf", "tp53", "brca", "pik3ca",
    "immunotherapy", "checkpoint", "pd-l1", "pd-1", "ctla4", "car-t",
    "adc", "antibody-drug", "degrader", "protac", "small molecule",
    "omics", "genomics", "proteomics", "transcriptomics", "rnaseq",
    "tcga", "depmap", "crispr", "clinical", "trial", "biomedical",
    "discovery", "screening", "efficacy", "toxicity", "apoptosis",
    "cell line", "in vivo", "in vitro", "preclinical",
    # Korean
    "암", "종양", "유방암", "폐암", "대장암", "위암", "간암", "췌장암", "혈액암",
    "표적", "치료", "치료제", "항체", "약물", "억제제", "수용체", "단백질",
    "유전자", "바이오마커", "전사", "발현", "돌연변이", "임상", "연구",
    "신약", "후보", "발굴", "스크리닝", "오믹스", "의약", "항암",
}

_OBVIOUS_NONTOPIC = re.compile(
    r"^(안녕|hello+|hi+|hey+|what\s*'?s\s*(your\s*name|the\s*weather|the\s*time)"
    r"|write\s+(a\s+|me\s+)?(poem|story|song|joke|essay)"
    r"|tell\s+me\s+(a\s+)?joke|how\s+are\s+you|좋은\s*(아침|저녁)|감사합니다|고마워"
    r"|날씨|음식|레시피|요리법|\d+\s*[\+\-\*/]\s*\d+|what\s+is\s+\d+|calculate)",
    re.IGNORECASE,
)

_INVALID_MSG = (
    "연구 관련 질문을 입력해 주세요.\n"
    "예: 'TNBC에서 HER2 표적 치료제 후보 발굴' 또는 "
    "'Identify therapeutic targets in triple-negative breast cancer'"
)


def _validate_research_question(question: str) -> tuple:
    """Return (is_valid: bool, error_message: str).
    Falls back to permissive (True) when the LLM is unavailable.
    """
    q_lower = question.lower()

    # 1. Obvious non-topic pattern → reject immediately
    if _OBVIOUS_NONTOPIC.match(question.strip()):
        return False, _INVALID_MSG

    # 2. Any recognised biomedical term → accept
    if any(term in q_lower for term in _BIOMEDICAL_TERMS):
        return True, ""

    # 3. Very short with no biomedical term → reject
    if len(question) < 20:
        return False, _INVALID_MSG

    # 4. LLM classification for ambiguous mid-length inputs
    try:
        from .llm_client import complete, LLMUnavailable  # noqa: F401
        resp = complete(
            "You are a classifier for an oncology drug-discovery research platform. "
            "Reply VALID if the user input is a biomedical or drug-discovery research question. "
            "Reply INVALID if it is unrelated (casual chat, math, creative writing, etc.). "
            "Reply with VALID or INVALID only — no other text.",
            {"question": question},
            role="discovery",
        )
        if resp.text.strip().upper().startswith("INVALID"):
            return False, _INVALID_MSG
    except Exception:
        pass  # LLM unavailable → allow through

    return True, ""

RESULT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed", "agent_results"
)

REPORT_RUNS_ROOT = Path(os.environ.get("OMICRAFT_RUNS_ROOT", Path.cwd() / "runs")).resolve()

app = FastAPI(title="OmiCraft Agent Dashboard")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

# 이벤트 큐 (스레드-안전)
event_queue: deque = deque()

# 파이프라인 콜백 레지스트리 — LangGraph state에 넣으면 msgpack 직렬화 실패하므로 분리 관리
from .pipeline_callbacks import pipeline_callbacks as _pipeline_callbacks
run_state = {"running": False, "result": None, "history": [], "run_id": 0,
             "research_question": "", "question": "", "task_label": "", "execution_kind": "pipeline"}


def restore_dashboard_state():
    path = os.environ.get('OMICRAFT_DASHBOARD_STATE')
    if path and Path(path).is_file():
        saved = json.loads(Path(path).read_text())
        if not saved.get('running'):
            run_state.update({key: saved[key] for key in ('run_id', 'started_at', 'result', 'report_run_id', 'erbb2_report_run_id', 'pipeline_report_run_id', 'workflow') if key in saved})
            report = saved.get('report_run_id') or ''
            kind = saved.get('execution_kind') or ('lab' if report.startswith('lab_') else
                                                  'design' if report.startswith('design_') else 'pipeline')
            question = saved.get('research_question', saved.get('question', '') if kind == 'pipeline' else '')
            run_state.update(research_question=question, question=question, execution_kind=kind,
                             task_label=saved.get('task_label', saved.get('question', '')))
            run_state['history'] = saved.get('events', [])
            state_path = (REPORT_RUNS_ROOT / report / 'pipeline_state.json').resolve()
            state = json.loads(state_path.read_text()) if report and state_path.is_relative_to(REPORT_RUNS_ROOT) and state_path.is_file() else {}
            hidden = state_hidden_genes(state)
            run_state['history'] = [{**event, **{key: [row for row in event[key] if visible_gene(row.get('gene'), hidden)]
                for key in ('evaluated_targets', 'advance_targets') if key in event}}
                for event in run_state['history']]


restore_dashboard_state()


run_lock = threading.Lock()


def begin_pipeline_run(task_label, *, execution_kind="pipeline", research_question=None):
    with run_lock:
        if run_state["running"]:
            return None
        cancel_event.clear()
        run_id = run_state.get("run_id", 0) + 1
        question = research_question if research_question is not None else (task_label if execution_kind == 'pipeline' else run_state.get('research_question', ''))
        run_state.update(running=True, cancelling=False, workflow=None, result=None, history=[], run_id=run_id,
                         question=question, research_question=question, task_label=task_label,
                         execution_kind=execution_kind, started_at=time.time(), report_run_id=None)
        event_queue.clear()
        return run_id


def _execute_pipeline_sync(question: str, options=None, *, reserved=False, record=None):
    """Run one graph and always publish a terminal event on failure."""
    if not reserved and begin_pipeline_run(question) is None:
        return
    start = time.time()
    errors = []
    final = {}
    observer = None
    pending_nodes = []
    thread_id = str(run_state.get("run_id", record.directory.name if record else "default"))
    try:
        graph = build_graph()
        initial_state = {
            "research_question": question, "evidence_cards": [], "messages": [],
            "errors": [], "critic_count": 0,
        }
        initial_state.update(options or {})
        from .workflow_monitor import WorkflowObserver
        monitor_id = record.directory.name if record else 'workflow_' + _uuid.uuid4().hex
        def publish_workflow(snapshot):
            snapshot['execution_id'] = run_state['run_id']
            with run_lock:
                run_state['workflow'] = snapshot
        try:
            observer = WorkflowObserver(graph, monitor_id, record.directory if record else REPORT_RUNS_ROOT / monitor_id, publish_workflow)
        except Exception:
            __import__('logging').getLogger(__name__).exception('Workflow observer unavailable')
        if record:
            initial_state.update(run_id=record.directory.name, _run_directory=str(record.directory))
        def design_progress(route, status, message, report_url):
            entry = {'type': 'node', 'node': 'therapeutic_design', 'elapsed': round(time.time()-start, 1),
                     'execution_status': status, 'messages': [message], 'report_url': report_url}
            with run_lock:
                run_state['history'].append(entry)
                event_queue.append(entry)
            if record:
                record.event(route, status, messages=[message], report_url=report_url)
        def progress(stage):
            entry = {"type": "progress", "node": stage, "elapsed": round(time.time()-start, 1),
                     "messages": [stage + " 실행 중"], "execution_status": "RUNNING"}
            with run_lock:
                run_state["history"].append(entry)
                event_queue.append(entry)
            if record:
                record.event(stage, "RUNNING")
        final = dict(initial_state)
        # 콜백은 state가 아닌 레지스트리에 저장 (LangGraph msgpack 직렬화 오류 방지)
        _pipeline_callbacks[thread_id] = {
            "progress": progress,
            "design_progress": design_progress,
        }
        graph_config = {"recursion_limit": 64, "configurable": {"thread_id": thread_id}, "callbacks": [observer] if observer else []}

        def _stream_graph(input_state, cfg):
            """그래프를 스트리밍하고 interrupt 시 all_eligible이면 자동 재개."""
            for event in graph.stream(input_state, config=cfg):
                for node_name, node_output in event.items():
                    if not isinstance(node_output, dict):
                        # LangGraph __interrupt__ 이벤트 등 비-dict 항목 skip
                        continue
                    yield node_name, node_output
            # Keep the request's design mode when nodes emit only state updates.
            # Resume failures must reach the terminal error handler.
            if hasattr(graph, "get_state"):
                gs = graph.get_state(cfg)
                if gs and gs.next and final.get("design_mode") == "all_eligible":
                    for event in graph.stream(None, config=cfg):
                        for node_name, node_output in event.items():
                            if isinstance(node_output, dict):
                                yield node_name, node_output
                    gs = graph.get_state(cfg)
                pending_nodes[:] = list(gs.next) if gs else []

        for node_name, node_output in _stream_graph(initial_state, graph_config):
            display_hidden = state_hidden_genes({**final, **node_output})
            entry = {
                "type": "node", "node": node_name,
                "elapsed": round(time.time() - start, 1),
                "messages": node_output.get("messages", []),
                "errors": node_output.get("errors", []),
                "execution_status": node_output.get("dossier_status") if node_name == "dossier_review" else node_output.get("execution_status"),
                "review_issues": node_output.get("review_issues", []),
                "gate_status": node_output.get("gate_status"),
                "verdict": node_output.get("critic_verdict"),
                "regression_target": node_output.get("critic_regression_target"),
                "has_qualification_counts": any(k in node_output for k in ("advance_targets", "hold_targets", "reject_targets")),
                "n_hold": len(node_output.get("hold_targets", [])),
                "n_reject": len(node_output.get("reject_targets", [])),
                "evaluated_targets": [
                    {"gene": t.get("gene_name", ""), "tier": t.get("tier", "UNKNOWN"),
                     "modality": t.get("modality", ""), "decision": decision}
                    for key, decision in (("advance_targets", "ADVANCE"), ("hold_targets", "HOLD"), ("reject_targets", "REJECT"))
                    for t in node_output.get(key, []) if visible_gene(t.get("gene_name"), display_hidden)
                ][:20],
                "n_advance": len(node_output.get("advance_targets", [])),
                "advance_targets": [
                    {"gene": t.get("gene_name", ""),
                     "tier": f"{t.get('bio_tier', '')}{t.get('dev_tier', '')}",
                     "modality": t.get("modality", "")}
                    for t in node_output.get("advance_targets", []) if visible_gene(t.get("gene_name"), display_hidden)
                ][:10],
            }
            errors.extend(str(error) for error in entry["errors"])
            with run_lock:
                run_state["history"].append(entry)
                event_queue.append(entry)
            prior_issues = final.get("review_issues", [])
            prior_reused = final.get("reused_stages", [])
            final.update(node_output)
            final["review_issues"] = list(dict.fromkeys(prior_issues + node_output.get("review_issues", [])))
            final["reused_stages"] = prior_reused + node_output.get("reused_stages", [])
            if record:
                from .upstream_reporting import sync_artifacts
                sync_artifacts(record, final.get("r_analysis_result", {}).get("execution_mode", "FRESH_ANALYSIS"))
                record.event(node_name, entry["execution_status"] or ("FAILED" if entry["errors"] else "COMPLETED"),
                             messages=entry["messages"], errors=entry["errors"])
                (record.directory / "pipeline_state.json").write_text(json.dumps(final, ensure_ascii=False, indent=2, default=str))
        terminal = {
            "type": "done", "total_time": round(time.time() - start, 1),
            "execution_status": (final.get("execution_status") if final.get("execution_status") in
                                 {"FAILED", "BLOCKED", "BLOCKED_INPUT", "BLOCKED_ENVIRONMENT", "NOT_RUN"}
                                 else "COMPLETED_WITH_ERRORS" if errors else "PARTIAL" if any(h.get("execution_status") == "PARTIAL" for h in run_state["history"]) else "COMPLETED"),
            "blockers": final.get("blockers", []),
            "errors": list(dict.fromkeys(errors)),
            "review_issues": list(dict.fromkeys(final.get("review_issues", []))),
            "design_outcome": final.get("design_outcome", {}),
            "reused_stages": final.get("reused_stages", []),
            "critic_count": final.get("critic_count", 0),
            "final_verdict": final.get("critic_verdict", ""),
            "demo_pipeline_status": final.get("demo_pipeline_status"),
            "scientific_validation_status": final.get("scientific_validation_status"),
            "report_run_id": final.get("run_id"),
            "node_path": " → ".join(h["node"] for h in run_state["history"] if h.get("type") == "node"),
        }
    except Exception as exc:
        terminal = {
            "type": "error", "execution_status": "FAILED",
            "message": str(exc), "error_type": type(exc).__name__,
            "total_time": round(time.time() - start, 1),
        }
    finally:
        _pipeline_callbacks.pop(thread_id, None)
    if record:
        from .upstream_reporting import finish_report
        terminal["report_run_id"] = record.directory.name
        try:
            finish_report(record, final, terminal)
            terminal["scientific_validation_status"] = record.data["validation_decision"]
        except Exception as exc:
            terminal["execution_status"] = "COMPLETED_WITH_ERRORS" if terminal["type"] == "done" else "FAILED"
            terminal.setdefault("errors", []).append("report: " + str(exc))
    if observer:
        try:
            observer.finish(terminal, pending_nodes)
        except Exception:
            __import__("logging").getLogger(__name__).exception("Workflow monitor finalization failed")
    with run_lock:
        run_state["history"].append(terminal)
        run_state["result"] = terminal
        run_state["running"] = False
        event_queue.append(terminal)


def run_pipeline_sync(question: str, options=None, *, reserved=False):
    from datetime import datetime, timezone

    from .configuration import configuration_scope, get_config
    from .run_records import RunRecord
    if not reserved and begin_pipeline_run(question) is None:
        return
    options = options or {}
    if options.get("execution_profile") == "competition_demo" or options.get("advance_targets"):
        return _execute_pipeline_sync(question, options, reserved=True)
    config = get_config().model_copy(deep=True)
    name = "upstream_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ_") + _uuid.uuid4().hex[:8]
    directory = Path(config.data.base_dir) / name
    config.data.base_dir = str(directory)
    try:
        record = RunRecord(directory, "UPSTREAM")
        record.data.update(question=question, execution_profile="upstream_analysis", running=True, execution_status="PARTIAL")
        with run_lock:
            run_state["report_run_id"] = name
            run_state["pipeline_report_run_id"] = name
        from .resource_usage import resource_scope
        with configuration_scope(config), resource_scope(directory):
            _execute_pipeline_sync(question, options, reserved=True, record=record)
    except Exception as exc:
        terminal = {"type": "error", "execution_status": "FAILED", "message": str(exc), "total_time": 0}
        with run_lock:
            run_state.update(running=False, result=terminal)
            run_state["history"].append(terminal)
            event_queue.append(terminal)


@app.websocket("/ws")
async def websocket_endpoint(ws: WebSocket):
    await ws.accept()
    # 기존 히스토리 전송
    for h in run_state["history"]:
        await ws.send_text(json.dumps(h, ensure_ascii=False, default=str))

    cursor = len(run_state["history"])
    run_id = run_state.get("run_id", 0)
    try:
        while True:
            # 새 run이 시작되면 cursor 리셋
            current_run_id = run_state.get("run_id", 0)
            if current_run_id != run_id:
                run_id = current_run_id
                cursor = 0

            # 새 이벤트 전송
            history = run_state["history"]
            while cursor < len(history):
                await ws.send_text(json.dumps(history[cursor], ensure_ascii=False, default=str))
                cursor += 1

            await asyncio.sleep(0.15)  # 150ms polling
    except WebSocketDisconnect:
        pass


@app.post("/api/run")
async def api_run(body: dict = None):
    from .research_scenarios import RESEARCH_SCENARIOS
    body = dict(body or {})
    scenario = body.get('scenario')
    if scenario is not None:
        if not isinstance(scenario, str) or scenario not in RESEARCH_SCENARIOS:
            raise HTTPException(422, '지원하지 않는 연구질문입니다.')
        if type(body.get('reuse_results', True)) is not bool:
            raise HTTPException(422, 'reuse_results must be boolean')
        body['question'] = RESEARCH_SCENARIOS[scenario]
        if scenario == 'erbb2_design':
            from .web_execution import start_erbb2_batch
            return JSONResponse(start_erbb2_batch({
                'parent_run_id': 'reference_erbb2', 'gene': 'ERBB2', 'modality': 'ALL',
                'input': {}, 'reuse_results': body.get('reuse_results', True),
            }, research_question=body['question']))
        body['design_mode'] = 'all_eligible'
    question = body.get("question", "")
    if not isinstance(question, str) or not question.strip() or len(question) > 10000:
        raise HTTPException(status_code=422, detail="Enter a question of 1–10000 characters")
    question = question.strip()
    is_valid, validation_error = (True, "") if scenario is not None else _validate_research_question(question)
    if not is_valid:
        raise HTTPException(status_code=422, detail=validation_error)
    design_mode = (body or {}).get('design_mode', 'none')
    if design_mode not in {'none', 'all_eligible'}:
        raise HTTPException(status_code=422, detail='Unknown design_mode')
    if type((body or {}).get('reuse_results', True)) is not bool:
        raise HTTPException(status_code=422, detail='reuse_results must be boolean')
    structural_backend = (body or {}).get('structural_backend', 'af3')
    if structural_backend not in {'af3', 'boltz2_legacy'}:
        raise HTTPException(status_code=422, detail='Unknown structural_backend: choose af3 or boltz2_legacy')
    run_id = begin_pipeline_run(question)
    if run_id is None:
        return JSONResponse({"status": "already_running", "run_id": run_state["run_id"]})
    option_fields = {
        "advance_targets", "user_selection_events", "execution_inputs", "modality_config",
        "output_dir", "dry_run", "execution_profile", "competition_demo", "design_mode",
        "structural_backend", "reuse_results",
    }
    options = {key: value for key, value in (body or {}).items() if key in option_fields}
    try:
        start_worker(run_pipeline_sync, question, options, reserved=True)
    except RuntimeError as exc:
        with run_lock:
            run_state["running"] = False
        raise HTTPException(status_code=503, detail="Could not start the pipeline worker") from exc
    return JSONResponse({"status": "started", "question": question, "run_id": run_id})


@app.post("/api/stop")
async def api_stop(body: dict):
    with run_lock:
        if type(body.get('run_id')) is not int or body['run_id'] != run_state['run_id']:
            raise HTTPException(409, '실행이 변경되었습니다. 현재 실행을 확인하세요.')
        if not run_state['running']:
            return {'status': 'not_running', 'run_id': run_state['run_id']}
        cancel_event.set()
        run_state['cancelling'] = True
        return {'status': 'stopping', 'run_id': run_state['run_id']}


@app.get("/api/events")
async def api_events(run_id: int = 0, cursor: int = 0):
    if run_id < 0 or cursor < 0:
        raise HTTPException(status_code=422, detail="Invalid event cursor")
    with run_lock:
        history = run_state["history"]
        if run_id != run_state["run_id"] or cursor > len(history):
            cursor = 0
        return {
            "run_id": run_state["run_id"], "running": run_state["running"],
            "cancelling": run_state.get("cancelling", False),
            "cursor": len(history), "events": list(history[cursor:]),
            "question": run_state.get("research_question", ""),
            "research_question": run_state.get("research_question", ""),
            "task_label": run_state.get("task_label", ""),
            "execution_kind": run_state.get("execution_kind", "pipeline"),
            "started_at": run_state.get("started_at"), "result": run_state["result"],
            "report_run_id": run_state.get("report_run_id"),
            "erbb2_report_run_id": run_state.get("erbb2_report_run_id"),
            "pipeline_report_run_id": run_state.get("pipeline_report_run_id"),
        }


@app.get("/health")
async def health():
    return {"status": "ok", "version": "8"}


@app.get('/api/resources')
async def resource_status():
    from starlette.concurrency import run_in_threadpool
    from .gpu_scheduler import available_devices
    try:
        available = await run_in_threadpool(available_devices)
        ready = {device.physical for device in available}
        return {'devices': [{'index': index, 'available': index in ready} for index in (4, 5, 6, 7)]}
    except ValueError:
        return {'devices': [{'index': index, 'available': None} for index in (4, 5, 6, 7)]}


@app.get("/api/status")
async def api_status():
    return JSONResponse({
        "running": run_state["running"],
        "n_events": len(run_state["history"]),
        "run_id": run_state["run_id"],
        "result": run_state["result"],
    })


@app.get("/api/dossier")
async def api_dossier(run_id: str | None = None):
    from .run_records import run_directory
    selected = run_id or run_state.get('pipeline_report_run_id')
    if not selected:
        raise HTTPException(404, 'No pipeline dossier')
    try:
        directory = run_directory(REPORT_RUNS_ROOT, run_id=selected)
        path = (directory / 'agent_results/final_dossier.json').resolve()
        if not path.is_relative_to(directory):
            raise ValueError('dossier_outside_run')
        if not path.is_file():
            raise HTTPException(404, 'Dossier has not been generated for this run')
        return JSONResponse(json.loads(path.read_text()))
    except (ValueError, OSError) as exc:
        raise HTTPException(400, str(exc)) from exc


@app.get("/api/omics_charts")
async def api_omics_charts():
    """DEG/Pathway/Cell-of-origin 차트 데이터 제공"""
    import pandas as pd
    charts = {}

    # Volcano data (top 200 genes for rendering)
    deg_path = os.path.join(RESULT_DIR, "..", "deg_results", "deg_tnbc_vs_nontnbc_full.csv")
    if os.path.exists(deg_path):
        df = pd.read_csv(deg_path, index_col=0)
        df = df.dropna(subset=["log2FoldChange", "padj"])
        df = df[df["padj"] > 0]
        df["nlog10"] = -(df["padj"].apply(lambda x: __import__('math').log10(max(x, 1e-300))))
        # Sample: top 100 up + top 100 down + random 200
        up = df[df["log2FoldChange"] > 1].nlargest(100, "nlog10")
        down = df[df["log2FoldChange"] < -1].nlargest(100, "nlog10")
        ns = df[(df["log2FoldChange"].abs() <= 1)].sample(min(200, len(df)), random_state=42)
        sample = pd.concat([up, down, ns])
        charts["volcano"] = [
            {"x": round(r["log2FoldChange"], 2), "y": round(min(r["nlog10"], 300), 1),
             "g": r.get("gene_name", ""), "s": "up" if r["log2FoldChange"] > 1 else "down" if r["log2FoldChange"] < -1 else "ns"}
            for _, r in sample.iterrows()
        ]

    # GSEA data
    gsea_path = os.path.join(RESULT_DIR, "..", "pathway_results", "gsea_tnbc_vs_nontnbc_combined.csv")
    if os.path.exists(gsea_path):
        gdf = pd.read_csv(gsea_path)
        sig = gdf[gdf["FDR q-val"] < 0.25].sort_values("NES", ascending=False)
        charts["gsea"] = [
            {"term": row["Term"][:40], "nes": round(row["NES"], 2), "fdr": round(row["FDR q-val"], 4), "lib": row.get("library", "")}
            for _, row in sig.head(15).iterrows()
        ]

    # PCA + Distance matrix — computed from count matrix + metadata
    import numpy as np
    counts_path = os.path.join(RESULT_DIR, "..", "tcga_brca_counts.csv")
    meta_path = os.path.join(RESULT_DIR, "..", "sample_metadata.csv")
    if os.path.exists(counts_path) and os.path.exists(meta_path):
        try:
            meta = pd.read_csv(meta_path)
            subtype_map = dict(zip(meta["sample_id"], meta["subtype"]))
            cdf = pd.read_csv(counts_path, index_col=0, nrows=1000)  # top 1000 genes for speed
            cdf = cdf.drop(columns=["gene_name"], errors="ignore")
            # Subsample 80 samples for rendering
            tnbc_cols = [c for c in cdf.columns if subtype_map.get(c) == "TNBC"][:40]
            nontnbc_cols = [c for c in cdf.columns if subtype_map.get(c) in ("Luminal_A", "Luminal_B", "HER2_enriched")][:40]
            sel_cols = tnbc_cols + nontnbc_cols
            if len(sel_cols) >= 20:
                mat = cdf[sel_cols].T
                mat = np.log2(mat + 1)
                # PCA
                from sklearn.decomposition import PCA
                pca = PCA(n_components=2, random_state=42)
                coords = pca.fit_transform(mat.values)
                labels = ["TNBC" if c in tnbc_cols else "nonTNBC" for c in sel_cols]
                charts["pca"] = [
                    {"x": round(float(coords[i, 0]), 2), "y": round(float(coords[i, 1]), 2), "g": labels[i]}
                    for i in range(len(sel_cols))
                ]
                charts["pca_var"] = [round(float(pca.explained_variance_ratio_[0]) * 100, 1),
                                     round(float(pca.explained_variance_ratio_[1]) * 100, 1)]
                # Distance matrix — use first 40 samples, compute pairwise correlation
                sub = mat.iloc[:40]
                corr = np.corrcoef(sub.values)
                dist = 1 - corr
                charts["distance"] = {
                    "matrix": [[round(float(dist[i][j]), 2) for j in range(len(dist))] for i in range(len(dist))],
                    "labels": labels[:40],
                    "n": len(dist),
                }
        except Exception as e:
            charts["pca_error"] = str(e)

    return JSONResponse(charts)


STRUCTURE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed", "structure_predictions"
)


@app.get("/api/structure/{gene}")
async def api_structure(gene: str):
    from .web_execution import structure_info
    try:
        return structure_info(gene)
    except ValueError as exc:
        raise HTTPException(404, str(exc)) from exc


@app.get("/api/structure/{gene}/plddt.png")
async def api_structure_plddt(gene: str):
    """pLDDT 플롯 이미지 반환"""
    from fastapi.responses import FileResponse
    gene_l = gene.lower()
    img_path = os.path.join(STRUCTURE_DIR, gene_l, f"{gene_l}_boltz2_plddt.png")
    if not os.path.exists(img_path):
        return JSONResponse({"error": "not found"}, status_code=404)
    return FileResponse(img_path, media_type="image/png")


@app.get("/structure/{gene}/3d", response_class=HTMLResponse)
async def structure_3d_raw(gene: str):
    from .web_execution import asset, viewer_document
    try:
        return viewer_document(asset(gene=gene.upper()))
    except ValueError as exc:
        raise HTTPException(404, str(exc)) from exc


@app.get("/structure/{gene}", response_class=HTMLResponse)
async def structure_viewer(gene: str):
    return await structure_3d_raw(gene)


@app.get("/dashboard", response_class=HTMLResponse)
@app.get("/", response_class=HTMLResponse)
async def dashboard():
    return dashboard_document()


@app.get("/scenario1", response_class=HTMLResponse)
async def scenario1():
    from .scenario1_ui import SCENARIO1_HTML
    return SCENARIO1_HTML


@app.get("/scenario2", response_class=HTMLResponse)
async def scenario2():
    from .scenario2_ui import SCENARIO2_HTML
    return SCENARIO2_HTML


@app.post("/api/run_s2")
async def api_run_s2(body: dict = None):
    """시나리오 2 파이프라인 실행"""
    if run_state["running"]:
        return JSONResponse({"status": "already_running"})
    question = (body or {}).get("question", "TNBC 표적 발굴 + 자원 인지형 스크리닝")

    def run_s2_sync(q):
        from .orchestrator import build_graph as build_graph_s2
        graph = build_graph_s2()
        initial = {
            "research_question": q,
            "evidence_cards": [], "messages": [], "errors": [], "critic_count": 0,
        }
        run_state["running"] = True
        run_state["history"] = []
        run_state["result"] = None
        run_state["run_id"] = run_state.get("run_id", 0) + 1
        start = time.time()
        final = {}
        for event in graph.stream(initial):
            for node_name, node_output in event.items():
                display_hidden = state_hidden_genes({**final, **node_output})
                elapsed = round(time.time() - start, 1)
                entry = {
                    "type": "node", "node": node_name, "elapsed": elapsed,
                    "messages": node_output.get("messages", []),
                    "errors": node_output.get("errors", []),
                    "review_issues": node_output.get("review_issues", []),
                    "gate_status": node_output.get("gate_status"),
                    "verdict": node_output.get("critic_verdict"),
                    "regression_target": node_output.get("critic_regression_target"),
                    "n_advance": len(node_output.get("advance_targets", [])),
                    "advance_targets": [
                        {"gene": t.get("gene_name", ""), "tier": f"{t.get('bio_tier','')}{t.get('dev_tier','')}",
                         "modality": t.get("modality", "")}
                        for t in node_output.get("advance_targets", []) if visible_gene(t.get("gene_name"), display_hidden)
                    ],
                }
                run_state["history"].append(entry)
                prior_issues = final.get("review_issues", [])
                final.update(node_output)
                final["review_issues"] = list(dict.fromkeys(prior_issues + node_output.get("review_issues", [])))
        total = round(time.time() - start, 1)
        done_msg = {
            "type": "done", "total_time": total,
            "critic_count": final.get("critic_count", 0),
            "final_verdict": final.get("critic_verdict", ""),
            "node_path": " → ".join(h["node"] for h in run_state["history"]),
        }
        run_state["history"].append(done_msg)
        run_state["result"] = done_msg
        run_state["running"] = False

    t = threading.Thread(target=run_s2_sync, args=(question,), daemon=True)
    t.start()
    return JSONResponse({"status": "started", "question": question})


@app.get("/api/screening")
async def api_screening():
    """스크리닝 결과 조회"""
    path = os.path.join(RESULT_DIR, "screening_output.json")
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return JSONResponse(json.load(f))
    return JSONResponse({"error": "not found"}, status_code=404)


def dashboard_document():
    root = Path(__file__).parent
    return (root / 'dashboard.html').read_text().replace(
        '__DASHBOARD_SCRIPT__', (root / 'assets/dashboard.js').read_text()).replace(
        '__SAVED_REPORT__', os.environ.get('OMICRAFT_SAVED_REPORT', ''))



@app.get("/api/reports")
async def api_reports():
    from .report_history import report_history
    return report_history(REPORT_RUNS_ROOT, os.environ.get(
        'OMICRAFT_SAVED_REPORT', ''))


@app.get("/api/report")
async def api_report(run_id: str | None = None, run_dir: str | None = None, artifact: int | None = None):
    if run_id is not None or run_dir is not None:
        from fastapi.responses import FileResponse

        from .run_records import artifact_path, read_report
        try:
            if artifact is not None:
                path = artifact_path(REPORT_RUNS_ROOT, run_id, run_dir, artifact)
                return FileResponse(path, filename=path.name)
            return JSONResponse(read_report(REPORT_RUNS_ROOT, run_id, run_dir))
        except (ValueError, OSError, KeyError, TypeError) as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
    """파이프라인 전체 리포트 데이터 — DEG / PCA / GSEA / 타겟 통합."""
    import math
    data: dict = {}

    # ── 1. DEG (Volcano) ──
    deg_path = os.path.join(RESULT_DIR, "..", "deg_results", "deg_tnbc_vs_nontnbc_full.csv")
    if os.path.exists(deg_path):
        try:
            import pandas as pd
            df = pd.read_csv(deg_path, index_col=0).dropna(subset=["log2FoldChange", "padj"])
            df = df[df["padj"] > 0]
            df["nl10"] = df["padj"].apply(lambda p: min(-math.log10(max(p, 1e-300)), 300))
            lfc_cut, fdr_cut = 1.0, 0.05
            is_up   = (df["log2FoldChange"] >  lfc_cut) & (df["padj"] < fdr_cut)
            is_down = (df["log2FoldChange"] < -lfc_cut) & (df["padj"] < fdr_cut)
            is_ns   = ~(is_up | is_down)
            data["deg_stats"] = {"n_total": len(df), "n_up": int(is_up.sum()),
                                  "n_down": int(is_down.sum()), "lfc_cut": lfc_cut,
                                  "fdr_cut": fdr_cut}
            gcol = "gene_name" if "gene_name" in df.columns else None
            def _samp(mask, label, n=150):
                sub = df[mask].sample(min(n, int(mask.sum())), random_state=42)
                return [{"x": round(float(r["log2FoldChange"]), 2),
                         "y": round(float(r["nl10"]), 1),
                         "g": str(r[gcol]) if gcol else str(idx),
                         "s": label}
                        for idx, r in sub.iterrows()]
            data["volcano"] = _samp(is_up,"up") + _samp(is_down,"dn") + _samp(is_ns,"ns",200)
            def _top(mask, n=8):
                sub = df[mask].nlargest(n, "nl10")
                return [{"gene": str(r[gcol]) if gcol else str(i),
                         "lfc":  round(float(r["log2FoldChange"]), 2),
                         "padj": f"{r['padj']:.1e}"}
                        for i, r in sub.iterrows()]
            data["top_up"] = _top(is_up)
            data["top_dn"] = _top(is_down)
        except Exception as e:
            data["deg_error"] = str(e)

    # ── 2. PCA ──
    counts_path = os.path.join(RESULT_DIR, "..", "tcga_brca_counts.csv")
    meta_path   = os.path.join(RESULT_DIR, "..", "sample_metadata.csv")
    if os.path.exists(counts_path) and os.path.exists(meta_path):
        try:
            import numpy as np
            import pandas as pd
            from sklearn.decomposition import PCA as _PCA
            meta = pd.read_csv(meta_path)
            smap = dict(zip(meta["sample_id"], meta["subtype"]))
            cdf  = pd.read_csv(counts_path, index_col=0, nrows=1000).drop(columns=["gene_name"], errors="ignore")
            tnbc_c = [c for c in cdf.columns if smap.get(c) == "TNBC"][:40]
            non_c  = [c for c in cdf.columns if smap.get(c) in ("Luminal_A","Luminal_B","HER2_enriched")][:40]
            sel = tnbc_c + non_c
            if len(sel) >= 20:
                mat    = np.log2(cdf[sel].T.values + 1)
                pca    = _PCA(n_components=2, random_state=42)
                coords = pca.fit_transform(mat)
                labels = ["TNBC"]*len(tnbc_c) + ["nonTNBC"]*len(non_c)
                data["pca"] = [{"x": round(float(coords[i,0]),2),
                                 "y": round(float(coords[i,1]),2),
                                 "g": labels[i]} for i in range(len(sel))]
                data["pca_var"] = [round(float(pca.explained_variance_ratio_[0])*100,1),
                                   round(float(pca.explained_variance_ratio_[1])*100,1)]
        except Exception as e:
            data["pca_error"] = str(e)

    # ── 3. GSEA ──
    gsea_path = os.path.join(RESULT_DIR, "..", "pathway_results", "gsea_tnbc_vs_nontnbc_combined.csv")
    if os.path.exists(gsea_path):
        try:
            import pandas as pd
            gdf = pd.read_csv(gsea_path)
            sig = gdf[gdf["FDR q-val"] < 0.25].sort_values("NES")
            up15  = sig[sig["NES"] > 0].nlargest(8, "NES")
            dn15  = sig[sig["NES"] < 0].nsmallest(8, "NES")
            rows  = pd.concat([dn15, up15]).sort_values("NES")
            data["gsea"] = [{"term": r["Term"][:55],
                              "nes":  round(float(r["NES"]), 2),
                              "fdr":  round(float(r["FDR q-val"]), 4)}
                             for _, r in rows.iterrows()]
        except Exception as e:
            data["gsea_error"] = str(e)

    # ── 4. ADVANCE targets (dossier + qualification_full.csv + cell_context.tsv) ──
    dossier_path = os.path.join(RESULT_DIR, "final_dossier.json")
    if os.path.exists(dossier_path):
        try:
            import csv as _csv
            with open(dossier_path, encoding="utf-8") as f:
                dossier = json.load(f)

            # Load enrichment: qualification_full.csv
            qual_map: dict = {}
            qual_csv = os.path.join(RESULT_DIR, "..", "qualification_results", "qualification_full.csv")
            if os.path.exists(qual_csv):
                with open(qual_csv, encoding="utf-8") as fq:
                    for row in _csv.DictReader(fq):
                        qual_map[row["gene_name"]] = row

            # Load enrichment: cell_context.tsv (best specificity per gene)
            ctx_map: dict = {}
            ctx_tsv = os.path.join(RESULT_DIR, "..", "cell_context", "cell_context.tsv")
            if os.path.exists(ctx_tsv):
                with open(ctx_tsv, encoding="utf-8") as fc:
                    for row in _csv.DictReader(fc, delimiter="\t"):
                        gene = row["gene"]
                        if gene not in ctx_map or float(row.get("specificity") or 0) > float(ctx_map[gene].get("specificity") or 0):
                            ctx_map[gene] = row

            advance = dossier.get("advance_targets", [])
            enriched = []
            for t in advance:
                gene = t.get("gene_name", "")
                if not visible_gene(gene):
                    continue
                q    = qual_map.get(gene, {})
                ctx  = ctx_map.get(gene, {})

                # subcellular location → list from main_location (semicolon-separated)
                main_loc = q.get("main_location") or ""
                loc_list = [location.strip() for location in main_loc.split(";") if location.strip()]

                enriched.append({
                    "gene":      gene,
                    "tier":      (q.get("bio_tier") or "") + (q.get("dev_tier") or "") or t.get("bio_tier", "") + t.get("dev_tier", ""),
                    "bio_conf":  q.get("bio_tier") or t.get("bio_confidence") or "",
                    "dev_tier":  q.get("dev_tier") or t.get("dev_tier") or "",
                    "modality":  q.get("modality") or t.get("modality") or "",
                    "location":  loc_list,
                    "subcell":   q.get("subcellular_location") or "",
                    "depmap":    q.get("depmap_status") or "",
                    "chronos":   q.get("depmap_pct_dep") or "",
                    "cell_ctx":  ctx.get("cell_class") or "",
                    "cell_spec": ctx.get("specificity") or "",
                    "direction": t.get("direction") or "",
                    "safety":    t.get("safety_verdict") or q.get("safety_verdict") or "",
                    "role":      t.get("intervention_role") or "",
                })
            data["advance_targets"] = enriched[:25]
            data["research_question"] = dossier.get("research_question", "")
            data["final_verdict"]     = dossier.get("critic_verdict", "")
        except Exception as e:
            data["dossier_error"] = str(e)

    data["provenance"] = {
        "kind": "legacy_tcga", "deg_source": deg_path, "gsea_source": gsea_path,
        "pca_source": counts_path, "metadata_source": meta_path,
        "volcano_sampling": "up <=150, down <=150, nonsignificant <=200; random_state=42",
        "pca_sampling": "first 1000 file rows, <=40 samples/group; log2(count+1), no library-size normalization",
        "pca_computed": "pca" in data,
        "deg_available": "deg_stats" in data, "gsea_available": "gsea" in data,
    }
    return JSONResponse(data)



RUN_REPORT_SCRIPT = (Path(__file__).parent / 'templates/run-report-script.html').read_text(encoding='utf-8')

REPORT_HTML = (Path(__file__).parent / 'templates/report.html').read_text(encoding='utf-8')


def report_document(html, selected_run=False):
    if selected_run:
        # The run view uses native DOM tables and does not need the TCGA chart library.
        html = re.sub(r'<script\b[^>]*\bsrc=[^>]*>\s*</script>', '', html, flags=re.IGNORECASE)
        html = re.sub(r'<link\b[^>]*https?://[^>]*>', '', html, flags=re.IGNORECASE)
        html = re.sub(r'@import\s+url\([^;]+;', '', html, flags=re.IGNORECASE)
    from .report_theme import style_report
    assets = Path(__file__).parent / 'assets'
    integrated = '<style>' + (assets / 'integrated-report.css').read_text() + '</style>'
    integrated += '<script>' + (assets / 'integrated-report.js').read_text() + '</script>'
    integrated += '<script>' + (assets / 'pipeline-design-report.js').read_text() + '</script>'
    return style_report(html.replace("</head>", integrated + RUN_REPORT_SCRIPT + "</head>"))


def pipeline_report_document(integrated=False):
    """Scenario 1 layout with recorded data; no demo data or external JS required."""
    base = Path(__file__).parent
    assets = base / 'assets'
    styles = ''.join('<style>' + (assets / name).read_text() + '</style>' for name in
                     ('integrated-report.css', 'pipeline-scenario.css'))
    styles = re.sub(r'@import\s+url\([^;]+;', '', styles)
    scripts = ''.join('<script>' + (assets / name).read_text() + '</script>' for name in
                      ('integrated-report.js', 'pipeline-design-report.js', 'pipeline-scenario.js'))
    html = (base / 'pipeline_scenario.html').read_text()
    if integrated:
        sections = [('question', '연구질문 · 입력 조건'), ('overview', '모달리티 실행 요약'),
                    ('comparison', '모달리티 비교'), ('design', '모달리티별 상세 결과'),
                    ('validation', '검증 범위 · 후속 실험'), ('files', '실행 기록')]
        nav = ''.join(f'<a class="nv" href="#s-{key}"><span class="mono">{i:02d}</span>{title}</a>'
                      for i, (key, title) in enumerate(sections, 1))
        content = ''.join(f'<section id="s-{key}"><div class="shead"><span class="num">{i:02d}</span>'
                          f'<h2>{title}</h2></div><div id="pr-{key}"></div></section>'
                          for i, (key, title) in enumerate(sections, 1))
        html = re.sub(r'<nav class="nav".*?</nav>', '<nav class="nav" aria-label="보고서 목차">' + nav + '</nav>', html)
        html = re.sub(r'<div class="content">.*?</div><footer', '<div class="content">' + content + '</div><footer', html)
        html = html.replace('renderPipelineScenario(d)', 'renderDesignScenario(d)')
        scripts += '<script>' + (assets / 'design-scenario.js').read_text() + '</script>'
    return html.replace('<!--PIPELINE_ASSETS-->', styles + scripts)


@app.get("/report", response_class=HTMLResponse)
async def report_page(run_id: str | None = None, run_dir: str | None = None):
    from .run_records import run_directory
    if run_id is None and run_dir is None:
        # The default report opens the latest recorded pipeline with analysis results.
        for path in sorted(Path(REPORT_RUNS_ROOT).glob('upstream_*/run_summary.json'), reverse=True):
            if not path.resolve().is_relative_to(Path(REPORT_RUNS_ROOT).resolve()):
                continue
            try:
                saved = json.loads(path.read_text())
                if saved.get('execution_profile') == 'upstream_analysis' and saved.get('summary', {}).get('candidates'):
                    from fastapi.responses import RedirectResponse
                    return RedirectResponse('/report?run_id=' + path.parent.name)
            except (OSError, ValueError, TypeError, AttributeError):
                continue
    else:
        try:
            directory = run_directory(REPORT_RUNS_ROOT, run_id, run_dir)
            summary_path = directory / 'run_summary.json'
            if summary_path.resolve().is_relative_to(Path(REPORT_RUNS_ROOT).resolve()):
                saved = json.loads(summary_path.read_text())
                if saved.get('execution_profile') == 'upstream_analysis':
                    return pipeline_report_document()
                if saved.get('execution_profile') == 'therapeutic_design' and saved.get('modality') == 'ALL':
                    return pipeline_report_document(integrated=True)
        except (OSError, ValueError, TypeError, AttributeError):
            pass
    return report_document(REPORT_HTML, run_id is not None or run_dir is not None)


@app.get("/demo", response_class=HTMLResponse)
async def demo_page(run_id: str | None = None, run_dir: str | None = None):
    html_path = os.path.join(os.path.dirname(__file__), "demo_report.html")
    if os.path.exists(html_path):
        return HTMLResponse(report_document(Path(html_path).read_text(), run_id is not None or run_dir is not None))
    return HTMLResponse("<body>demo_report.html not found</body>", status_code=404)


# ─── Structure Lab (ADC / SM / Binder / Degrader) ──────────────────────────────────────

LAB_OUTPUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "dataset", "tcga_brca", "processed", "lab_docking"
)
PDB_CACHE_DIR = os.path.join(LAB_OUTPUT_DIR, "pdb_cache")

_lab_executor = ThreadPoolExecutor(max_workers=2)
_lab_sessions: dict = {}  # sid -> {"status": "running"|"done"|"error", "result": dict}


def _get_lab():
    """lab_docking 모듈 지연 임포트."""
    _dir = os.path.dirname(__file__)
    if _dir not in sys.path:
        sys.path.insert(0, _dir)
    from . import lab_docking as _ld
    return _ld


def web_request(function, *args):
    try:
        return function(*args)
    except (ValueError, KeyError, TypeError, FileNotFoundError) as exc:
        raise HTTPException(422, str(exc)) from exc


@app.get("/assets/3Dmol-min.js")
def viewer_library():
    from fastapi.responses import FileResponse
    return FileResponse(Path(__file__).parent / 'assets/3Dmol-min.js', media_type='application/javascript')


@app.get("/api/structures")
def api_structures():
    from .web_execution import catalog
    return catalog()


@app.post("/api/structures")
def api_upload_structure(body: dict):
    from .web_execution import upload
    return web_request(upload, body)


@app.get("/api/structures/{identifier}/viewer", response_class=HTMLResponse)
def api_asset_viewer(identifier: str):
    from .web_execution import asset, viewer_document
    return web_request(lambda: viewer_document(asset(identifier)))


@app.get("/api/design/candidates")
def api_design_candidates():
    from .web_execution import design_candidates
    return design_candidates()


@app.post("/api/design")
def api_design(body: dict):
    from .web_execution import start_design
    return web_request(start_design, body)


@app.post("/api/lab/{mode}")
def api_lab(mode: str, body: dict):
    from .web_execution import start_lab
    if mode not in {'adc', 'sm', 'protac'}:
        raise HTTPException(404, 'Unknown Lab mode')
    return web_request(start_lab, mode, body)


@app.get("/api/lab/status/{sid}")
def api_lab_status(sid: str):
    from .web_execution import lab_status
    return web_request(lab_status, sid)


@app.get("/api/lab/viewer/{sid}", response_class=HTMLResponse)
def api_lab_viewer(sid: str):
    from .web_execution import lab_directory, lab_status
    status = web_request(lab_status, sid)
    if status['status'] != 'done':
        raise HTTPException(409, status.get('error', 'Still running'))
    return (lab_directory(sid) / 'viewer.html').read_text()


LAB_HTML = (Path(__file__).parent / 'templates/lab.html').read_text(encoding='utf-8')


@app.get("/lab", response_class=HTMLResponse)
async def lab_page():
    page = LAB_HTML.replace('<div>\n      <h2>Status</h2>', Path(__file__).with_name('lab_design.html').read_text(encoding='utf-8') + '<div>\n      <h2>Status</h2>')
    return page


if __name__ == "__main__":
    import uvicorn
    print("=" * 60)
    print("OmiCraft Agent Dashboard")
    print("http://localhost:8000")
    print("=" * 60)
    uvicorn.run(app, host="0.0.0.0", port=8000, log_level="info")


def write_competition_report(data, directory):
    """Export the existing report renderer with local artifact links and embedded data."""
    import copy

    directory = Path(directory).resolve()
    payload = copy.deepcopy(data)
    from .report_figures import label_figure
    for artifact in payload.get("artifacts", []):
        if isinstance(artifact, dict):
            label_figure(artifact)
    if payload.get('modality') == 'ALL' and payload.get('execution_profile') == 'therapeutic_design':
        from .integrated_report import assemble_report
        payload = assemble_report(payload, directory.parent)
    elif payload.get('modality') == 'DEGRADER':
        from .integrated_report import assemble_degrader_report
        payload = assemble_degrader_report(payload, directory.parent)
    elif payload.get('modality') == 'DE_NOVO_BINDER':
        from .integrated_report import assemble_binder_report
        payload = assemble_binder_report(payload, directory.parent)
    elif payload.get('execution_profile') == 'upstream_analysis':
        from .pipeline_report import assemble_pipeline_report
        payload = assemble_pipeline_report(payload, directory.parent)
    payload['kind'] = 'therapeutic'
    design = payload.get('design_report', {})
    reference = design.get('reference') or {}
    panels = (payload.get('integrated_report', {}).get('panels', [])
              + design.get('panels', []) + design.get('followups', []) + reference.get('panels', [])
              + ([payload['binder_report']] if payload.get('binder_report') else [])
              + ([payload['degrader_report']] if payload.get('degrader_report') else []))
    for panel in panels:
        for item in panel['artifacts'] + panel['structures'] + panel['figures']:
            path = Path(item['path']).resolve()
            item.update(available=path.is_file(), href=os.path.relpath(path, directory))
        if panel.get('run_id'):
            panel['report_url'] = '../' + panel['run_id'] + '/report.html'
    for candidate in payload.get('integrated_report', {}).get('passing_candidates', []):
        candidate['report_url'] = '../' + candidate['run_id'] + '/report.html'
    if reference:
        reference['report_url'] = '../' + reference['run_id'] + '/report.html'
    if payload.get('integrated_report'):
        payload['integrated_report']['previous_report_url'] = '../visual_diagnostic_20260922T152713Z/report.html'
    if panels:
        import shutil
        viewer = Path(__file__).parent / 'assets/3Dmol-min.js'
        if viewer.is_file():
            shutil.copy2(viewer, directory / '3Dmol-min.js')
    for item in payload['artifacts']:
        path = Path(item['path']).resolve()
        item.update(available=path.is_file(), href=os.path.relpath(path, directory))
    serialized = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    html = (pipeline_report_document() if payload.get('execution_profile') == 'upstream_analysis'
            else pipeline_report_document(integrated=True) if payload.get('integrated_report')
            else report_document(REPORT_HTML, selected_run=True))
    html = html.replace("await fetch('/api/report' + location.search).then(r=>r.json())", 'window.__COMPETITION_RUN__')
    html = html.replace('</head>', '<script>window.__COMPETITION_RUN__=' + serialized + ';</script></head>')
    path = directory / 'report.html'
    path.write_text(html, encoding='utf-8')
    return path


# Independent read-only workflow monitor endpoints.
from .workflow_api import register as _register_workflow
_register_workflow(app, sys.modules[__name__])

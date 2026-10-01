"""Workflow topology and runtime delivery, separate from scientific APIs."""
import asyncio
import copy
import json
import re

from fastapi import HTTPException, Request, WebSocket, WebSocketDisconnect
from fastapi.responses import RedirectResponse, StreamingResponse


def register(app, server):
    def snapshot(identifier=None):
        with server.run_lock:
            current = server.run_state.get('workflow')
            if current and (not identifier or current['run_id']==identifier):
                return copy.deepcopy(current)
        if identifier:
            if not re.fullmatch(r'[A-Za-z0-9_-]{1,120}', identifier):
                raise HTTPException(422, 'Invalid workflow run ID')
            path = server.REPORT_RUNS_ROOT / identifier / 'workflow.json'
            if not path.is_file() or not path.resolve().is_relative_to(server.REPORT_RUNS_ROOT):
                raise HTTPException(404, 'Workflow run not found')
            return json.loads(path.read_text())
        return {'run_id': None, 'version': 0, 'status': 'IDLE', 'nodes': {}, 'events': []}

    @app.get('/workflow')
    async def page():
        return RedirectResponse('/dashboard', status_code=302)

    @app.get('/api/workflow/graph')
    async def definition():
        from .workflow_monitor import graph_definition
        return graph_definition(server.build_graph())

    @app.get('/api/workflow/runs')
    async def runs():
        rows = []
        for path in sorted(server.REPORT_RUNS_ROOT.glob('*/workflow.json'), key=lambda p:p.stat().st_mtime, reverse=True)[:100]:
            if not path.resolve().is_relative_to(server.REPORT_RUNS_ROOT):
                continue
            try:
                data=json.loads(path.read_text())
                rows.append({k:data.get(k) for k in ('run_id','status','started_at','ended_at')})
            except (OSError, ValueError):
                continue
        return rows

    @app.get('/api/workflow/state')
    async def state(run_id: str = None):
        return snapshot(run_id)

    @app.get('/api/workflow/events')
    async def events(request: Request, run_id: str = None):
        snapshot(run_id)  # validate before starting response
        async def stream():
            previous = None
            heartbeat = 0
            while not await request.is_disconnected():
                data = snapshot(run_id)
                key=(data['run_id'],data['version'])
                if key != previous:
                    yield 'event: snapshot\ndata: '+json.dumps(data,ensure_ascii=False)+'\n\n'
                    previous = key
                heartbeat += 1
                if heartbeat % 50 == 0:
                    yield ': heartbeat\n\n'
                await asyncio.sleep(.15)
        return StreamingResponse(stream(),media_type='text/event-stream',headers={'Cache-Control':'no-cache, no-transform','X-Accel-Buffering':'no'})

    @app.websocket('/ws/workflow')
    async def websocket(ws: WebSocket):
        await ws.accept()
        identifier=ws.query_params.get('run_id')
        previous=None
        ticks=0
        try:
            while True:
                data=snapshot(identifier); key=(data['run_id'],data['version'])
                ticks += 1
                if key != previous or ticks % 50 == 0:
                    await ws.send_json(data);previous=key
                await asyncio.sleep(.15)
        except (WebSocketDisconnect, RuntimeError):
            pass
        except HTTPException:
            await ws.close(code=1008)

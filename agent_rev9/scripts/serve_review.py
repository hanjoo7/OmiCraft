"""Serve the dashboard and one approved report with one active execution at a time."""
import json
import os
import re
import threading
import socket
import select
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlencode, urlsplit

RUN = Path(os.environ.get('OMICRAFT_PUBLIC_RUN', Path.cwd() / 'runs/visual_diagnostic_20260922T152713Z')).resolve()
UPSTREAM = os.environ.get('OMICRAFT_UPSTREAM', 'http://127.0.0.1:8088')
ARTIFACT_COUNT = (len(json.loads((RUN / 'run_summary.json').read_text()).get('artifacts', []))
                  if (RUN / 'run_summary.json').is_file() else 0)
POLICY_PATH = RUN / 'cloudflare_access/execution_policy.json'
RUN_REQUEST_LOCK = threading.Lock()
DASHBOARD_PATH = '/omicraft-version-9'
EXECUTION_ROUTES = {'/api/run', '/api/design', '/api/lab/adc', '/api/lab/sm', '/api/lab/protac'}


def approved_report_count(name):
    if name == RUN.name:
        return ARTIFACT_COUNT
    if not re.fullmatch(r"(?:upstream|design|lab)_\d{8}T\d{6}Z_[0-9a-f]{8}", name):
        return None
    summary = RUN.parent / name / "run_summary.json"
    if not summary.is_file() or not summary.resolve().is_relative_to(RUN.parent.resolve()):
        return None
    data = json.loads(summary.read_text())
    if data.get("execution_profile") not in {"upstream_analysis", "therapeutic_design", "structure_lab"}:
        return None
    return len(data.get("artifacts", []))

def default_pipeline_report():
    """Use a recorded pipeline with results for the report landing page."""
    for path in sorted(RUN.parent.glob('upstream_*/run_summary.json'), reverse=True):
        try:
            if approved_report_count(path.parent.name) is None:
                continue
            data = json.loads(path.read_text())
            if data.get('execution_profile') == 'upstream_analysis' and data.get('summary', {}).get('candidates'):
                return path.parent.name
        except (OSError, ValueError, TypeError, AttributeError):
            continue
    return RUN.name


def execution_policy():
    return json.loads(POLICY_PATH.read_text()) if POLICY_PATH.is_file() else {'enabled': False}



def dashboard_view(content):
    return content


class Handler(BaseHTTPRequestHandler):
    server_version = 'OmiCraftReview'

    def do_GET(self):
        self.respond(False)

    def do_HEAD(self):
        self.respond(True)

    def json_response(self, status, value, head=False):
        content = json.dumps(value).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(content)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        if not head:
            self.wfile.write(content)

    def execution_enabled(self):
        if not execution_policy().get('enabled'):
            self.json_response(503, {'detail': 'Execution access is not enabled'})
            return False
        return True

    def do_POST(self):
        request = urlsplit(self.path)
        routes = {'/api/stop', '/api/run', '/api/design', '/api/structures', '/api/lab/adc', '/api/lab/sm', '/api/lab/protac'}
        if request.query or request.path not in routes:
            self.json_response(404, {'detail': 'Endpoint not exposed'})
            return
        if not self.execution_enabled():
            return
        try:
            size = int(self.headers.get('Content-Length', '0'))
            limit = 6000000 if request.path == '/api/structures' else 65536
            if not 0 < size <= limit:
                self.json_response(413, {'detail': 'Request body too large or empty'})
                return
            body = json.loads(self.rfile.read(size))
            policy = execution_policy()
            if not isinstance(body, dict):
                raise ValueError('Expected object')
            if request.path == '/api/run':
                if not ('question' in body or 'scenario' in body) or set(body) - {'scenario', 'question', 'design_mode', 'structural_backend', 'reuse_results'}:
                    raise ValueError('Public pipeline accepts scenario or question, design_mode, structural_backend and reuse_results')
                if type(body.get('reuse_results', True)) is not bool:
                    raise ValueError('reuse_results must be boolean')
                if body.get('design_mode', 'none') not in {'none', 'all_eligible'}:
                    raise ValueError('Unknown design_mode')
                if body.get('structural_backend', 'af3') not in {'af3', 'boltz2_legacy'}:
                    raise ValueError('Unknown structural_backend')
                if 'scenario' in body:
                    if body['scenario'] not in ('tnbc_discovery', 'erbb2_design'):
                        raise ValueError('Unknown research scenario')
                else:
                    question = body['question']
                    if not isinstance(question, str) or not question.strip() or len(question) > policy['max_question_characters']:
                        raise ValueError('Enter a question of 1–10000 characters')
                    body['question'] = question.strip()
            with RUN_REQUEST_LOCK:
                if request.path in EXECUTION_ROUTES:
                    with urllib.request.urlopen(UPSTREAM + '/api/status', timeout=5) as response:
                        status = json.load(response)
                    if status['running']:
                        self.json_response(409, {'detail': 'Another task is running. Wait for it to finish.',
                                                 'status': 'already_running', 'run_id': status.get('run_id')})
                        return
                upstream = urllib.request.Request(UPSTREAM + request.path, data=json.dumps(body).encode(),
                    headers={'Content-Type': 'application/json'}, method='POST')
                with urllib.request.urlopen(upstream, timeout=30) as response:
                    data = json.load(response)
                    self.json_response(response.status, data)
        except urllib.error.HTTPError as exc:
            try:
                data = json.load(exc)
            except (ValueError, OSError):
                data = {'detail': 'Request rejected by upstream'}
            self.json_response(exc.code, data)
        except (ValueError, TypeError) as exc:
            self.json_response(422, {'detail': str(exc)})
        except (OSError, TimeoutError):
            self.json_response(502, {'detail': 'Pipeline server unavailable'})

    def respond(self, head):
        request = urlsplit(self.path)
        if request.path in {'/workflow', '/workflow/'}:
            self.send_response(302)
            self.send_header('Location', DASHBOARD_PATH)
            self.send_header('Cache-Control', 'no-store')
            self.send_header('Content-Length', '0')
            self.end_headers()
            return
        if request.path == '/ws/workflow':
            params=parse_qs(request.query, keep_blank_values=True)
            if (set(params)-{'run_id'} or any(len(v)!=1 or not re.fullmatch(r'[A-Za-z0-9_-]{1,120}',v[0]) for v in params.values())
                or self.headers.get('Upgrade','').lower()!='websocket'
                or self.headers.get('Sec-WebSocket-Version')!='13'
                or not re.fullmatch(r'[A-Za-z0-9+/=]{20,32}',self.headers.get('Sec-WebSocket-Key',''))):
                self.json_response(400, {'detail':'Invalid workflow WebSocket request'},head);return
            target=urlsplit(UPSTREAM)
            self.close_connection=True
            try:
                with socket.create_connection((target.hostname,target.port or 80),timeout=10) as upstream:
                    headers=(f'GET {self.path} HTTP/1.1\r\nHost: {target.netloc}\r\n'
                             'Upgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Version: 13\r\n'
                             f'Sec-WebSocket-Key: {self.headers["Sec-WebSocket-Key"]}\r\n\r\n')
                    upstream.sendall(headers.encode('ascii'))
                    response=b''
                    while b'\r\n\r\n' not in response and len(response)<65536:
                        part=upstream.recv(4096)
                        if not part:break
                        response+=part
                    self.connection.sendall(response)
                    if not response.startswith(b'HTTP/1.1 101 '):return
                    self.log_request(101)
                    upstream.settimeout(None)
                    while True:
                        readable,_,_=select.select([self.connection,upstream],[],[],20)
                        for source in readable:
                            chunk=source.recv(65536)
                            if not chunk:return
                            (upstream if source is self.connection else self.connection).sendall(chunk)
            except (OSError,ValueError):
                pass
            return
        if request.path == '/api/workflow/events':
            params = parse_qs(request.query, keep_blank_values=True)
            if set(params) - {'run_id'} or any(len(v)!=1 or not re.fullmatch(r'[A-Za-z0-9_-]{1,120}',v[0]) for v in params.values()):
                self.json_response(422, {'detail': 'Invalid workflow run'}, head);return
            try:
                with urllib.request.urlopen(UPSTREAM + self.path, timeout=30) as response:
                    self.send_response(200)
                    self.send_header('Content-Type','text/event-stream')
                    self.send_header('Cache-Control','no-cache, no-transform')
                    self.send_header('X-Accel-Buffering','no')
                    self.end_headers()
                    if head:return
                    while True:
                        line=response.readline()
                        if not line:break
                        self.wfile.write(line)
                        self.wfile.flush()
            except urllib.error.HTTPError as exc:
                self.json_response(exc.code, {'detail':'Workflow run unavailable'})
            except (OSError,TimeoutError):
                pass
            return
        if request.path in {'/api/workflow/graph','/api/workflow/runs','/api/workflow/state'}:
            params=parse_qs(request.query, keep_blank_values=True)
            if (set(params)-{'run_id'} or (request.path!='/api/workflow/state' and params)
                or any(len(v)!=1 or not re.fullmatch(r'[A-Za-z0-9_-]{1,120}',v[0]) for v in params.values())):
                self.json_response(422, {'detail':'Invalid workflow query'}, head);return
            try:
                with urllib.request.urlopen(UPSTREAM+self.path,timeout=30) as response:
                    content=response.read();self.send_response(response.status)
                    self.send_header('Content-Type',response.headers.get('Content-Type','application/json'))
                    self.send_header('Content-Length',str(len(content)));self.send_header('Cache-Control','no-store');self.end_headers()
                    if not head:self.wfile.write(content)
            except urllib.error.HTTPError as exc:self.json_response(exc.code, {'detail':'Workflow unavailable'}, head)
            except OSError:self.json_response(502, {'detail':'Workflow server unavailable'}, head)
            return
        dashboard = request.path in {'/', '/dashboard', '/omicraft-version-7-2', '/omicraft-version-7-2/', DASHBOARD_PATH, DASHBOARD_PATH + '/'} and not request.query
        if request.path in {'/report', '/demo'} and not request.query:
            self.send_response(302)
            self.send_header('Location', request.path+'?'+urlencode({'run_id':default_pipeline_report() if request.path == '/report' else RUN.name}))
            self.end_headers()
            return
        if request.path == '/favicon.ico':
            self.send_response(204)
            self.end_headers()
            return
        if request.path in {'/api/events', '/api/status'}:
            if not self.execution_enabled():
                return
            params = parse_qs(request.query, keep_blank_values=True)
            if (set(params) - {'run_id', 'cursor'} or request.path == '/api/status' and params
                    or any(len(v) != 1 or not v[0].isdecimal() for v in params.values())):
                self.json_response(422, {'detail': 'Invalid event cursor'}, head)
                return
            route = request.path + ('?' + urlencode({k: v[0] for k, v in params.items()}) if params else '')
        elif request.path == '/lab':
            params = parse_qs(request.query, keep_blank_values=True)
            if set(params) - {'parent', 'gene', 'modality'} or any(len(v) != 1 or len(v[0]) > 100 for v in params.values()):
                self.json_response(422, {'detail': 'Invalid Lab query'}, head)
                return
            route = '/lab'
        elif not request.query and (request.path in {'/api/reports', '/api/structures', '/api/design/candidates', '/api/resources', '/assets/3Dmol-min.js'}
                or re.fullmatch(r'/structure/[A-Za-z0-9_-]{1,40}(?:/3d)?', request.path)
                or re.fullmatch(r'/api/structure/[A-Za-z0-9_-]{1,40}', request.path)
                or re.fullmatch(r'/api/structures/[A-Za-z0-9_-]{1,80}/viewer', request.path)
                or re.fullmatch(r'/api/lab/(?:status|viewer)/[0-9a-f]{32}', request.path)):
            route = request.path
        elif dashboard:
            route = '/'
        elif request.path == '/health' and not request.query:
            route = '/health'
        else:
            params = parse_qs(request.query, keep_blank_values=True)
            selectors = set(params) & {'run_id','run_dir'}
            valid = (request.path in {'/report','/demo','/api/report'} and len(selectors)==1
                and not (set(params)-{'run_id','run_dir','artifact'})
                and all(len(v)==1 for v in params.values())
                and approved_report_count(params[next(iter(selectors))][0]) is not None)
            if 'artifact' in params:
                value = params['artifact'][0]
                valid = valid and request.path == '/api/report' and value.isdecimal() and 0 <= int(value) < (approved_report_count(params[next(iter(selectors))][0]) or 0)
            if not valid:
                self.send_error(404, 'Report resource not found')
                return
            route = request.path+'?'+urlencode({k:v[0] for k,v in params.items()})
        try:
            with urllib.request.urlopen(UPSTREAM+route, timeout=30) as response:
                content = response.read()
                if dashboard:
                    content = dashboard_view(content)
                self.send_response(response.status)
                for name in ['Content-Type','Content-Disposition']:
                    if response.headers.get(name):self.send_header(name,response.headers[name])
                self.send_header('Content-Length',str(len(content)))
                self.send_header('X-Content-Type-Options','nosniff')
                self.send_header('Cache-Control','no-store')
                self.end_headers()
                if not head:self.wfile.write(content)
        except urllib.error.HTTPError as exc:
            self.send_error(exc.code)
        except (OSError,TimeoutError,ValueError):
            self.send_error(502,'Review server unavailable')

if __name__=='__main__':
    ThreadingHTTPServer(('127.0.0.1',int(os.environ.get('OMICRAFT_PROXY_PORT', '8089'))),Handler).serve_forever()

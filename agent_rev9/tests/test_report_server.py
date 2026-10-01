"""Opt-in HTTP and browser smoke; all artifacts live under tmp_path."""

import base64
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import pytest

from ..modality_dispatch import result
from ..run_records import RunRecord


class Firefox:
    def __init__(self, port):
        self.socket = socket.create_connection(("127.0.0.1", port), timeout=20)
        self.socket.settimeout(20)
        self.sequence = 0
        self.receive()
        self.call("WebDriver:NewSession", {"capabilities": {"alwaysMatch": {"pageLoadStrategy": "none"}}})

    def receive(self):
        length = b""
        while not length.endswith(b":"):
            chunk = self.socket.recv(1)
            if not chunk:
                raise RuntimeError("Firefox connection closed")
            length += chunk
        length = int(length[:-1])
        payload = b""
        while len(payload) < length:
            chunk = self.socket.recv(length - len(payload))
            if not chunk:
                raise RuntimeError("Firefox response incomplete")
            payload += chunk
        return json.loads(payload)

    def call(self, command, params=None):
        self.sequence += 1
        payload = json.dumps([0, self.sequence, command, params or {}]).encode()
        self.socket.sendall(str(len(payload)).encode() + b":" + payload)
        response = self.receive()
        if response[2]:
            raise RuntimeError(response[2])
        return response[3]

    def script(self, script):
        return self.call(
            "WebDriver:ExecuteScript",
            {"script": script, "args": [], "newSandbox": True, "sandbox": "default"},
        )["value"]

    def wait(self, script):
        deadline = time.monotonic() + 15
        while time.monotonic() < deadline:
            value = self.script(script)
            if value:
                return value
            time.sleep(0.1)
        raise AssertionError("Browser condition timed out: " + script)

    def screenshot(self, path):
        value = self.call("WebDriver:TakeScreenshot", {"id": None, "full": False, "scroll": False})[
            "value"
        ]
        path.write_bytes(base64.b64decode(value))

    def close(self):
        self.call("WebDriver:DeleteSession")
        self.socket.close()



def free_port():
    with socket.socket() as sock:
        sock.bind(('127.0.0.1', 0))
        return sock.getsockname()[1]


@pytest.mark.external_tool
def test_report_server_http_and_browser(tmp_path):
    if not shutil.which('firefox'):
        pytest.skip('Firefox is not installed')
    source = Path(__file__).resolve().parents[2]
    record = RunRecord(tmp_path / 'runs/check', 'ADC')
    artifact = tmp_path / 'runs/original.csv'
    artifact.write_text('metric,value\ncontact,1\n')
    record.finish({**result('ADC', 'COMPLETED', 'cdr_mapping',
                           artifacts=[{'path': str(artifact), 'label': 'Source table'}], is_mock=True),
                   'validation_decision': 'HOLD', 'critic': {'decision': 'HOLD'}})
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONPATH': str(source),
           'OMICRAFT_RUNS_ROOT': str(tmp_path / 'runs')}
    port = free_port()
    profile = tmp_path / 'firefox'
    profile.mkdir()
    browser = None
    firefox = None
    with (tmp_path / 'server.log').open('w') as log, (tmp_path / 'browser.log').open('w') as browser_log:
        process = subprocess.Popen([sys.executable, '-m', 'uvicorn', 'agent_rev9.server:app',
            '--host', '127.0.0.1', '--port', str(port)], env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            base = f'http://127.0.0.1:{port}'
            deadline = time.monotonic() + 20
            while True:
                try:
                    with urllib.request.urlopen(base + '/api/report?run_id=check', timeout=2) as response:
                        data = json.load(response)
                    break
                except OSError:
                    if time.monotonic() > deadline:
                        pytest.fail('server failed to start: ' + (tmp_path / 'server.log').read_text())
                    time.sleep(.1)
            assert data['execution_status'] == 'COMPLETED' and data['validation_decision'] == 'HOLD'
            for path in ['/api/report', '/report', '/demo', '/scenario1', '/scenario2', data['artifacts'][0]['href']]:
                with urllib.request.urlopen(base + path, timeout=10) as response:
                    assert response.status == 200
            browser_port = free_port()
            (profile / 'user.js').write_text(f'user_pref("marionette.port",{browser_port});\n')
            firefox = subprocess.Popen(['firefox', '--headless', '--no-remote', '--marionette',
                '--profile', str(profile)], stdout=browser_log, stderr=subprocess.STDOUT)
            deadline = time.monotonic() + 25
            while browser is None:
                try:
                    browser = Firefox(browser_port)
                except OSError:
                    if time.monotonic() > deadline:
                        raise
                    time.sleep(.1)
            browser.call('WebDriver:SetTimeouts', {'pageLoad': 25000})
            for page in ['report', 'demo']:
                browser.call('WebDriver:Navigate', {'url': base + '/' + page + '?run_id=check'})
                browser.wait("return document.body.textContent.includes('Therapeutic run: check')")
                browser.wait("return document.body.textContent.includes('Validation: HOLD')")
                browser.wait("return [...document.querySelectorAll('a')].some(a=>a.textContent==='Source table' && a.href.includes('artifact=0'))")
                browser.screenshot(tmp_path / (page + '.png'))
            (tmp_path / 'verification.json').write_text(json.dumps({'http': 'passed', 'report': 'passed',
                'demo': 'passed', 'artifact': 'passed', 'is_mock': True}))
        finally:
            if browser:
                try:
                    browser.call('Marionette:Quit', {'flags': ['eForceQuit']})
                except (OSError, RuntimeError):
                    pass
            for child in [firefox, process]:
                if child and child.poll() is None:
                    child.terminate()
                    child.wait(timeout=15)

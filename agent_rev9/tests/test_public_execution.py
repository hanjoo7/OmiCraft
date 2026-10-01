import pytest
import concurrent.futures
import importlib.util
import io
import json
import tempfile
import threading
import unittest
from http.client import HTTPConnection
from pathlib import Path
from unittest.mock import patch


@pytest.mark.local_io
class PublicExecutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()
        run = Path(cls.folder.name)
        (run / 'run_summary.json').write_text('{"artifacts": []}')
        (run / 'cloudflare_access').mkdir()
        (run / 'cloudflare_access/execution_policy.json').write_text(
            '{"enabled": true, "max_question_characters": 10000}')
        source = Path(__file__).resolve().parents[1] / 'scripts/serve_review.py'
        spec = importlib.util.spec_from_file_location('public_execution_under_test', source)
        cls.proxy = importlib.util.module_from_spec(spec)
        with patch.dict('os.environ', OMICRAFT_PUBLIC_RUN=str(run)):
            spec.loader.exec_module(cls.proxy)

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def setUp(self):
        self.running = False
        self.starts = 0
        self.forwarded = []
        self.unavailable = False
        def upstream(request, **kwargs):
            method = request.get_method() if hasattr(request, 'get_method') else 'GET'
            if method == 'GET':
                if self.unavailable:
                    raise OSError('unavailable')
                data = {'running': self.running, 'run_id': self.starts}
            else:
                self.forwarded.append(json.loads(request.data))
                self.starts += 1
                self.running = True
                data = {'status': 'started', 'run_id': self.starts}
            response = io.BytesIO(json.dumps(data).encode())
            response.status = 200
            return response
        self.mock = patch.object(self.proxy.urllib.request, 'urlopen', upstream)
        self.mock.start()
        self.log = patch.object(self.proxy.Handler, 'log_message', lambda *args: None)
        self.log.start()
        self.http = self.proxy.ThreadingHTTPServer(('127.0.0.1', 0), self.proxy.Handler)
        self.thread = threading.Thread(target=self.http.serve_forever, kwargs={'poll_interval': .01})
        self.thread.start()

    def tearDown(self):
        self.http.shutdown()
        self.thread.join()
        self.http.server_close()
        self.mock.stop()
        self.log.stop()

    def request(self, path='/api/run', method='POST', body=None):
        connection = HTTPConnection(*self.http.server_address, timeout=5)
        connection.request(method, path, json.dumps({'question': 'ERBB2'} if body is None else body),
                           {'Content-Type': 'application/json'})
        response = connection.getresponse()
        status, body = response.status, response.read()
        connection.close()
        return status, body

    def test_dashboard_pipeline_options_are_forwarded(self):
        for reuse in [True, False]:
            with self.subTest(reuse_results=reuse):
                self.running = False
                body = {'question': 'TNBC therapeutic targets', 'design_mode': 'all_eligible',
                        'structural_backend': 'af3', 'reuse_results': reuse}
                status, response = self.request(body=body)
                self.assertEqual(status, 200, response)
                self.assertEqual(json.loads(response)['status'], 'started')
                self.assertEqual(self.forwarded[-1], body)
        self.assertEqual(self.starts, 2)

    def test_fixed_scenarios_are_forwarded(self):
        for scenario in ('tnbc_discovery', 'erbb2_design'):
            self.running = False
            body = {'scenario': scenario, 'reuse_results': True, 'structural_backend': 'af3', 'design_mode': 'all_eligible'}
            status, response = self.request(body=body)
            self.assertEqual(status, 200, response)
            self.assertEqual(self.forwarded[-1], body)
        self.running = False
        self.assertEqual(self.request(body={'scenario': 'unknown'})[0], 422)
        self.assertEqual(self.request(body={'scenario': 'erbb2_design', 'dry_run': True})[0], 422)
        self.assertEqual(self.starts, 2)

    def test_reuse_option_requires_a_boolean(self):
        for value in ['true', 'false', 0, 1, None, [], {}]:
            with self.subTest(value=value):
                status, response = self.request(body={'question': 'TNBC', 'reuse_results': value})
                self.assertEqual(status, 422)
                self.assertEqual(json.loads(response)['detail'], 'reuse_results must be boolean')
        self.assertEqual(self.starts, 0)

    def test_private_pipeline_options_remain_rejected(self):
        self.assertEqual(self.request(body={'question': 'TNBC', 'output_dir': '/tmp/run'})[0], 422)
        self.assertEqual(self.starts, 0)

    def test_finished_jobs_have_no_hourly_quota(self):
        for _ in range(5):
            self.running = False
            self.assertEqual(self.request()[0], 200)
        self.assertEqual(self.starts, 5)

    def test_all_execution_buttons_are_blocked_only_while_busy(self):
        self.running = True
        for route in self.proxy.EXECUTION_ROUTES:
            status, body = self.request(route)
            self.assertEqual(status, 409, route)
            self.assertEqual(json.loads(body)['status'], 'already_running')
        self.assertEqual(self.starts, 0)
        self.running = False
        self.assertEqual(self.request('/api/design')[0], 200)

    def test_simultaneous_requests_start_one_job(self):
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            statuses = [r[0] for r in pool.map(lambda _: self.request(), range(4))]
        self.assertEqual(sorted(statuses), [200, 409, 409, 409])
        self.assertEqual(self.starts, 1)

    def test_unavailable_status_does_not_launch_a_job(self):
        self.unavailable = True
        self.assertEqual(self.request()[0], 502)
        self.assertEqual(self.starts, 0)

    def test_stop_is_forwarded_while_busy(self):
        self.running = True
        self.assertEqual(self.request('/api/stop')[0], 200)

    def test_neon_page_is_not_public(self):
        self.assertEqual(self.request('/omicraft-version-7-1/neon', 'GET')[0], 404)


if __name__ == '__main__':
    unittest.main()

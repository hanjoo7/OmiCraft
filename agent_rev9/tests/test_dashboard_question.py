import asyncio
import json
import tempfile
import unittest
from collections import deque
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from .. import server, web_execution as web


class DashboardQuestionTests(unittest.TestCase):
    def setUp(self):
        self.question = 'TNBC · ERBB2 치료 표적을 분석해줘'
        self.saved_state, self.saved_queue = server.run_state, server.event_queue
        server.run_state = dict(running=False, run_id=0, result=None, history=[],
                                research_question=self.question, question=self.question)
        server.event_queue = deque()
        self.directory = tempfile.TemporaryDirectory()
        self.config = SimpleNamespace(model_copy=lambda **kw: SimpleNamespace())

    def tearDown(self):
        server.run_state, server.event_queue = self.saved_state, self.saved_queue
        self.directory.cleanup()

    def test_lab_start_preserves_research_question(self):
        with patch.object(web, 'asset', return_value={'gene': 'ESR1'}), \
             patch.object(web, 'lab_directory', return_value=Path(self.directory.name) / 'session'), \
             patch.object(web, 'get_config', return_value=self.config), \
             patch.object(server, 'start_worker'):
            response = web.start_lab('sm', {'gene': 'ESR1', 'ligand': 'tamoxifen'})
        self.assertEqual(response['status'], 'running')
        data = asyncio.run(server.api_events())
        self.assertEqual(data['research_question'], self.question)
        self.assertEqual(data['question'], self.question)
        self.assertEqual(data['task_label'], 'ESR1 Structure Lab docking')
        self.assertEqual(data['execution_kind'], 'lab')

    def test_erbb2_design_preserves_research_question(self):
        with patch.object(web, 'get_config', return_value=self.config), \
             patch.object(server, 'start_worker'):
            response = web.start_design({'parent_run_id': 'reference_erbb2', 'gene': 'ERBB2',
                                         'modality': 'ALL', 'input': {}})
        self.assertEqual(response['status'], 'started')
        self.assertEqual(server.run_state['research_question'], self.question)
        self.assertEqual(server.run_state['execution_kind'], 'design')

    def test_pipeline_updates_question_but_busy_request_does_not(self):
        self.assertEqual(server.begin_pipeline_run('New question'), 1)
        self.assertEqual(server.run_state['research_question'], 'New question')
        self.assertIsNone(server.begin_pipeline_run('Other question'))
        self.assertEqual(server.run_state['research_question'], 'New question')

    def test_restart_restores_question_separately_from_task(self):
        snapshot = Path(self.directory.name) / 'state.json'
        snapshot.write_text(json.dumps(dict(running=False, run_id=4, research_question=self.question,
            question='ESR1 Structure Lab docking', task_label='ESR1 Structure Lab docking',
            execution_kind='lab', events=[])))
        with patch.dict('os.environ', OMICRAFT_DASHBOARD_STATE=str(snapshot)):
            server.restore_dashboard_state()
        self.assertEqual(server.run_state['research_question'], self.question)
        self.assertEqual(server.run_state['task_label'], 'ESR1 Structure Lab docking')

    def test_legacy_lab_snapshot_does_not_become_a_question(self):
        snapshot = Path(self.directory.name) / 'state.json'
        snapshot.write_text(json.dumps(dict(running=False, run_id=4, question='ESR1 Structure Lab docking',
                                           report_run_id='lab_20260924T000000Z_12345678', events=[])))
        with patch.dict('os.environ', OMICRAFT_DASHBOARD_STATE=str(snapshot)):
            server.restore_dashboard_state()
        self.assertEqual(server.run_state['research_question'], '')
        self.assertEqual(server.run_state['execution_kind'], 'lab')

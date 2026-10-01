import json
import subprocess
import time
import threading
from collections import deque

import pytest
from fastapi.testclient import TestClient

from .. import server
from .. import web_execution as web
from ..configuration import get_config
from ..run_records import read_report


@pytest.fixture
def client(tmp_path, monkeypatch):
    config = get_config().model_copy(deep=True)
    config.data.base_dir = str(tmp_path)
    monkeypatch.setattr(web, 'get_config', lambda: config)
    monkeypatch.setattr(server, 'REPORT_RUNS_ROOT', tmp_path)
    monkeypatch.setattr(server, 'run_state', dict(running=False, run_id=0, result=None, history=[]))
    monkeypatch.setattr(server, 'event_queue', deque())
    def start_local(function, *args, **kwargs):
        # Test design behavior with mocked tools; process control is tested separately.
        if function.__name__ == 'run_bounded_design':
            function = web.run_design
        threading.Thread(target=function, args=args, kwargs=kwargs, daemon=True).start()
    monkeypatch.setattr(server, 'start_worker', start_local)
    with TestClient(server.app) as client:
        yield client


def parent(tmp_path, tier='2B', veto=False):
    folder = tmp_path / 'upstream_20260923T000000Z_12345678'
    folder.mkdir()
    target = dict(best_tier=tier, tier_results=[dict(modality='ADC', Tier=tier)],
                  safety_veto={'ADC': {'has_veto': veto}})
    (folder / 'pipeline_state.json').write_text(json.dumps({'qualification_results': {'ERBB2': target}}))
    return folder.name


def test_lab_keeps_direct_buttons_and_has_no_confirmation(client):
    page = client.get('/lab').text
    assert all(name in page for name in ['btnADC', 'btnSM', 'btnPROTAC', 'runDesign', 'designAsset'])
    assert 'confirm(' not in page
    assert '예상 시간' not in page
    assert 'Boltz2 구조' not in page


@pytest.mark.parametrize('tier,veto', [('HOLD', False), ('REJECT', False), ('3', False), ('2B', True)])
def test_design_rejects_unqualified_route(client, tmp_path, tier, veto):
    name = parent(tmp_path, tier, veto)
    r = client.post('/api/design', json=dict(parent_run_id=name, gene='ERBB2', modality='ADC'))
    assert r.status_code == 422
    assert not server.run_state['running']


def test_design_never_accepts_runtime_paths_or_other_gene_reference(client):
    for body in [dict(parent_run_id='reference_erbb2', gene='FOXA1', modality='ADC'),
                 dict(parent_run_id='reference_erbb2', gene='ERBB2', modality='ADC', input={'command': 'whoami'}),
                 dict(parent_run_id='../outside', gene='ERBB2', modality='ADC')]:
        assert client.post('/api/design', json=body).status_code == 422


def test_gpu_selection_stays_within_authorized_devices(monkeypatch):
    from .. import gpu_scheduler
    def query(command, **kwargs):
        if '--query-compute-apps=gpu_uuid,pid' in command:
            return subprocess.CompletedProcess(command, 0, stdout='GPU-5, 12345')
        return subprocess.CompletedProcess(command, 0, stdout='0, GPU-0, 195000, 0\n4, GPU-4, 5000, 0\n5, GPU-5, 130000, 0\n6, GPU-6, 150000, 99\n7, GPU-7, 140000, 0')
    monkeypatch.delenv('CUDA_VISIBLE_DEVICES', raising=False)
    monkeypatch.setattr(gpu_scheduler.subprocess, 'run', query)
    assert web.choose_gpu() == 7
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '7,4')
    assert web.choose_gpu() == 0
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '4')
    with pytest.raises(ValueError):
        web.choose_gpu()


def test_busy_run_is_not_replaced(client):
    server.run_state['running'] = True
    for path,body in [('/api/design',dict(parent_run_id='reference_erbb2', gene='ERBB2', modality='ADC')),
                      ('/api/lab/sm',dict(gene='ESR1', ligand='tamoxifen'))]:
        assert client.post(path, json=body).status_code == 409
    assert server.run_state['run_id'] == 0


def test_worker_failure_is_persisted_and_releases_lock(client, monkeypatch, tmp_path):
    from .. import modality_dispatch
    monkeypatch.setattr(web, 'design_input', lambda b: ({}, {}, None))
    monkeypatch.setattr(modality_dispatch, 'readiness', lambda *a: {'blockers': []})
    def fail(*args):
        raise RuntimeError('backend failed')
    monkeypatch.setattr(modality_dispatch, 'execute', fail)
    r = client.post('/api/design', json=dict(parent_run_id='reference_erbb2', gene='ERBB2', modality='ADC')).json()
    deadline = time.monotonic() + 3
    while server.run_state['running'] and time.monotonic() < deadline:
        time.sleep(.01)
    assert not server.run_state['running']
    report = read_report(tmp_path, run_id=r['report_run_id'])
    assert report['execution_status'] == 'FAILED'
    assert report['blockers'] == ['backend failed']
    assert report['selection_event']['tier_at_selection'] == 'UNKNOWN'


def test_structure_upload_is_bounded_and_gene_checked(client, tmp_path):
    # A minimal real-format CA record is enough for registry validation.
    content = 'ATOM      1  CA  ALA A   1      11.000  12.000  13.000  1.00 20.00           C  \nEND\n'
    response = client.post('/api/structures', json=dict(gene='FOXA1', format='pdb', content=content))
    assert response.status_code == 200
    row = response.json()
    assert 'path' not in row
    assert web.asset(row['id'], 'FOXA1')['gene'] == 'FOXA1'
    with pytest.raises(ValueError):
        web.asset(row['id'], 'ERBB2')
    assert client.post('/api/structures', json=dict(gene='../FOXA1', format='pdb', content=content)).status_code == 422
    assert client.post('/api/structures', json=dict(gene='FOXA1', format='pdb', content='no atoms')).status_code == 422


def test_lab_rejects_mismatched_reference_and_session_paths(client):
    assert client.post('/api/lab/sm', json=dict(gene='ERBB2', pdb_id='3ERT', ligand='tak285')).status_code == 422
    assert client.get('/api/lab/status/not-a-session').status_code == 422


def test_execution_completion_includes_admet_and_comparison(tmp_path, monkeypatch):
    from .. import modality_dispatch, small_molecule_workflow
    summary = tmp_path / 'candidate_summary.csv'
    summary.write_text('candidate_id\nx\n')
    candidate = dict(candidate_id='x', preparation={'receptor': {'status': 'success'}, 'ligand': {'status': 'success'}},
                     docking={'status': 'success'}, pose_qc=[{'pose_qc_status': 'fail'}],
                     af3_ligand={'status': 'success'}, consensus={'comparisons': [{}]},
                     validation={'screening_call': 'INDETERMINATE'}, critic={'verdict': 'NEEDS_REVIEW'},
                     admet={'status': 'skipped_structural_screening'})
    def run(*args, **kwargs):
        return dict(candidates=[candidate], candidate_summary_path=str(summary))
    monkeypatch.setattr(small_molecule_workflow, 'run_small_molecule_workflow', run)
    config = {'output_dir': str(tmp_path / 'run')}
    data = {'candidates': [], 'context': {}}
    result = modality_dispatch.execute_small_molecule(data, config)
    assert result['execution_status'] == 'PARTIAL'
    assert result['blockers'] == ['admet_incomplete:x']
    candidate['admet']['status'] = 'success'
    result = modality_dispatch.execute_small_molecule(data, config)
    assert result['execution_status'] == 'COMPLETED'
    assert candidate['validation']['screening_call'] == 'INDETERMINATE'

    candidate['admet']['status'] = 'not_run'
    config['small_molecule'] = {'run_admet': False}
    result = modality_dispatch.execute_small_molecule(data, config)
    assert result['execution_status'] == 'COMPLETED'
    stage = next(row for row in result['summary']['execution_checklist'] if row['stage'] == 'admet')
    assert stage['execution_status'] == 'NOT_REQUESTED'


def test_binder_sampling_budget_is_bounded_and_preserved(monkeypatch, tmp_path):
    structure=tmp_path/'target.pdb'
    structure.write_text('ATOM      1  CA  ALA A   1      11.000  12.000  13.000  1.00 20.00           C  \nEND\n')
    monkeypatch.setattr(web,'candidate_data',lambda *a: {'tier':'2A'})
    monkeypatch.setattr(web,'asset',lambda *a: {'path':str(structure),'id':'fixture','gene':'MMP7','kind':'protein'})
    body={'parent_run_id':'parent','gene':'MMP7','modality':'DE_NOVO_BINDER',
          'input':{'asset_id':'fixture','target_chain':'A','num_designs':4,'num_sequences':2,'design_seed':73}}
    data,cfg,_=web.design_input(body)
    assert data['protein_design_input']['num_designs']==4
    assert data['protein_design_input']['seed']==73
    assert cfg['screening_options']['max_rfd3_candidates']==4
    assert cfg['screening_options']['max_mpnn_candidates']==8
    assert cfg['screening_options']['max_af3_candidates']==8
    for key,value in [('num_designs',0),('num_designs',9),('num_designs',True),('design_seed',-1),('num_sequences',9)]:
        invalid={**body,'input':{**body['input'],key:value}}
        with pytest.raises(ValueError):web.design_input(invalid)


def test_resume_reuses_verified_sequences_and_rejects_wrong_target(monkeypatch,tmp_path):
    from Bio import SeqIO
    from ..protein_mpnn_tool import ProteinMPNNTool
    cfg=get_config().model_copy(deep=True);cfg.data.base_dir=str(tmp_path)
    monkeypatch.setattr(web,'get_config',lambda: cfg)
    identifier='design_20261001T000000Z_12345678'
    folder=tmp_path/identifier;target=folder/'execution/binder/screening_test/target_000'
    mpnn=target/'target_000_candidate_0000/protein_mpnn';mpnn.mkdir(parents=True)
    backbone=target/'rfd3';backbone.mkdir()
    (backbone/'rfd3_input_protein_binder_0_model_0.cif.gz').write_bytes(b'fixture')
    fasta=mpnn/'sample.fa';fasta.write_text('>backbone, sample=1\nAAAAA\n')
    sequence={'candidate_id':'candidate_0_seq_1','sequence':'AAAAA','score':1.0,'status':'MPNN_COMPLETED','sequence_path':str(fasta),'is_mock':False}
    (mpnn/'protein_mpnn_run.json').write_text(json.dumps({'status':'success','is_mock':False,'sequences':[sequence]}))
    body={'gene':'MMP7','modality':'DE_NOVO_BINDER','parent_run_id':'parent'}
    (folder/'worker_input.json').write_text(json.dumps({'body':body,'data':{'protein_design_input':{'binder_chain':'A','design_length':5}},'config':{}}))
    data,_,_=web.resume_binder_input(body,identifier,{'Tier':'2A'})
    assert data['protein_design_result']['candidates'][0]['protein_mpnn']['sequence']=='AAAAA'
    assert data['protein_design_result']['source_run_id']==identifier
    with pytest.raises(ValueError):web.resume_binder_input({**body,'gene':'ERBB2'},identifier,{})
    with pytest.raises(ValueError):web.resume_binder_input({**body,'parent_run_id':'other'},identifier,{})
    from ..result_reuse import resume_partial
    saved=json.loads((folder/'worker_input.json').read_text())
    for key,name in [('body','request'),('data','input'),('config','config')]:
        (folder/(name+'.json')).write_text(json.dumps(saved[key]))
    (folder/'worker_input.json').unlink()
    (folder/'run_summary.json').write_text(json.dumps({'execution_status':'FAILED','running':False}))
    resumed=resume_partial(tmp_path,body,saved['data'],saved['config'],{'Tier':'2A'})
    assert resumed[1]['candidate_count']==1
    assert resumed[1]['source_run_id']==identifier
    assert resume_partial(tmp_path,{**body,'reuse_results':False},saved['data'],saved['config'],{}) is None
    fasta.write_text('>backbone, sample=1\nGGGGG\n')
    with pytest.raises(ValueError):web.resume_binder_input(body,identifier,{})


def test_viewer_cif_long_ligand_name_preserves_pdb_coordinates(tmp_path):
    from Bio.PDB import Atom, Chain, Model, Structure, Residue, MMCIFIO
    import numpy as np
    from ..web_execution import pdb_text
    structure = Structure.Structure('complex')
    model = Model.Model(0)
    chain = Chain.Chain('B')
    ligand = Residue.Residue(('H_LIG_B', 1, ' '), 'LIG_B', ' ')
    ligand.add(Atom.Atom('C1', np.array([8.455, -2.032, 15.042]), 20.06, 1., ' ', ' C1 ', 1, element='C'))
    chain.add(ligand); model.add(chain); structure.add(model)
    path = tmp_path / 'ligand.cif'
    writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(path))
    original = path.read_bytes()
    line = next(line for line in pdb_text({'path': str(path)}).splitlines() if line.startswith('HETATM'))
    assert line[17:20] == 'LIG' and line[21] == 'B'
    assert [float(line[a:b]) for a,b in [(30,38),(38,46),(46,54)]] == [8.455,-2.032,15.042]
    assert path.read_bytes() == original

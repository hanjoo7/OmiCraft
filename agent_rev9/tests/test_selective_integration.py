import ast
import gzip
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Residue
from fastapi.testclient import TestClient
from rdkit import Chem

from .. import modality_dispatch as dispatch
from .. import server
from ..antibody_mapping import contact_pairs, map_sequence
from ..antibody_numbering import AnarciiBackend, cdr_region, normalize
from ..configuration import ProtenixConfig
from ..docking_backend import (
    DockingRequest,
    NativeDockingBackend,
    resolve_executable,
    tool_readiness,
)
from ..protein_design_route import _read_af3_json
from ..routes import MODALITY_ROUTES, normalize_route_state, validate_handoff_routes
from ..run_records import RunRecord
from ..small_molecule_config import DockingConfig
from ..user_selection import check_design_gate, create_selection_event


@pytest.mark.parametrize('modality,route', MODALITY_ROUTES.items())
def test_routes_and_dispatch(modality, route, tmp_path, monkeypatch):
    assert normalize_route_state({'route': modality})['route'] == route
    validate_handoff_routes([{'modality': modality, 'route': route}], route)
    monkeypatch.setattr(dispatch, 'readiness', lambda *a: {'tool_status': 'ready', 'blockers': []})
    called = []
    for name in ('small_molecule', 'binder', 'adc', 'degrader'):
        def executor(data, config, name=name):
            called.append(name)
            return dispatch.result(modality, 'COMPLETED', 'test', is_mock=True)
        monkeypatch.setattr(dispatch, 'execute_' + name, executor)
    output = dispatch.execute({'modality': modality}, {'output_dir': tmp_path / 'run', 'dry_run': False})
    assert called == [{'SMALL_MOLECULE': 'small_molecule', 'DE_NOVO_BINDER': 'binder',
                       'ADC': 'adc', 'DEGRADER': 'degrader'}[modality]]
    assert output['validation_decision'] == 'HOLD'
    assert {p.name for p in (tmp_path / 'run').iterdir()} == {'events.jsonl', 'run_summary.json'}


@pytest.mark.parametrize('modality,wrong', [('ADC', 'protein_binder_rfd3'), ('DEGRADER', 'small_molecule')])
def test_wrong_mapping_rejected(modality, wrong):
    with pytest.raises(ValueError, match='conflicts'):
        validate_handoff_routes([{'modality': modality, 'route': wrong}])


def test_user_selection_requires_known_tier_and_user():
    event = create_selection_event('ERBB2', 'ADC', 'ERBB2:ADC', '1A', rationale='reviewed')
    assert {'candidate', 'modality', 'timestamp', 'reason'} <= event.keys()
    assert check_design_gate('ERBB2', 'ADC', 'UNKNOWN', [event])['gate_status'] == 'BLOCKED_UNKNOWN_TIER'
    for actor in ('auto', '', 'model', None):
        assert not check_design_gate('ERBB2', 'ADC', '1A', [{**event, 'selected_by': actor}])['allowed']
    assert check_design_gate('ERBB2', 'ADC', '1A', [event])['allowed']


def test_gnina_priority_and_no_backend_fallback(tmp_path, monkeypatch):
    from .. import docking_backend
    monkeypatch.setattr(docking_backend.shutil, "which", lambda value, path=None: str(Path(value) if Path(value).is_absolute() else Path(path or "") / value) if (Path(value) if Path(value).is_absolute() else Path(path or "") / value).is_file() else None)
    paths = []
    for name in ('cli', 'config', 'environment', 'gnina', 'vina'):
        p = tmp_path / name
        p.write_text('#!/bin/sh\nexit 0\n')
        p.chmod(0o755)
        paths.append(str(p))
    cli, cfg, env_path, gnina, vina = paths
    env = {'GNINA_EXECUTABLE': env_path, 'PATH': str(tmp_path), 'OMICRAFT_VINA': vina}
    assert resolve_executable(cli=cli, configured=cfg, environ=env)['path'] == cli
    assert resolve_executable(configured=cfg, environ=env)['path'] == cfg
    assert resolve_executable(environ=env)['path'] == env_path
    env.pop('GNINA_EXECUTABLE')
    assert resolve_executable(environ=env)['path'] == gnina
    assert resolve_executable('vina', configured=cfg, environ=env)['path'] == vina
    assert resolve_executable(configured='/missing/gnina', environ=env)['blocker'] == 'gnina_executable_missing'
    assert resolve_executable(environ={'PATH': ''})['path'] is None


def test_gnina_version_failure_and_autobox(tmp_path, monkeypatch):
    from .. import docking_backend
    monkeypatch.setattr(docking_backend.shutil, "which", lambda value, **kw: str(value) if value and Path(value).is_file() else None)
    binary = tmp_path / 'gnina'
    binary.touch()
    binary.chmod(0o755)
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=127, stdout='', stderr='library missing'))
    assert tool_readiness(configured=str(binary))['blocker'] == 'gnina_version_failed'
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stdout='gnina v1', stderr=''))
    receptor, ligand = tmp_path / 'r.pdb', tmp_path / 'l.sdf'
    receptor.touch()
    ligand.touch()
    backend = NativeDockingBackend(DockingConfig(backend='gnina', executable=str(binary), cnn_scoring='none', device=5))
    job = backend.prepare(DockingRequest('x', str(receptor), str(ligand), (), (), 'CC', str(tmp_path / 'out'), autobox_ligand_path=str(ligand)))
    assert '--autobox_ligand' in job.command and '--center_x' not in job.command
    assert backend.gpu_device == 5 and backend.version == 'gnina v1'


def test_score_missing_is_failure(tmp_path):
    mol = Chem.MolFromSmiles('CC')
    mol.SetProp('CNNscore', '0.9')
    with Chem.SDWriter(str(tmp_path / 'poses.sdf')) as writer:
        writer.write(mol)
    row = NativeDockingBackend(DockingConfig(backend='gnina')).parse(tmp_path, 'CC')[0]
    assert row['docking_status'] == 'failed'
    assert row['error_message'] == 'docking_score_missing'
    assert 'CNNscore' in row['raw_score_fields']


def test_numbering_unavailable_blocks_mapping(tmp_path, monkeypatch):
    monkeypatch.setattr(AnarciiBackend, 'availability', lambda s: {'status': 'not_installed', 'error': 'missing'})
    output = AnarciiBackend().number('H', 'ACD')
    assert output.numbering_status == 'backend_unavailable' and not output.numbering
    structure = tmp_path / 'complex.cif'
    structure.touch()
    run = dispatch.execute({'modality': 'ADC', 'input': {'structure_path': str(structure), 'antigen_chain': 'C'}},
                           {'output_dir': tmp_path / 'run', 'dry_run': False})
    assert run['execution_status'] == 'BLOCKED'
    assert run['summary']['cdr_mapping_status'] == 'not_run'
    assert run['validation_decision'] == 'NOT_EVALUATED'


@pytest.mark.parametrize('kind', ['H', 'K', 'L'])
def test_heavy_light_numbering_and_error(kind):
    raw = {'chain_type': kind, 'score': 99, 'query_start': 0, 'query_end': 2,
           'scheme': 'imgt', 'numbering': [((27, ''), 'A'), ((28, 'A'), 'C'), ((29, ''), 'D')]}
    out = normalize('x', 'ACD', raw, 'ANARCII', '2.0.8')
    assert out.numbering_status == 'success' and out.chain_type == kind
    assert out.score == 99 and out.backend_version == '2.0.8'
    assert out.numbering[1]['imgt_insertion_code'] == 'A'
    assert cdr_region(out.numbering[1]['imgt_position']) == 'CDR1'
    failed = normalize('x', 'ACD', {**raw, 'error': 'cannot number'}, 'ANARCII', '2.0.8')
    assert failed.numbering == [] and failed.error_message == 'cannot number'


def residues(cid='H', offset=0):
    chain = Chain.Chain(cid)
    for idx, (aa, code) in enumerate([('ALA', ''), ('CYS', 'A'), ('ASP', '')]):
        residue = Residue.Residue((' ', idx+10, code or ' '), aa, '')
        residue.add(Atom.Atom('CA', np.array([idx*4., offset, 0.]), 80., 1., ' ', ' CA ', idx, element='C'))
        chain.add(residue)
    return list(chain)


def test_mapping_insertion_missing_and_contacts():
    observed = residues()
    rows, qc = map_sequence('ACD', observed)
    assert qc['identity'] == qc['coverage'] == 1
    assert rows[1]['insertion_code'] == 'A'
    _, partial = map_sequence('ACED', observed, min_coverage=.9)
    assert partial['status'] == 'failed' and partial['missing_sequence_indices']
    contacts = contact_pairs(observed, residues('C', 3.0))
    assert len(contacts) == 3 and sum(c['severe_clash_count'] for c in contacts) == 0
    assert sum(c['severe_clash_count'] for c in contact_pairs(observed, residues('C', .5))) == 3


def test_gzip_and_invalid_af3_json(tmp_path):
    value = {'atom_plddts': [80, 90]}
    plain, compressed = tmp_path / 'out.json', tmp_path / 'out.json.gz'
    plain.write_text(json.dumps(value))
    compressed.write_bytes(gzip.compress(plain.read_bytes()))
    assert _read_af3_json(plain) == _read_af3_json(compressed) == value
    compressed.write_bytes(b'invalid gzip')
    with pytest.raises((OSError, ValueError)):
        _read_af3_json(compressed)
    compressed.write_bytes(gzip.compress(b'invalid json'))
    with pytest.raises(ValueError):
        _read_af3_json(compressed)


def test_protenix_cannot_be_enabled_and_not_called(state, fake_models, monkeypatch):
    from pydantic import ValidationError

    from ..critic_agent import _apply_deterministic_rules
    from ..screening_agent import screening_node
    from ..screening_model_tools import ProtenixTool
    monkeypatch.setattr(ProtenixTool, 'run', lambda *a, **k: pytest.fail('Protenix must not run'))
    output = screening_node(state)
    rows = output['screening_result']['candidates']
    assert rows and all(c['cross_validation_status'] == 'not_used' for c in rows)
    assert all(c['cross_validation_required'] is False for c in rows)
    assert output['screening_result']['protenix_tool_calls'] == 0
    assert _apply_deterministic_rules(output)['verdict'] == 'APPROVE'
    with pytest.raises(ValidationError):
        ProtenixConfig(enabled=True)


def test_existing_report_pages_and_run_artifact(tmp_path, monkeypatch):
    monkeypatch.setattr(server, 'REPORT_RUNS_ROOT', tmp_path)
    record = RunRecord(tmp_path / 'run1', 'ADC')
    raw = tmp_path / 'original.csv'
    raw.write_text('a,b\n1,2\n')
    record.finish(dispatch.result('ADC', 'COMPLETED', 'mapping', artifacts=[{'path': str(raw)}]))
    with TestClient(server.app) as client:
        response = client.get('/api/report?run_id=run1')
        assert response.status_code == 200
        data = response.json()
        assert data['execution_status'] == 'COMPLETED'
        assert client.get(data['artifacts'][0]['href']).text == raw.read_text()
        assert client.get('/api/report?run_dir=../outside').status_code == 400
        assert client.get('/api/report?run_id=missing').status_code == 400
        assert client.get('/api/report').status_code == 200
        for page in ('/report', '/demo'):
            html = client.get(page + '?run_id=run1')
            assert html.status_code == 200 and 'renderTherapeuticRun' in html.text
            assert 'Chart' in html.text


def test_source_has_no_other_project_imports():
    root = Path(server.__file__).parent
    for path in root.rglob('*.py'):
        if 'tests' in path.parts or '__pycache__' in path.parts:
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            names = [a.name for a in node.names] if isinstance(node, ast.Import) else [node.module or ''] if isinstance(node, ast.ImportFrom) and not node.level else []
            assert not any(n.startswith(('agent_rev5', 'agent_rev4', 'omicraft.agent')) for n in names), path


def test_degrader_component_assembly_is_independent(tmp_path):
    output = dispatch.execute({'modality': 'DEGRADER', 'input': {
        'warhead_smiles': 'CC[*:1]', 'e3_ligand_smiles': 'CC[*:2]',
        'linker_smiles': '[*:3]CC[*:4]', 'warhead_map': 1, 'e3_map': 2, 'linker_maps': [3, 4], 'expected_full_smiles': 'CCCCCC'}},
        {'output_dir': tmp_path / 'run', 'dry_run': False})
    assert output['execution_status'] == 'PARTIAL', output
    assert output['validation_decision'] == 'HOLD'
    assert Path(output['artifacts'][0]['path']).is_file()
    assert output['summary']['assembly_status'] == 'success'


def test_adc_numbering_to_structure_connection(tmp_path):
    from Bio.PDB import PDBIO, Model, Structure

    from ..antibody_mapping import analyze

    structure = Structure.Structure('fixture')
    model = Model.Model(0)
    structure.add(model)
    for chain, offset in [('H', 0), ('L', 6), ('C', 3)]:
        model.add(residues(chain, offset)[0].parent)
    writer = PDBIO()
    writer.set_structure(structure)
    path = tmp_path / 'complex.pdb'
    writer.save(str(path))

    class NumberingFixture:
        def number(self, cid, sequence):
            return normalize(cid, sequence, {'chain_type': 'H' if cid == 'H' else 'K',
                'scheme': 'imgt', 'score': 98, 'query_start': 0, 'query_end': len(sequence)-1,
                'numbering': [((27+i, ''), aa) for i, aa in enumerate(sequence)]}, 'fixture', 'test')

    data = {'structure_path': str(path), 'antigen_chain': 'C', 'antibody_sequences': {'H': 'ACD', 'L': 'ACD'}}
    output = analyze(data, {}, backend=NumberingFixture())
    assert output['cdr_mapping_status'] == 'success'
    assert output['cdr_contact_residues'] == 6
    assert output['severe_clashes'] == 0
    assert output['mapping'][1]['insertion_code'] == 'A'
    data['antibody_sequences']['H'] = 'ACED'
    failed = analyze(data, {}, backend=NumberingFixture())
    assert failed['cdr_mapping_status'] == 'not_run'
    assert failed['contacts'] == []
    assert failed['mapping_qc']['H']['missing_sequence_indices'] == [3]


def test_1n8z_contact_regression():
    # Two observed atoms from 1N8Z; this fragment is a distance regression only.
    chains = []
    for cid, number, aa, atom_name, element, xyz in [
        ('A', 92, 'TYR', 'OH', 'O', (45.008, 94.273, 92.612)),
        ('C', 569, 'LYS', 'NZ', 'N', (43.324, 95.132, 94.886)),
    ]:
        chain = Chain.Chain(cid)
        residue = Residue.Residue((' ', number, ' '), aa, '')
        residue.add(Atom.Atom(atom_name, np.array(xyz), 0., 1., ' ', atom_name, 1, element=element))
        chain.add(residue)
        chains.append(list(chain))
    row = contact_pairs(*chains)[0]
    assert row['minimum_distance'] == pytest.approx(2.957164, abs=1e-5)
    assert row['heavy_atom_contact_count'] == 1 and row['severe_clash_count'] == 0


def test_agent_selection_dispatch_validation_critic(tmp_path, monkeypatch):
    from ..orchestrator import run_screening

    monkeypatch.setattr(dispatch, 'readiness', lambda *a: {'tool_status': 'ready', 'blockers': []})
    calls = []
    def adc(data, cfg):
        calls.append(data)
        return dispatch.result('ADC', 'COMPLETED', 'cdr_mapping', summary={'severe_clashes': 0})
    monkeypatch.setattr(dispatch, 'execute_adc', adc)
    selection = create_selection_event('ERBB2', 'ADC', 'ERBB2:ADC', '1A', rationale='explicit choice')
    state = {'advance_targets': [{'gene_name': 'ERBB2', 'modality': 'SMALL_MOLECULE', 'tier': '1A'}],
             'user_selection_events': [selection], 'execution_inputs': {'ADC': {'fixture': True}},
             'output_dir': str(tmp_path / 'run'), 'dry_run': False}
    output = run_screening(state)
    assert len(output['modality_assessments']) == 4 and len(calls) == 1
    item = output['modality_results'][0]
    assert item['modality'] == 'ADC' and item['validation_decision'] == 'HOLD'
    assert item['critic']['decision'] == 'HOLD'
    assert item['selection_event'] == selection
    stages = [s['stage'] for s in item['stages']]
    assert stages.index('validation') < stages.index('critic')


def test_safety_veto_blocks_selected_executor(tmp_path, monkeypatch):
    from ..orchestrator import run_screening

    monkeypatch.setattr(dispatch, 'readiness', lambda *a: {'tool_status': 'ready', 'blockers': []})
    monkeypatch.setattr(dispatch, 'execute_adc', lambda *a: pytest.fail('safety veto bypass'))
    selection = create_selection_event('ERBB2', 'ADC', 'ERBB2:ADC', '1A')
    output = run_screening({'advance_targets': [{'gene_name': 'ERBB2', 'tier': '1A', 'safety_verdict': 'VETO'}],
        'user_selection_events': [selection], 'execution_inputs': {'ADC': {}},
        'output_dir': tmp_path / 'run', 'dry_run': False})
    assert not output['modality_results']
    assert all(not a['gate']['allowed'] for a in output['modality_assessments'])


@pytest.mark.external_tool
def test_installed_gnina_version():
    ready = tool_readiness()
    if ready['blocker'] == 'gnina_executable_missing':
        pytest.skip('GNINA_EXECUTABLE or PATH is not configured')
    assert ready['tool_status'] == 'ready', ready


@pytest.mark.external_tool
def test_installed_anarcii_version():
    ready = AnarciiBackend().availability()
    if ready['status'] == 'not_installed':
        pytest.skip('ANARCII is not configured')
    assert ready['status'] == 'ready' and ready['version']


def test_binder_configuration_is_scoped(tmp_path, monkeypatch):
    from .. import screening_agent
    from ..configuration import get_config

    original = get_config().rfd3.executable
    seen = []
    def screen(state):
        seen.append((get_config().rfd3.executable, state['screening_options']['max_af3_candidates']))
        return {'screening_result': {'status': 'COMPLETED', 'candidates': []}}
    monkeypatch.setattr(screening_agent, 'screening_node', screen)
    dispatch.execute_binder({}, {'output_dir': tmp_path, 'rfd3': {'executable': '/configured/rfd3'},
                                 'screening_options': {'max_af3_candidates': 2}})
    assert seen == [('/configured/rfd3', 2)]
    assert get_config().rfd3.executable == original


def test_report_manifest_validation_and_path_boundaries(tmp_path, monkeypatch):
    monkeypatch.setattr(server, 'REPORT_RUNS_ROOT', tmp_path / 'runs')
    root = tmp_path / 'runs'
    root.mkdir()
    outside = tmp_path / 'outside'
    outside.mkdir()
    (root / 'escape').symlink_to(outside, target_is_directory=True)
    (root / 'empty').mkdir()
    malformed = root / 'malformed'
    malformed.mkdir()
    client = TestClient(server.app)
    for directory in ['../outside', str(outside), 'escape', 'missing', 'empty', 'empty/../empty']:
        assert client.get('/api/report', params={'run_dir': directory}).status_code == 400
    for content in ['{', '[]', '{}', '{"execution_status":"COMPLETED"}']:
        (malformed / 'run_summary.json').write_text(content)
        assert client.get('/api/report?run_id=malformed').status_code == 400
    record = RunRecord(root / 'valid', 'ADC')
    record.finish(dispatch.result('ADC', 'PARTIAL', 'mapping'))
    assert client.get('/api/report', params={'run_dir': str(root / 'valid')}).status_code == 200
    assert client.get('/api/report', params={'run_id': 'valid', 'run_dir': 'valid'}).status_code == 400


def test_af3_unique_sample_gzip_and_missing_metrics(tmp_path):
    from ..protein_design_route import AlphaFold3BinderAdapter

    sample = tmp_path / 'job/seed-11_sample-0'
    sample.mkdir(parents=True)
    for directory, prefix in [(sample, 'job_seed-11_sample-0_'), (sample.parent, 'job_')]:
        (directory / (prefix + 'model.cif')).touch()
        payload = {'ranking_score': .8, 'iptm': .7, 'has_clash': False}
        summary = directory / (prefix + 'summary_confidences.json')
        summary.write_bytes(json.dumps(payload).encode())
    summary = sample / 'job_seed-11_sample-0_summary_confidences.json'
    compressed = summary.with_suffix('.json.gz')
    compressed.write_bytes(gzip.compress(summary.read_bytes()))
    summary.unlink()
    (sample / 'job_seed-11_sample-0_confidences.json.gz').write_bytes(gzip.compress(b'{"atom_plddts":[]}'))
    parsed = AlphaFold3BinderAdapter.parse_output(tmp_path)
    assert parsed['sample_count'] == 1
    assert parsed['samples'][0]['seed'] == 11 and parsed['samples'][0]['sample_id'] == 0
    assert parsed['samples'][0]['ptm'] is None
    assert parsed['plddt_summary']['mean'] is None
    assert parsed['iptm'] == .7


def test_binder_sample_budget_reaches_official_command(tmp_path):
    from ..protein_design_route import AlphaFold3BinderAdapter

    adapter = AlphaFold3BinderAdapter(python_path='python', script_path='run_alphafold.py',
        model_dir='models', database_dir='db', num_samples=1)
    out = adapter.run(target_sequence='ACDE', binder_sequence='ACDE', output_dir=str(tmp_path), dry_run=True)
    assert '--num_diffusion_samples=1' in out.command
    with pytest.raises(ValueError):
        AlphaFold3BinderAdapter(python_path='', script_path='', model_dir='', database_dir='', num_samples=0)


def test_receptor_altloc_consistent_across_formats(tmp_path):
    from Bio.PDB import MMCIFIO, PDBParser

    from ..ligand_geometry import protein_atoms

    pdb = tmp_path / 'r.pdb'
    pdb.write_text('ATOM      1  CA ASER A 974      18.487  30.545  49.650  0.60 10.00           C  \n'
                   'ATOM      2  CA BSER A 974      18.498  30.530  49.652  0.40 10.00           C  \nEND\n')
    structure = PDBParser(QUIET=True).get_structure('r', pdb)
    cif = tmp_path / 'r.cif'
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(cif))
    assert protein_atoms(cif, 'A')[0]['ca'] == pytest.approx(protein_atoms(pdb, 'A')[0]['ca'])


@pytest.mark.parametrize('modality', list(MODALITY_ROUTES))
def test_all_agent_choices_dispatch_only_selected_executor(tmp_path, monkeypatch, modality):
    from ..orchestrator import run_screening

    names = {'SMALL_MOLECULE': 'small_molecule', 'DE_NOVO_BINDER': 'binder', 'ADC': 'adc', 'DEGRADER': 'degrader'}
    monkeypatch.setattr(dispatch, 'readiness', lambda *a: {'tool_status': 'ready', 'blockers': []})
    calls = []
    for key, name in names.items():
        def spy(data, config, key=key):
            calls.append(key)
            return dispatch.result(key, 'NOT_RUN', 'routing_only')
        monkeypatch.setattr(dispatch, 'execute_' + name, spy)
    event = create_selection_event('ERBB2', modality, 'ERBB2:' + modality, '1A')
    state = {'advance_targets': [{'gene_name': 'ERBB2', 'tier': '1A'}],
             'execution_inputs': {key: {'routing_only': True} for key in names},
             'user_selection_events': [event], 'output_dir': tmp_path / 'selected', 'dry_run': False}
    output = run_screening(state)
    assert calls == [modality]
    item = output['modality_results'][0]
    assert item['validation_decision'] == item['critic']['decision'] == 'NOT_EVALUATED'
    assert [s['stage'] for s in item['stages']][-3:] == ['validation', 'critic', 'routing_only']
    for case, target, selection in [
        ('auto', {'gene_name': 'ERBB2', 'tier': '1A'}, {**event, 'selected_by': 'auto'}),
        ('unknown', {'gene_name': 'ERBB2', 'tier': 'UNKNOWN'}, event),
        ('veto', {'gene_name': 'ERBB2', 'tier': '1A', 'safety_verdict': 'VETO'}, event),
    ]:
        calls.clear()
        output = run_screening({**state, 'advance_targets': [target], 'user_selection_events': [selection],
                                'output_dir': tmp_path / case})
        assert not calls and not output['modality_results']

"""Structure Lab assets and persisted design jobs for the existing dashboard."""

import copy
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import nullcontext
import gzip
import html
import io
import json
import re
import shutil
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from fastapi import HTTPException

from .configuration import get_config
from .runtime_paths import config_path, load_json, workspace
from .gpu_scheduler import available_devices, reserve_gpu
from .run_records import RunRecord, register_demo_artifact, run_directory
from .small_molecule_io import write_json

PACKAGE = Path(__file__).resolve().parent
WORKSPACE = workspace()
MODELS = config_path('workbench_models', 'OMICRAFT_MODELS_CONFIG')
ASSETS = config_path('structure_assets', 'OMICRAFT_ASSETS_CONFIG')
REFERENCE = Path(os.environ.get('OMICRAFT_REFERENCE_RUN', WORKSPACE / 'runs/real_validation_20260922T133906Z'))


def storage():
    return Path(get_config().data.base_dir).resolve() / 'structure_lab'


def assets():
    rows = load_json(ASSETS)['assets'] if ASSETS.is_file() else []
    for path in (storage() / 'uploads').glob('*/asset.json'):
        rows.append(json.loads(path.read_text()))
    return [row for row in rows if Path(row['path']).is_file()]


def asset(identifier=None, gene=None):
    rows = assets()
    found = next((row for row in rows if row['id'] == identifier), None) if identifier else next(
        (row for row in rows if row['gene'] == gene and row['kind'] == 'protein' and row['role'] == 'target'), None)
    if not found:
        raise ValueError('구조 입력이 없습니다. Structure Lab에서 PDB/mmCIF 파일을 등록하세요.')
    if gene and found['gene'] != gene:
        raise ValueError('선택한 구조의 표적과 설계 표적이 다릅니다.')
    return found


def catalog():
    return [{k: v for k, v in row.items() if k != 'path'} for row in assets()]


def upload(body):
    from .structure_io import protein_residues
    gene = body.get('gene', '')
    if not isinstance(gene, str):
        raise ValueError('Invalid gene')
    gene = gene.upper()
    kind = body.get('format')
    text = body.get('content')
    if not re.fullmatch(r'[A-Z0-9_-]{1,40}', gene) or kind not in {'pdb', 'cif'}:
        raise ValueError('유전자명과 pdb/cif 형식을 확인하세요.')
    if not isinstance(text, str) or not 0 < len(text.encode()) <= 5_000_000:
        raise ValueError('구조 파일은 5 MB 이하여야 합니다.')
    identifier = 'uploaded_' + uuid.uuid4().hex
    folder = storage() / 'uploads' / identifier
    folder.mkdir(parents=True)
    path = folder / ('structure.' + kind)
    try:
        path.write_text(text)
        residues = protein_residues(str(path))
        if not residues or len(residues) > 10000:
            raise ValueError('표준 단백질 잔기를 가진 구조가 필요합니다 (최대 10,000 잔기).')
        row = dict(id=identifier, gene=gene, path=str(path), label=gene + ' · uploaded structure',
                   kind='protein', role='target', source='User supplied; target identity not independently verified',
                   chains=sorted({r[0] for r in residues}))
        write_json(folder / 'asset.json', row)
        return {k: v for k, v in row.items() if k != 'path'}
    except Exception:
        shutil.rmtree(folder)
        raise


def _fetch_alphafold_structure(gene: str) -> dict:
    """AlphaFold DB에서 예측 구조를 다운로드하여 asset으로 등록."""
    import urllib.request

    # UniProt ID 조회 (reviewed/canonical 우선)
    for query_suffix in ["+AND+reviewed:true", ""]:
        uniprot_url = (f"https://rest.uniprot.org/uniprotkb/search?"
                       f"query=gene_exact:{gene}+AND+organism_id:9606{query_suffix}&format=json&size=1")
        try:
            with urllib.request.urlopen(uniprot_url, timeout=15) as resp:
                data = json.loads(resp.read())
            if data.get("results"):
                accession = data["results"][0]["primaryAccession"]
                break
        except Exception:
            continue
    else:
        raise ValueError(f"{gene}: UniProt 매핑 실패. Structure Lab에서 PDB/mmCIF 파일을 직접 등록하세요.")

    # AlphaFold API로 최신 CIF URL 조회
    folder = storage() / 'uploads' / f'af_{gene.lower()}'
    folder.mkdir(parents=True, exist_ok=True)
    cif_path = folder / f'AF-{accession}.cif'
    if not cif_path.is_file():
        try:
            af_api = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
            with urllib.request.urlopen(af_api, timeout=15) as resp:
                af_data = json.loads(resp.read())
            cif_url = af_data[0]["cifUrl"]
            urllib.request.urlretrieve(cif_url, str(cif_path))
        except Exception:
            raise ValueError(f"{gene}: AlphaFold 구조 다운로드 실패 (accession={accession}). PDB/mmCIF를 직접 등록하세요.")

    row = dict(id=f'af_{gene.lower()}', gene=gene, path=str(cif_path),
               label=f'{gene} · AlphaFold predicted structure',
               kind='protein', role='target',
               source=f'AlphaFold DB (UniProt {accession})',
               source_url=f'https://alphafold.ebi.ac.uk/entry/{accession}',
               confidence_type='pLDDT', chains=['A'])
    write_json(folder / 'asset.json', row)
    return row


def pdb_text(row, chain=None):
    from Bio.PDB import PDBIO, MMCIFParser, PDBParser, Select
    path = Path(row['path'])
    opener = gzip.open if path.suffix == '.gz' else open
    parser = MMCIFParser(QUIET=True) if '.cif' in path.suffixes else PDBParser(QUIET=True)
    with opener(path, 'rt') as stream:
        model = parser.get_structure('target', stream)[0]
    if chain and chain not in model:
        raise ValueError('구조에 없는 chain입니다: ' + chain)
    class Selection(Select):
        def accept_chain(self, item):
            return chain is None or item.id == chain
    # mmCIF permits longer component IDs (e.g. AF3 LIG_B). PDB's residue
    # name occupies three columns; longer IDs shift every coordinate field.
    for residue in model.get_residues():
        if len(residue.resname) > 3:
            residue.resname = residue.resname[:3]
    output = io.StringIO()
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(output, Selection())
    return output.getvalue()


def structure_info(gene):
    row = asset(gene=gene.upper())
    from .structure_io import protein_residues
    residues = protein_residues(row['path'])
    return {**{k: v for k, v in row.items() if k != 'path'}, 'residues': len(residues),
            'chains': sorted({r[0] for r in residues}), 'viewer_url': '/structure/' + gene.upper() + '/3d',
            'confidence_status': 'NOT_APPLICABLE' if row['source'].startswith('PDB ') else 'NOT_LOADED'}


def viewer_document(row, *, surface=False, second=None, ligand=None, title=None):
    models = [(Path(row['path']).read_text(encoding='utf-8'), 'sdf')] if row['kind'] == 'ligand' else [(pdb_text(row), 'pdb')]
    if second:
        lines = []
        for line in pdb_text(second).splitlines(True):
            if line.startswith(('ATOM', 'HETATM')):
                line = line[:30] + f'{float(line[30:38]) + 70:8.3f}' + line[38:]
            lines.append(line)
        models.append((''.join(lines), 'pdb'))
    if ligand:
        models.append((Path(ligand).read_text(encoding='utf-8'), 'sdf'))
    data = json.dumps(models).replace('<', '\\u003c')
    label = title or row['label']
    note = ('Independent structures placed side by side; not a predicted ternary complex.' if second else row['source'])
    confidence = row.get('confidence_type') == 'pLDDT' and second is None
    controls = '<label>Color: <select id="colorBy"><option value="chain">Chain</option>'
    if confidence:
        controls += '<option value="plddt" selected>Confidence (pLDDT)</option>'
    controls += '</select></label>'
    if row.get('source_url', '').startswith('https://'):
        controls += ' · <a target="_blank" rel="noopener" href="' + html.escape(row['source_url'], quote=True) + '">Source</a>'
    if confidence:
        controls += '<span id="confidenceLegend">pLDDT: <i style="color:#0053d6">≥90</i> · <i style="color:#65cbf3">70–90</i> · <i style="color:#ffdb13">50–70</i> · <i style="color:#ff7d45">&lt;50</i></span>'
    return '''<!doctype html><html><head><meta charset="utf-8"><script src="/assets/3Dmol-min.js" onerror="this.onerror=null;this.src='https://3dmol.org/build/3Dmol-min.js'"></script>
<style>body{margin:0;background:#0d1117;color:#eee;font:14px system-ui;display:flex;flex-direction:column;height:100vh}header{padding:12px;line-height:1.6}header a{color:#65cbf3}#confidenceLegend{margin-left:16px}#confidenceLegend i{font-style:normal}#viewer{position:relative;flex:1;min-height:0;width:100%}</style></head>
<body><header>''' + html.escape(label) + '<br><small>' + html.escape(note) + '</small><br>' + controls + '''</header><div id="viewer"></div>
<script>''' + (PACKAGE / 'assets/structure-fallback.js').read_text(encoding='utf-8') + '''</script><script>window.viewerReady=false;window.viewerError=null;const models=''' + data + ''';
try{const viewer=$3Dmol.createViewer('viewer',{backgroundColor:'#0d1117'});
models.forEach(([text,format])=>viewer.addModel(text,format));
const colorBy=document.getElementById('colorBy');
function applyColor(){
 const confidence=colorBy.value==='plddt';
 const colorfunc=atom=>!Number.isFinite(atom.b)?'#999999':atom.b>=90?'#0053d6':atom.b>=70?'#65cbf3':atom.b>=50?'#ffdb13':'#ff7d45';
 models.forEach((model,i)=>{
  const color=confidence&&i===0?{colorfunc}:{colorscheme:'chain'};
  viewer.setStyle({model:i},{cartoon:{...color},stick:{radius:0.13,...color}});
 });
 const legend=document.getElementById('confidenceLegend');if(legend)legend.hidden=!confidence;
 viewer.render();
}
colorBy.addEventListener('change',applyColor);applyColor();
''' + ("viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.45,color:'#69b8df'},{model:0});" if surface else '') + '''
viewer.zoomTo();viewer.render();window.viewerMode='webgl';window.viewerReady=true;window.structureViewer=viewer;
}catch(e){document.getElementById('colorBy').disabled=true;const legend=document.getElementById('confidenceLegend');if(legend)legend.hidden=true;try{showStructureFallback(models,e);}catch(f){window.viewerError=String(f);document.querySelector('header').textContent='3D viewer error: '+f.message;}}</script></body></html>'''


def model_config(modality):
    config = copy.deepcopy(load_json(MODELS)[modality])
    cache_dir = str(Path(get_config().data.base_dir).resolve() / 'cache/af3_target_msa')
    if modality == 'DE_NOVO_BINDER':
        config.setdefault('af3_binder', {}).setdefault('target_msa_cache_dir', cache_dir)
    elif modality == 'DEGRADER':
        config.setdefault('af3', {}).setdefault('target_msa_cache_dir', cache_dir)
    elif modality == 'SMALL_MOLECULE':
        config.setdefault('small_molecule', {}).setdefault('af3', {}).setdefault('target_msa_cache_dir', cache_dir)
    return config


def candidate_data(parent_id, gene, modality):
    from .routes import MODALITY_ROUTES
    if modality not in MODALITY_ROUTES:
        raise ValueError('지원하지 않는 모달리티입니다.')
    if parent_id == 'reference_erbb2':
        if gene != 'ERBB2':
            raise ValueError('기존 입력은 ERBB2에만 해당합니다.')
        return None
    folder = run_directory(Path(get_config().data.base_dir), run_id=parent_id)
    state = json.loads((folder / 'pipeline_state.json').read_text())
    target = (state.get('qualification_results') or {}).get(gene)
    if not target:
        raise ValueError('해당 실행의 후보가 아닙니다.')
    route = next((r for r in target['tier_results'] if r['modality'] == modality), None)
    if not route or route['Tier'] not in {'1A', '1B', '2A', '2B'}:
        raise ValueError('이 후보의 해당 모달리티는 현재 설계 가능한 Tier가 아닙니다.')
    if target.get('safety_veto', {}).get(modality, {}).get('has_veto'):
        raise ValueError('해당 모달리티의 safety veto가 적용됐습니다.')
    return target


def choose_gpu():
    ready = available_devices()
    if ready:
        return ready[0].logical
    raise ValueError('GPU 4–7 중 사용 가능한 GPU가 없습니다.')


def resume_binder_input(body, identifier, qualification):
    """Resume verified local backbone/FASTA artifacts; AF3 outputs go to a new run."""
    from .protein_mpnn_tool import ProteinMPNNTool
    folder = run_directory(Path(get_config().data.base_dir), run_id=identifier)
    if (folder / 'worker_input.json').is_file():
        job = json.loads((folder / 'worker_input.json').read_text())
    else:
        job = {key: json.loads((folder / filename).read_text()) for key, filename in
               [('body', 'request.json'), ('data', 'input.json'), ('config', 'config.json')]}
    original = job['body']
    if (body['gene'], body['modality'], body['parent_run_id']) != (original['gene'], original['modality'], original['parent_run_id']):
        raise ValueError('Resume source must match the target, modality and parent analysis')
    data, config = copy.deepcopy(job['data']), copy.deepcopy(job['config'])
    if body['modality'] != 'DE_NOVO_BINDER':
        raise ValueError('Only binder candidate validation can be resumed')
    previous = []
    design = data['protein_design_input']
    for record in sorted((folder / 'execution/binder').rglob('protein_mpnn_run.json')):
        saved = json.loads(record.read_text())
        if saved.get('status') != 'success' or saved.get('is_mock'):
            continue
        backbone_id = record.parents[1].name
        match = re.fullmatch(r'target_\d+_candidate_(\d+)', backbone_id)
        if not match:
            continue
        backbones = list((record.parents[2] / 'rfd3').glob('*_model_' + str(int(match[1])) + '.cif.gz'))
        if len(backbones) != 1 or not backbones[0].resolve().is_relative_to(folder.resolve()):
            continue
        for sequence in saved.get('sequences', []):
            if sequence.get('is_mock') or not Path(sequence.get('sequence_path', '')).resolve().is_relative_to(folder.resolve()):
                continue
            if not ProteinMPNNTool.reusable(sequence, design['design_length']):
                continue
            previous.append({'gene_name': body['gene'], 'candidate_id': sequence['candidate_id'],
                             'rfd3_candidate_id': backbone_id, 'rfd3_structure_path': str(backbones[0]),
                             'binder_chain': design['binder_chain'], 'is_mock': False,
                             'protein_mpnn': {**sequence, 'provenance': saved.get('provenance', {})}})
    if not previous:
        raise ValueError('No completed, non-mock backbone/sequence artifacts to resume')
    data['protein_design_result'] = {'rfd3': {'status': 'success', 'metadata': {'is_mock': False}},
                                     'candidates': previous, 'is_mock': False, 'source_run_id': identifier}
    data['qualification'] = copy.deepcopy(qualification)
    return data, config, qualification


def design_input(body):
    gene, modality = body['gene'].upper(), body['modality']
    qualification = candidate_data(body['parent_run_id'], gene, modality)
    config = model_config(modality)
    supplied = body.get('input', {})
    if isinstance(supplied, dict) and 'resume_run_id' in supplied:
        if set(supplied) != {'resume_run_id'} or not isinstance(supplied['resume_run_id'], str):
            raise ValueError('Resume accepts only a source run ID')
        return resume_binder_input(body, supplied['resume_run_id'], qualification)
    allowed = {'asset_id', 'target_chain', 'antigen_chain', 'hotspot_residues', 'design_length',
               'ligand_smiles', 'docking_box', 'warhead_smiles', 'e3_ligand_smiles', 'linker_smiles',
               'warhead_map', 'e3_map', 'linker_maps',
               'expected_full_smiles', 'sequence_redesign', 'binder_chain', 'num_sequences', 'num_designs', 'design_seed'}
    if not isinstance(supplied, dict) or set(supplied) - allowed:
        raise ValueError('지원하지 않는 설계 입력입니다. 실행 경로·명령어는 서버 설정을 사용합니다.')
    if 'sequence_redesign' in supplied and type(supplied['sequence_redesign']) is not bool:
        raise ValueError('sequence_redesign must be boolean')
    if any(key in supplied for key in ('sequence_redesign', 'binder_chain', 'num_sequences')) and modality != 'DE_NOVO_BINDER':
        raise ValueError('Binder redesign options require DE_NOVO_BINDER')
    for key, lower, upper in [('num_designs', 1, 8), ('design_seed', 0, 2147483647)]:
        if key in supplied:
            if modality != 'DE_NOVO_BINDER' or type(supplied[key]) is not int or not lower <= supplied[key] <= upper:
                raise ValueError(f'{key} requires a binder integer in {lower}–{upper}')
            if key == 'design_seed':
                config['rfd3']['seed'] = supplied[key]
                config['protein_mpnn']['seed'] = supplied[key]
    if 'num_sequences' in supplied:
        count = supplied['num_sequences']
        if type(count) is not int or not 1 <= count <= 8:
            raise ValueError('num_sequences must be 1–8')
        config.setdefault('protein_mpnn', {})['num_sequences_per_backbone'] = count
        for key in ('max_mpnn_candidates', 'max_af3_candidates'):
            config.setdefault('screening_options', {})[key] = count
    if supplied.get('sequence_redesign'):
        row = asset(supplied.get('asset_id'), gene)
        from .binder_analysis import protein_chains
        chains = protein_chains(row['path'])
        target_chain, binder_chain = supplied.get('target_chain', 'A'), supplied.get('binder_chain', 'B')
        if target_chain == binder_chain or target_chain not in chains or binder_chain not in chains:
            raise ValueError('Select distinct target and binder protein chains in the supplied complex')
        if len(chains[target_chain]['sequence']) > 1200 or not 5 <= len(chains[binder_chain]['sequence']) <= 150:
            raise ValueError('Redesign supports targets up to 1200 residues and binders of 5–150 residues')
        config['af3_binder']['target_msa_path'] = ''
        return {'advance_targets':[{'gene_name':gene,'modality':modality}], 'protein_design_input':{
            **config['rfd3'], 'enabled':True, 'dry_run':False, 'target_structure':row['path'],
            'backbone_paths':[row['path']], 'target_chain':target_chain, 'binder_chain':binder_chain,
            'target_sequence':chains[target_chain]['sequence'], 'design_length':len(chains[binder_chain]['sequence']),
            'num_designs':1, 'hotspot_residues':supplied.get('hotspot_residues',[])}}, config, qualification
    if body['parent_run_id'] == 'reference_erbb2' and not supplied:
        folder = {'SMALL_MOLECULE': 'small_molecule', 'DE_NOVO_BINDER': 'binder', 'ADC': 'adc'}.get(modality)
        if folder:
            data = json.loads((REFERENCE / folder / 'input.json').read_text())
            if modality == 'SMALL_MOLECULE':
                config['small_molecule']['admet_policy'] = 'valid_input'
            if modality == 'DE_NOVO_BINDER':
                data['protein_design_input'].pop('output_dir', None)
                target = PACKAGE / 'examples/binder/target_3RCD_chain_A.pdb'
                from .structure_io import protein_residues
                mapping = json.loads((PACKAGE / 'examples/binder/reference_numbering.json').read_text())
                observed = ''.join(aa for chain, _, aa in protein_residues(str(target)) if chain == 'A')
                if observed != mapping['target_sequence'] or observed != data['protein_design_input']['target_sequence']:
                    raise ValueError('ERBB2 reference sequence differs from verified residue mapping')
                data['protein_design_input']['hotspot_residues'] = [mapping['residue_mapping'][h]
                    for h in data['protein_design_input']['hotspot_residues']]
                data['protein_design_input']['target_structure'] = str(target)
                data['reference_numbering'] = mapping
        else:
            from .modality_dispatch import select_degrader_reference
            data = select_degrader_reference({'target': gene})
    elif modality == 'DEGRADER':
        data = {k: v for k, v in supplied.items() if k not in {'asset_id', 'target_chain'}}
        data['input_mode'] = 'component_assembly'
    else:
        row = asset(supplied.get('asset_id'), gene)
        chain = supplied.get('target_chain', 'A')
        if not isinstance(chain, str) or not re.fullmatch(r'[A-Za-z0-9]{1,4}', chain):
            raise ValueError('Invalid target chain')
        if modality == 'ADC':
            from .structure_io import protein_residues
            antigen = supplied.get('antigen_chain', chain)
            if antigen not in {r[0] for r in protein_residues(row['path'])}:
                raise ValueError('Antigen chain is absent (mmCIF uses label_asym_id)')
            data = {'structure_path': row['path'], 'antigen_chain': antigen}
        else:
            from Bio.Data.PDBData import protein_letters_3to1
            from Bio.PDB import PDBParser
            protein = PDBParser(QUIET=True).get_structure('target', io.StringIO(pdb_text(row, chain)))[0]
            sequence = ''.join(protein_letters_3to1.get(residue.resname, '') for residue in protein[chain])
            if not sequence or len(sequence) > 1200:
                raise ValueError('선택 chain은 1–1200개의 표준 아미노산이어야 합니다. 큰 표적은 domain 구조를 사용하세요.')
        if modality == 'DE_NOVO_BINDER':
            length = supplied.get('design_length', 100)
            if type(length) is not int or not 30 <= length <= 150:
                raise ValueError('Binder 길이는 30–150 범위여야 합니다.')
            config['af3_binder']['target_msa_path'] = ''
            data = {'advance_targets': [{'gene_name': gene, 'modality': modality}], 'protein_design_input': {
                **config['rfd3'], 'enabled': True, 'dry_run': False, 'target_structure': row['path'],
                'target_chain': chain, 'target_sequence': sequence, 'hotspot_residues': supplied.get('hotspot_residues', []),
                'num_designs': supplied.get('num_designs', 1), 'design_length': length}}
            if 'num_designs' in supplied:
                backbones = supplied['num_designs']
                sequences = supplied.get('num_sequences', config['protein_mpnn'].get('num_sequences_per_backbone', 1))
                options = config.setdefault('screening_options', {})
                options.update(max_rfd3_candidates=backbones, max_mpnn_candidates=backbones * sequences,
                               max_af3_candidates=backbones * sequences)
        elif modality == 'SMALL_MOLECULE':
            smiles = supplied.get('ligand_smiles')
            box = supplied.get('docking_box')
            if not smiles or not isinstance(box, dict) or set(box) != {'center', 'size'}:
                raise ValueError('Small molecule은 ligand_smiles와 docking_box(center/size)가 필요합니다.')
            for key, values in box.items():
                if not isinstance(values, list) or len(values) != 3 or any(type(v) not in {int, float} or not -1000 < v < 1000 for v in values):
                    raise ValueError('docking_box는 유한한 숫자 3개씩 필요합니다.')
                if key == 'size' and any(not 1 <= v <= 40 for v in values):
                    raise ValueError('도킹 box 크기는 각 축 1–40 Å입니다.')
            prepared = storage() / 'prepared' / (row['id'] + '_' + re.sub('[^A-Za-z0-9]', '', chain) + '.pdb')
            prepared.parent.mkdir(parents=True, exist_ok=True)
            prepared.write_text(''.join(line for line in pdb_text(row, chain).splitlines(True) if line.startswith(('ATOM', 'TER', 'END'))))
            config['small_molecule']['af3']['target_msa_path'] = ''
            data = {'candidates': [{'candidate_id': gene + '_ligand', 'smiles': smiles}], 'context': {
                'gene_name': gene, 'target_sequence': sequence, 'prepared_receptor_path': str(prepared),
                'target_chain': chain, 'docking_box': box}}
    if qualification is not None:
        data['qualification'] = copy.deepcopy(qualification)
    return data, config, qualification


def start_design(body):
    from . import server
    if isinstance(body, dict) and body.get('modality') == 'ALL':
        return start_erbb2_batch(body)
    if not isinstance(body, dict) or set(body) - {'parent_run_id', 'gene', 'modality', 'input', 'reuse_results'}:
        raise ValueError('Invalid design request')
    if not all(isinstance(body.get(k), str) for k in ['parent_run_id', 'gene', 'modality']):
        raise ValueError('parent_run_id, gene, modality가 필요합니다.')
    if server.run_state['running']:
        raise HTTPException(409, '다른 분석이 실행 중입니다.')
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,40}', body['gene']):
        raise ValueError('Invalid gene')
    if type(body.get('reuse_results', True)) is not bool:
        raise ValueError('reuse_results must be boolean')
    body = {**body, 'gene': body['gene'].upper()}
    data, config, qualification = design_input(body)
    from .modality_dispatch import readiness
    ready = readiness(body['modality'], data, config)
    if ready['blockers']:
        raise ValueError('; '.join(ready['blockers']))
    question = body['gene'] + ' · ' + body['modality'] + ' design'
    if server.begin_pipeline_run(question, execution_kind='design') is None:
        raise HTTPException(409, '다른 분석이 실행 중입니다.')
    identifier = 'design_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_') + uuid.uuid4().hex[:8]
    with server.run_lock:
        server.run_state['report_run_id'] = identifier
    cfg = get_config().model_copy(deep=True)
    from .design_execution import run_bounded_design
    try:
        server.start_worker(run_bounded_design, identifier, body, data, config, qualification, cfg, finalize=True)
    except RuntimeError:
        with server.run_lock:
            server.run_state['running'] = False
        raise
    return {'status': 'started', 'report_run_id': identifier, 'report_url': '/report?run_id=' + identifier}


from .resource_usage import measured_design


@measured_design
def run_design(identifier, body, data, config, qualification, app_config, *, finalize=True):
    from . import server
    from .configuration import configuration_scope
    from .modality_dispatch import execute
    from .user_selection import create_selection_event
    start = time.monotonic()
    directory = Path(app_config.data.base_dir) / identifier
    terminal = {'type': 'done', 'report_run_id': identifier, 'execution_status': 'FAILED', 'errors': []}
    def publish(kind, node, status, messages):
        event = {'type': kind, 'node': node, 'execution_status': status,
                 'elapsed': round(time.monotonic() - start, 1), 'messages': messages}
        if config.get('_event_path'):
            with Path(config['_event_path']).open('a') as stream:
                stream.write(json.dumps(event, ensure_ascii=False)+'\n')
        else:
            with server.run_lock:
                server.run_state['history'].append(event)
                server.event_queue.append(event)
    from .result_reuse import find_design, reused_design, remember_design, reuse_scope, resume_partial, route_key
    reuse_config = copy.deepcopy(config)
    reuse_data = copy.deepcopy(data)
    record = None
    try:
        record = RunRecord(directory, body['modality'])
        record.data.update(execution_profile='therapeutic_design', parent_run_id=body['parent_run_id'],
                           gene=body['gene'], batch_run_id=body.get('batch_run_id'), running=True, execution_status='PARTIAL')
        selection = create_selection_event(body['gene'], body['modality'], route_id=body['gene'] + ':' + body['modality'],
                                           tier=next(r['Tier'] for r in qualification['tier_results'] if r['modality'] == body['modality']) if qualification else 'UNKNOWN',
                                           selected_by='user', rationale=('Run Pipeline request: execute all eligible routes' if body.get('selection_scope') == 'all_eligible' else 'Direct dashboard design request'))
        if body.get('selection_scope'):
            selection['selection_scope'] = body['selection_scope']
        record.data['selection_event'] = selection
        record.event('therapeutic_design', 'RUNNING')
        publish('progress', 'therapeutic_design', 'RUNNING', [body['gene'] + ':' + body['modality'] + ' · starting'])
        write_json(directory / 'request.json', body)
        record.data['input_origin'] = 'prior_ERBB2_reference' if qualification is None else 'qualified_upstream_candidate'
        write_json(directory / 'input.json', data)
        write_json(directory / 'reuse_request.json', {'key': route_key(body, reuse_data, reuse_config)})
        def progress(stage, status, candidate_id):
            record.event(stage, status, candidate_id=candidate_id)
            publish('progress' if status == 'RUNNING' else 'node', body['modality'].lower() + '/' + stage,
                    status, [str(candidate_id or body['gene']) + ': ' + status])
        def snapshot(output):
            record.data.update(summary=copy.deepcopy(output.get('summary', {})),
                               execution_status=output.get('execution_status', 'PARTIAL'),
                               execution_scope=output.get('execution_scope'),
                               missing_steps=output.get('missing_steps', []),
                               blockers=output.get('blockers', []), report_status='LIVE', running=True)
            known = {a['path']: a for a in record.data['artifacts']}
            for artifact in output.get('artifacts', []):
                path = Path(artifact['path']).resolve()
                if not path.is_file() or not path.is_relative_to((directory / 'execution').resolve()):
                    continue
                if artifact['path'] in known:
                    known[artifact['path']].update(artifact)
                else:
                    record.data['artifacts'].append(dict(artifact))
            write_json(directory / 'run_summary.json', record.data)
            server.write_competition_report(record.data, directory)
        config.update(output_dir=str(directory / 'execution'), dry_run=False,
                      _progress_callback=progress, _snapshot_callback=snapshot)
        cached = find_design(app_config.data.base_dir, body, reuse_data, reuse_config)
        resumed = None if cached else resume_partial(app_config.data.base_dir, body, reuse_data, reuse_config, qualification)
        if resumed:
            data, record.data['resume'] = resumed
            publish('node', 'result_reuse', 'COMPLETED', ['완료된 후보 재사용 · ' + record.data['resume']['source_run_id']])
        needs_gpu = not cached and (body['modality'] in {'SMALL_MOLECULE', 'DE_NOVO_BINDER'} or (
            body['modality'] == 'DEGRADER' and config.get('run_ternary', False)))
        if needs_gpu:
            publish('progress', 'therapeutic_design', 'QUEUED',
                    [body['gene'] + ':' + body['modality'] + ' · waiting for GPU 4–7'])
        with reserve_gpu(timeout=config.get('gpu_wait_timeout_seconds', 120)) if needs_gpu else nullcontext(None) as gpu:
            if gpu is not None:
                config.setdefault('screening_options', {})['gpu_device'] = gpu.logical
                if body['modality'] == 'SMALL_MOLECULE':
                    config.setdefault('small_molecule', {}).setdefault('docking', {})['device'] = gpu.logical
                record.data['gpu_device'] = gpu.physical
                publish('progress', 'therapeutic_design', 'RUNNING',
                        [body['gene'] + ':' + body['modality'] + f' · GPU {gpu.physical}'])
            write_json(directory / 'config.json', {k: v for k, v in config.items() if not k.startswith('_')})
            with configuration_scope(app_config), reuse_scope(app_config.data.base_dir, body.get('reuse_results', True), progress):
                if cached:
                    from .resource_usage import record_usage
                    record_usage('cache_hit', tool=body['modality'], source=cached['source_run_id'])
                    result = reused_design(cached, reuse_config)
                    record.data['reuse'] = result['reuse']
                    record.event('result_reuse', 'COMPLETED', **result['reuse'])
                    publish('node', 'result_reuse', 'COMPLETED',
                            ['기존 결과 재사용 · ' + cached['source_run_id']])
                else:
                    result = execute({'modality': body['modality'], 'input': data, 'selection_event': selection}, config)
        if app_config.agent_logging and not cached:
            try:
                from .agent_narration import narrate
                note_config = app_config.model_copy(deep=True)
                note_config.data.base_dir = str(directory)
                latest_stages = {row['stage']: row for row in result.get('stages', [])}
                evidence = {'gene': body['gene'], 'modality': body['modality'],
                            'execution_status': result.get('execution_status'),
                            'validation_decision': result.get('validation_decision'),
                            'critic': {k: result.get('critic', {}).get(k) for k in ('decision', 'reasons')},
                            'validation': result.get('validation', {}),
                            'stages': [{k: row.get(k) for k in ('stage', 'execution_status')} for row in latest_stages.values()],
                            'blockers': result.get('blockers', []), 'is_mock': result.get('is_mock', False)}
                with configuration_scope(note_config):
                    note = narrate({}, 'design', evidence=evidence)
                result.setdefault('summary', {})['agent_notes'] = note.get('agent_notes', {})
                publish('node', 'design', result['execution_status'], note.get('messages', []))
            except Exception as exc:
                result.setdefault("summary", {})["narration_warning"] = str(exc)
        record.data['stages'].extend(result.get('stages', []))
        known_artifacts = {a['path']: a for a in record.data['artifacts']}
        for artifact in result.get('artifacts', []):
            if artifact['path'] in known_artifacts:
                known_artifacts[artifact['path']].update(artifact)
            else:
                record.data['artifacts'].append(artifact)
        result['artifacts'] = record.data['artifacts']
        record.data['running'] = False
        record.finish({**result, 'run_id': identifier, 'stages': record.data['stages'],
                       'execution_profile': 'therapeutic_design', 'parent_run_id': body['parent_run_id'],
                       'execution_mode': 'REUSED_RESULT' if cached else 'RESUMED_FROM_ARTIFACTS' if resumed else 'FRESH_ANALYSIS', 'running': False,
                       'workflow_status': 'FINISHED', 'report_status': 'READY'})
        try:
            from .design_figures import write_design_figures
            record.data['artifacts'].extend(write_design_figures(record.data, directory))
        except Exception as exc:
            record.data.setdefault('report_warnings', []).append('Figure generation: ' + str(exc))
        server.write_competition_report(record.data, directory)
        register_demo_artifact(record, directory / 'report.html', body['modality'], 'report', artifact_type='report', execution_mode='REPORT_RENDER')
        write_json(directory / 'run_summary.json', record.data)
        remember_design(app_config.data.base_dir, body, reuse_data, reuse_config, record.data)
        publish('node', 'therapeutic_design', result['execution_status'], result.get('blockers', []))
        publish('node', 'critic_node', 'COMPLETED', [result.get('critic', {}).get('decision', 'NOT_EVALUATED')])
        terminal.update(execution_status=result['execution_status'], scientific_validation_status=result['validation_decision'],
                        errors=result.get('blockers', []), execution_mode=record.data['execution_mode'],
                        reuse=record.data.get('reuse'), node_path='therapeutic_design → critic_node')
    except Exception as exc:
        terminal.update(type='error', message=str(exc), errors=[str(exc)])
        if record:
            record.data.update(running=False, execution_status='FAILED', blockers=[str(exc)],
                               workflow_status='FINISHED', report_status='READY')
            record.event('therapeutic_design', 'FAILED', message=str(exc))
            server.write_competition_report(record.data, directory)
            terminal.update(type='done', report_status='READY', workflow_status='FINISHED')
    finally:
        terminal['total_time'] = round(time.monotonic() - start, 1)
        if finalize:
            with server.run_lock:
                server.run_state.update(running=False, result=terminal)
                server.run_state['history'].append(terminal)
                server.event_queue.append(terminal)
    return terminal



ERBB2_MODALITIES = ('SMALL_MOLECULE', 'DE_NOVO_BINDER', 'ADC', 'DEGRADER')


def start_erbb2_batch(body, *, research_question=None):
    from . import server
    if (set(body) - {'parent_run_id', 'gene', 'modality', 'input', 'reuse_results'}
            or body.get('parent_run_id') != 'reference_erbb2'
            or body.get('gene') != 'ERBB2' or body.get('input', {}) != {}):
        raise ValueError('전체 실행은 저장된 ERBB2 입력을 사용합니다.')
    if type(body.get('reuse_results', True)) is not bool:
        raise ValueError('reuse_results must be boolean')
    cfg = get_config().model_copy(deep=True)
    if server.begin_pipeline_run('ERBB2 · 네 모달리티 전체 실행', execution_kind='design', research_question=research_question) is None:
        raise HTTPException(409, '다른 분석이 실행 중입니다.')
    identifier = 'design_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_') + uuid.uuid4().hex[:8]
    with server.run_lock:
        server.run_state['report_run_id'] = identifier
    try:
        server.start_worker(run_erbb2_batch, identifier, body, cfg)
    except RuntimeError:
        with server.run_lock:
            server.run_state['running'] = False
        raise
    return {'status': 'started', 'report_run_id': identifier, 'report_url': '/report?run_id=' + identifier}


def run_erbb2_batch(identifier, body, app_config):
    from . import server
    from .configuration import configuration_scope
    from .modality_dispatch import readiness
    from .design_execution import run_bounded_design
    started = time.monotonic()
    directory = Path(app_config.data.base_dir) / identifier
    terminal = {'type': 'done', 'report_run_id': identifier, 'execution_status': 'FAILED'}
    record = None
    try:
        record = RunRecord(directory, 'ALL')
        rows = [{'gene':'ERBB2', 'modality':modality, 'run_id': 'design_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_') + uuid.uuid4().hex[:8],
                 'execution_status':'RUNNING', 'running':True, 'execution_mode':'FRESH_ANALYSIS',
                 'validation_status':'NOT_EVALUATED', 'critic_decision':'NOT_EVALUATED', 'blocker':''}
                for modality in ERBB2_MODALITIES]
        record.data.update(execution_profile='therapeutic_design', parent_run_id='reference_erbb2',
                           gene='ERBB2', running=True, execution_status='PARTIAL',
                           execution_mode='FRESH_ANALYSIS', summary={'modalities': rows})
        write_json(directory / 'request.json', body)
        record.event('modality_batch', 'RUNNING')
        server.write_competition_report(record.data, directory)
        with server.run_lock:
            server.run_state['erbb2_report_run_id'] = identifier
        def execute_modality(modality):
            row = dict(next(row for row in rows if row['modality'] == modality))
            row.update(execution_status='FAILED', is_mock=False)
            child_body = {**body, 'modality': modality, 'batch_run_id': identifier}
            child_id = row['run_id']
            row.update(run_id=child_id, report_url='/report?run_id=' + child_id)
            try:
                with configuration_scope(app_config):
                    data, config, qualification = design_input(child_body)
                    ready = readiness(modality, data, config)
                if ready['blockers']:
                    raise ValueError('; '.join(ready['blockers']))
                result = run_bounded_design(child_id, child_body, data, config, qualification, app_config, finalize=False)
                row.update(run_id=child_id, report_url='/report?run_id=' + child_id,
                           execution_mode=result.get('execution_mode', 'FRESH_ANALYSIS'), execution_status=result['execution_status'],
                           validation_status=result.get('scientific_validation_status', 'NOT_EVALUATED'),
                           blocker='; '.join(result.get('errors', [])))
                child_summary = Path(app_config.data.base_dir) / child_id / 'run_summary.json'
                if child_summary.is_file():
                    child = json.loads(child_summary.read_text())
                    row['critic_decision'] = child.get('critic', {}).get('decision', 'NOT_EVALUATED')
                    row['gpu_device'] = child.get('gpu_device')
                    row['missing_steps'] = child.get('missing_steps', [])
                    row['execution_scope'] = child.get('execution_scope')
            except (ValueError, FileNotFoundError) as exc:
                row.update(execution_status='BLOCKED', blocker=str(exc))
            except Exception as exc:
                row.update(execution_status='FAILED', blocker=str(exc))
            row['running'] = False
            child_path = Path(app_config.data.base_dir) / child_id
            if not (child_path / 'run_summary.json').exists():
                from .automatic_design import _blocked_report
                _blocked_report(child_path, 'reference_erbb2', 'ERBB2', modality,
                                row['execution_status'], row['blocker'], None)
            return row

        with ThreadPoolExecutor(max_workers=4, thread_name_prefix='erbb2-design') as pool:
            pending = {}
            for index, modality in enumerate(ERBB2_MODALITIES):
                record.event(modality, 'RUNNING')
                pending[pool.submit(execute_modality, modality)] = index
            for future in as_completed(pending):
                row = future.result()
                rows[pending[future]] = row
                modality = row['modality']
                record.event(modality, row['execution_status'], message=row['blocker'])
                server.write_competition_report(record.data, directory)
                event = {'type': 'node', 'node': 'therapeutic_design', 'elapsed': round(time.monotonic() - started, 1),
                         'execution_status': row['execution_status'],
                         'report_url': row['report_url'],
                         'messages': ['ERBB2:' + modality + ': ' + row['execution_status'], row['blocker']]}
                with server.run_lock:
                    server.run_state['history'].append(event)
                    server.event_queue.append(event)
        status = 'COMPLETED' if all(r['execution_status'] == 'COMPLETED' for r in rows) else 'PARTIAL'
        record.finish({'stage': 'modality_batch', 'execution_status': status, 'running': False,
                       'workflow_status':'FINISHED', 'report_status':'READY',
                       'validation_decision': 'NOT_EVALUATED',
                       'blockers': [r['modality'] + ': ' + r['blocker'] for r in rows if r['blocker']]})
        server.write_competition_report(record.data, directory)
        terminal.update(execution_status=status, scientific_validation_status='NOT_EVALUATED', erbb2_batch=True,
                        workflow_status='FINISHED', report_status='READY',
                        message='네 모달리티 실행 종료. 과학적 판정은 각 보고서를 확인하세요.')
        with server.run_lock:
            server.run_state['erbb2_report_run_id'] = identifier
    except Exception as exc:
        terminal.update(type='error', message=str(exc))
        if record:
            record.data.update(running=False, execution_status='FAILED', blockers=[str(exc)],
                               workflow_status='FINISHED', report_status='READY')
            record.event('modality_batch', 'FAILED', message=str(exc))
            server.write_competition_report(record.data, directory)
            terminal.update(type='done', workflow_status='FINISHED', report_status='READY')
    finally:
        terminal['total_time'] = round(time.monotonic() - started, 1)
        with server.run_lock:
            server.run_state.update(running=False, result=terminal)
            server.run_state['history'].append(terminal)
            server.event_queue.append(terminal)


def design_candidates():
    from .routes import MODALITY_ROUTES
    from .candidate_display import state_hidden_genes, visible_gene
    rows = []
    root = Path(get_config().data.base_dir)
    for path in sorted(root.glob('upstream_*/pipeline_state.json'), reverse=True)[:10]:
        state = json.loads(path.read_text())
        hidden = state_hidden_genes(state)
        for gene, target in (state.get('qualification_results') or {}).items():
            if not visible_gene(gene, hidden):
                continue
            for route in target.get('tier_results', []):
                modality = route['modality']
                if modality in MODALITY_ROUTES and route['Tier'] in {'1A', '1B', '2A', '2B'} and not target.get('safety_veto', {}).get(modality, {}).get('has_veto'):
                    rows.append(dict(parent_run_id=path.parent.name, gene=gene, modality=modality,
                                     tier=route['Tier'], origin='upstream_analysis'))
    rows.extend(dict(parent_run_id='reference_erbb2', gene='ERBB2', modality=modality,
                     tier='UNKNOWN', origin='prior_reference_inputs') for modality in MODALITY_ROUTES)
    return rows


def lab_directory(sid):
    if not re.fullmatch(r'[0-9a-f]{32}', sid):
        raise ValueError('Invalid session')
    return storage() / 'sessions' / sid


def lab_status(sid):
    path = lab_directory(sid) / 'status.json'
    if not path.is_file():
        raise HTTPException(404, 'Session not found')
    return json.loads(path.read_text(encoding='utf-8'))


def start_lab(mode, body):
    from . import server
    allowed = {'gene', 'pdb_id', 'target_pdb', 'e3_pdb', 'ligand'}
    if not isinstance(body, dict) or set(body) - allowed:
        raise ValueError('Invalid Structure Lab input')
    gene = body.get('gene', 'ESR1')
    if not isinstance(gene, str):
        raise ValueError('Invalid gene')
    gene = gene.upper()

    # 등록된 구조 검색 → 없으면 AlphaFold에서 자동 다운로드
    _KNOWN_ASSETS = {'ESR1': 'esr1_3ert', 'ERBB2': 'erbb2_3rcd'}
    _KNOWN_PDBS = {'ESR1': '3ERT', 'ERBB2': '3RCD'}
    if gene in _KNOWN_ASSETS:
        target = asset(_KNOWN_ASSETS[gene], gene)
        pdb_id = _KNOWN_PDBS[gene]
    else:
        try:
            target = asset(gene=gene)
        except ValueError:
            # AlphaFold 구조 자동 다운로드
            target = _fetch_alphafold_structure(gene)
        pdb_id = body.get('pdb_id', body.get('target_pdb', 'AF-' + gene))
    if gene in _KNOWN_PDBS and any(body.get(k, _KNOWN_PDBS[gene]) != _KNOWN_PDBS[gene] for k in ['pdb_id', 'target_pdb']):
        raise ValueError('PDB does not match selected gene')
    if body.get('e3_pdb', '4TZ4') not in {'4TZ4', '1VCB'}:
        raise ValueError('Supported E3 structures: CRBN 4TZ4 / VHL 1VCB')
    sid = uuid.uuid4().hex
    directory = lab_directory(sid)
    if mode == 'sm':
        from .lab_docking import LIGANDS
        if gene in _KNOWN_PDBS:
            allowed_ligs = set(LIGANDS) if gene == 'ESR1' else {'tak285'}
            if body.get('ligand', 'tamoxifen') not in allowed_ligs:
                raise ValueError('Unsupported ligand for the selected Lab reference')
        if server.begin_pipeline_run(gene + ' Structure Lab docking', execution_kind='lab') is None:
            raise HTTPException(409, '다른 분석이 실행 중입니다.')
        directory.mkdir(parents=True)
        write_json(directory / 'status.json', dict(session=sid, status='running'))
        app_config = get_config().model_copy(deep=True)
        try:
            server.start_worker(run_lab_docking, sid, target, body, app_config)
        except Exception:
            with server.run_lock:
                server.run_state['running'] = False
            raise
        return dict(session=sid, status='running', poll_url='/api/lab/status/' + sid)
    directory.mkdir(parents=True)
    second = asset('crbn_4tz4' if body.get('e3_pdb', '4TZ4') == '4TZ4' else 'vhl_1vcb') if mode == 'protac' else None
    content = viewer_document(target, surface=mode == 'adc', second=second)
    (directory / 'viewer.html').write_text(content, encoding='utf-8')
    message = 'Independent target/E3 structures; not a ternary prediction.' if second else 'Protein surface view; ADC eligibility and conjugation are not evaluated here.'
    result = dict(session=sid, status='done', viewer_url='/api/lab/viewer/' + sid, message=message,
                  scientific_validation_status='NOT_EVALUATED', source=target['source'])
    write_json(directory / 'status.json', result)
    return result


def run_lab_docking(sid, target, body, app_config):
    from Bio.PDB import PDBParser

    from . import server
    from .configuration import configuration_scope
    from .docking_backend import DockingRequest, docking_backend, run_docking
    from .lab_docking import LIGANDS
    from .ligand_preparation import prepare_ligand, validate_candidate
    from .small_molecule_config import DockingConfig
    start = time.monotonic()
    with configuration_scope(app_config):
        directory = lab_directory(sid)
    identifier = 'lab_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ_') + sid[:8]
    folder = Path(app_config.data.base_dir) / identifier
    with server.run_lock:
        server.run_state['report_run_id'] = identifier
        event = dict(type='progress', node='therapeutic_design', elapsed=0, messages=['CPU docking'])
        server.run_state['history'].append(event)
    terminal = dict(type='done', execution_status='COMPLETED', scientific_validation_status='NOT_EVALUATED', report_run_id=identifier)
    record = None
    try:
        record = RunRecord(folder, 'SMALL_MOLECULE')
        record.data.update(execution_profile='structure_lab', running=True, execution_status='PARTIAL')
        record.event('docking', 'RUNNING')
        gene = target['gene']
        ligand_name = body.get('ligand', 'tamoxifen')
        if gene == 'ERBB2':
            reference = json.loads((REFERENCE / 'small_molecule/input.json').read_text(encoding='utf-8'))
            candidate = reference['candidates'][0]
            context = reference['context']
            receptor_text = Path(context['prepared_receptor_path']).read_text(encoding='utf-8')
            center, size = context['docking_box']['center'], context['docking_box']['size']
        else:
            candidate = {'candidate_id': ligand_name, 'smiles': LIGANDS[ligand_name]}
            raw = pdb_text(target, 'A')
            protein = PDBParser(QUIET=True).get_structure('target', io.StringIO(raw))[0]
            ligands = [res for res in protein['A'] if res.id[0].startswith('H_') and len(list(res.get_atoms())) > 12]
            if not ligands:
                raise ValueError('No crystallographic ligand available to define the docking box')
            native = max(ligands, key=lambda res: len(list(res.get_atoms())))
            atoms = [a.coord for a in native if a.element != 'H']
            center = [float(sum(a[i] for a in atoms) / len(atoms)) for i in range(3)]
            size = [22., 22., 22.]
            receptor_text = ''.join(line for line in raw.splitlines(True) if line.startswith(('ATOM', 'TER', 'END')))
        receptor = folder / 'receptor.pdb'
        receptor.write_text(receptor_text, encoding='utf-8')
        validation = validate_candidate(candidate)
        ligand = prepare_ligand(candidate, validation, folder / 'ligand_preparation')
        if ligand['status'] != 'success':
            raise ValueError(ligand.get('error_message') or 'Ligand preparation rejected: ' + ', '.join(validation.get('warning', [])))
        options = model_config('SMALL_MOLECULE')['small_molecule']['docking']
        options.update(use_gpu=False, cnn_scoring='none', timeout_seconds=600)
        config = DockingConfig.model_validate(options)
        request = DockingRequest(candidate['candidate_id'], str(receptor), ligand['prepared_ligand_path'],
                                 tuple(center), tuple(size), candidate['smiles'], str(folder / 'docking'))
        result = run_docking(docking_backend(config), request)
        if result['status'] != 'success':
            raise ValueError(result.get('error_message') or result.get('reason') or result['status'])
        poses = [pose for pose in result['poses'] if pose['docking_status'] == 'success']
        pose = min(poses, key=lambda p: p['docking_score'])
        pose_path = next(pose[k] for k in ['pose_path', 'docked_pose_path', 'ligand_pose_path'] if pose.get(k))
        viewer = viewer_document({**target, 'path': str(receptor)}, ligand=pose_path, title=gene + ' · ' + ligand_name + ' docking')
        (directory / 'viewer.html').write_text(viewer, encoding='utf-8')
        (folder / 'viewer.html').write_text(viewer, encoding='utf-8')
        warning = 'GNINA with Vina scoring (CPU); protonation as supplied. No AF3, ADMET, affinity or efficacy validation.'
        result.update(message=warning, binding_site_center=center, binding_site_size=size, receptor_preparation='Chain selection and heteroatom removal; no pH assignment or missing-atom repair')
        write_json(folder / 'docking_result.json', result)
        record.finish(dict(execution_status='COMPLETED', validation_decision='NOT_EVALUATED', stage='docking',
                           summary=result, blockers=[], is_mock=False, running=False))
        for path in folder.rglob('*'):
            if path.is_file() and path.name not in {'run_summary.json', 'events.jsonl'}:
                register_demo_artifact(record, path, 'SMALL_MOLECULE', 'docking', execution_mode='FRESH_ANALYSIS')
        server.write_competition_report(record.data, folder)
        register_demo_artifact(record, folder / 'report.html', 'SMALL_MOLECULE', 'report', execution_mode='REPORT_RENDER')
        write_json(folder / 'run_summary.json', record.data)
        status = dict(session=sid, status='done', viewer_url='/api/lab/viewer/' + sid, report_url='/report?run_id=' + identifier,
                      best_affinity_kcal_mol=pose['docking_score'], message=warning, scientific_validation_status='NOT_EVALUATED')
    except Exception as exc:
        status = dict(session=sid, status='error', error=str(exc), report_url='/report?run_id=' + identifier)
        terminal.update(type='error', execution_status='FAILED', message=str(exc))
        if record:
            record.finish(dict(execution_status='FAILED', stage='docking', blockers=[str(exc)], running=False))
    finally:
        try:
            write_json(directory / 'status.json', status)
        except OSError as exc:
            terminal.update(type='error', execution_status='FAILED', message=str(exc))
        terminal['total_time'] = round(time.monotonic() - start, 1)
        with server.run_lock:
            server.run_state.update(running=False, result=terminal)
            server.run_state['history'].append(terminal)
            server.event_queue.append(terminal)

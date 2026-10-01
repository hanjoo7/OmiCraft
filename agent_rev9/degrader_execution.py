"""Validate and assemble explicitly supplied degrader components."""

from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

VERSION = rdBase.rdkitVersion


def component_qc(component_id, component_type, smiles, allow_dummy=False):
    row = dict(
        component_id=component_id,
        component_type=component_type,
        name=component_id,
        raw_smiles=smiles,
        canonical_isomeric_smiles=None,
        formal_charge=None,
        stereo_status="not_assessed",
        protonation_status="preserved_input_not_enumerated",
        molecular_weight=None,
        source="user_input",
        version=VERSION,
        execution_status="not_run",
        warnings=[],
        error_message=None,
    )
    if not smiles:
        row["error_message"] = "MISSING_COMPONENT"
        return row
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        row.update(execution_status="invalid_input", error_message="SMILES_PARSE_OR_VALENCE_FAILED")
        return row
    maps = [a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()]
    if len(maps) != len(set(maps)):
        row.update(execution_status="invalid_input", error_message="DUPLICATE_ATOM_MAP")
        return row
    supported = {1, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 19, 20, 35, 53} | (
        {0} if allow_dummy else set()
    )
    if any(a.GetAtomicNum() not in supported for a in mol.GetAtoms()):
        row.update(execution_status="invalid_input", error_message="UNSUPPORTED_ATOM")
        return row
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    fragments = len(Chem.GetMolFrags(mol))
    row.update(
        canonical_isomeric_smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
        formal_charge=Chem.GetFormalCharge(mol),
        molecular_weight=Descriptors.MolWt(mol),
        stereo_status="unspecified_centers"
        if any(c == "?" for _, c in centers)
        else "specified"
        if centers
        else "no_tetrahedral_centers",
        fragment_count=fragments,
        execution_status="success",
        explicit_hydrogen_count=sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()),
    )
    if fragments > 1:
        row["warnings"].append("MULTI_FRAGMENT_OR_SALT_PRESERVED")
    if row["stereo_status"] == "unspecified_centers":
        row["warnings"].append("STEREOCHEMISTRY_UNSPECIFIED")
    return row


def mapped_atom(mol, number):
    if number is None:
        raise ValueError("MISSING_ATTACHMENT_ATOM_MAP")
    matches = [a for a in mol.GetAtoms() if a.GetAtomMapNum() == number]
    if len(matches) != 1:
        raise ValueError("MISSING_OR_DUPLICATE_ATTACHMENT_ATOM_MAP")
    atom = matches[0]
    if atom.GetAtomicNum() == 0:
        if atom.GetDegree() != 1 or atom.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE:
            raise ValueError("ATTACHMENT_DUMMY_MUST_HAVE_ONE_SINGLE_BOND")
    elif atom.GetTotalNumHs() < 1:
        raise ValueError("ATTACHMENT_VALENCE_UNAVAILABLE")
    return atom


def assemble(warhead, e3_ligand, linker, warhead_map, e3_map, linker_maps):
    if not linker_maps or len(set(linker_maps)) != 2:
        raise ValueError("MISSING_OR_DUPLICATE_LINKER_ATTACHMENT_MAP")
    inputs = [(warhead, [warhead_map]), (e3_ligand, [e3_map]), (linker, list(linker_maps))]
    molecules, ports, provenance = [], [], {}
    for component, (smiles, maps) in enumerate(inputs):
        qc = component_qc(str(component), "assembly_component", smiles, allow_dummy=True)
        if qc["execution_status"] != "success" or qc["fragment_count"] != 1:
            raise ValueError(qc["error_message"] or "DISCONNECTED_COMPONENT")
        mol = Chem.MolFromSmiles(smiles)
        selected = [mapped_atom(mol, n).GetIdx() for n in maps]
        port = []
        for atom in mol.GetAtoms():
            original = atom.GetAtomMapNum()
            new = (component + 1) * 10000 + atom.GetIdx() + 1
            provenance[str(new)] = {
                "component": ["warhead", "e3_ligand", "linker"][component],
                "original_atom_map": original,
                "original_atom_index": atom.GetIdx(),
                "isotope": atom.GetIsotope(),
                "formal_charge": atom.GetFormalCharge(),
                "chiral_tag": int(atom.GetChiralTag()),
            }
            atom.SetAtomMapNum(new)
        for index in selected:
            atom = mol.GetAtomWithIdx(index)
            endpoint = atom.GetNeighbors()[0] if atom.GetAtomicNum() == 0 else atom
            port.append(
                (
                    endpoint.GetAtomMapNum(),
                    atom.GetAtomMapNum() if atom.GetAtomicNum() == 0 else None,
                )
            )
        molecules.append(mol)
        ports.append(port)
    combined = Chem.RWMol(
        Chem.CombineMols(Chem.CombineMols(molecules[0], molecules[1]), molecules[2])
    )
    joins = [(ports[0][0], ports[2][0]), (ports[1][0], ports[2][1])]
    for left, right in joins:
        by_map = {a.GetAtomMapNum(): a for a in combined.GetAtoms()}
        for number, dummy in (left, right):
            atom = by_map[number]
            if dummy is None and atom.GetNumExplicitHs():
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
        combined.AddBond(by_map[left[0]].GetIdx(), by_map[right[0]].GetIdx(), Chem.BondType.SINGLE)
    removed = {dummy for pair in joins for _, dummy in pair if dummy is not None}
    for index in sorted(
        [a.GetIdx() for a in combined.GetAtoms() if a.GetAtomMapNum() in removed], reverse=True
    ):
        combined.RemoveAtom(index)
    mol = combined.GetMol()
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("ASSEMBLY_DISCONNECTED_OR_UNRESOLVED_DUMMY")
    for atom in mol.GetAtoms():
        before = provenance[str(atom.GetAtomMapNum())]
        if (int(atom.GetChiralTag()), atom.GetFormalCharge(), atom.GetIsotope()) != (
            before["chiral_tag"],
            before["formal_charge"],
            before["isotope"],
        ):
            raise ValueError("ASSEMBLY_STEREO_CHARGE_OR_ISOTOPE_CHANGED")
    mapped = Chem.MolToSmiles(mol, isomericSmiles=True)
    plain = Chem.Mol(mol)
    for atom in plain.GetAtoms():
        atom.SetAtomMapNum(0)
    return dict(
        assembly_status="success",
        canonical_isomeric_smiles=Chem.MolToSmiles(plain, isomericSmiles=True),
        mapped_smiles=mapped,
        atom_map_provenance=provenance,
        endpoint_maps=[ports[0][0][0], ports[1][0][0]],
        formal_charge=Chem.GetFormalCharge(mol),
        molecular_weight=Descriptors.MolWt(mol),
        stereochemistry_status="preserved",
        component_smiles={"warhead": warhead, "e3_ligand": e3_ligand, "linker": linker},
    )


from .resource_usage import measured_tool


@measured_tool
def execute_assembly(input_data, config):
    from .modality_dispatch import result

    required = COMPONENT_FIELDS + ('expected_full_smiles',)
    missing = [k for k in required if input_data.get(k) is None]
    if missing:
        return result('DEGRADER', 'BLOCKED', 'input_validation', blockers=['missing:' + k for k in missing])
    product = assemble(*(input_data[k] for k in COMPONENT_FIELDS))
    expected = Chem.MolFromSmiles(input_data['expected_full_smiles'])
    actual = Chem.MolFromSmiles(product['canonical_isomeric_smiles'])
    if expected is None or not same_molecule(actual, expected):
        return result('DEGRADER', 'FAILED', 'assembly_identity',
                      summary={'assembly_status': 'FAILED', 'expected_identity_match': False},
                      blockers=['ASSEMBLED_MOLECULE_DIFFERS_FROM_EXPECTED'])
    product['expected_identity_match'] = True
    Path(config['output_dir']).mkdir(parents=True, exist_ok=True)
    heavy = Chem.MolFromSmiles(product['mapped_smiles'])
    mol = Chem.AddHs(heavy)
    product.update(atom_count_before_hydrogens=heavy.GetNumAtoms(),
                   atom_count_after_hydrogens=mol.GetNumAtoms(), fragment_count=len(Chem.GetMolFrags(heavy)),
                   sanitization_status='success', components=[component_qc(name, name, smiles, allow_dummy=True)
                   for name, smiles in product['component_smiles'].items()])
    params = AllChem.ETKDGv3()
    params.randomSeed = int(config.get('seed', 11))
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError('degrader_embedding_failed')
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError('degrader_mmff_parameters_missing')
    converged = AllChem.MMFFOptimizeMolecule(mol, maxIters=int(config.get('max_iterations', 500))) == 0
    output = Path(config['output_dir']) / 'degrader.sdf'
    with Chem.SDWriter(str(output)) as writer:
        writer.write(mol)
    by_map = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms()}
    a, b = (mol.GetConformer().GetAtomPosition(by_map[m]) for m in product['endpoint_maps'])
    product.update(endpoint_distance=float((a-b).Length()), mmff_converged=converged,
                   protonation_status='preserve_input_not_ph_predicted', backend='RDKit',
                   backend_version=VERSION, scope='user_components_assembly_and_geometry')
    import numpy as np

    topology = Chem.GetDistanceMatrix(mol)
    positions = mol.GetConformer().GetPositions()
    distances = np.linalg.norm(positions[:, None, :] - positions[None, :, :], axis=2)
    heavy_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    nonbonded = [float(distances[i, j]) for i in heavy_indices for j in heavy_indices if i < j and topology[i, j] > 2]
    product.update(nonbonded_heavy_atom_clashes=sum(d < 1.5 for d in nonbonded),
                   minimum_nonbonded_heavy_atom_distance=min(nonbonded, default=None),
                   clash_distance_angstrom=1.5, geometry_scope='free_conformer_not_ternary_complex')
    blockers = ['ternary_complex_not_evaluated', 'degradation_efficacy_not_evaluated']
    if not converged:
        blockers.append('mmff_not_converged')
    return result('DEGRADER', 'PARTIAL', 'geometry', summary=product,
                  artifacts=[{'path': str(output.resolve()), 'label': 'Assembled degrader'}],
                  blockers=blockers)


COMPONENT_FIELDS = ('warhead_smiles', 'e3_ligand_smiles', 'linker_smiles',
                    'warhead_map', 'e3_map', 'linker_maps')


def same_molecule(left, right):
    if left is None or right is None:
        return False
    left, right = Chem.Mol(left), Chem.Mol(right)
    for mol in (left, right):
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    return (left.GetNumAtoms() == right.GetNumAtoms()
            and left.GetNumBonds() == right.GetNumBonds()
            and left.HasSubstructMatch(right, useChirality=True)
            and right.HasSubstructMatch(left, useChirality=True)
            and Chem.MolToSmiles(left, isomericSmiles=True) == Chem.MolToSmiles(right, isomericSmiles=True)
            and Chem.MolToInchiKey(left) == Chem.MolToInchiKey(right))


def input_blockers(data):
    mode = data.get('input_mode', 'component_assembly')
    if mode == 'component_assembly':
        return ['missing:' + k for k in COMPONENT_FIELDS + ('expected_full_smiles',)
                if data.get(k) is None or data.get(k) == '']
    if mode != 'published_full_molecule':
        return ['unknown_degrader_input_mode']
    fields = ('candidate_id', 'full_degrader_smiles', 'candidate_origin', 'target', 'e3_ligase',
              'primary_doi', 'primary_pmid', 'source_url', 'source_sha256',
              'source_accessed_at', 'structure_source', 'curation_status')
    errors = ['missing:' + k for k in fields if not data.get(k)]
    if data.get('candidate_origin') != 'published_reference' or data.get('novel_design') is not False:
        errors.append('published_reference_identity_required')
    if data.get('curation_status') != 'VERIFIED':
        errors.append('CURATION_BLOCKED')
    checksum = data.get('source_sha256', '')
    if len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum.lower()):
        errors.append('source_checksum_required')
    if data.get('is_mock'):
        errors.append('mock_reference_not_allowed')
    if not data.get('structure_crosscheck'):
        errors.append('structure_crosscheck_required')
    return errors


def full_molecule_qc(data):
    mol = Chem.MolFromSmiles(data['full_degrader_smiles'])
    if mol is None:
        raise ValueError('FULL_SMILES_PARSE_OR_SANITIZATION_FAILED')
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError('FULL_MOLECULE_MUST_BE_CONNECTED_WITHOUT_DUMMY_ATOMS')
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if any(c == '?' for _, c in centers):
        raise ValueError('REFERENCE_STEREOCHEMISTRY_UNSPECIFIED')
    row = dict(original_smiles=data['full_degrader_smiles'], sanitization_status='success',
               canonical_smiles=Chem.MolToSmiles(mol, isomericSmiles=False),
               isomeric_smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
               inchi=Chem.MolToInchi(mol), inchikey=Chem.MolToInchiKey(mol),
               molecular_formula=rdMolDescriptors.CalcMolFormula(mol), molecular_weight=Descriptors.MolWt(mol),
               exact_mass=Descriptors.ExactMolWt(mol), formal_charge=Chem.GetFormalCharge(mol),
               fragment_count=len(Chem.GetMolFrags(mol)), stereocenter_count=len(centers),
               stereocenters=centers, hbd=Descriptors.NumHDonors(mol), hba=Descriptors.NumHAcceptors(mol),
               rotatable_bonds=Descriptors.NumRotatableBonds(mol), tpsa=Descriptors.TPSA(mol),
               logp=Descriptors.MolLogP(mol), ring_count=Descriptors.RingCount(mol), backend_version=VERSION)
    cross = data['structure_crosscheck']
    for key in ('inchikey', 'molecular_formula', 'stereocenter_count'):
        if key not in cross or cross[key] != row[key]:
            raise ValueError('SOURCE_STRUCTURE_MISMATCH:' + key)
    if abs(float(cross['molecular_weight']) - row['molecular_weight']) > 0.1:
        raise ValueError('SOURCE_STRUCTURE_MISMATCH:molecular_weight')
    if not same_molecule(mol, Chem.MolFromSmiles(cross['isomeric_smiles'])):
        raise ValueError('SOURCE_STRUCTURE_MISMATCH:isomeric_smiles')
    beyond = row['molecular_weight'] > 500 or row['logp'] > 5 or row['hbd'] > 5 or row['hba'] > 10
    row.update(full_molecule_validation='COMPLETED', property_profile='beyond_rule_of_five' if beyond else 'within_rule_of_five',
               development_warning=beyond, warning_message=['BEYOND_RULE_OF_FIVE_NOT_AN_EFFICACY_FAILURE'] if beyond else [])
    return mol, row


def ternary_readiness(data, config):
    from .af3_ligand_tool import AF3LigandTool
    from .small_molecule_config import LigandAF3Config
    from .small_molecule_io import sha256

    model = LigandAF3Config(**config.get('af3', {}))
    dependency = AF3LigandTool(model).dependency_status()
    checks = []
    for name in ('ERBB2', 'VHL', 'ElonginB', 'ElonginC'):
        entry = config.get('ternary_sequences', {}).get(name, {})
        path = Path(entry.get('path') or '')
        available = path.is_file() and bool(entry.get('source') and entry.get('sequence'))
        if available:
            from Bio.PDB import MMCIFParser
            from Bio.SeqUtils import seq1
            try:
                chain = MMCIFParser(QUIET=True).get_structure(name, str(path))[0][entry['source_chain']]
                observed = ''.join(seq1(res.resname) for res in chain if res.id[0] == ' ')
                available = observed == entry['sequence']
            except (OSError, ValueError, KeyError):
                available = False
        checks.append({'entity': name, 'available': available, 'chain': entry.get('chain'),
                       'source': entry.get('source'), 'path': str(path) if available else None,
                       'sha256': sha256(path) if available else None,
                       'sequence_length': len(entry.get('sequence', '')) if available else None})
    blockers = []
    if not all(c['available'] for c in checks if c['entity'] in {'ERBB2', 'VHL'}):
        blockers.append('BLOCKED_SEQUENCE_MISSING')
    if not model.model_dir or not any(Path(model.model_dir).glob('af3.bin*')):
        blockers.append('BLOCKED_WEIGHTS_MISSING')
    if not model.database_dir or not Path(model.database_dir).is_dir() or not any(Path(model.database_dir).iterdir()):
        blockers.append('BLOCKED_DATABASE_MISSING')
    if dependency != 'available' and dependency not in {'weights_missing', 'database_missing'}:
        blockers.append('BLOCKED_AF3_RUNTIME_MISSING')
    return dict(status=blockers[0] if blockers else 'READY', blockers=blockers, sequences=checks,
                adapter='AF3TernaryTool.prepare_input', dependency_status=dependency,
                weights_present='BLOCKED_WEIGHTS_MISSING' not in blockers,
                database_directory_present='BLOCKED_DATABASE_MISSING' not in blockers,
                database_verification_scope='existing local directory inventory; data pipeline not rerun',
                adapter_supported=True, input_json_generation=False,
                reason='Local target and E3 sequence provenance checked; AF3 target/E3/ligand adapter configured.',
                required_entities=['ERBB2 kinase domain', 'VHL (Elongin B/C optional)', 'SJF1528'],
                ligand_smiles_valid=Chem.MolFromSmiles(data['full_degrader_smiles']) is not None,
                chain_definition_status='READY' if not blockers else 'BLOCKED',
                prediction_status='NOT_RUN', physical_gpu_id=None, allowed_physical_gpus=[7, 4, 5, 6],
                model_calls=0, retries=0, degradation_efficacy='NOT_EVALUATED')


def prepare_reference(data, config):
    from .modality_dispatch import result
    from .small_molecule_io import write_json

    directory = Path(config['output_dir'])
    directory.mkdir(parents=True, exist_ok=True)
    mol, qc = full_molecule_qc(data)
    provenance = {k: data.get(k) for k in ('candidate_id', 'display_name', 'candidate_origin', 'novel_design',
        'primary_doi', 'primary_pmid', 'source_url', 'source_accessed_at', 'source_sha256', 'structure_source',
        'curation_status', 'curation_message', 'structure_crosscheck', 'source_records')}
    provenance.update(original_smiles=data['full_degrader_smiles'], source_file_sha256=data['source_sha256'])
    # Component ports are independent evidence; a whole-molecule record does not establish them.
    components = dict(component_mapping_status='NOT_AVAILABLE', assembly_status='NOT_APPLICABLE_REFERENCE_INPUT',
                      assembly_performed=False, decomposition='NOT VERIFIED',
                      reason='Verified full molecule supplied. No curated component-port evidence is available.',
                      fields={k: data.get(k) for k in COMPONENT_FIELDS})
    hydrogens = Chem.AddHs(mol)
    qc.update(hydrogen_count_before=sum(a.GetAtomicNum() == 1 for a in mol.GetAtoms()),
              hydrogen_count_after=sum(a.GetAtomicNum() == 1 for a in hydrogens.GetAtoms()),
              conformer_generation_status='NOT_RUN', optimization_status='NOT_RUN', force_field=None,
              energy_if_available=None, protonation_status='preserved_input_not_ph_predicted')
    params = AllChem.ETKDGv3()
    params.randomSeed = int(config.get('seed', 11))
    params.numThreads = 1
    if AllChem.EmbedMolecule(hydrogens, params) != 0:
        qc['conformer_generation_status'] = 'FAILED'
    else:
        qc['conformer_generation_status'] = 'COMPLETED'
        if AllChem.MMFFHasAllMoleculeParams(hydrogens):
            qc['force_field'] = 'MMFF94'
            force = AllChem.MMFFGetMoleculeForceField(hydrogens, AllChem.MMFFGetMoleculeProperties(hydrogens))
        elif AllChem.UFFHasAllMoleculeParams(hydrogens):
            qc['force_field'] = 'UFF'
            qc['warning_message'].append('UFF_FALLBACK_MMFF_UNSUPPORTED')
            force = AllChem.UFFGetMoleculeForceField(hydrogens)
        else:
            force = None
            qc['warning_message'].append('NO_SUPPORTED_FORCE_FIELD')
        if force is not None:
            code = force.Minimize(maxIts=int(config.get('max_iterations', 2000)))
            qc['optimization_status'] = 'CONVERGED' if code == 0 else 'NOT_CONVERGED'
            qc['energy_if_available'] = float(force.CalcEnergy())
            qc['energy_unit'] = 'kcal/mol; force-field energy, not binding affinity'
            if code != 0:
                qc['warning_message'].append('OPTIMIZATION_NOT_CONVERGED')
        hydrogens.SetProp('_Name', data['candidate_id'])
        with Chem.SDWriter(str(directory / 'conformer.sdf')) as writer:
            writer.write(hydrogens)
        restored = next(iter(Chem.SDMolSupplier(str(directory / 'conformer.sdf'), removeHs=True)))
        qc['conformer_identity_preserved'] = same_molecule(mol, restored)
        if not qc['conformer_identity_preserved']:
            qc['conformer_generation_status'] = 'FAILED'
    readiness = ternary_readiness(data, config)
    prediction = {'status': 'not_run', 'samples': []}
    callback = config.get('_progress_callback', lambda *args: None)
    snapshot = config.get('_snapshot_callback', lambda *args: None)
    for name, value in [('input', data), ('provenance', provenance), ('rdkit_validation', qc),
                        ('component_status', components), ('ternary_readiness', readiness)]:
        write_json(directory / (name + '.json'), value)
    import csv
    with (directory / 'molecular_properties.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['property', 'value'])
        writer.writerows((k, v) for k, v in qc.items() if isinstance(v, (str, int, float, bool)))
    figure_warnings = []
    try:
        reference_figures(mol, hydrogens if qc['conformer_generation_status'] == 'COMPLETED' else None, directory)
    except Exception as exc:
        figure_warnings.append('Figure generation: ' + str(exc))

    def current_output():
        from .degrader_status import ternary_execution_status
        predicted = prediction['status'] == 'success'
        status = 'COMPLETED' if qc['conformer_generation_status'] == 'COMPLETED' and predicted else 'PARTIAL'
        missing = ([] if predicted else ['ternary_prediction']) + ['experimental_degradation']
        summary = dict(candidate_id=data['candidate_id'], candidate_origin='published_reference', novel_design=False,
                       input_mode='published_full_molecule', target_context='ERBB2/EGFR', e3_ligase=data['e3_ligase'],
                       full_molecule_validation=qc['full_molecule_validation'], **components, qc=qc, provenance=provenance,
                       ternary_readiness=readiness, ternary_prediction=prediction,
                       ternary_prediction_performed=bool(prediction.get('command') or prediction.get('inference_started')),
                       erbb2_specific_novel_design='NOT_CLAIMED', degradation_efficacy='NOT_EVALUATED',
                       demo_route_status='COMPLETE' if status == 'COMPLETED' else 'PARTIAL', modality_execution_status=status,
                       execution_scope='published_reference_chemistry_and_ternary_structure',
                       report_warnings=figure_warnings)
        summary['ternary_status'] = ternary_execution_status(summary)
        summary['execution_checklist'] = [
            {'stage': 'chemical_preparation', 'execution_status': qc['conformer_generation_status']},
            {'stage': 'ternary_prediction', 'execution_status': summary['ternary_status']},
            {'stage': 'experimental_degradation', 'execution_status': 'NOT_EVALUATED'}]
        artifacts = [{'path': str(p.resolve()), 'label': str(p.relative_to(directory))} for p in sorted(directory.rglob('*'))
                     if p.is_file() and p.name not in {'events.jsonl', 'run_summary.json', 'live_result.json'}]
        blockers = readiness['blockers'] + ['degradation_efficacy_not_evaluated']
        if not predicted and prediction['status'] != 'running':
            blockers.append('ternary_prediction_incomplete: ' + (prediction.get('error_message') or prediction['status']))
        output = result('DEGRADER', status, 'reference_evaluation', summary=summary, artifacts=artifacts, blockers=blockers)
        output.update(missing_steps=missing, execution_scope=summary['execution_scope'])
        return output

    def checkpoint():
        write_json(directory / 'ternary_prediction.json', prediction)
        write_json(directory / 'ternary_readiness.json', readiness)
        output = current_output()
        write_json(directory / 'live_result.json', output)
        snapshot(output)
        return output

    checkpoint()
    callback('chemical_preparation', qc['conformer_generation_status'], data['candidate_id'])
    if not readiness['blockers']:
        from .af3_ternary_tool import AF3TernaryTool
        from .small_molecule_config import LigandAF3Config
        entries = config['ternary_sequences']
        model = LigandAF3Config(**config.get('af3', {}))
        context = {'target_sequence': entries['ERBB2']['sequence'],
                   'partner_proteins': [entries[name] for name in ('VHL', 'ElonginB', 'ElonginC')
                       if any(row['entity'] == name and row['available'] for row in readiness['sequences'])]}
        candidate = {'candidate_id': data['candidate_id'], 'canonical_smiles': qc['isomeric_smiles']}
        tool = AF3TernaryTool(model, gpu_device=config.get('screening_options', {}).get('gpu_device', 0))
        enabled = bool(config.get('run_ternary')) and not config.get('dry_run', False)
        def phase_progress(stage, status, details):
            prediction.update(status='running', stage_details=details,
                              protein_msa_cache=details.get('chains', []),
                              inference_started=bool(details.get('inference_calls')))
            checkpoint()
            callback(stage, status, data['candidate_id'])
        tool.progress_callback = phase_progress
        if enabled:
            prediction['status'] = 'running'
            checkpoint()
            callback('ternary_prediction', 'RUNNING', data['candidate_id'])
        try:
            prediction = tool.run(candidate=candidate, context=context, output_dir=directory / 'ternary', dry_run=not enabled)
        except Exception as exc:
            prediction = {'status':'backend_failed', 'error_message':str(exc), 'samples':[],
                          'stage_details':prediction.get('stage_details', {}), 'phase_status':'FAILED'}
        readiness.update(input_json_generation=bool(prediction.get('input_json')),
                         prediction_status=prediction['status'], model_calls=int(bool(prediction.get('command'))),
                         data_pipeline_calls=prediction.get('stage_details', {}).get('data_pipeline_calls', 0))
        output = checkpoint()
        if enabled:
            callback('ternary_prediction', output['summary']['ternary_status'], data['candidate_id'])
    return checkpoint()


def reference_figures(mol, conformer, directory):
    from rdkit.Chem import Draw
    Draw.MolToFile(mol, str(directory / 'degrader_2d.png'), size=(1600, 850),
                   legend='SJF1528 | Published HER2/EGFR degrader | Full molecule; component maps unverified')
    if conformer is None:
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    heavy = Chem.RemoveHs(conformer)
    xyz = heavy.GetConformer().GetPositions()
    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111, projection='3d')
    for bond in heavy.GetBonds():
        pair = xyz[[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]]
        ax.plot(*pair.T, color='#9baaba', lw=2)
    colors = {'C': '#465b71', 'N': '#397fd0', 'O': '#d95454', 'S': '#b79524', 'F': '#59a34f', 'Cl': '#318a56'}
    ax.scatter(*xyz.T, c=[colors.get(a.GetSymbol(), '#888888') for a in heavy.GetAtoms()], s=45, depthshade=True)
    span = max(xyz.max(0) - xyz.min(0)) / 2 + 1
    center = xyz.mean(0)
    for setter, value in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), center):
        setter(value - span, value + span)
    ax.set_box_aspect((1, 1, 1))
    ax.set_axis_off()
    ax.set_title('SJF1528 · RDKit free conformer\nElement colors; no component assignment or ternary pose', fontsize=13)
    fig.savefig(directory / 'degrader_3d.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


def execute(input_data, config):
    from .modality_dispatch import result
    errors = input_blockers(input_data)
    if errors:
        return result('DEGRADER', 'BLOCKED', 'input_validation', blockers=errors)
    if input_data.get('input_mode', 'component_assembly') == 'component_assembly':
        return execute_assembly(input_data, config)
    try:
        return prepare_reference(input_data, config)
    except (ValueError, RuntimeError, KeyError) as exc:
        return result('DEGRADER', 'FAILED', 'reference_evaluation',
                      summary={'full_molecule_validation': 'FAILED', 'original_smiles': input_data['full_degrader_smiles'],
                               'assembly_status': 'NOT_APPLICABLE_REFERENCE_INPUT'}, blockers=[str(exc)])

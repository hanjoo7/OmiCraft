"""IMGT mapping and contacts in an observed antibody–antigen complex."""

from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser
from Bio.SeqUtils import seq1

from .antibody_numbering import AnarciiBackend, LegacyAnarciBackend, cdr_region


def residues(chain):
    return [r for r in chain if r.id[0] == " " and seq1(r.resname) != "X"]


def map_sequence(sequence, observed, *, min_identity=0.90, min_coverage=0.80):
    coordinate_sequence = "".join(seq1(r.resname) for r in observed)
    if not sequence or not coordinate_sequence:
        return [], {"status": "failed", "reason": "empty_sequence"}
    aligner = PairwiseAligner(
        mode="global", match_score=2, mismatch_score=-4, open_gap_score=-6, extend_gap_score=-0.5
    )
    alignment = aligner.align(sequence, coordinate_sequence)[0]
    pairs = [(i, j) for a, b in zip(*alignment.aligned) for i, j in zip(range(*a), range(*b))]
    identity = (
        sum(sequence[i] == coordinate_sequence[j] for i, j in pairs) / len(pairs) if pairs else 0
    )
    coverage = len(pairs) / len(sequence)
    duplicate = len({r.id for r in observed}) != len(observed)
    missing = sorted(set(range(1, len(sequence) + 1)) - {i + 1 for i, _ in pairs})
    inserted = sorted(set(range(1, len(observed) + 1)) - {j + 1 for _, j in pairs})
    qc = {
        "status": "success"
        if identity >= min_identity and coverage >= min_coverage and not duplicate
        else "failed",
        "identity": identity,
        "coverage": coverage,
        "missing_sequence_indices": missing,
        "inserted_coordinate_indices": inserted,
        "duplicate_mapping": duplicate,
        "method": "Bio.Align global sequence alignment",
        "min_identity": min_identity,
        "min_coverage": min_coverage,
    }
    rows = [
        {
            "original_sequence_index": i + 1,
            "coordinate_sequence_index": j + 1,
            "chain_id": observed[j].parent.id,
            "structure_residue_id": observed[j].id[1],
            "insertion_code": observed[j].id[2].strip(),
            "residue_name": observed[j].resname,
            "amino_acid": coordinate_sequence[j],
            "sequence_match": sequence[i] == coordinate_sequence[j],
        }
        for i, j in pairs
    ]
    return rows, qc


def contact_pairs(antibody, antigen, cutoff=4.5, clash_cutoff=1.5):
    rows = []
    antigen_atoms = [(r, a) for r in antigen for a in r if a.element not in ("H", "D")]
    xyz = np.array([a.coord for _, a in antigen_atoms])
    for residue in antibody:
        atoms = [a for a in residue if a.element not in ("H", "D")]
        if not atoms or not len(xyz):
            continue
        distances = np.linalg.norm(np.array([a.coord for a in atoms])[:, None] - xyz[None], axis=2)
        groups = {}
        for i, j in zip(*np.where(distances <= cutoff)):
            other, _ = antigen_atoms[j]
            key = (other.parent.id, other.id)
            record = groups.setdefault(
                key,
                {
                    "antibody_chain": residue.parent.id,
                    "antibody_residue_id": residue.id[1],
                    "antibody_insertion_code": residue.id[2].strip(),
                    "antigen_chain": other.parent.id,
                    "antigen_residue_id": other.id[1],
                    "antigen_insertion_code": other.id[2].strip(),
                    "antigen_residue_name": other.resname,
                    "heavy_atom_contact_count": 0,
                    "minimum_distance": float("inf"),
                    "severe_clash_count": 0,
                },
            )
            record["heavy_atom_contact_count"] += 1
            record["minimum_distance"] = min(record["minimum_distance"], float(distances[i, j]))
            record["severe_clash_count"] += int(distances[i, j] < clash_cutoff)
        rows.extend(groups.values())
    return rows


from .resource_usage import measured_tool


@measured_tool
def analyze(input_data, config, backend=None):
    source = Path(input_data['structure_path'])
    parser = MMCIFParser(QUIET=True) if source.suffix.lower() in {'.cif', '.mmcif'} else PDBParser(QUIET=True)
    model = parser.get_structure('complex', str(source))[0]
    antigen_id = input_data['antigen_chain']
    if antigen_id not in model:
        raise ValueError('antigen_chain_missing')
    observed = {c.id: residues(c) for c in model if c.id != antigen_id and residues(c)}
    sequences = dict(input_data.get('antibody_sequences') or {})
    sequence_source = 'user_sequence'
    if not sequences and source.suffix.lower() in {'.cif', '.mmcif'}:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        meta = MMCIF2Dict(str(source))
        sequences = {cid.strip(): seq.replace('\n', '').replace(' ', '')
                     for ids, seq in zip(meta.get('_entity_poly.pdbx_strand_id', []),
                                         meta.get('_entity_poly.pdbx_seq_one_letter_code_can', []))
                     for cid in ids.split(',') if cid.strip() in observed}
        sequence_source = 'mmcif_polymer_sequence'
    if not sequences:
        raise ValueError('antibody_sequence_missing: supply full sequences keyed by chain ID')
    name = config.get('numbering_backend', 'anarcii').lower()
    if name not in {'anarcii', 'anarci'}:
        raise ValueError('unknown_numbering_backend')
    backend = backend or (LegacyAnarciBackend() if name == 'anarci'
                          else AnarciiBackend(config.get('numbering_python_path')))
    numbered = {cid: backend.number(cid, seq) for cid, seq in sequences.items()}
    result = {'numbering': {cid: v.to_dict() for cid, v in numbered.items()},
              'numbering_status': 'success', 'cdr_mapping_status': 'not_run',
              'sequence_source': sequence_source, 'mapping': [], 'mapping_qc': {},
              'contacts': [], 'blockers': [], 'is_mock': False,
              'structure_path': str(source.resolve())}
    if any(v.numbering_status == 'backend_unavailable' for v in numbered.values()):
        result.update(numbering_status='backend_unavailable', blockers=['numbering_backend_unavailable'])
        return result
    valid = {cid: v for cid, v in numbered.items() if v.numbering_status == 'success'}
    roles = [v.chain_type for v in valid.values()]
    if roles.count('H') != 1 or sum(x in {'K', 'L'} for x in roles) != 1:
        result.update(numbering_status='failed', blockers=['heavy_light_pair_unresolved'])
        return result
    for cid, value in valid.items():
        rows, qc = map_sequence(sequences[cid], observed.get(cid, []),
                                min_identity=config.get('min_identity', .90),
                                min_coverage=config.get('min_coverage', .80))
        result['mapping_qc'][cid] = qc
        positions = {p['sequence_index']: p for p in value.numbering}
        for row in rows:
            annotation = positions.get(row['original_sequence_index'])
            row.update(chain_role='heavy' if value.chain_type == 'H' else 'light',
                       imgt_position=annotation['imgt_position'] if annotation else None,
                       imgt_insertion_code=annotation['imgt_insertion_code'] if annotation else None,
                       region=cdr_region(annotation['imgt_position']) if annotation else 'unassigned')
        result['mapping'].extend(rows)
    if any(q['status'] != 'success' for q in result['mapping_qc'].values()):
        result['blockers'] = ['sequence_structure_mapping_failed']
        return result
    cutoff, clash = config.get('contact_cutoff', 4.5), config.get('clash_cutoff', 1.5)
    if not 0 < clash < cutoff:
        raise ValueError('invalid_contact_cutoffs')
    for cid in valid:
        lookup = {(r['structure_residue_id'], r['insertion_code']): r
                  for r in result['mapping'] if r['chain_id'] == cid}
        for contact in contact_pairs(observed[cid], residues(model[antigen_id]), cutoff, clash):
            row = lookup.get((contact['antibody_residue_id'], contact['antibody_insertion_code']))
            contact.update(region=row['region'] if row else 'unassigned')
            result['contacts'].append(contact)
    cdr = [r for r in result['contacts'] if r['region'].startswith('CDR')]
    result.update(cdr_mapping_status='success', contact_cutoff=cutoff, clash_cutoff=clash,
                  cdr_heavy_atom_contacts=sum(r['heavy_atom_contact_count'] for r in cdr),
                  cdr_contact_residues=len({(r['antibody_chain'], r['antibody_residue_id'],
                                            r['antibody_insertion_code']) for r in cdr}),
                  severe_clashes=sum(r['severe_clash_count'] for r in result['contacts']))
    return result


def conjugation_feasibility(analysis):
    """Rank mapped Lys/Cys for review without selecting a conjugation site."""
    from Bio.PDB.SASA import ShrakeRupley

    path = Path(analysis['structure_path'])
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {'.cif', '.mmcif'} else PDBParser(QUIET=True)
    model = parser.get_structure('adc', str(path))[0]
    ShrakeRupley(n_points=100).compute(model, level='R')
    mapped = {(r['chain_id'], r['structure_residue_id'], r['insertion_code']): r for r in analysis['mapping']}
    interface = {(r['antibody_chain'], r['antibody_residue_id'], r['antibody_insertion_code']) for r in analysis['contacts']}
    atoms, cdr_atoms, interface_atoms, sulfurs = [], [], [], []
    for chain in model:
        for residue in residues(chain):
            key = (chain.id, residue.id[1], residue.id[2].strip())
            heavy = [a for a in residue if a.element not in {'H', 'D'}]
            atoms.extend(heavy)
            if mapped.get(key, {}).get('region', '').startswith('CDR'):
                cdr_atoms.extend(heavy)
            if key in interface:
                interface_atoms.extend(heavy)
            if residue.resname == 'CYS' and 'SG' in residue:
                sulfurs.append(residue['SG'])
    rows = []
    for key, annotation in mapped.items():
        if annotation['residue_name'] not in {'LYS', 'CYS'}:
            continue
        residue = model[key[0]][(' ', key[1], key[2] or ' ')]
        atom_name = 'NZ' if residue.resname == 'LYS' else 'SG'
        if atom_name not in residue:
            continue
        anchor = residue[atom_name]
        def distance(targets):
            return min((float(np.linalg.norm(anchor.coord-a.coord)) for a in targets), default=None)
        cdr = annotation['region'].startswith('CDR')
        at_interface = key in interface
        disulfide = atom_name == 'SG' and any(a is not anchor and np.linalg.norm(anchor.coord-a.coord) <= 2.3 for a in sulfurs)
        reasons = (["inside_CDR"] if cdr else []) + (["antigen_interface"] if at_interface else []) + (["possible_disulfide"] if disulfide else [])
        rows.append({'chain': key[0], 'residue_id': key[1], 'insertion_code': key[2],
            'residue': residue.resname, 'anchor_atom': atom_name, 'region': annotation['region'],
            'inside_cdr': cdr, 'at_antigen_interface': at_interface, 'cdr_distance_angstrom': distance(cdr_atoms),
            'interface_distance_angstrom': distance(interface_atoms), 'sasa_angstrom2': float(residue.sasa),
            'exposure_metric': 'ShrakeRupley_complex_residue_SASA', 'possible_disulfide': disulfide,
            'exclusion_reasons': reasons, 'experimental_accessibility': 'NOT_EVALUATED', 'selected': False})
    rows.sort(key=lambda r: (bool(r['exclusion_reasons']), -r['sasa_angstrom2'], -(r['cdr_distance_angstrom'] or 0)))
    for rank, row in enumerate(rows, 1):
        row['review_rank'] = rank
    return {'status': 'COMPLETED', 'candidates': rows, 'site_selection': 'NOT_PERFORMED',
            'method': 'mapped Lys NZ/Cys SG; SASA and distances in observed complex',
            'policy': {'sasa_sphere_points': 100, 'possible_disulfide_distance_angstrom': 2.3,
                       'threshold_status': 'uncalibrated_demo_heuristics'},
            'experimental_accessibility': 'NOT_EVALUATED'}

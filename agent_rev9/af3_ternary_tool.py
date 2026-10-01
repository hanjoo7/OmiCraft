"""AF3 target–E3–ligand input using the shared runtime and output parser."""

import json
from pathlib import Path

from .af3_ligand_tool import AF3LigandTool, confidence_chain_ids
from .protein_design_route import AlphaFold3BinderAdapter, _read_af3_json
from .small_molecule_io import finite, write_json


from .resource_usage import measured_tool

class AF3TernaryTool(AF3LigandTool):
    @measured_tool
    def run(self, *, candidate, context, output_dir, dry_run=False, is_mock=False):
        if dry_run or is_mock:
            return super().run(candidate=candidate, context=context, output_dir=output_dir,
                               dry_run=dry_run, is_mock=is_mock)
        from .af3_staged import StagedAF3Runner
        return StagedAF3Runner(self, output_dir, getattr(self, 'progress_callback', None)).run(candidate, context)

    def prepare_input(self, candidate, context, directory):
        self.ligand_smiles = candidate['canonical_smiles']
        path = super().prepare_input(candidate, context, directory)
        payload = json.loads(path.read_text())
        entities = context.get('partner_proteins', [])
        if not entities:
            raise ValueError('At least one E3 protein is required')
        chains = {self.config.protein_chain_id, self.config.ligand_chain_id}
        if len(chains) != 2:
            raise ValueError('Duplicate target/ligand chain ID')
        for entity in entities:
            chain = entity['chain']
            if not isinstance(chain, str) or not chain.isalpha() or not chain.isupper() or chain in chains:
                raise ValueError('Invalid or duplicate partner chain ID')
            chains.add(chain)
            adapter = AlphaFold3BinderAdapter(
                python_path=self.config.python_path, script_path=self.config.script_path,
                model_dir=self.config.model_dir, database_dir=self.config.database_dir,
                target_msa_cache_dir=self.config.target_msa_cache_dir,
                use_templates=self.config.use_templates)
            protein = adapter.input_payload(target_sequence=entity['sequence'], binder_sequence='A',
                                            seeds=self.config.seeds)['sequences'][0]['protein']
            protein['id'] = chain
            payload['sequences'].insert(-1, {'protein': protein})
        self.expected_chain_ids = chains
        payload['name'] = payload['name'].replace('ligand_', 'ternary_', 1)
        return write_json(Path(directory) / 'af3_ternary_input.json', payload)

    def _sample(self, directory, candidate_id, seed, sample):
        row = super()._sample(directory, candidate_id, seed, sample)
        if row.get('summary_path'):
            data = _read_af3_json(Path(row['summary_path']))
            try:
                details = _read_af3_json(Path(row['confidences_path'])) if row.get('confidences_path') else {}
                ids = confidence_chain_ids(data, details)
            except ValueError as exc:
                row.update(af3_status='backend_failed', error_message=str(exc))
                return row
            expected = getattr(self, 'expected_chain_ids', None)
            if expected is None or set(ids) != expected:
                row.update(af3_status='backend_failed', error_message='TERNARY_CHAIN_IDS_MISSING_OR_MISMATCH')
            elif row.get('model_path'):
                from .screening_gates import cif_atoms
                observed = {atom.label_chain_id for atom in cif_atoms(row['model_path'])}
                if not expected.issubset(observed):
                    row.update(af3_status='backend_failed', error_message='TERNARY_STRUCTURE_CHAIN_MISSING')
            pairs = []
            for i, left in enumerate(ids):
                for j, right in enumerate(ids):
                    if i >= j:
                        continue
                    pair = {'chain_a': left, 'chain_b': right}
                    for key in ('chain_pair_iptm', 'chain_pair_pae_min'):
                        try:
                            pair[key] = finite(data[key][i][j])
                        except (KeyError, TypeError, IndexError):
                            pair[key] = None
                    pairs.append(pair)
            row['chain_pair_metrics'] = pairs
        if row.get('model_path') and row.get('af3_status') == 'success':
            row['ternary_geometry'] = ternary_geometry(row['model_path'], self.config.ligand_chain_id)
            if getattr(self, 'ligand_smiles', None):
                try:
                    from .ligand_geometry import ligand_from_cif, chirality_valid
                    row['chirality_valid'] = chirality_valid(ligand_from_cif(
                        row['model_path'], self.config.ligand_chain_id, self.ligand_smiles), self.ligand_smiles)
                except (ValueError, RuntimeError, OSError) as exc:
                    row.setdefault('warning', []).append('TERNARY_LIGAND_GEOMETRY:' + str(exc))
        row['interpretation_scope'] = 'complex_structure_confidence_not_degradation_efficacy'
        return row


def ternary_geometry(path, ligand_chain):
    """Observed heavy-atom contacts; not an affinity or degradation prediction."""
    import numpy as np
    from .screening_gates import cif_atoms
    try:
        atoms = [a for a in cif_atoms(path) if not a.is_hydrogen]
        ligand = np.asarray([a.coord for a in atoms if a.label_chain_id == ligand_chain])
        if not ligand.size or not np.isfinite(ligand).all():
            raise ValueError('missing_or_nonfinite_ligand')
        contacts, clashes = {}, 0
        for chain in sorted({a.label_chain_id for a in atoms} - {ligand_chain}):
            xyz = np.asarray([a.coord for a in atoms if a.label_chain_id == chain])
            if not np.isfinite(xyz).all():
                raise ValueError('nonfinite_partner_coordinates')
            distances = np.linalg.norm(xyz[:, None, :] - ligand[None, :, :], axis=2)
            contacts[chain] = int((distances <= 4.5).sum())
            clashes += int((distances < 1.5).sum())
        return {'status': 'COMPLETED', 'ligand_contacts_by_chain': contacts, 'severe_clash_count': clashes,
                'contact_cutoff_angstrom': 4.5, 'clash_cutoff_angstrom': 1.5}
    except (OSError, ValueError, KeyError) as exc:
        return {'status': 'FAILED', 'error': str(exc)}

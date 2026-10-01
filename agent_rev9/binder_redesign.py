"""Validate supplied target–binder backbones before ProteinMPNN redesign."""
from pathlib import Path

from .binder_analysis import protein_chains
from .protein_mpnn_tool import ProteinMPNNTool
from .small_molecule_io import write_json


def import_backbones(settings, directory, target_id, limit):
    directory=Path(directory)
    if settings.target_chain==settings.binder_chain:
        raise ValueError('redesign_requires_distinct_target_and_binder_chains')
    candidates=[]
    for index,source in enumerate(settings.backbone_paths[:limit]):
        chains=protein_chains(source)
        if chains.get(settings.target_chain,{}).get('sequence')!=settings.target_sequence:
            raise ValueError('redesign_target_sequence_mismatch')
        sequence=chains.get(settings.binder_chain,{}).get('sequence','')
        if not sequence or len(sequence)!=settings.design_length:
            raise ValueError('redesign_binder_sequence_or_length_mismatch')
        folder=directory/f'backbone_{index:03d}';folder.mkdir(parents=True,exist_ok=True)
        normalized,mapping=ProteinMPNNTool.prepare_backbone(source,settings.binder_chain,folder)
        write_json(folder/'mapping.json',mapping)
        binder_chain=mapping['binder_pdb_chain']
        candidates.append({'candidate_id':f'{target_id}_redesign_{index:04d}',
            'rfd3_sequence':sequence,'rfd3_structure_path':str(normalized),
            'target_chain':mapping['chain_mapping'][settings.target_chain], 'binder_chain':binder_chain,
            'generated':False,'is_mock':False,'status':'BACKBONE_PROVIDED','failure_reasons':[],
            'metadata':{'generated_by':'USER_BACKBONE','source_structure':str(Path(source).resolve())},
            'rfd3':{'status':'not_run','reason':'sequence_redesign_uses_supplied_backbone','structure_path':str(normalized)}})
    return candidates

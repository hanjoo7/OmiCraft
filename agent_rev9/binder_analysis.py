"""Backbone agreement, confidence filtering and sequence diversity for AF3 binders."""
from collections import Counter, defaultdict
import gzip
import math
from pathlib import Path
import shutil
import zipfile

import numpy as np

from .small_molecule_io import write_csv, write_json


def protein_chains(path):
    from Bio.Data.PDBData import protein_letters_3to1
    from Bio.PDB import MMCIFParser, PDBParser
    path = Path(path)
    parser = (MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False)
              if '.cif' in path.suffixes else PDBParser(QUIET=True))
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as handle:
        model = next(parser.get_structure('protein', handle).get_models())
    result = {}
    for chain in model:
        residues = [r for r in chain if r.resname in protein_letters_3to1 and r.id[0] == ' ']
        if not residues:
            continue
        if any('CA' not in r for r in residues):
            raise ValueError('missing_CA:' + chain.id)
        coords = np.asarray([r['CA'].coord for r in residues], dtype=float)
        if not np.isfinite(coords).all():
            raise ValueError('nonfinite_coordinates:' + chain.id)
        result[chain.id] = {'sequence': ''.join(protein_letters_3to1[r.resname] for r in residues),
                            'ca': coords, 'residue_ids': [str(r.id[1])+r.id[2].strip() for r in residues]}
    return result


def fit_transform(mobile, reference):
    mobile, reference = np.asarray(mobile, float), np.asarray(reference, float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3 or len(mobile) < 3:
        raise ValueError('alignment_requires_three_matching_CA_atoms')
    if not np.isfinite(mobile).all() or not np.isfinite(reference).all():
        raise ValueError('nonfinite_alignment_coordinates')
    x, y = mobile-mobile.mean(0), reference-reference.mean(0)
    if min(np.linalg.matrix_rank(x),np.linalg.matrix_rank(y)) < 2:
        raise ValueError('degenerate_alignment_coordinates')
    u, _, vt = np.linalg.svd(x.T @ y)
    rotation = u @ np.diag([1.,1.,-1. if np.linalg.det(u @ vt)<0 else 1.]) @ vt
    return rotation, reference.mean(0)-mobile.mean(0) @ rotation


def rmsd(x, y):
    return float(np.sqrt(np.mean(np.sum((x-y)**2, axis=1))))


def backbone_metrics(candidate, context):
    ref = protein_chains(candidate['rfd3_structure_path'])
    pred = protein_chains(candidate['af3_metrics']['model_path'])
    binder = candidate.get('binder_chain', 'A')
    matches = [key for key,value in ref.items() if key != binder and value['sequence'] == context['target_sequence']]
    if len(matches) != 1 or binder not in ref:
        raise ValueError('reference_target_chain_not_uniquely_mapped')
    if pred.get('A',{}).get('sequence') != context['target_sequence'] or pred.get('B',{}).get('sequence') != candidate['sequence']:
        raise ValueError('AF3_sequence_mismatch')
    target_ref, binder_ref = ref[matches[0]]['ca'], ref[binder]['ca']
    target_pred, binder_pred = pred['A']['ca'], pred['B']['ca']
    if len(binder_ref) != len(binder_pred):
        raise ValueError('binder_length_mismatch')
    rotation, shift = fit_transform(target_pred,target_ref)
    br, bt = fit_transform(binder_pred,binder_ref)
    all_pred, all_ref = np.concatenate([target_pred,binder_pred]), np.concatenate([target_ref,binder_ref])
    cr, ct = fit_transform(all_pred,all_ref)
    return {'binder_rmsd':rmsd(binder_pred@rotation+shift,binder_ref),
            'binder_internal_rmsd':rmsd(binder_pred@br+bt,binder_ref),
            'complex_rmsd':rmsd(all_pred@cr+ct,all_ref),
            'target_rmsd':rmsd(target_pred@rotation+shift,target_ref),
            'reference_target_chain':matches[0], 'reference_binder_chain':binder,
            'rmsd_units':'angstrom', 'rmsd_scope':'binder_CA_after_exact_target_sequence_CA_alignment'}


def confidence_metrics(candidate, context):
    from .protein_design_route import _read_af3_json
    metrics = candidate.get('af3_metrics', {})
    full = _read_af3_json(Path(metrics['confidences_path']))
    chains = np.asarray(full.get('token_chain_ids', []))
    pae = np.asarray(full.get('pae', []),dtype=float)
    target, binder = np.flatnonzero(chains=='A'), np.flatnonzero(chains=='B')
    if (len(target)!=len(context['target_sequence']) or len(binder)!=len(candidate['sequence'])
            or pae.shape!=(len(chains),len(chains)) or not np.isfinite(pae).all() or np.any(pae<0)):
        raise ValueError('AF3_PAE_token_mapping_invalid')
    forward, reverse = float(pae[np.ix_(target,binder)].mean()), float(pae[np.ix_(binder,target)].mean())
    plddt = metrics.get('plddt_summary',{}).get('by_chain_mean',{}).get('B')
    if not isinstance(plddt,(int,float)) or not math.isfinite(plddt) or not 0<=plddt<=100:
        raise ValueError('binder_plddt_unavailable')
    return {'binder_plddt':plddt, 'plddt_scale':'0-100', 'plddt_scope':'AF3_binder_atom_mean',
            'pae_mean':float(pae.mean()), 'i_pae':(forward+reverse)/2,
            'i_pae_target_to_binder':forward, 'i_pae_binder_to_target':reverse,
            'pae_units':'angstrom', 'iptm':metrics.get('iptm'), 'ptm':metrics.get('ptm')}


def candidate_metrics(candidate, context, config):
    row = {'candidate_id':candidate['candidate_id'], 'gene':candidate.get('gene_name',''),
           'sequence':candidate.get('sequence',''), 'length':len(candidate.get('sequence','')),
           'original_verdict':candidate.get('final_verdict','NOT_EVALUATED'),
           'is_mock':bool(candidate.get('is_mock')), 'errors':[]}
    if row['is_mock'] or not candidate.get('af3_metrics',{}).get('model_path'):
        row.update(analysis_status='NOT_EVALUATED',quality_status='NOT_EVALUATED')
        return row
    for function in (backbone_metrics,confidence_metrics):
        try:
            row.update(function(candidate,context))
        except Exception as exc:
            row['errors'].append(function.__name__+': '+str(exc))
    failures = []
    if not config.min_length<=row['length']<=config.max_length:
        failures.append('LENGTH_OUTSIDE_RANGE')
    for key,threshold,minimum in [('binder_plddt',config.plddt_min,True),('i_pae',config.ipae_max,False),('binder_rmsd',config.rmsd_max,False)]:
        value = row.get(key)
        if value is not None and (value<threshold if minimum else value>threshold):
            failures.append(key+'_OUTSIDE_RANGE')
    row.update(analysis_status='PARTIAL' if row['errors'] else 'COMPLETED',filter_reasons=failures,
               quality_status='FAIL' if failures else 'NOT_EVALUATED' if row['errors'] else 'PASS')
    return row


def rank_rows(rows, config):
    descending = config.rank_by in {'binder_plddt','iptm','ptm'}
    def key(row):
        value=row.get(config.rank_by)
        valid=isinstance(value,(int,float)) and math.isfinite(value)
        return (row['quality_status']!='PASS',not valid,(-value if descending else value) if valid else 0,row['candidate_id'])
    ordered=sorted(rows,key=key)
    seen={}
    for index,row in enumerate(ordered,1):
        identity=(row['gene'],row['sequence'])
        row.update(rank=index,duplicate_of=seen.get(identity),ranking_metric=config.rank_by)
        if row['sequence']:seen.setdefault(identity,row['candidate_id'])
    return ordered


def sequence_variation(candidates, rows, directory):
    by_id={c['candidate_id']:c for c in candidates}
    groups=defaultdict(list)
    for row in rows:
        if row['is_mock'] or not row['sequence']:
            continue
        candidate=by_id[row['candidate_id']]
        group=(row['gene'],candidate.get('rfd3_candidate_id') or candidate.get('rfd3_structure_path'),row['length'])
        if not any(item['sequence'] == row['sequence'] for item in groups[group]):
            groups[group].append(row)
    profile=[]
    for number,((gene,backbone,length),group) in enumerate(groups.items(),1):
        reference=group[0]['sequence']
        alignment=directory/f'alignment_{number:03d}.fasta'
        alignment.write_text(''.join('>'+r['candidate_id']+'\n'+r['sequence']+'\n' for r in group))
        for position in range(length):
            counts=Counter(r['sequence'][position] for r in group)
            consensus,count=counts.most_common(1)[0]
            profile.append({'group':number,'gene':gene,'backbone_id':backbone,'position':position+1,
                            'reference':reference[position],'consensus':consensus,'sequence_count':len(group),
                            'mutation_frequency':1-counts[reference[position]]/len(group),
                            'entropy_bits':-sum((n/len(group))*math.log2(n/len(group)) for n in counts.values()),
                            'amino_acid_counts':dict(counts)})
    write_csv(directory/'sequence_variation.csv',profile)
    return profile


def write_plots(rows, profile, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    valid=[r for r in rows if r.get('binder_plddt') is not None and r.get('i_pae') is not None]
    if valid:
        fig,ax=plt.subplots(figsize=(6,4))
        for status,color in [('PASS','#138c84'),('FAIL','#b76755'),('NOT_EVALUATED','#748397')]:
            subset=[r for r in valid if r['quality_status']==status]
            if subset:ax.scatter([r['binder_plddt'] for r in subset],[r['i_pae'] for r in subset],c=color,label=status)
        ax.set(xlabel='Binder pLDDT (0–100)',ylabel='Bidirectional iPAE (Å)',xlim=(0,100),title='Binder quality analysis · separate from primary verdict')
        ax.legend();fig.savefig(directory/'quality_map.png',dpi=140,bbox_inches='tight');plt.close(fig)
    for group in sorted({r['group'] for r in profile})[:6]:
        values=[r for r in profile if r['group']==group]
        fig,ax=plt.subplots(figsize=(8,3))
        ax.bar([r['position'] for r in values],[r['mutation_frequency'] for r in values],color='#138c84')
        ax.set(xlabel='Designed backbone position',ylabel='Mutation frequency',ylim=(0,1),
               title=f"Within-backbone sequence variation · {values[0]['sequence_count']} unique sequences")
        fig.savefig(directory/f'sequence_variation_{group:03d}.png',dpi=140,bbox_inches='tight');plt.close(fig)


def analyze_binders(screening, directory, config):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    managed_names = {'ranking.csv', 'quality_pass.csv', 'sequences.fasta', 'sequence_variation.csv',
                     'quality_map.png', 'best.pdb', 'best.cif', 'best.json', 'analysis.json', 'binder_results.zip'}
    for path in directory.iterdir():
        if path.name in managed_names or (path.name.startswith('alignment_') and path.suffix == '.fasta') or (path.name.startswith('sequence_variation_') and path.suffix == '.png'):
            if path.is_file():
                path.unlink()
    candidates=screening.get('candidates',[])
    contexts={t['target_id']:t.get('design_conditions',{}) for t in screening.get('targets',[])}
    rows=[]
    for candidate in candidates:
        row=candidate_metrics(candidate,contexts.get(candidate.get('target_id'),{}),config)
        rows.append(row);candidate['design_quality']=row
    rows=rank_rows(rows,config)
    profile=sequence_variation(candidates,rows,directory)
    write_csv(directory/'ranking.csv',rows)
    passed=[r for r in rows if r['quality_status']=='PASS' and not r['is_mock'] and not r.get('duplicate_of')]
    write_csv(directory/'quality_pass.csv',passed)
    (directory/'sequences.fasta').write_text(''.join('>'+r['candidate_id']+' | '+r['original_verdict']+'\n'+r['sequence']+'\n' for r in rows if r['sequence'] and not r['is_mock'] and not r.get('duplicate_of')))
    summary={'status':'COMPLETED','predictor':'AlphaFold3','settings':config.model_dump(),'ranking':rows,
             'quality_pass_count':len(passed),'unique_sequence_count':len({(r['gene'],r['sequence']) for r in rows if r['sequence'] and not r['is_mock']}),
             'candidate_count':len(rows),'evaluated_candidate_count':sum(r['analysis_status']=='COMPLETED' for r in rows),
             'sequence_alignment_scope':'same_backbone_equal_length_positional_alignment',
             'reference_sequence':'highest-ranked unique sequence within each backbone',
             'changes_primary_verdict':False,'best_candidate_id':None,'best_structure_path':None,
             'competition_status':'NOT_REQUESTED','warnings':[]}
    by_id={c['candidate_id']:c for c in candidates}
    best=next((r for r in rows if r['quality_status']=='PASS' and not r['is_mock'] and r['original_verdict']=='PASS'),None)
    if best:
        candidate=by_id[best['candidate_id']];source=Path(candidate['af3_metrics']['model_path'])
        destination=directory/('best.pdb' if '.pdb' in source.suffixes else 'best.cif')
        if source.suffix=='.gz':
            with gzip.open(source,'rb') as src,destination.open('wb') as dest:shutil.copyfileobj(src,dest)
        else:shutil.copy2(source,destination)
        summary.update(best_candidate_id=best['candidate_id'],best_structure_path=str(destination.resolve()))
        write_json(directory/'best.json',{'candidate':best,'source_model':str(source),'selection':'primary_PASS_and_quality_PASS'})
    else:
        summary['warnings'].append('No non-mock candidate passed both primary and additional quality filters; best structure not exported')
    try:write_plots(rows,profile,directory)
    except Exception as exc:summary['warnings'].append('plots: '+str(exc))
    archive=directory/'binder_results.zip'
    summary['archive_path']=str(archive.resolve())
    write_json(directory/'analysis.json',summary)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(directory.iterdir()):
            if path.is_file() and path!=archive and (path.name in managed_names or path.name.startswith(('alignment_', 'sequence_variation_'))):
                bundle.write(path,path.name)
    summary['artifacts']=[{'path':str(p.resolve()),'label':p.name,'artifact_type':'plot' if p.suffix=='.png' else 'structure' if p.suffix in {'.pdb','.cif'} else 'archive' if p.suffix=='.zip' else 'table'} for p in sorted(directory.iterdir()) if p.is_file() and (p.name in managed_names or p.name.startswith(('alignment_', 'sequence_variation_')))]
    return summary

"""Figures derived only from stored model outputs and analysis tables."""
from pathlib import Path


def write_design_figures(data, directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from .small_molecule_io import write_csv
    from .protein_design_route import _read_af3_json

    root = Path(directory)/'figures'
    root.mkdir(exist_ok=True)
    artifacts = []
    summary = data.get('summary', {})
    modality = data.get('modality')
    def save(fig, name, label):
        path = root/(name+'.png')
        fig.savefig(path,dpi=140,bbox_inches='tight',facecolor='white')
        plt.close(fig)
        artifacts.append({'path':str(path.resolve()),'label':label,'artifact_type':'plot',
                          'execution_mode':'DERIVED_FROM_RECORDED_OUTPUTS'})
    def source(rows,name):
        path = write_csv(root/(name+'.csv'),rows)
        artifacts.append({'path':str(path.resolve()),'label':name+' source table','artifact_type':'table'})
    with plt.rc_context({'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#092b50',
                         'text.color':'#092b50','font.size':10}):
        if modality == 'DE_NOVO_BINDER':
            candidates = summary.get('candidates', [])
            for index,c in enumerate(candidates[:5]):
                m = c.get('af3_metrics', {})
                rows = [{'metric':key,'value':m.get(key)} for key in ('iptm','ptm') if isinstance(m.get(key),(int,float))]
                if rows:
                    fig,ax=plt.subplots(figsize=(6,3.5))
                    bars=ax.bar([r['metric'] for r in rows],[r['value'] for r in rows],color=['#138c84','#426991'])
                    ax.set_ylim(0,1);ax.bar_label(bars,fmt='%.2f');ax.set_title(c['candidate_id']+' · '+c.get('final_verdict',''))
                    save(fig,f'binder_{index}_confidence','Binder AF3 confidence');source(rows,f'binder_{index}_confidence')
                path=m.get('confidences_path')
                if path and Path(path).is_file():
                    confidence=_read_af3_json(Path(path))
                    matrix=np.asarray(confidence.get('pae',[]),dtype=float)
                    if matrix.ndim==2 and matrix.shape[0]==matrix.shape[1] and np.isfinite(matrix).all():
                        fig,ax=plt.subplots(figsize=(6,5));image=ax.imshow(matrix,vmin=0,vmax=40,cmap='viridis_r')
                        fig.colorbar(image,ax=ax,label='PAE (Å)');ax.set(xlabel='Token index',ylabel='Token index',title='AF3 predicted aligned error')
                        chains=confidence.get('token_chain_ids',[])
                        for k in range(1,len(chains)):
                            if chains[k]!=chains[k-1]:ax.axhline(k-.5,color='white',lw=.7);ax.axvline(k-.5,color='white',lw=.7)
                        save(fig,f'binder_{index}_pae','AF3 PAE · recorded confidence output')
        elif modality == 'SMALL_MOLECULE':
            rows=[]
            for c in summary.get('candidates',[]):
                for pose in c.get('docking',{}).get('poses',[]):
                    if isinstance(pose.get('docking_score'),(int,float)):
                        rows.append({'pose':c.get('candidate_id','')+' / '+str(pose.get('pose_rank')), 'score':pose['docking_score']})
            if rows:
                fig,ax=plt.subplots(figsize=(8,4));ax.bar([r['pose'] for r in rows],[r['score'] for r in rows],color='#138c84')
                ax.set(ylabel='Docking score (kcal/mol)',title='Recorded docking poses');ax.tick_params(axis='x',rotation=45)
                save(fig,'docking_scores','Docking scores');source(rows,'docking_scores')
        elif modality == 'ADC':
            rows=summary.get('conjugation',{}).get('candidates',[])[:15]
            if rows:
                fig,ax=plt.subplots(figsize=(9,4));ax.bar([str(r['chain'])+str(r['residue_id']) for r in rows],
                    [r['sasa_angstrom2'] for r in rows],color=['#c7d5df' if r['exclusion_reasons'] else '#138c84' for r in rows])
                ax.set(ylabel='Residue SASA (Å²)',title='Conjugation review candidates · grey = exclusion flag')
                ax.tick_params(axis='x',rotation=45);save(fig,'conjugation_sasa','Conjugation candidates · surface accessibility');source(rows,'conjugation_sasa')
    return artifacts

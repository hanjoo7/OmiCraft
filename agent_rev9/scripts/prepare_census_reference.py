"""Download a pinned, donor-stratified breast cancer Census reference."""
import argparse
import hashlib
import json
from pathlib import Path


def prepare(directory, version='2025-11-08', disease='triple-negative breast carcinoma', max_cells=20000):
    import cellxgene_census as census_api
    import numpy as np
    from scipy import sparse
    from scipy.io import mmwrite

    out = Path(directory).resolve()
    out.mkdir(parents=True, exist_ok=True)
    query = f"tissue_general == 'breast' and disease == {disease!r} and is_primary_data == True"
    with census_api.open_soma(census_version=version) as census:
        exp = census['census_data']['homo_sapiens']
        columns = ['soma_joinid', 'cell_type', 'disease', 'donor_id', 'dataset_id', 'tissue']
        obs = exp.obs.read(value_filter=query, column_names=columns).concat().to_pandas()
        obs = obs[~obs.donor_id.astype(str).str.lower().isin(['unknown', 'na', 'nan', ''])].copy()
        obs['donor_id'] = obs.dataset_id.astype(str) + ':' + obs.donor_id.astype(str)
        obs['priority'] = obs.soma_joinid.map(lambda i: hashlib.sha256(f'11:{i}'.encode()).hexdigest())
        # Round-robin sampling preserves rare donor/cell-type strata before filling common ones.
        obs = obs.sort_values('priority')
        obs['stratum_rank'] = obs.groupby(['donor_id','cell_type'], observed=True).cumcount()
        selected = obs.sort_values(['stratum_rank','priority']).head(max_cells).sort_values('soma_joinid')
        if selected.empty:
            raise ValueError('No eligible cells for ' + query)
        print('Selected',len(selected),'cells',selected.donor_id.nunique(),'donors',flush=True)
        adata = census_api.get_anndata(census, organism='Homo sapiens', obs_coords=selected.soma_joinid.to_numpy(),
                                      obs_column_names=columns[1:], var_column_names=['feature_id','feature_name'])
        if not sparse.issparse(adata.X) or np.any(adata.X.data < 0):
            raise ValueError('Expected nonnegative sparse raw counts')
        datasets_info = census["census_info"]["datasets"].read().concat().to_pandas()
        datasets_info = datasets_info[datasets_info.dataset_id.isin(selected.dataset_id.unique())]
        all_vars = exp.ms["RNA"].var.read(column_names=["soma_joinid","feature_id","feature_name"]).concat().to_pandas()
        presence = census_api.get_presence_matrix(census, organism="Homo sapiens")
        coverage = np.asarray(presence[datasets_info.soma_joinid.to_numpy(), :].sum(axis=0)).ravel()
        measured = set(all_vars.loc[coverage[all_vars.soma_joinid.to_numpy()] == len(datasets_info), "feature_name"])
        var = adata.var.reset_index(drop=True)
        keep = ~var.feature_name.duplicated() & var.feature_name.notna()
        expr = adata.X[:, keep.to_numpy()].T.tocsr()
        var = var.loc[keep].copy()
        cells = adata.obs.copy()
        cells['cell_id'] = [str(i) for i in cells.index]
        cells['donor_id'] = cells.dataset_id.astype(str) + ':' + cells.donor_id.astype(str)
        # Original ontology annotations are retained; epithelial is not relabeled malignant.
        cells['cell_class'] = cells.cell_type.astype(str)
        mmwrite(out/'expression.mtx', expr)
        var.to_csv(out/'genes.tsv',sep='\t',index=False)
        cells.to_csv(out/'cells.tsv',sep='\t',index=False)
        datasets = sorted(cells.dataset_id.astype(str).unique())
        provenance = {'status':'COMPLETED', 'census_version':version, 'query':query, 'disease_filter_applied':True,
            'requested_disease':disease,'datasets':datasets,'n_cells':len(cells),'n_donors':cells.donor_id.nunique(),
            'n_genes':len(var),'measured_genes':sorted(measured),'dataset_metadata':datasets_info.to_dict('records'),'sampling':'SHA256(seed=11, soma_joinid); round-robin donor/cell_type strata',
            'max_cells':max_cells,'expression':'Census raw nonnegative counts',
            'limitations':['Ontology cell type does not establish malignant status',
                            'A sampled reference is not an independent validation cohort for bulk DE']}
        provenance['files'] = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                               [out/'expression.mtx',out/'genes.tsv',out/'cells.tsv']}
        (out/'provenance.json').write_text(json.dumps(provenance,indent=2,default=str))
        print(json.dumps({k:v for k,v in provenance.items() if k not in {'files','measured_genes','dataset_metadata'}}),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory',default='datasets/upstream_evidence/census')
    parser.add_argument('--version',default='2025-11-08')
    parser.add_argument('--disease',default='triple-negative breast carcinoma')
    parser.add_argument('--max-cells',type=int,default=20000)
    args=parser.parse_args()
    prepare(args.directory,args.version,args.disease,args.max_cells)

"""Fetch official breast dependency and normal protein expression snapshots."""
import argparse
import concurrent.futures
import json
import zipfile
from pathlib import Path

from .prepare_tcga_brca_inputs import digest, download


def prepare(directory):
    root=Path(directory).resolve()
    root.mkdir(parents=True,exist_ok=True)
    files=[]
    releases={27993248:('CRISPRGeneEffect.csv','Model.csv'),6025238:('D2_combined_gene_dep_scores.csv',)}
    for article, names in releases.items():
        metadata=root/f'figshare_{article}.json'
        download(f'https://api.figshare.com/v2/articles/{article}',metadata)
        record=json.loads(metadata.read_text())
        selected=[{**item,'article_id':article,'article_title':record['title']} for item in record['files'] if item['name'] in names]
        if len(selected)!=len(names):
            raise ValueError(f'Missing files in official article {article}')
        files.extend(selected)
    def fetch(item):
        path=root/item['name']
        download(item['download_url'],path,item.get('computed_md5'))
        return {**item,'path':str(path),'sha256':digest(path)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records=list(pool.map(fetch,files))
    (root/'depmap_manifest.json').write_text(json.dumps(records,indent=2))
    url='https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip'
    archive=root/'normal_ihc_data.tsv.zip'
    download(url,archive)
    with zipfile.ZipFile(archive) as zipped:
        with zipped.open('normal_ihc_data.tsv') as source:
            (root/'normal_ihc_data.tsv').write_bytes(source.read())
    (root/'hpa_manifest.json').write_text(json.dumps({'version':'25.1','url':url,'sha256':digest(archive)},indent=2))
    return records


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory',default='datasets/upstream_evidence')
    args=parser.parse_args()
    for item in prepare(args.directory):
        print(item['name'],item['sha256'])

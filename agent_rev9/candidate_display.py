"""rev7_2 candidate visibility rules; never change saved scientific decisions."""

import csv
import os
from functools import lru_cache
from pathlib import Path


def symbol(value):
    return str(value or '').strip().upper()


@lru_cache(maxsize=64)
def _read_hidden(path, mtime, size, kind):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        rows = csv.DictReader(stream, delimiter='\t' if kind == 'volcano' else ',')
        return frozenset(symbol(row.get('gene') or row.get('gene_name') or
                                row.get('gene_symbol') or row.get('gene_label'))
                         for row in rows if kind == 'essential' or
                         (kind == 'contrast' and row.get('c1_direction') == 'Down') or
                         (kind == 'volcano' and row.get('direction') == 'Down_in_TNBC')) - {''}


def hidden_genes(volcano_path=None):
    """Read explicit exclusion lists and, optionally, the current run's DE file."""
    workspace = Path(__file__).resolve().parents[2]
    roots = [workspace / 'datasets', workspace / 'dataset', workspace.parent / 'dataset']
    if os.environ.get('OMICRAFT_DATASET_ROOT'):
        roots = [Path(os.environ['OMICRAFT_DATASET_ROOT'])]
    files = [(root / 'depmap/common_essentials_computed.csv', 'essential') for root in roots]
    files += [(root / 'tcga_brca/processed/deg_results/dual_contrast_classification.csv', 'contrast') for root in roots]
    if volcano_path:
        files.append((Path(volcano_path), 'volcano'))
    hidden = {'CDKN2A'}
    for path, kind in files:
        try:
            stat = path.stat()
            hidden.update(_read_hidden(str(path), stat.st_mtime_ns, stat.st_size, kind))
        except FileNotFoundError:
            continue
    return hidden


def state_hidden_genes(state):
    output = (state.get('r_analysis_result') or {}).get('output_dir')
    return hidden_genes(Path(output) / '04_de_plots/volcano_data.tsv' if output else None)


def visible_gene(gene, hidden=None):
    return symbol(gene) not in (hidden_genes() if hidden is None else hidden)

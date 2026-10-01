"""Register a canonical human AlphaFold DB structure for an eligible Binder target."""
import hashlib
import re
from urllib.parse import urlsplit
import urllib.request

from .small_molecule_io import now, write_json
from .structure_io import protein_residues
from .target_qualification import _http_get_json, query_uniprot


def ensure_target_structure(gene, uniprot=None):
    from .web_execution import assets, storage
    rows = [row for row in assets() if row['gene'] == gene and row['kind'] == 'protein'
            and row.get('source_type') == 'alphafold_db'
            and row.get('uniprot_id') and '-' not in row['uniprot_id']]
    if rows:
        return sorted(rows, key=lambda row: row['id'])[0]
    if not re.fullmatch(r'[A-Za-z0-9_-]{1,40}', gene):
        raise ValueError('invalid_target_gene')
    uniprot = uniprot or query_uniprot(gene)
    accession = uniprot.get('uniprot_id') or ''
    if uniprot.get('status') != 'FOUND' or not re.fullmatch(r'[A-Z0-9]{6,10}', accession):
        raise ValueError('canonical_uniprot_accession_unavailable')
    records = _http_get_json('https://alphafold.ebi.ac.uk/api/prediction/' + accession)
    matches = [row for row in records if isinstance(row,dict)
               and row.get('uniprotAccession') == accession and str(row.get('taxId')) == '9606'
               and row.get('gene', '').upper() == gene.upper() and not row.get('isComplex')
               and row.get('entryId') == 'AF-' + accession + '-F1'] if isinstance(records,list) else []
    if not matches:
        raise ValueError('canonical_human_alphafold_structure_unavailable')
    model = max(matches, key=lambda row: int(row.get('latestVersion',0)))
    sequence = model.get('uniprotSequence') or ''
    if not 1 <= len(sequence) <= 1200:
        raise ValueError('target_requires_explicit_domain_structure: canonical sequence exceeds 1200 residues or is missing')
    if model.get('uniprotStart') != 1 or model.get('uniprotEnd') != len(sequence):
        raise ValueError('partial_alphafold_structure_requires_explicit_domain_selection')
    url = model.get('pdbUrl') or ''
    address = urlsplit(url)
    if address.scheme != 'https' or address.netloc != 'alphafold.ebi.ac.uk' or not address.path.startswith('/files/'):
        raise ValueError('invalid_alphafold_download_url')
    with urllib.request.urlopen(url, timeout=20) as response:
        content = response.read(5_000_001)
    if len(content) > 5_000_000:
        raise ValueError('alphafold_structure_too_large')
    identifier = 'alphafold_' + gene.lower() + '_' + accession.lower()
    directory = storage() / 'uploads' / identifier
    directory.mkdir(parents=True,exist_ok=True)
    path = directory / 'structure.pdb'
    path.write_bytes(content)
    try:
        residues = protein_residues(str(path))
        chains = {chain for chain,_,_ in residues}
        if len(chains) != 1 or ''.join(aa for _,_,aa in residues) != sequence:
            raise ValueError('alphafold_structure_sequence_mismatch')
    except Exception:
        path.unlink(missing_ok=True)
        raise
    row = {'id':identifier, 'gene':gene, 'kind':'protein', 'role':'target',
           'label':gene + ' · AlphaFold DB', 'path':str(path.resolve()),
           'source':'AlphaFold DB prediction; downloaded for Binder design, not an AF3 run',
           'source_type':'alphafold_db', 'source_url':'https://alphafold.ebi.ac.uk/entry/' + accession,
           'download_url':url, 'uniprot_id':accession, 'confidence_type':'pLDDT',
           'sequence_start':1, 'sequence_end':len(sequence), 'model_version':model.get('latestVersion'),
           'mean_plddt':model.get('globalMetricValue'), 'registered_at':now(),
           'sha256':hashlib.sha256(content).hexdigest(), 'chains':sorted(chains)}
    write_json(directory/'source.json',model)
    write_json(directory/'asset.json',row)
    return row

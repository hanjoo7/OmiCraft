"""Prepare a reproducible receptor-defined TCGA-BRCA cohort from public sources."""

import argparse
import concurrent.futures
import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SOURCES = {
    'xena_phenotype.tsv': 'https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix',
    'msigdb_hallmark.gmt': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.symbols.gmt',
    'msigdb_gobp.gmt': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c5.go.bp.v2025.1.Hs.symbols.gmt',
}


def digest(path, algorithm='sha256'):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def download(url, path, expected_md5=None):
    path = Path(path)
    if path.is_file() and (expected_md5 is None or digest(path, 'md5') == expected_md5):
        return
    temporary = path.with_suffix(path.suffix + '.part')
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'OmiCraft-public-data/1.0'})
            with urllib.request.urlopen(request, timeout=180) as response, temporary.open('wb') as output:
                while chunk := response.read(1024 * 1024):
                    output.write(chunk)
            if expected_md5 and digest(temporary, 'md5') != expected_md5:
                raise ValueError('GDC MD5 mismatch: ' + path.name)
            temporary.replace(path)
            return
        except (OSError, ValueError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def receptor_group(row):
    er = row.get('breast_carcinoma_estrogen_receptor_status', '').strip().lower()
    pr = row.get('breast_carcinoma_progesterone_receptor_status', '').strip().lower()
    ihc = row.get('lab_proc_her2_neu_immunohistochemistry_receptor_status', '').strip().lower()
    ish = row.get('lab_procedure_her2_neu_in_situ_hybrid_outcome_type', '').strip().lower()
    her2 = ihc
    if ihc not in {'positive', 'negative'} and ish in {'positive', 'negative'}:
        her2 = ish
    elif ihc in {'positive', 'negative'} and ish in {'positive', 'negative'} and ihc != ish:
        her2 = 'discordant'
    values = [er, pr, her2]
    group = 'TNBC' if values == ['negative'] * 3 else 'Non_TNBC' if 'positive' in values else None
    return group, er, pr, her2


def write_table(path, rows):
    if not rows:
        raise ValueError('Cannot write an empty cohort')
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def prepare(directory, workers=4):
    directory = Path(directory).resolve()
    raw = directory / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    star = raw / 'star_counts'
    star.mkdir(exist_ok=True)
    for name, url in SOURCES.items():
        download(url, raw / name)
    filters = {'op': 'and', 'content': [
        {'op': '=', 'content': {'field': key, 'value': value}} for key, value in
        [('cases.project.project_id', 'TCGA-BRCA'), ('analysis.workflow_type', 'STAR - Counts'), ('access', 'open')]]}
    query = {'filters': json.dumps(filters), 'fields': 'file_id,file_name,file_size,md5sum,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type',
             'size': 2000, 'format': 'JSON'}
    manifest_url = 'https://api.gdc.cancer.gov/files?' + urllib.parse.urlencode(query)
    download(manifest_url, raw / 'gdc_star_manifest.json')
    manifest = json.loads((raw / 'gdc_star_manifest.json').read_text())
    if manifest['data']['pagination']['total'] > len(manifest['data']['hits']):
        raise ValueError('GDC manifest was truncated')
    with (raw / 'xena_phenotype.tsv').open() as handle:
        clinical = {row['sampleID']: row for row in csv.DictReader(handle, delimiter='\t')}
    patients, excluded = {}, Counter()
    for item in sorted(manifest['data']['hits'], key=lambda row: row['file_id']):
        for case in item.get('cases', []):
            for sample in case.get('samples', []):
                sample_id = sample['submitter_id']
                if sample.get('sample_type') != 'Primary Tumor':
                    excluded['not_primary_tumor'] += 1
                    continue
                row = clinical.get(sample_id) or clinical.get(sample_id[:15])
                if row is None:
                    excluded['clinical_sample_not_matched'] += 1
                    continue
                group, er, pr, her2 = receptor_group(row)
                if group is None:
                    excluded['receptor_status_undetermined'] += 1
                    continue
                patient = case['submitter_id']
                if patient in patients:
                    excluded['duplicate_patient_aliquot'] += 1
                    continue
                vital = row.get('vital_status', '').upper()
                event = 1 if vital == 'DECEASED' else 0 if vital == 'LIVING' else ''
                duration = row.get('days_to_death' if event == 1 else 'days_to_last_followup', '') if event != '' else ''
                patients[patient] = {
                    'sample_id': sample_id, 'patient_id': patient, 'group': group,
                    'tss': patient.split('-')[1], 'er_status': er, 'pr_status': pr, 'her2_status': her2,
                    'os_event': event, 'os_time_days': duration, 'file_id': item['file_id'],
                    'source_md5': item['md5sum'], 'counts_path': str(star / (item['file_id'] + '.tsv')),
                }
    rows = [patients[key] for key in sorted(patients)]
    groups = Counter(row['group'] for row in rows)
    if min(groups.get('TNBC', 0), groups.get('Non_TNBC', 0)) < 10:
        raise ValueError('Insufficient receptor-defined comparison groups: ' + str(groups))
    write_table(directory / 'cohort.tsv', rows)
    (directory / 'sample_manifest.json').write_text(json.dumps({'groups': groups, 'has_normal_tissue': False,
        'cohort': 'TCGA-BRCA primary tumors; one STAR file per patient; clinical ER/PR/HER2 definition',
        'excluded': excluded}, indent=2))
    print('Selected cohort:', dict(groups), 'excluded:', dict(excluded), flush=True)

    def fetch(row):
        url = 'https://api.gdc.cancer.gov/data/' + row['file_id']
        download(url, row['counts_path'], row['source_md5'])
        return {'file_id': row['file_id'], 'sample_id': row['sample_id'], 'patient_id': row['patient_id'],
                'url': url, 'path': row['counts_path'], 'md5': row['source_md5'],
                'sha256': digest(row['counts_path'])}

    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = [pool.submit(fetch, row) for row in rows]
        for job in concurrent.futures.as_completed(jobs):
            records.append(job.result())
            if len(records) % 25 == 0 or len(records) == len(rows):
                print('Verified STAR files:', len(records), '/', len(rows), flush=True)
                (directory / 'download_progress.json').write_text(json.dumps({'complete': len(records), 'total': len(rows)}))
    sources = [{'path': str(raw / name), 'url': url, 'sha256': digest(raw / name)} for name, url in SOURCES.items()]
    sources.append({'path': str(raw / 'gdc_star_manifest.json'), 'url': manifest_url, 'sha256': digest(raw / 'gdc_star_manifest.json')})
    provenance = {'created_at': datetime.now(timezone.utc).isoformat(), 'sources': sources,
        'counts': sorted(records, key=lambda row: row['patient_id']), 'groups': groups, 'excluded': excluded,
        'count_measure': 'GDC STAR unstranded raw integer counts; not TPM or log-transformed expression',
        'selection': 'Primary tumors; first eligible file in UUID lexical order per patient; no outcome-based selection',
        'tnbc_rule': 'ER negative AND PR negative AND HER2 negative; unresolved IHC uses ISH; discordant HER2 is unknown',
        'non_tnbc_rule': 'At least one measured receptor is positive; others may be unknown',
        'normal_contrast': 'not included; no normal contrast or replication evidence inferred'}
    (directory / 'source_manifest.json').write_text(json.dumps(provenance, indent=2))
    return provenance


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, default=Path('datasets/tcga_brca'))
    parser.add_argument('--workers', type=int, choices=range(1, 9), default=4)
    args = parser.parse_args()
    prepare(args.directory, args.workers)

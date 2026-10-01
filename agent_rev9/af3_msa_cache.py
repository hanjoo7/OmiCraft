"""Persistent target MSAs, keyed by sequence and local AF3 search provenance."""

import sys as _sys
if _sys.platform == "win32":
    import types as _types
    _sys.modules.setdefault("fcntl", _types.ModuleType("fcntl"))

import fcntl
import hashlib
import json
import os
import tempfile
import time
from threading import Event
from contextlib import contextmanager
from pathlib import Path


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _query(a3m):
    lines = a3m.splitlines()
    if not lines or not lines[0].startswith('>'):
        return None
    query = []
    for line in lines[1:]:
        if line.startswith('>'):
            break
        query.append(''.join(c for c in line.strip() if c.isupper()))
    return ''.join(query)


class TargetMSACache:
    def __init__(self, directory, script_path, database_dir):
        self.directory = Path(directory).expanduser().resolve() if directory else None
        self.script = Path(script_path).expanduser().resolve()
        self.database = Path(database_dir).expanduser().resolve()

    def provenance(self):
        code = {}
        paths = [self.script, *[
            self.script.parent / 'src/alphafold3/data' / name
            for name in ('pipeline.py', 'msa.py', 'tools/jackhmmer.py')]]
        for path in paths:
            code[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        databases = {}
        # These are the protein MSA databases selected by run_alphafold.py defaults.
        for name in ('uniref90_2022_05.fa', 'mgy_clusters_2022_05.fa',
                     'uniprot_all_2021_04.fa', 'bfd-first_non_consensus_sequences.fasta'):
            path = self.database / name
            stat = path.stat() if path.is_file() else None
            databases[str(path)] = [stat.st_size, stat.st_mtime_ns] if stat else None
        return {'schema': 1, 'code': code, 'databases': databases}

    def key(self, sequence):
        return _digest({'sequence': sequence, 'provenance': self.provenance()})

    def path(self, sequence):
        return self.directory / (self.key(sequence) + '.json')

    @staticmethod
    def valid(protein, sequence):
        return (protein.get('sequence') == sequence
                and isinstance(protein.get('unpairedMsa'), str)
                and _query(protein['unpairedMsa']) == sequence
                and isinstance(protein.get('pairedMsa'), str)
                and (not protein['pairedMsa'] or _query(protein['pairedMsa']) == sequence))

    def load(self, sequence):
        if self.directory is None:
            return None
        try:
            value = json.loads(self.path(sequence).read_text())
            protein = value['protein']
            if (value['provenance'] == self.provenance()
                    and value['sha256'] == _digest(protein) and self.valid(protein, sequence)):
                return protein
        except (OSError, ValueError, KeyError, TypeError, AttributeError):
            pass
        return None

    @contextmanager
    def lock(self, sequence, *, enabled=True, timeout=None):
        if self.directory is None or not enabled:
            yield
            return
        self.directory.mkdir(parents=True, exist_ok=True)
        with self.path(sequence).with_suffix('.lock').open('a') as handle:
            if timeout is None:
                fcntl.flock(handle, fcntl.LOCK_EX)
            else:
                deadline = time.monotonic() + timeout
                while True:
                    try:
                        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        break
                    except BlockingIOError:
                        remaining = deadline - time.monotonic()
                        if remaining <= 0:
                            raise TimeoutError('MSA cache lock wait exceeded data pipeline budget')
                        Event().wait(min(.2, remaining))
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)

    def capture(self, predictions, sequence, chain_id):
        if self.directory is None:
            return False
        for source in sorted(Path(predictions).rglob('*_data.json')):
            try:
                data = json.loads(source.read_text())
                if data.get('is_mock') or data.get('dialect') != 'alphafold3':
                    continue
                protein = next((item['protein'] for item in data.get('sequences', [])
                                if item.get('protein', {}).get('id') == chain_id), {})
                if not self.valid(protein, sequence):
                    continue
                protein = {k: protein[k] for k in ('sequence', 'unpairedMsa', 'pairedMsa')}
                value = {'protein': protein, 'provenance': self.provenance(),
                         'sha256': _digest(protein), 'source': str(source.resolve())}
                self.directory.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile(mode='w', dir=self.directory, delete=False) as handle:
                    temporary = Path(handle.name)
                    json.dump(value, handle)
                try:
                    os.replace(temporary, self.path(sequence))
                finally:
                    temporary.unlink(missing_ok=True)
                return True
            except (OSError, ValueError, KeyError, TypeError, AttributeError):
                continue
        return False

    @contextmanager
    def execution(self, sequence, directory, predictions, chain_id, *, explicit=False, dry_run=False):
        enabled = self.directory is not None and not explicit and not dry_run
        # Serialize only cold searches for one target; cache hits can infer concurrently.
        cold = enabled and self.load(sequence) is None
        with self.lock(sequence, enabled=cold):
            hit = enabled and self.load(sequence) is not None
            audit = {'status': 'EXPLICIT' if explicit else 'DISABLED' if not enabled else 'HIT' if hit else 'MISS',
                     'stored': False}
            if enabled:
                audit.update(key=self.key(sequence), path=str(self.path(sequence)))
            provenance = self.provenance() if enabled else None
            fresh = not Path(predictions).exists() or not any(Path(predictions).iterdir())
            try:
                yield audit
            finally:
                if enabled and not hit and fresh and self.provenance() == provenance:
                    audit['stored'] = self.capture(predictions, sequence, chain_id)
                if not dry_run:
                    path = Path(directory) / 'target_msa_cache.json'
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(json.dumps(audit, indent=2))

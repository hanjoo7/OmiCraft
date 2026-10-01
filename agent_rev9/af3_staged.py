"""Run protein data pipelines separately, then infer the assembled AF3 complex."""
import copy
import json
import os
from pathlib import Path
import subprocess
import time

from .af3_msa_cache import TargetMSACache
from .model_runtime import gpu_environment
from .small_molecule_io import now, redact, unavailable, write_json


class StagedAF3Runner:
    def __init__(self, tool, directory, callback=None):
        self.tool = tool
        self.config = tool.config
        self.directory = Path(directory)
        self.callback = callback or (lambda *args: None)
        self.cache = TargetMSACache(self.config.target_msa_cache_dir,
                                   self.config.script_path, self.config.database_dir)
        self.audit = {'phase': 'data_pipeline', 'status': 'RUNNING', 'chains': [],
                      'data_pipeline_calls': 0, 'inference_calls': 0,
                      'started_at': now(), 'commands': []}

    def publish(self, stage, status):
        self.audit.update(phase=stage, status=status, updated_at=now())
        write_json(self.directory / 'stage_status.json', self.audit)
        self.callback(stage, status, copy.deepcopy(self.audit))

    def command(self, input_path, output_dir, *, inference):
        c = self.config
        return [c.python_path, c.script_path, f'--json_path={input_path.resolve()}',
                f'--model_dir={Path(c.model_dir).resolve()}',
                f'--db_dir={Path(c.database_dir).resolve()}',
                f'--output_dir={output_dir.resolve()}', '--gpu_device=0',
                f'--num_diffusion_samples={c.num_samples}',
                f'--run_data_pipeline={"false" if inference else "true"}',
                f'--run_inference={"true" if inference else "false"}']

    def execute(self, command, log_path, timeout, *, inference):
        import psutil

        env = gpu_environment(self.tool.gpu_device, environment=self.config.environment)
        if not inference:
            env.update(CUDA_VISIBLE_DEVICES='', JAX_PLATFORMS='cpu')
        if self.config.hmmer_bin_dir:
            env['PATH'] = str(Path(self.config.hmmer_bin_dir).resolve()) + os.pathsep + env.get('PATH', '')
        self.audit['commands'].append(redact(command))
        with log_path.open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=env)
            try:
                code = process.wait(timeout=max(.01, timeout))
                if code:
                    raise RuntimeError(f'AF3 exited {code}; see {log_path.name}')
            finally:
                if process.poll() is None:
                    # Keep the worker process group intact so its supervisor can also stop it.
                    try:
                        parent = psutil.Process(process.pid)
                        parent.suspend()
                        children = parent.children(recursive=True)
                        for child in children:
                            try:
                                child.kill()
                            except psutil.NoSuchProcess:
                                pass
                        parent.kill()
                        psutil.wait_procs(children, timeout=2)
                    except psutil.NoSuchProcess:
                        pass
                    process.wait(timeout=5)

    def proteins(self, payload, deadline):
        prepared = copy.deepcopy(payload)
        for item in prepared['sequences']:
            if 'protein' not in item:
                continue
            protein = item['protein']
            chain, sequence = protein['id'], protein['sequence']
            row = {'chain': chain, 'sequence_length': len(sequence), 'status': 'RUNNING',
                   'cache_status': 'DISABLED' if self.cache.directory is None else 'MISS', 'stored': False}
            self.audit['chains'].append(row)
            self.publish('msa_' + chain, 'RUNNING')
            try:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError('Data pipeline budget exhausted')
                with self.cache.lock(sequence, timeout=remaining):
                    cached = self.cache.load(sequence)
                    has_msa = all(isinstance(protein.get(k), str) for k in ('unpairedMsa', 'pairedMsa'))
                    if not has_msa and cached:
                        protein.update({k: cached[k] for k in ('unpairedMsa', 'pairedMsa')})
                        has_msa = True
                    explicit = chain == self.config.protein_chain_id and getattr(self, 'explicit_target', False)
                    row['cache_status'] = 'PROVIDED' if explicit else 'HIT' if cached else 'PROVIDED' if has_msa else row['cache_status']
                    if not has_msa or 'templates' not in protein:
                        folder = self.directory / 'data_pipeline' / chain
                        folder.mkdir(parents=True, exist_ok=True)
                        request = {**prepared, 'name': 'msa_' + chain.lower(),
                                   'sequences': [{'protein': {**protein, 'id': 'A'}}]}
                        input_path = write_json(folder / 'input.json', request)
                        output = folder / 'output'
                        command = self.command(input_path, output, inference=False)
                        remaining = deadline - time.monotonic()
                        if remaining <= 0:
                            raise TimeoutError('Data pipeline budget exhausted')
                        self.audit['data_pipeline_calls'] += 1
                        provenance = self.cache.provenance()
                        (folder / 'msa.log').touch()
                        row['log_path'] = str((folder / 'msa.log').resolve())
                        self.publish('msa_' + chain, 'RUNNING')
                        self.execute(command, folder / 'msa.log', remaining, inference=False)
                        files = list(output.rglob('*_data.json'))
                        if len(files) != 1:
                            raise ValueError('Expected one completed AF3 protein data file')
                        data = json.loads(files[0].read_text())
                        observed = data.get('sequences', [])
                        if data.get('is_mock') or data.get('dialect') != 'alphafold3' or len(observed) != 1:
                            raise ValueError('Invalid completed AF3 protein data file')
                        complete = observed[0].get('protein', {})
                        if not TargetMSACache.valid(complete, sequence):
                            raise ValueError('Processed MSA does not match the requested protein')
                        if 'templates' not in complete:
                            raise ValueError('Processed protein is missing its template field')
                        if self.cache.provenance() == provenance:
                            row['stored'] = self.cache.capture(output, sequence, 'A')
                        else:
                            row['cache_status'] = 'PROVENANCE_CHANGED'
                        protein.update({k: complete[k] for k in ('unpairedMsa', 'pairedMsa', 'templates')})
                        row['data_path'] = str(files[0].resolve())
                    row['status'] = 'COMPLETED'
                    if self.cache.directory:
                        row['cache_path'] = str(self.cache.path(sequence))
                self.publish('msa_' + chain, 'COMPLETED')
            except BaseException as exc:
                row['status'] = 'TIMED_OUT' if isinstance(exc, (TimeoutError, subprocess.TimeoutExpired)) else 'FAILED'
                row['error'] = str(exc)
                self.publish('msa_' + chain, row['status'])
                raise
        return prepared

    def run(self, candidate, context):
        self.directory.mkdir(parents=True, exist_ok=True)
        c = self.config
        input_path = None
        command = None
        self.explicit_target = bool(context.get('target_msa_path') or c.target_msa_path)
        try:
            input_path = self.tool.prepare_input(candidate, context, self.directory)
            dependency = self.tool.dependency_status()
            if dependency != 'available':
                return unavailable(dependency, 'AF3 local runtime is not ready', samples=[], input_json=str(input_path))
            predictions = self.directory / 'predictions'
            if predictions.exists() and any(predictions.iterdir()):
                raise ValueError('AF3 prediction directory must be fresh')
            if (self.directory / 'data_pipeline').exists():
                raise ValueError('AF3 data pipeline directory must be fresh')
            self.publish('data_pipeline', 'RUNNING')
            payload = self.proteins(json.loads(input_path.read_text()),
                                    time.monotonic() + c.data_pipeline_timeout_seconds)
            prepared = write_json(self.directory / 'af3_prepared_input.json', payload)
            self.publish('data_pipeline', 'COMPLETED')
            command = self.command(prepared, predictions, inference=True)
            self.audit['inference_calls'] = 1
            (self.directory / 'af3.log').touch()
            self.publish('af3_inference', 'RUNNING')
            self.execute(command, self.directory / 'af3.log',
                         c.inference_timeout_seconds or c.timeout_seconds, inference=True)
            result = self.tool.parse_output(predictions, candidate['candidate_id'])
            if result.get('is_mock'):
                raise ValueError('Synthetic AF3 outputs in a real backend run')
            self.publish('af3_inference', 'COMPLETED' if result['status'] == 'success' else 'FAILED')
        except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.SubprocessError) as exc:
            timed_out = isinstance(exc, (TimeoutError, subprocess.TimeoutExpired))
            stage = 'af3_inference' if self.audit['inference_calls'] else 'data_pipeline'
            self.publish(stage, 'TIMED_OUT' if timed_out else 'FAILED')
            result = unavailable('backend_failed', str(exc), samples=[],
                                 error_code=stage.upper() + ('_TIMEOUT' if timed_out else '_FAILED'))
        result.update(input_json=str(input_path) if input_path else None, command=redact(command) if command else None,
                      started_at=self.audit['started_at'], finished_at=now(), backend='alphafold3',
                      seeds=list(c.seeds), expected_sample_count=len(c.seeds)*c.num_samples,
                      gpu_device=self.tool.gpu_device, phase_status=self.audit['status'],
                      stage_details=self.audit, protein_msa_cache=self.audit['chains'],
                      target_msa_cache=self.audit['chains'][0] if self.audit['chains'] else {},
                      inference_started=bool(self.audit['inference_calls']))
        return result

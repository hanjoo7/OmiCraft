"""Private JSON worker for cancellable dashboard jobs; no HTTP entry point."""
import importlib
import importlib.util
import json
from pathlib import Path
import sys
import threading


def main():
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('_execution_launch', root / 'launch.py')
    launch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(launch)
    launch.load_package()
    from agent_rev9 import server
    from agent_rev9.configuration import OmiCraftConfig, configuration_scope
    source = Path(sys.argv[1])
    job = json.loads(source.read_text())
    server.run_state.update(job['state'])
    done = threading.Event()
    snapshot = source.parent / 'state.json'
    def flush():
        with server.run_lock:
            content = json.dumps(server.run_state, ensure_ascii=False, default=str)
        temporary = snapshot.with_suffix('.tmp')
        temporary.write_text(content)
        temporary.replace(snapshot)
    def publish():
        while not done.wait(.1):
            flush()
    reporter = threading.Thread(target=publish, daemon=True)
    reporter.start()
    try:
        module = importlib.import_module('agent_rev9.' + job['module'])
        args = [OmiCraftConfig(**item['__omicraft_config__'])
                if isinstance(item, dict) and '__omicraft_config__' in item else item for item in job['args']]
        with configuration_scope(OmiCraftConfig(**job['config'])):
            getattr(module, job['function'])(*args, **job['kwargs'])
    finally:
        done.set()
        reporter.join()
        flush()


if __name__ == '__main__':
    main()

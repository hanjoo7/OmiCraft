"""Isolated design worker; invoked by the web execution supervisor."""
import importlib.util
import json
from pathlib import Path
import sys


def main():
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('_launch',root/'launch.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.load_package()
    from agent_rev9.configuration import OmiCraftConfig
    from agent_rev9.small_molecule_io import write_json
    from agent_rev9.web_execution import run_design
    job_path=Path(sys.argv[1])
    job=json.loads(job_path.read_text())
    job['config']['_event_path']=str(job_path.parent/'worker_events.jsonl')
    result=run_design(job['identifier'],job['body'],job['data'],job['config'],job['qualification'],
                      OmiCraftConfig(**job['app_config']),finalize=False)
    write_json(job_path.parent/'worker_result.json',result)


if __name__=='__main__':
    main()

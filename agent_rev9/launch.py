"""Launch the package stored in the agent_rev9 directory."""
import argparse
import importlib.util
import os
from pathlib import Path
import runpy
import sys


def load_package():
    root = Path(__file__).resolve().parent
    name = 'agent_rev9'
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, root / '__init__.py',
                                                     submodule_search_locations=[str(root)])
        package = importlib.util.module_from_spec(spec)
        sys.modules[name] = package
        spec.loader.exec_module(package)
    return root


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('service', choices=('backend', 'proxy'))
    parser.add_argument('--port', type=int)
    parser.add_argument('--host', default='127.0.0.1')
    parser.add_argument('--config', help='Pipeline configuration JSON')
    args = parser.parse_args(argv)
    root = load_package()
    if args.config:
        os.environ['OMICRAFT_CONFIG'] = args.config
    elif (root / 'configs/upstream.local.json').is_file():
        os.environ.setdefault('OMICRAFT_CONFIG', str(root / 'configs/upstream.local.json'))
    if args.service == 'backend':
        import uvicorn
        uvicorn.run('agent_rev9.server:app', host=args.host, port=args.port or 8088)
    else:
        if args.port:
            os.environ['OMICRAFT_PROXY_PORT'] = str(args.port)
        runpy.run_module('agent_rev9.scripts.serve_review', run_name='__main__')


if __name__ == '__main__':
    main()

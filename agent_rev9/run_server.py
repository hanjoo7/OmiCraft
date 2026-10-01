"""Convenience wrapper for the backend launcher."""
import os
import sys

if __package__:
    from .launch import main
else:
    from launch import main

if __name__ == '__main__':
    main(['backend', '--port', os.environ.get('OMICRAFT_PORT', '8000'), *sys.argv[1:]])

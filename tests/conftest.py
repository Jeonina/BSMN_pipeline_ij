"""
conftest.py — pytest session setup.
Adds project root to sys.path so tests can import scripts.* modules.
"""
import sys
from pathlib import Path

# Make project root importable (e.g. scripts.auto_params)
PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

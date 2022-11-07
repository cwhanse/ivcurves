import sys
from pathlib import Path


ROOT_DIR = (Path(__file__).parent / '..' / '..' / '..').resolve()
# a str must be inserted into sys.path in order to import modules from the
# GITHUB_WORKFLOW_UTILS_PATH directory. Path objects do not work.
GITHUB_WORKFLOW_UTILS_PATH = str(ROOT_DIR / '.github' / 'workflows' / 'utils')
sys.path.insert(0, GITHUB_WORKFLOW_UTILS_PATH)


import record_scores




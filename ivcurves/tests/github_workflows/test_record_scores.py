import sys
from pathlib import Path


# this must be a str, not Path, to be inserted into sys.path
GITHUB_WORKFLOW_UTILS_PATH = str(
    Path(__file__).parent /
    '..' / '..' / '.github' / 'workflows' / 'utils'
)
sys.path.insert(0, GITHUB_WORKFLOW_UTILS_PATH)

import record_scores


import sys
from pathlib import Path
import pytest


# this must be a str, not Path, to be inserted into sys.path
BASE_DIR = (Path(__file__).parent / '..' / '..' / '..').resolve()
GITHUB_WORKFLOW_UTILS_PATH = str(BASE_DIR / '.github' / 'workflows' / 'utils')
sys.path.insert(0, GITHUB_WORKFLOW_UTILS_PATH)

import print_json_as_env


@pytest.mark.parametrize('pr_config, bad_key', [
    ({'RUN_SCORER': False,
      'REQUIREMENTS': str(BASE_DIR / 'ivcurves' / 'requirements.txt'),
      'SUBMISSION_MAIN': str(BASE_DIR / 'ivcurves' / 'precise.py')
    }, None),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(BASE_DIR / 'ivcurves' / 'requirements.txt'),
      'SUBMISSION_MAIN': str(BASE_DIR / 'ivcurves' / 'precise.py')
    }, 'RUN_SCORER'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(BASE_DIR / 'ivcurves' / 'does_not_exist.txt'),
      'SUBMISSION_MAIN': str(BASE_DIR / 'ivcurves' / 'precise.py')
    }, 'REQUIREMENTS'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(BASE_DIR / 'ivcurves' / 'requirements.txt'),
      'SUBMISSION_MAIN': str(BASE_DIR / 'ivcurves' / 'does_not_exist.py')
    }, 'SUBMISSION_MAIN'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(BASE_DIR / 'ivcurves' / 'requirements.txt'),
    }, 'SUBMISSION_MAIN'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(BASE_DIR / 'ivcurves' / 'requirements.txt'),
      'SUBMISSION_MAIN': str(BASE_DIR / 'ivcurves' / 'precise.py'),
      'EXTRA_KEY': 'str'
    }, 'EXTRA_KEY')
])
def test_validate_pr_config(pr_config, bad_key):
    try:
        print_json_as_env.validate_pr_config(pr_config)
        # no error thrown
        assert pr_config['REQUIREMENTS'].exists()
        assert pr_config['SUBMISSION_MAIN'].exists()
    except ValueError as e:
        if bad_key is None:
            assert False, f'No error should have been raised: {str(e)}'
        else:
            assert bad_key in str(e), f"'{bad_key}' did not cause the validation error."


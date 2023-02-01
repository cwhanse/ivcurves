import sys
from pathlib import Path
import pytest


ROOT_DIR = (Path(__file__).parent / '..' / '..' / '..').resolve()
# a str must be inserted into sys.path in order to import modules from the
# GITHUB_WORKFLOW_UTILS_PATH directory. Path objects do not work.
GITHUB_WORKFLOW_UTILS_PATH = str(ROOT_DIR / '.github' / 'workflows' / 'utils')
sys.path.insert(0, GITHUB_WORKFLOW_UTILS_PATH)


import print_json_as_env


@pytest.mark.parametrize('pr_config, bad_key', [
    ({'RUN_SCORER': False,
      'REQUIREMENTS': str(ROOT_DIR / 'ivcurves' / 'precise.py'),
      'SUBMISSION_MAIN': str(ROOT_DIR / 'ivcurves' / 'precise.py')
    }, None),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(ROOT_DIR / 'ivcurves' / 'precise.py'),
      'SUBMISSION_MAIN': str(ROOT_DIR / 'ivcurves' / 'precise.py')
    }, 'RUN_SCORER'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(ROOT_DIR / 'ivcurves' / 'does_not_exist.txt'),
      'SUBMISSION_MAIN': str(ROOT_DIR / 'ivcurves' / 'precise.py')
    }, 'REQUIREMENTS'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(ROOT_DIR / 'ivcurves' / 'precise.py'),
      'SUBMISSION_MAIN': str(ROOT_DIR / 'ivcurves' / 'does_not_exist.py')
    }, 'SUBMISSION_MAIN'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(ROOT_DIR / 'ivcurves' / 'precise.py'),
    }, 'SUBMISSION_MAIN'),
    ({'RUN_SCORER': True,
      'REQUIREMENTS': str(ROOT_DIR / 'ivcurves' / 'precise.py'),
      'SUBMISSION_MAIN': str(ROOT_DIR / 'ivcurves' / 'precise.py'),
      'EXTRA_KEY': 'str'
    }, 'EXTRA_KEY')
])
def test_validate_pr_config(pr_config, bad_key):
    try:
        print_json_as_env.validate_pr_config(ROOT_DIR, pr_config)
        # no error was thrown
        assert pr_config['REQUIREMENTS'].exists()
        assert pr_config['SUBMISSION_MAIN'].exists()
    except ValueError as e:
        if bad_key is None:
            assert False, f'No error should have been raised: {str(e)}'
        else:
            # bad_key will appear in the ValueError message if bad_key caused it
            assert bad_key in str(e), f"'{bad_key}' did not cause the validation error."


def test_format_bool_variables():
    validated_dict = dict(is_true=True, is_false=False)
    print_json_as_env.format_variable_values(validated_dict)

    assert validated_dict['is_true'] == 'true'
    assert validated_dict['is_false'] == 'false'


def test_format_path_variables():
    key = 'REQUIREMENTS'
    directory = 'ivcurves'
    file = 'requirements.txt'
    validated_dict = {key: ROOT_DIR / directory / file}
    print_json_as_env.format_variable_values(validated_dict, split_path_variables=True)

    assert validated_dict[key] == str(ROOT_DIR / directory / file)
    assert validated_dict[f'{key}_FILENAME'] == str(file)
    assert validated_dict[f'{key}_PATH'] == str(ROOT_DIR / directory)

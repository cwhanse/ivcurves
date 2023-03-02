import sys
from pathlib import Path
import pytest
import ivcurves.utils as utils


# a str must be inserted into sys.path in order to import modules from the
# GITHUB_WORKFLOW_UTILS_PATH directory. Path objects do not work.
GITHUB_WORKFLOW_UTILS_PATH = str(utils.REPO_ROOT_DIR / '.github' / 'workflows' / 'utils')
sys.path.insert(0, GITHUB_WORKFLOW_UTILS_PATH)


import print_json_as_env


@pytest.mark.parametrize('pr_config, bad_key', [
    (dict(
        RUN_SCORER=False,
        REQUIREMENTS=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        SUBMISSION_MAIN=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        TEST_SETS_FOR_SCORING=''
    ), None),
    (dict(
        RUN_SCORER=True,
        REQUIREMENTS=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'does_not_exist.txt'),
        SUBMISSION_MAIN=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        TEST_SETS_FOR_SCORING=''
    ), 'REQUIREMENTS'),
    (dict(
        RUN_SCORER=True,
        REQUIREMENTS=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        SUBMISSION_MAIN=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'does_not_exist.py'),
        TEST_SETS_FOR_SCORING=''
    ), 'SUBMISSION_MAIN'),
    (dict(
        RUN_SCORER=True,
        REQUIREMENTS=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        TEST_SETS_FOR_SCORING=''
    ), 'SUBMISSION_MAIN'),
    (dict(
        RUN_SCORER=True,
        REQUIREMENTS=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        SUBMISSION_MAIN=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        TEST_SETS_FOR_SCORING='caseX'
    ), 'TEST_SETS_FOR_SCORING'),
    (dict(
        RUN_SCORER=True,
        REQUIREMENTS=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        SUBMISSION_MAIN=str(utils.REPO_ROOT_DIR / 'ivcurves' / 'precise.py'),
        TEST_SETS_FOR_SCORING='',
        EXTRA_KEY='str'
    ), 'EXTRA_KEY')
])
def test_validate_pr_config(pr_config, bad_key):
    try:
        print_json_as_env.validate_pr_config(utils.REPO_ROOT_DIR, pr_config)
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
    validated_dict = {key: utils.REPO_ROOT_DIR / directory / file}
    print_json_as_env.format_variable_values(validated_dict, split_path_variables=True)

    assert validated_dict[key] == str(utils.REPO_ROOT_DIR / directory / file)
    assert validated_dict[f'{key}_FILENAME'] == str(file)
    assert validated_dict[f'{key}_PATH'] == str(utils.REPO_ROOT_DIR / directory)

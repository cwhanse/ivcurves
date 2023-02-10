import pytest
import jschon
import sys
import ivcurves.utils as utils


# a str must be inserted into sys.path in order to import modules from the
# GITHUB_WORKFLOW_UTILS_PATH directory. Path objects do not work.
GITHUB_WORKFLOW_UTILS_PATH = str(utils.REPO_ROOT_DIR / '.github' / 'workflows' / 'utils')
sys.path.insert(0, GITHUB_WORKFLOW_UTILS_PATH)


import record_scores


def test_scores_database_pass_jsonschema_validation(scores_database_json, scores_database_jsonschema_validator):
    result = scores_database_jsonschema_validator.evaluate(jschon.JSON(scores_database_json))
    validation_messages = result.output('basic')
    assert validation_messages['valid']


@pytest.mark.parametrize('overall_scores, bad_test_set', [
    ({'case1': '0', 'case2': '0'}, None),
    ({'case1': '0'}, 'case2'),
    ({'case1': 'a', 'case2': '0'}, 'case1')
])
def test_validate_overall_scores(overall_scores, bad_test_set):
    is_valid, msg = record_scores.validate_overall_scores(overall_scores)

    if bad_test_set is None:
        assert is_valid
        assert msg == ''
    else:
        assert not is_valid
        assert bad_test_set in msg, f"'{bad_test_set}' did not cause the validation error."

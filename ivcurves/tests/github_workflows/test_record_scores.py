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


def params_test_validate_overall_scores():
    valid_test_set_filenames = utils.get_filenames_in_directory(utils.TEST_SETS_DIR)

    overall_scores = {name: '0' for name in valid_test_set_filenames}

    test_params = []

    overall_scores1 = overall_scores.copy()
    test_params.append((overall_scores1, None))  # is valid

    overall_scores2 = overall_scores.copy()
    test_set_missing = list(overall_scores2.keys())[0]
    del overall_scores2[test_set_missing]  # remove a test set
    test_params.append((overall_scores2, test_set_missing))  # missing test set

    overall_scores3 = overall_scores.copy()
    test_set_NaN_score = list(overall_scores3.keys())[0]
    overall_scores3[test_set_NaN_score] = 'a'
    test_params.append((overall_scores3, test_set_NaN_score))

    return test_params


@pytest.mark.parametrize('overall_scores, bad_test_set',
                         params_test_validate_overall_scores())
def test_validate_overall_scores(overall_scores, bad_test_set):
    is_valid, msg = record_scores.validate_overall_scores(overall_scores)

    if bad_test_set is None:
        assert is_valid
        assert msg == ''
    else:
        assert not is_valid
        assert bad_test_set in msg, f"'{bad_test_set}' did not cause the validation error."

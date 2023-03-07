import json
import jschon


def test_test_sets_pass_jsonschema_validation(test_set_noisy_csv_info_and_json,
                                              ivcurve_jsonschema_validator):
    _, _, test_set_json = test_set_noisy_csv_info_and_json
    result = ivcurve_jsonschema_validator.evaluate(jschon.JSON(test_set_json))
    validation_messages = result.output('basic')
    assert validation_messages['valid'], json.dumps(validation_messages['errors'], indent=2)

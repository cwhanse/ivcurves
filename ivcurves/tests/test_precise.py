import json
import pytest
import jschon
from conftest import mp # same instance of mpmath's mp imported in ivcurves/conftest
import ivcurves.utils as utils

import ivcurves.precise as precise


def test_test_set_json_schema_validation():
    with open('ivcurves/ivcurve_jsonschema.json') as f:
        schema = json.load(f)
    schema = jschon.JSONSchema(schema)
    case_parameter_sets = utils.read_iv_curve_parameter_sets
    assert 

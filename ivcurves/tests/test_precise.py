import json
import pytest
import jschon
from conftest import mp # same instance of mpmath's mp imported in ivcurves/conftest
import ivcurves.utils as utils

import ivcurves.precise as precise


def test_test_sets_pass_jsonschema_validation(test_set_json,
                                              ivcurve_jsonschema_validator):
    result = ivcurve_jsonschema_validator.evaluate(jschon.JSON(test_set_json))
    assert result.output('basic')['valid']


def test_diff_lhs_rhs(test_set_csv_info, constants):
    _, parameter_set = test_set_csv_info
    vth = constants['vth']
    for il, io, rs, rsh, n, ns in parameter_set.values():
        v, i = 0, il # I_sc
        val = precise.diff_lhs_rhs(v, i, il, io, rs, rsh, n, vth, ns)
        assert val == 0


def test_test_sets_ivcurve_current_precision(test_set_as_pandas_df, constants):
    vth, atol = constants['vth'], constants['atol']
    df = test_set_as_pandas_df
    for _, row in df.iterrows():
        il, io, rs, rsh, n, ns = row[utils.IV_PARAMETER_NAMES]
        diff = lambda v, i: precise.diff_lhs_rhs(v, i, il, io, rs, rsh, n, vth, ns)
        for v, i in zip(row['Voltages'], row['Currents']):
            assert abs(diff(v, i)) < atol


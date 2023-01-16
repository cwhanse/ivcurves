import json
import pytest
import jschon
from conftest import mp # same instance of mpmath's mp imported in ivcurves/conftest
import ivcurves.utils as utils

import ivcurves.precise as precise


def test_test_sets_pass_jsonschema_validation(test_set_json,
                                              ivcurve_jsonschema_validator):
    result = ivcurve_jsonschema_validator.evaluate(jschon.JSON(test_set_json))
    validation_messages = result.output('basic')
    assert validation_messages['valid']


def test_test_sets_precision(test_set_as_pandas_df, constants):
    vth, atol = constants['vth'], constants['atol']
    df = test_set_as_pandas_df
    zero = mp.mpmathify(0)

    for _, row in df.iterrows():
        il, io, rs, rsh, n, ns = row[utils.IV_PARAMETER_NAMES]
        abs_diff = lambda v, i: abs(precise.diff_lhs_rhs(v, i, il, io, rs, rsh, n, vth, ns))

        assert len(row['Voltages']) == len(row['Currents'])
        for v, i in zip(row['Voltages'], row['Currents']):
            assert abs_diff(v, i) < atol

        # diode_voltage is optional so it may be an empty array
        if len(row['diode_voltage']) > 0:
            assert len(row['diode_voltage']) == len(row['Currents'])
            for vd, i in zip(row['diode_voltage'], row['Currents']):
                assert abs_diff(vd - rs * i, i) < atol

        assert abs_diff(row['v_oc'], zero) < atol
        assert abs_diff(zero, row['i_sc']) < atol

        # these values are computed again for comparison with the ones saved to
        # the test sets. It is assumed they are correct.
        v_mp, i_mp, p_mp = precise.max_power_pt_finder(il, io, rs, rsh, n, vth,
                                                       ns, atol)
        v_oc = precise.lambert_v_from_i(zero, il, io, rs, rsh, n, vth, ns)
        i_x = precise.lambert_i_from_v(v_oc / 2, il, io, rs, rsh, n, vth, ns)
        i_xx = precise.lambert_i_from_v((v_oc + v_mp) / 2, il, io, rs, rsh, n, vth, ns)

        assert abs_diff(row['v_mp'], i_mp) < atol
        assert abs_diff(v_mp, row['i_mp']) < atol
        assert abs(row['p_mp'] - p_mp) < atol
        assert abs(row['i_x'] - i_x) < atol
        assert abs(row['i_xx'] - i_xx) < atol

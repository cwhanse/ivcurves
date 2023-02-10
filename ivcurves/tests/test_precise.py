import pytest
import jschon
from conftest import mp  # same instance of mpmath's mp imported in ivcurves/conftest
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

    for _, row in df.iterrows():
        il, io, rs, rsh, n, ns = row[utils.IV_PARAMETER_NAMES]
        params = il, io, rs, rsh, n, vth, ns
        abs_diff = lambda v, i: abs(precise.diff_lhs_rhs(v, i, *params))

        assert len(row['Voltages']) == len(row['Currents'])
        for v, i in zip(row['Voltages'], row['Currents']):
            assert abs_diff(v, i) < atol

        # diode_voltage is optional so it may be an empty array
        if len(row['diode_voltage']) > 0:
            assert len(row['diode_voltage']) == len(row['Currents'])
            for vd, i in zip(row['diode_voltage'], row['Currents']):
                assert abs_diff(vd - rs * i, i) < atol

        assert abs_diff(row['v_oc'], 0) < atol
        assert abs_diff(0, row['i_sc']) < atol

        # these values are computed again for comparison with the ones saved to
        # the test sets. It is assumed they are correct.
        v_mp, i_mp, p_mp = precise.max_power_pt_finder(*params, atol=atol)
        v_oc = precise.lambert_v_from_i(0, *params)
        i_x = precise.lambert_i_from_v(v_oc / 2, *params)
        i_xx = precise.lambert_i_from_v((v_oc + v_mp) / 2, *params)

        assert abs_diff(row['v_mp'], i_mp) < atol
        assert abs_diff(v_mp, row['i_mp']) < atol
        assert abs(row['p_mp'] - p_mp) < atol
        assert abs(row['i_x'] - i_x) < atol
        assert abs(row['i_xx'] - i_xx) < atol


rho = (1/2) * (3 - mp.sqrt(5))


@pytest.mark.parametrize('func, interval, maximizer', [
    (mp.sin, [0, mp.pi], mp.pi / 2),
    (lambda x: mp.exp(-x) * mp.sin(x), [0, mp.pi], mp.pi / 4),
    (lambda x: -x**2, [-1, 1], 0),
    (lambda x: (x - rho) * (x - (1 - rho)), [0, 1], 1)  # implementation detail, 0 is correct too
])
def test_golden_search(func, interval, maximizer, constants):
    (a, b), atol, max_iter = interval, constants['atol'], 10_000
    assert abs(precise.golden_search(a, b, func, atol, max_iter) - maximizer) < atol

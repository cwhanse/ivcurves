import pytest
import numpy as np

from ivcurves.utils import mp
import ivcurves.compare_curves as compare_curves


@pytest.mark.parametrize('a, b', [
    (1, 1),
    (1, 2),
    (2, 1),
    # (1, 0.00001),
    (0.00001, 1)
])
def test__find_x_intersection_ellipses(a, b, constants):
    r = 1
    a, b, r = map(mp.mpmathify, (a, b, r))
    ellipse_x = lambda y: a * (r**2 - (y / b)**2)**(1 / 2)
    y_r = a * b * r / mp.sqrt(a**2 + b**2)  # where y = x on ellipse
    x_r = y_r
    atol = constants['atol']
    x_r_calc = compare_curves._find_x_intersection(
        (0, r), lambda _, y: ellipse_x(y), (2 * x_r, 2 * y_r), atol=atol
    )
    assert abs(x_r - x_r_calc) < atol


def test_score_curve(constants):
    vth, atol = constants['vth'], constants['atol']
    iv_known = list(map(mp.mpmathify,
                        [6.0, 3.8500023e-06, 1.6816000000000002,
                         8832.800000000005, 1.4200000000000004, 72]))
    iv_fitted = list(map(mp.mpmathify,
                         [4.2, 6.500087e-07, 0.9453, 17881.40000000001,
                          1.6300000000000006, 72]))

    assert compare_curves.score_curve(iv_known, iv_fitted, vth, 10, atol) - \
        mp.mpmathify('3.6019766551633213') < atol


def test_score_parameters():
    params_known = [6., 1e-9, 1.68, 1000., 0.5]
    params_fit = [6. - 1.2, 1e-9 + 2e-10, 1.68*0.8, 1200., 0.5]
    score = compare_curves.score_parameters(params_known, params_fit)
    assert np.isclose(score, 0.8)

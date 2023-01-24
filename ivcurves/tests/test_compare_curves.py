import pytest
import numpy as np
from conftest import mp # same instance of mpmath's mp imported in ivcurves/conftest
import ivcurves.utils as utils

import ivcurves.compare_curves as compare_curves


@pytest.mark.parametrize('a, b', [
    (1, 1),
    (1, 2),
    (2, 1),
    # (1, 0.00001),
    (0.00001, 1)
])
def test_find_x_intersection_ellipses(a, b, constants):
    r = 1
    a, b, r = map(mp.mpmathify, (a, b, r))
    ellipse_x = lambda y: a * (r**2 - (y / b)**2)**(1 / 2)
    ys_r = np.fromiter(map(mp.mpmathify, np.linspace(0, float(b), num=100)), dtype=mp.mpf)
    xs_r = ellipse_x(ys_r)
    y_r = a * b * r / mp.sqrt(a**2 + b**2)  # where y = x on ellipse
    x_r = y_r
    atol = constants['atol']
    x_r_calc = compare_curves.find_x_intersection(
        lambda _, y: ellipse_x(y),
        known_xs=xs_r, known_ys=ys_r,
        xp=2 * x_r, yp=2 * y_r,
        num_segments=99, atol=atol
    )
    assert abs(x_r - x_r_calc) < atol

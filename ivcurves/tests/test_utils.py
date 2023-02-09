import pytest

from ivcurves.utils import mp
import ivcurves.utils as utils


@pytest.mark.parametrize('mp_float, expected', [
    (mp.mpmathify('0.5'), 0),
    (mp.mpmathify('123.45'), 3),
    (mp.mpmathify('-123.45'), 3)
])
def test_mp_num_digits_left_of_decimal(mp_float, expected):
    assert utils.mp_num_digits_left_of_decimal(mp_float) == expected

import json
from pathlib import Path
import pytest
from mpmath import mp
import jschon
import numpy as np
import pandas as pd

# Functions from these are the target of the unit tests, so their functions'
# output should not be expected to be correct when used here. Functions used
# from ivcurves should be very simple or also tested elsewhere, and make the
# test suite easier to maintain due to their use here.
import ivcurves.utils as utils
import ivcurves.precise as precise


TEST_DIR = Path(__file__).parent
IVCURVES_DIR = TEST_DIR / '..'
mp.dps = 40


def test_set_to_pandas_df(test_set_parameter_sets, test_set_json):
    """
    Creates a pandas DataFrame from IV Curve JSON.
    All of the numerical values stored as strings in the JSON as parsed to
    mp.mpf.

    Parameters
    ----------
    test_set_parameter_sets : dict
        A dict mapping test set indices to their IV curve parameters.
        See :func:`utils.read_iv_curve_parameter_sets`.

    test_set_json : dict
        IV Curve JSON dict.

    Returns
    -------
        A pandas DataFrame.
    """
    params = pd.DataFrame([[Index] + params for Index, params in test_set_parameter_sets.items()],
                          columns=['Index'] + utils.IV_PARAMETER_NAMES)
    curves = pd.DataFrame(test_set_json['IV Curves'])

    # set up index
    params['Index'] = params['Index'].astype(int)
    params = params.set_index('Index')
    curves['Index'] = curves['Index'].astype(int)
    curves = curves.set_index('Index')

    # convert Voltages and Currents from string array to mp.mpf array
    arrays = ['Voltages', 'Currents']
    curves[arrays] = curves[arrays].applymap(
        lambda a: np.fromiter(map(mp.mpmathify, a), dtype=mp.mpf)
    )

    # convert from string to mp.mpf
    scalars = ['v_oc', 'i_sc', 'v_mp', 'i_mp', 'p_mp', 'Temperature']
    curves[scalars] = curves[scalars].applymap(mp.mpmathify)

    joined = params.merge(curves, on='Index', how='inner',
                          suffixes=(None, '_drop'), validate='one_to_one')
    joined = joined[(c for c in joined.columns if not c.endswith('_drop'))]

    return joined


@pytest.fixture()
def constants():
    """
    Commonly used constants of the ivcurves scripts.
    These are copied from ivcurves/utils.py. If they change there, they should
    also change here.
    """
    num_pts = 100
    precision = 16
    atol = mp.mpmathify(1e-16)

    # Boltzmann's const (J/K), electron charge (C), temp (K)
    k, q, temp_cell = map(mp.mpmathify, [1.380649e-23, 1.60217663e-19, 298.15])
    vth = (k * temp_cell) / q

    return {'k': k, 'q': q, 'temp_cell': temp_cell, 'vth': vth, 'atol': atol,
            'precision': precision, 'num_pts': num_pts}


@pytest.fixture()
def ivcurve_jsonschema():
    with open(f'{IVCURVES_DIR}/ivcurve_jsonschema.json') as f:
        return json.load(f)


@pytest.fixture()
def ivcurve_jsonschema_validator(ivcurve_jsonschema):
    jschon.create_catalog('2020-12') # identify the JSON Schema version used
    schema_validator = jschon.JSONSchema(ivcurve_jsonschema)
    return schema_validator


@pytest.fixture(scope='function', params=list(map(
    lambda name: f'{utils.TEST_SETS_DIR}/{name}',
    utils.get_filenames_in_directory(utils.TEST_SETS_DIR)
)))
def test_set_csv_info(request):
    filename = request.param
    return filename, utils.read_iv_curve_parameter_sets(filename)


@pytest.fixture()
def test_set_csv_info_and_json(test_set_csv_info, constants):
    vth, temp_cell, atol, num_pts = (constants['vth'], constants['temp_cell'],
                                     constants['atol'], constants['num_pts'])
    filename, parameter_set = test_set_csv_info
    test_set_json = precise.build_test_set_json(filename, parameter_set, vth,
                                                temp_cell, atol, num_pts)
    return test_set_csv_info, test_set_json


@pytest.fixture()
def test_set_json(test_set_csv_info_and_json):
    _, test_set_json = test_set_csv_info_and_json
    return test_set_json


@pytest.fixture()
def test_set_as_pandas_df(test_set_csv_info_and_json):
    (_, parameter_sets), test_set_json = test_set_csv_info_and_json
    return test_set_to_pandas_df(parameter_sets, test_set_json)


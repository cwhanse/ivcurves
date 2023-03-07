import pytest
import jschon
import numpy as np
import pandas as pd
import scipy

# Functions from these are the target of the unit tests, so their functions'
# output should not be expected to be correct when used here. Functions used
# from ivcurves should be very simple or also tested elsewhere, and make the
# test suite easier to maintain due to their use here.
from ivcurves.utils import mp
import ivcurves.utils as utils
import ivcurves.precise as precise
import ivcurves.build_case3 as build_case3


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
    curves_metadata = test_set_json
    curves = pd.DataFrame(curves_metadata['IV Curves'])
    curves['cells_in_series'] = curves_metadata['cells_in_series']
    joined = params.merge(curves, on='Index', how='inner',
                          suffixes=(None, '_drop'), validate='one_to_one')
    joined = joined[(c for c in joined.columns if not c.endswith('_drop'))]

    # parse strings to np.float64
    is_array = ['Currents', 'Voltages', 'diode_voltage']
    joined[is_array] = joined[is_array].applymap(
        lambda a: np.fromiter(map(mp.mpmathify, a), dtype=mp.mpf)
    )
    is_number = ['v_oc', 'i_sc', 'v_mp', 'i_mp', 'p_mp', 'i_x', 'i_xx',
                 'Temperature']
    joined[is_number] = joined[is_number].applymap(mp.mpmathify)

    joined['Boltzmann'] = scipy.constants.Boltzmann
    joined['Elementary Charge'] = scipy.constants.elementary_charge
    joined['Vth'] = (
        joined['Boltzmann'] * joined['Temperature']
        / joined['Elementary Charge']
    )

    return joined


@pytest.fixture()
def constants():
    return utils.constants()


@pytest.fixture()
def ivcurve_jsonschema():
    return utils.load_json(utils.IVCURVES_DIR / 'ivcurve_jsonschema.json')


@pytest.fixture()
def ivcurve_jsonschema_validator(ivcurve_jsonschema):
    jschon.create_catalog('2020-12')  # identify the JSON Schema version used
    schema_validator = jschon.JSONSchema(ivcurve_jsonschema)
    return schema_validator


@pytest.fixture(scope='function', params=list(map(
    lambda name: utils.TEST_SETS_DIR / name,
    utils.get_filenames_in_directory(utils.TEST_SETS_DIR)
)))
def test_set_csv_info(request):
    filename = request.param
    return filename, utils.read_iv_curve_parameter_sets(filename)


@pytest.fixture()
def test_set_csv_info_and_json(test_set_csv_info, constants):
    vth, temp_cell, atol, num_pts = (constants['vth'], constants['temp_cell'],
                                     constants['atol'], constants['num_pts'])
    _, parameter_set = test_set_csv_info
    test_set_json = precise.build_precise_json(
        parameter_set, vth, temp_cell, atol, num_pts
    )
    return test_set_csv_info, test_set_json


@pytest.fixture(scope='function', params=build_case3.case3_output())
def test_set_noisy_csv_info_and_json(request):
    # unpack to indentify to tuple items
    filename, params_csv_df, test_set_json = request.param
    return filename, params_csv_df, test_set_json


@pytest.fixture()
def test_set_json(test_set_csv_info_and_json):
    _, test_set_json = test_set_csv_info_and_json
    return test_set_json


@pytest.fixture()
def test_set_as_pandas_df(test_set_csv_info_and_json):
    (_, parameter_sets), test_set_json = test_set_csv_info_and_json
    return test_set_to_pandas_df(parameter_sets, test_set_json)


@pytest.fixture()
def scores_database_json():
    return utils.load_json(utils.DOCS_DIR / 'scores_database.json')


@pytest.fixture()
def scores_database_jsonschema():
    return utils.load_json(utils.DOCS_DIR / 'scores_database_jsonschema.json')


@pytest.fixture()
def scores_database_jsonschema_validator(scores_database_jsonschema):
    jschon.create_catalog('2020-12')  # identify the JSON Schema version used
    schema_validator = jschon.JSONSchema(scores_database_jsonschema)
    return schema_validator

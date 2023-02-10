# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 11:17:56 2023

@author: cliff
"""

# these modules are part of the Python standard library
import os
import json

# these modules are installed
import numpy as np
import pandas as pd

# local imports from ivcurves

from ./precise import write_test_set_json


def json_file_to_df(filepath):
    """
    Creates a pandas DataFrame from an IV Curve JSON file.
    All of the numerical values stored as strings in the JSON as parsed to
    np.float64.

    Parameters
    ----------
    filepath : pathlib.Path
        A Path object pointing to the JSON file.

    Returns
    -------
        A pandas DataFrame.
    """
    curves_metadata = pd.read_json(filepath)
    curves = pd.DataFrame(curves_metadata['IV Curves'].values.tolist())
    curves['cells_in_series'] = curves_metadata['cells_in_series']

    # set up index
    curves['Index'] = curves['Index'].astype(int)
    curves = curves.set_index('Index')

    # convert Voltages, Currents, and diode_voltage from string arrays to
    # float arrays. This truncates from precise values to 64-bit precision.
    is_array = ['Voltages', 'Currents', 'diode_voltage']
    curves[is_array] = curves[is_array].applymap(
        lambda a: np.asarray(a, dtype=np.float64)
    )

    # convert from string to float
    is_number = ['v_oc', 'i_sc', 'v_mp', 'i_mp', 'p_mp', 'Temperature']
    curves[is_number] = curves[is_number].applymap(np.float64)

    return curves


dirn = 'D:\\PVmodeling\\ivcurves\\ivcurves\\test_sets'
filen = 'case3t.json'

with open(os.path.join(dirn, filen), 'r') as infile:
    df = json_file_to_df(infile)


# two base curves: one at high irradiance, high efficiency, one at low/low
base_hi = df.loc[19]
base_lo = df.loc[14]

N = len(base_hi['Currents'])  # base_lo['Currents'] will have same shape
iters = 50

# statistics for error distribution on current values: +/-1.5%
# we will assume a triangular distribution in order to have upper and lower
# bounds, and a central mode at 0.
r_hi = np.random.RandomState(seed=11724).triangular(-0.015, 0., 0.015,
                                                    size=(iters, N))
r_low = np.random.RandomState(seed=11724).triangular(-0.015, 0., 0.015,
                                                     size=(iters, N))




# -*- coding: utf-8 -*-
"""
Uses pvlib.ivtools.sde.fit_sandia_simple to estimate the single diode equation
parameters for each IV curve. The fitting method is documented in [1]_

References
----------
.. [1] C. B. Jones, C. W. Hansen, "Single Diode Parameter Extraction from
   In-Field Photovoltaic I-V Curves on a Single Board Computer", 46th IEEE
   Photovoltaic Specialist Conference, Chicago, IL, 2019

"""

import json
import pathlib
import numpy as np
import pandas as pd
import pvlib


def get_test_set_filepaths():
    """
    Returns a sorted list of pathlib.Path objects pointing to the JSON test set
    files. pathlib.Path objects can be passed directly to Python's ``open``
    function to open the JSON file.

    Returns
    -------
        A list.
    """
    path_to_test_sets = pathlib.Path.cwd() / '..' / '..' / 'ivcurves'/'test_sets'
    file_entries = list({path_to_test_sets / f'{entry.stem}.json'
                         for entry in path_to_test_sets.iterdir()
                         if entry.is_file()})
    file_entries.sort()
    return file_entries


def get_test_set_name(filepath):
    """
    Gets a test set filename from a filepath.

    Parameters
    ----------
    filepath : pathlib.Path
        A filepath pointing to a JSON test set file.

    Returns
    -------
        The test set name given a pathlib.Path object pointing to a JSON
        test set file.
    """
    return filepath.stem


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


def json_file_to_dict(filepath):
    """
    Returns a Python dict of the contents of a JSON file.

    Parameters
    ----------
    filepath : pathlib.Path
        The filepath pointing to a JSON file.

    Returns
    -------
        A Python dict.
    """
    with open(filepath, 'r') as file:
        return json.load(file)


if __name__ == "__main__":

    test_files = get_test_set_filepaths()
    
    # set up dataframe to hold results
    cols = ['photocurrent', 'saturation_current', 'resistance_series',
            'resistance_shunt', 'n', 'cells_in_series']
    
    # physical constants
    k = 1.38066E-23  # Boltzman J/K
    q = 1.60218E-19  # elementary charge (Coulomb)

    # fit each IV curve in data1 and data2
    
    for filen in test_files:
        # extract case name and make output file name
        casename = pathlib.Path(filen).name.strip('.json')
        outname = casename + '.csv'
        # read IV curves for this case
        data = json_file_to_df(filen)
        # set up Dataframe to contain each curve's diode equation parameters
        results = pd.DataFrame(index=data.index, columns=cols)
        # process each IV curve
        for d in data.index:
            il, io, rs, rsh, nNsVth = pvlib.ivtools.sde.fit_sandia_simple(
                voltage=data.loc[d, 'Voltages'],
                current=data.loc[d, 'Currents'],
                v_oc=data.loc[d, 'v_oc'],
                i_sc=data.loc[d, 'i_sc'],
                v_mp_i_mp=(data.loc[d, 'v_mp'], data.loc[d, 'v_mp'])
                )
            vth = k / q * data.loc[d, 'Temperature']
            n = nNsVth / data.loc[d, 'cells_in_series'] / vth

            results.loc[d, 'photocurrent'] = il
            results.loc[d, 'saturation_current'] = io
            results.loc[d, 'resistance_series'] = rs
            results.loc[d, 'resistance_shunt'] = rsh
            results.loc[d, 'n'] = n
            results.loc[d, 'cells_in_series'] = data.loc[d, 'cells_in_series']

        # for case3, save average rather than the per curve fits
        if casename in ['case3a', 'case3b', 'case3c', 'case3d']:
            results = results.mean().to_frame().T
            results['Index'] = 1
            results.index = results['Index']
    
        outfilen = pathlib.Path.cwd() / outname
        with open(outfilen, 'w') as outfile:
            results.to_csv(outfile)

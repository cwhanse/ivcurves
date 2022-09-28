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
from pathlib import Path
import numpy as np
import pandas as pd
import pvlib


def get_test_set_filenames():
    path_to_test_sets = Path('../../test_sets')
    return [f'{path_to_test_sets}/{entry.stem}.json'
            for entry in path_to_test_sets.iterdir()
            if entry.is_file()]


def read_json(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)


def json_to_df(filepath):
    data = read_json(filepath)
    df = pd.DataFrame.from_dict(data['IV Curves'], dtype=float)
    # set up index
    df['Index'] = df['Index'].astype(int)
    df = df.set_index('Index')
    
    # convert Voltages and Currents from string to float
    for s in ['Voltages', 'Currents']:
        df[s] = df[s].astype(object)
    for idx in df.index:
        for s in ['Voltages', 'Currents']:
            df[s][idx] = [np.asarray(df[s][idx], dtype=np.float)]
    return df


test_files = get_test_set_filenames()

# set up dataframe to hold results
cols = ['photocurrent', 'saturation_current', 'resistance_series',
        'resistance_shunt', 'n', 'cells_in_series']

# physical constants
k = 1.38066E-23  # Boltzman J/K
q = 1.60218E-19  # elementary charge (Coulomb)

# fit each IV curve in data1 and data2

for filen in test_files:
    # extract case name and make output file name
    casename = Path(filen).name.strip('.json')
    outname = casename + '.csv'
    # read IV curves for this case
    data = json_to_df(filen)
    # set up Dataframe to contain each curve's diode equation parameters
    results = pd.DataFrame(index=data.index, columns=cols)
    # process each IV curve
    for d in data.index:
        il, io, rs, rsh, nNsVth = pvlib.ivtools.sde.fit_sandia_simple(
            voltage=data.loc[d, 'Voltages'][0],
            current=data.loc[d, 'Currents'][0],
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


    outfilen = Path.cwd() / outname
    with open(outfilen, 'w') as outfile:
        results.to_csv(outfile)

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:33:11 2022

@author: cliff
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd
import pvlib


def json_file_to_dict(filepath):
    with open(filepath, 'r') as file:
        return json.load(file)


def dict_to_df(data):
    df = pd.DataFrame.from_dict(data, dtype=float)
    # set up index
    df['Index'] = df['Index'].astype(int)
    df = df.set_index('Index')
    
    # convert Voltages and Currents from string to float
    df['Voltages'] = df['Voltages'].astype(object)
    df['Currents'] = df['Currents'].astype(object)
    for idx in df.index:
        df['Voltages'][idx] = [np.asarray(df['Voltages'][idx],
                                          dtype=np.float)]
        df['Currents'][idx] = [np.asarray(df['Currents'][idx],
                                          dtype=np.float)]
    return df


here = Path.cwd()
filepath1 = here.parent.parent / 'test_sets' / 'case1.json'
filepath2 = here.parent.parent / 'test_sets' / 'case2.json'

case1 = json_file_to_dict(filepath1)
data1 = dict_to_df(case1['IV Curves'])
case2 = json_file_to_dict(filepath1)
data2 = dict_to_df(case2['IV Curves'])

cols = ['photocurrent', 'saturation_current', 'resistance_series',
        'resistance_shunt', 'n', 'cells_in_series']
results1 = pd.DataFrame(index=data1.index, columns=cols)
results2 = pd.DataFrame(index=data2.index, columns=cols)

# names for output files
outfiles = ['case1.csv', 'case2.csv']
k = 1.38066E-23
q = 1.60218E-19

# fit each IV curve in data1 and data2

for data, results, outn in zip([data1, data2], [results1, results2], outfiles):
    for d in data.index:
        il, io, rs, rsh, nNsVth = pvlib.ivtools.sde.fit_sandia_simple(
            voltage=data.loc[d, 'Voltages'][0],
            current=data.loc[d, 'Currents'][0],
            v_oc=data.loc[d, 'v_oc'],
            i_sc=data.loc[d, 'i_sc'],
            v_mp_i_mp=(data.loc[d, 'v_mp'], data1.loc[d, 'v_mp'])
            )
        vth = k / q * data.loc[d, 'Temperature']
        n = nNsVth / data.loc[d, 'cells_in_series'] / vth
        results.loc[d, 'photocurrent'] = il
        results.loc[d, 'saturation_current'] = io
        results.loc[d, 'resistance_series'] = rs
        results.loc[d, 'resistance_shunt'] = rsh
        results.loc[d, 'n'] = n
        results.loc[d, 'cells_in_series'] = data1.loc[d, 'cells_in_series']


    outfilen = here / outn
    with open(outfilen, 'w') as outfile:
        results1.to_csv(outfile)



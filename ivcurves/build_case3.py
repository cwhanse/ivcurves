# these modules are installed
import numpy as np
import pandas as pd

# local imports from ivcurves
from ivcurves.utils import TEST_SETS_DIR, save_json


# build in a list of seeds for the random number generator to make the
# result reproducible
seeds = np.array(
    [8466243, 7830653, 8545626, 7577970, 8801099, 2805418,
     9245802, 7810002, 7399955, 5811689, 6807587, 6453042,
     188361, 1585792, 46235, 3378609, 1284417, 2086151,
     9303880, 4277404, 8845973, 4413848, 6124386, 4314400,
     9917234, 9871633, 4133262, 2165557, 6749672, 1707332,
     5747407, 1332340, 2586889, 3470624, 9119847, 940712,
     745427, 5591200, 5919274, 3983076, 3237398, 1186979,
     1489190, 6546041, 4183251, 9823942, 855015, 6887809,
     2629362, 500492, 9299569, 4298767, 9411704, 7233997,
     8574591, 7415107, 3302935, 6049057, 6670139, 734617,
     3831627, 8665502, 5069815, 7054591, 4441412, 6441311,
     852088, 9202321, 9115242, 5838268, 9259171, 5501799,
     5141066, 4551869, 4020087, 3057684, 7438913, 7022044,
     9607465, 5709633, 6020128, 4430908, 7251659, 2985635,
     361155, 1944900, 1868949, 1053817, 2787375, 9353124,
     3549314, 5628355, 8406221, 1654505, 379335, 5635444,
     3310782, 4739149, 4867500, 4148131, 3072551, 5271572,
     4762255, 1012816, 4149552, 1740142, 2709454, 7229607,
     8045390, 5758685, 9225621, 7488482, 7217453, 2097634,
     1256789, 7792249, 2981597, 6450989, 5079555, 2610293,
     2612778, 585579, 8328074, 6919901, 6879644, 7115600,
     8169834, 731479, 4997481, 5240554, 6564742, 8011437,
     4923227, 8847485, 4053790, 498749, 9926682, 806283,
     4564828, 7029113, 318340, 760405, 4784796, 9705194,
     2425602, 5990078, 4389102, 1566746, 7602949, 4189657,
     6128628, 3622223, 7767334, 9339279, 5025177, 9207281,
     9120821, 2017036, 1974756, 4250032, 622857, 4017069,
     9542110, 4906580, 7796260, 1420480, 3857621, 7550123,
     3393958, 1765965, 5453168, 9749620, 1546224, 5444781,
     9786495, 1042603, 4019419, 6615556, 9792843, 9156351,
     4815128, 7869654, 3616106, 7957564, 2177906, 5031677,
     299727, 5935315, 8504223, 6319429, 2133969, 6111877,
     1600009, 9173877, 6848020, 2215352, 342433, 2262764,
     4214296, 6954867]
     )


def corr_normal_ran(g, tau):
    ''' Modifies a sequence of length n of normally-distributed random numbers
    to have autocorrelation C(m, m+k) = exp(-k/tau).

    Algorithm from [1]_.

    Parameters
    ----------
    g : numeric, 1d
        Vector of uncorrelated, normally-distributed random numbers
    tau : float
        parameter defining autocorrelation

    References
    ----------
    [1] M. Deserno. How to generate exponentially correlated Gaussian random
    numbers. https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
    '''

    # define constants needed to impose correlation
    f = np.exp(-1./tau)
    sf = np.sqrt(1. - f**2.)

    # walk through and impose correlation
    def fun(x, g, f, sf):
        return f * x + sf * g

    x = np.zeros_like(g)
    x[0] = g[0]

    for i in np.arange(1, len(x)):
        x[i] = fun(x[i-1], g[i], f, sf)

    return x


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


def apply_error(currents, rans):
    iters = rans.shape[0]
    rcur = np.tile(currents, (iters, 1))
    return rcur * (1. + rans)


def _nparray_to_str(x):
    ''' Converts a numpy array to a list of strings
    '''
    return [str(v) for v in x]


def _npint_to_int(x):
    ''' Converts a numpy int64 or int32 to python int
    '''
    return int(x)


def case3_output():
    # two base curves: one at high irradiance, high efficiency, one at low/low
    curve_list = {'case1': [19, 14],
                  'case2': [19, 14]}

    output_suffix = {'case1': {19: 'a', 14: 'b'},
                     'case2': {19: 'c', 14: 'd'}}

    iters = 50
    # statistics for error distribution on current values
    # we will assume a normal distribution with zero mean and 0.1% scale.
    # later we'll autocorrelate these errors to have autocorrelation function
    # exp(-1/tau)
    mu = 0.
    sig = 0.001
    tau = 60.

    # statistics for error distribution in voltage
    # we will assume a uniform distribution with zero mean and 0.01% width
    # voltage errors are assumed to be uncorrelated
    w = (-0.0005, 0.0005)

    seed_idx = 0

    outputs = []
    for case in ['case1', 'case2']:
        filen = case + '.json'
        with open(TEST_SETS_DIR / filen, 'r') as infile:
            curves = json_file_to_df(infile)

        filen = case + '.csv'
        with open(TEST_SETS_DIR / filen, 'r') as infile:
            params = pd.read_csv(infile, index_col=0)

        for idx in curve_list[case]:
            base = curves.loc[idx]
            N = len(base['Currents'])

            outdf = pd.DataFrame(
                index=np.arange(1, iters + 1),
                columns=['Index', *curves.columns],
                data=dict(Index=range(1, iters+1)),
                dtype=object
            )

            # random values to modify current
            r = np.random.RandomState(seed=seeds[seed_idx]).normal(
                mu, sig, size=(iters, N))
            seed_idx += 1
            # random values to modify voltage
            rv = np.random.RandomState(seed=seeds[seed_idx]).uniform(
                w[0], w[1], size=(iters, N))
            seed_idx += 1

            for i in outdf.index:
                # add random noise to current
                g = corr_normal_ran(r[i-1, :], tau)
                cur = curves.loc[idx, 'Currents'] * (1 + g)
                # add random noise to voltage
                vol = curves.loc[idx, 'Voltages'] * (1 + rv[i-1, :])
                diode_voltage = vol + params.loc[idx, 'resistance_series']*cur
                outdf.loc[i, 'i_sc'] = cur[0]
                p = cur * vol
                midx = np.argmax(p)
                outdf.loc[i, 'i_mp'] = cur[midx]
                outdf.loc[i, 'v_mp'] = vol[midx]
                outdf.loc[i, 'v_oc'] = np.max(vol)
                outdf.loc[i, 'p_mp'] = p[midx]
                # Convert numpy arrays to lists of str
                outdf.loc[i, 'Currents'] = _nparray_to_str(cur)
                outdf.loc[i, 'Voltages'] = _nparray_to_str(vol)
                outdf.loc[i, 'diode_voltage'] = _nparray_to_str(diode_voltage)

            outdf['Temperature'] = curves.loc[idx, 'Temperature']

            outdf['Sweep direction'] = ''
            outdf['Datetime'] = '1970-01-01T00:00:00Z'
            outdf = outdf.drop(['cells_in_series'], axis=1)

            # write json output
            output = {'Manufacturer': '', 'Model': '', 'Serial Number': '',
                      'Module ID': '',  'Description': '', 'Material': '',
                      'cells_in_series': int(params.loc[idx, 'cells_in_series']),
                      'IV Curves': outdf.replace(float('NaN'), None).to_dict(orient='records')}
            outfilen = 'case3' + output_suffix[case][idx]

            outparams = params[params.index == idx].copy()

            outputs.append((outfilen, outparams, output))

    return outputs


if __name__ == '__main__':
    outputs = case3_output()
    for outfilen, csv_df, output in outputs:
        save_json(output, TEST_SETS_DIR / f'{outfilen}.json')

        # write corresponding csv file
        new_idx = pd.Index(np.arange(1, len(csv_df) + 1), name='Index')
        csv_df.index = new_idx
        with open(TEST_SETS_DIR / f'{outfilen}.csv', 'w') as outfile:
            csv_df.to_csv(outfile, lineterminator='\n')


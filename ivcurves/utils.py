import csv
import pathlib
import scipy
from mpmath import mp


TEST_SETS_DIR = f'{pathlib.Path(__file__).parent}/../test_sets'
IV_PARAMETER_NAMES = ['photocurrent', 'saturation_current',
                      'resistance_series', 'resistance_shunt', 'n',
                      'cells_in_series']


def set_globals():
    r"""
    Sets library parameters that must be the same whenever the libraries are
    imported.

    ivcurves scripts should import these libraries from this script's
    namespace to use these library parameter settings.

    The following are set:

    - ``mpmath``: The precision of calculations (``mp.dps``) is set to 40
      decimal places.
    """
    mp.dps = 40 # set precision, 16*2 rounded up


def constants():
    r"""
    Commonly used constants of the ivcurves scripts.
    """
    num_pts = 100
    precision = 16
    atol = mp.mpmathify(1e-16)

    # Boltzmann's const (J/K), electron charge (C), temp (K)
    k, q, temp_cell = map(mp.mpmathify, [
        scipy.constants.Boltzmann,
        scipy.constants.elementary_charge,
        298.15
    ])
    vth = (k * temp_cell) / q

    return {'k': k, 'q': q, 'temp_cell': temp_cell, 'vth': vth, 'atol': atol,
            'precision': precision, 'num_pts': num_pts}


def mp_num_digits_left_of_decimal(num_mpf):
    r"""
    Finds the number of digits to the left of an mpmath float's decimal point.
    If the mpmath float is strictly between -1 and 1, the number of digits
    is zero.

    Parameters
    ----------
    num_mpf : numeric
        The mpmath float. [-]

    Returns
    -------
    int
        The number of digits to the left of the decimal point of ``num_mpf``,
        ignoring leading zeros.
    """
    if abs(num_mpf) < 1:
        # ignore leading zero
        return 0
    else:
        precision = constants()['precision']
        # force mpf to string in decimal format, no scientific notation
        # mpf string will have precision*2 significant digits
        # all leading zeros are stripped
        res = mp.nstr(num_mpf, n=precision*2, min_fixed=-mp.inf,
                      max_fixed=mp.inf).find('.')
        if num_mpf < 0:
            return res - 1 # ignore negative sign '-'
        else:
            return res


def mp_nstr_precision_func(num_mpf):
    r"""
    Converts an mpmath float to a string with at least 16 significant digits
    after the decimal place.

    Parameters
    ----------
    num_mpf : numeric
        The mpmath float. [-]

    Returns
    -------
    str
        A string representation of ``num_mpf`` with at least 16 significant
        digits after the decimal place.
    """
    precision = constants()['precision']
    ldigits = mp_num_digits_left_of_decimal(num_mpf)

    # Stringifying to 16 significant digits truncates or rounds the mp.mpf
    # value. The lost precision can cause error greater than 1e-16 when
    # calculating the difference between the left and right side of the single
    # diode equation. The number of required sigfigs in the string is increased
    # by 3 to ensure the JSON test set values are still precise enough when
    # converted into mp.mpf values. The required sigfigs may need to be
    # increased when new JSON test sets are added.
    sigfigs = ldigits + precision + 3
    return mp.nstr(num_mpf, n=sigfigs, strip_zeros=False)


def read_iv_curve_parameter_sets(filename):
    r"""
    Returns a dictionary of indices to a list of these values:
    photocurrent, saturation_current, resistance_series,
    resistance_shunt, n, and cells_in_series.
    The indices and values are read from the CSV file at ``filename``.

    Parameters
    ----------
    filename : str
        The path to a CSV file with these column names:
        Index, photocurrent, saturation_current, resistance_series,
        resistance_shunt, n, and cells_in_series.
        The path must exclude the file extension.

    Returns
    -------
    dict
    """
    with open(f'{filename}.csv', newline='') as file:
        reader = csv.DictReader(file, delimiter=',')
        mapping = {}
        for row in reader:
            mapping[int(row['Index'])] = [mp.mpmathify(row[col])
                                          for col in IV_PARAMETER_NAMES]
        return mapping


def make_iv_curve_name(test_set_name, index):
    r"""
    Returns a unique name for an IV curve created from parameters of
    a test set's test case. The unique name is of the form
    ``'{test_set_name}_case_{index}'``.

    Parameters
    ----------
    test_set_name : str
        The name of the test set that contains the test case of the IV curve's
        parameters.

    index : int
        The Index of the test case of the IV curve's parameters.

    Returns
    -------
    str
        A unique name for the IV curve.
    """
    return f'{test_set_name}_case_{index}'


def get_filenames_in_directory(directory_path):
    """
    Returns a set of entries in the directory ``directory_path``.
    The filenames do not have file extensions.

    Returns
    -------
    set
        A set of filenames without file extensions.
    """
    return {entry.stem for entry in pathlib.Path(directory_path).iterdir()}


set_globals()


import os
import csv
from mpmath import mp


TESTS_DIR = 'tests'


def set_globals():
    mp.dps = 40 # set precision, 16*2 rounded up


def constants():
    num_pts = 100
    precision = 16
    atol = mp.mpmathify(1e-16)

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = map(mp.mpmathify, [1.380649e-23, 1.60217663e-19, 298.15])
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
        The number of digits to the left of the decimal point of `num_mpf`,
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
        return mp.nstr(num_mpf, n=precision*2, min_fixed=-mp.inf,
                       max_fixed=mp.inf).find('.')


def mp_nstr_precision_func(num_mpf):
    r"""
    Converts an mpmath float to a string with 16 significant digits
    after the decimal place.

    Parameters
    ----------
    num_mpf : numeric
        The mpmath float. [-]

    Returns
    -------
    string
        A string representation of `num_mpf` with 16 siginigcant digits
        after the decimal place.
    """
    precision = constants()['precision']
    ldigits = mp_num_digits_left_of_decimal(num_mpf)
    return mp.nstr(num_mpf, n=ldigits+precision, strip_zeros=False)


def read_case_parameters(filename):
    r"""
    Returns a 2-D list of entries in a CSV file with these column names:
    Index, photocurrent, saturation_current, resistance_series,
    resistance_shunt, n, and cells_in_series.

    Parameters
    ----------
    filename : string
        The path to the CSV file excluding the file extension.

    Returns
    -------
    2-D list
        The 2-D list stores every value in the CSV file as a mpmath float,
        except values in the Index column which are stored as integers.
    """
    with open(f'{filename}.csv', newline='') as file:
        reader = csv.DictReader(file, delimiter=',')
        parameter_columns = ['photocurrent', 'saturation_current',
                             'resistance_series', 'resistance_shunt', 'n',
                             'cells_in_series']
        rows = []
        for row in reader:
            rows.append(
                [int(row['Index'])]
                + [mp.mpmathify(row[col], strings=True)
                    for col in parameter_columns]
            )
        return rows


def case_filenames():
    """
    Returns a sorted list of filenames in the directory `TESTS_DIR`.
    The filenames do not have file extensions.

    Returns
    -------
    list
        A sorted list of filenames without file extensions.
    """
    res = set()
    for _, _, filenames in os.walk(TESTS_DIR):
        for filename in filenames:
            res.add(filename.split(".")[0])
    res = list(res)
    res.sort()
    return res


def diff_lhs_rhs(v, i, il, io, rs, rsh, n, vth, ns):
    r"""
    Calculates the difference between the left hand side and right hand side of
    the single diode equation.

    Parameters
    ----------
    v : numeric
        Voltage [V]

    i : numeric
        Current [A]

    il : numeric
        Light-generated current :math:`I_L` (photocurrent) [A]

    io : numeric
        Diode saturation :math:`I_0` current under desired IV curve conditions.
        [A]

    rs : numeric
        Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

    rsh : numeric
        Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
        [ohm]

    n : numeric
        Diode ideality factor :math:`n`

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as
        :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K),
        :math:`T_c` is the temperature of the p-n junction in Kelvin, and
        :math:`q` is the charge of an electron (coulombs). 

    ns : numeric
        Number of cells in series :math:`N_s`

    Returns
    -------
    mpmath float
        Difference between the left hand and right hand sides of the single
        diode equation.
    """
    return (il - io*mp.expm1((v + i*rs)/(n*ns*vth)) - (v + i*rs)/(rsh) - i)


set_globals()


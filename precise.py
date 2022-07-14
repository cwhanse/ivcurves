import pvlib
import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp
import itertools
import json
import csv
import datetime
import max_power

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


def get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts):
    r"""
    Calculates precise solutions to the single diode equation for the given
    parameters, with an error of at most `atol`.

    Parameters
    ----------
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

    atol : numeric
        The error of each of the returned solution pairs is at most `atol`.

        Each returned voltage, current pair is a solution to the single diode
        equation. Because we are working with inexact numbers, these pairs are
        rarely exact solutions to this equation. Here, an error of at most
        `atol` means that for a given :math:`(V, I)` pair,

        .. math:: 

         atol > I_L - I_0 \left[ \exp \left( \frac{V+I R_s}{n N_s V_{th}}
         \right) - 1 \right] - \frac{V + I R_s}{R_{sh}} - I. 

    num_pts : int
        Number of points calculated on IV curve.

    Returns
    -------
    (vv, precise_i) : tuple of numpy arrays
        `vv` and `precise_i` are numpy arrays of mpmath floats of
        length `num_pts`.

    Notes
    -----
    Uses pvlib.pvsystem.singlediode to generate solution pairs, then uses
    max_power's `lambert_i_from_v` to sharpen the precision of the solutions
    if necessary.
    """
    # convert mpf to np.float64
    parameters_npfloat64 = map(lambda x: np.float64(x), [il, io, rs, rsh, n*vth*ns])
    res = pvlib.pvsystem.singlediode(*parameters_npfloat64, ivcurve_pnts=num_pts)

    # convert np.float64 to mpf
    vv = np.fromiter(map(mp.mpmathify, res['v']), dtype=mp.mpf)
    ii = np.fromiter(map(mp.mpmathify, res['i']), dtype=mp.mpf)

    # allocate array for new, more precise i's
    precise_i = np.zeros(ii.shape[0], dtype=mp.mpf)
    i_precise_enough = lambda c: abs(diff_lhs_rhs(v, c, il, io, rs, rsh, n, vth, ns)) < atol
    for idx, (v, i) in enumerate(zip(vv, ii, strict=True)):
        # check if i val already precise enough
        if i_precise_enough(i): 
            new_i = i
        else:
            new_i = max_power.lambert_i_from_v(v, il, io, rs, rsh, n, vth, ns)
            assert i_precise_enough(new_i), f'Index: {idx}'

        # updating array of i's
        precise_i[idx] = new_i

    # find precise v_oc and set as last coordinate
    precise_i[-1] = mp.mpmathify(0)
    vv[-1] = max_power.lambert_v_from_i(precise_i[-1], il, io, rs, rsh, n, vth,
                                        ns)

    assert vv[0] == 0, f'Must be zero: vv[0] = {vv[0]}'
    assert precise_i[-1] == 0, f'Must be zero: precise_i[-1] = {precise_i[-1]}'

    return vv, precise_i


def plotter(il, io, rs, rsh, n, vth, ns, atol, num_pts, case_title):
    # plot a single IV curve 
    plt.xlabel('Voltage')
    plt.ylabel('Current')

    plt.title(case_title)

    v_vals, i_vals = get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
    return plt.plot(v_vals, i_vals)


def plot_case_iv_curves(case_title, case_parameter_sets, vth, atol, num_pts):
    plt.style.use('seaborn-darkgrid')
    plot = plt.plot()
    for _, il, io, rs, rsh, n, ns in case_parameter_sets:
        plot += plotter(il, io, rs, rsh, n, vth, ns, atol, num_pts, case_title)
    plt.show()


def read_case_parameters(filename):
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


def write_case_tests(case_filename, case_parameter_sets, vth, temp_cell, atol,
                     num_pts):
    case_test_suite = {'Manufacturer': '', 'Sandia ID': '', 'Material': '',
                       'IV Curves': []}
    for test_idx, il, io, rs, rsh, n, ns in case_parameter_sets:
        vv, ii = get_precise_i(il, io, rs, rsh, n, vth, ns, atol,
                               num_pts)
        v_oc = vv.max()
        i_sc = ii.max()
        v_mp, i_mp, p_mp = max_power.max_power_pt_finder(il, io, rs, rsh, n,
                                                         vth, ns, atol)
        precision = 16
        # find max digits to the left of decimal
        all_mpf = itertools.chain(vv, ii, [v_oc, i_sc, v_mp, i_mp, p_mp])
        # force mpf to string in decimal format, no exponent
        mpf_to_str = lambda x: mp.nstr(x, n=precision*2, min_fixed=-mp.inf,
                                       max_fixed=mp.inf)
        ldigits = max(mpf_to_str(abs(num)).find('.') for num in all_mpf)

        nstr16 = lambda x: mp.nstr(x, n=ldigits+precision)
        vv_str_list = [nstr16(x) for x in vv]
        ii_str_list = [nstr16(x) for x in ii]

        case_test_suite['IV Curves'].append({
            'Index': test_idx, 'Voltages': vv_str_list,
            'Currents': ii_str_list, 'v_oc': nstr16(v_oc),
            'i_sc': nstr16(i_sc), 'v_mp': nstr16(v_mp),
            'i_mp': nstr16(i_mp), 'p_mp': nstr16(p_mp),
            'Temperature': mp.nstr(temp_cell, n=5), 'Irradiance': None,
            'Sweep direction': None, 'Datetime': None
        })

    with open(f'{case_filename}.json', 'w') as file:
        json.dump(case_test_suite, file, indent=2)


if __name__ == "__main__":
    mp.dps = 40 # set precision, 16*2 rounded up
    num_pts = 100
    atol = mp.mpmathify(1e-16)

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = map(mp.mpmathify, [1.380649e-23, 1.60217663e-19, 298.15])
    vth = (k * temp_cell) / q

    case_number = 1
    case_filename = f'tests/case{case_number}'
    case_title = f'Case {case_number}'
    case_parameter_sets = read_case_parameters(case_filename)

    write_case_tests(case_filename, case_parameter_sets, vth, temp_cell, atol,
                     num_pts)
    plot_case_iv_curves(case_title, case_parameter_sets, vth, atol, num_pts)


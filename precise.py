import pvlib
import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp
from itertools import product
import json
import csv
import datetime


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
    (vv, prec_i) : tuple of numpy arrays
        `vv` is a numpy array of float64 and `prec_i` is a numpy array of
        mpmath floats. Each array has `num_pts` entries.

    Notes
    -----
    Uses pvlib.pvsystem.singlediode to generate solution pairs, then uses
    mpmath's findroot to sharpen the precision of the solutions if necessary.
    """
    res = pvlib.pvsystem.singlediode(il, io, rs, rsh, n*vth*ns, num_pts)
    vv = res['v']
    ii = res['i']
    # initialize 'empty' array for new, more precise i's
    prec_i = np.array([0 for _ in range(ii.shape[0])], dtype=mp.mpf) 

    assert len(vv) == len(ii)
    for idx in range(len(vv)):
        i = mp.mpf(str(ii[idx]))
        # check if i val already precise enough
        if abs(diff_lhs_rhs(vv[idx], i, il, io, rs, rsh, n, vth, ns)) < atol: 
            new_i = i

        else:
            # findroot takes the function whose roots we want to find, a good
            # guess for where the root is, and a tolerance
            # default solver uses secant method
            new_i = mp.findroot(lambda x: diff_lhs_rhs(vv[idx], x, il, io, rs, rsh, n, vth, ns), i, tol=atol**2)
            # setting tol=atol**2 because findroot checks func(zero)**2 < tol

        # check that mp.findroot did what we wanted 
        assert abs(diff_lhs_rhs(vv[idx], new_i, il, io, rs, rsh, n, vth, ns)) < atol

        # updating array of i's
        prec_i[idx] = new_i

    # vv is a np.array of float64 and prec_i is a np.array of mp.mpf 
    return vv, prec_i


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
    for il, io, rs, rsh, n, ns in case_parameter_sets:
        plot += plotter(il, io, rs, rsh, n, vth, ns, atol, num_pts, case_title)
    plt.show()


def read_case_parameters(filename):
    with open(f'{filename}.csv', newline='') as file:
        reader = csv.DictReader(file, delimiter=',')
        columns = ['IL', 'IO', 'RS', 'RSH', 'N', 'ns']
        rows = []
        for row in reader:
            rows.append([float(row[col]) for col in columns])
        return rows


def write_case_tests(case_filename, case_parameter_sets, vth, atol, num_pts):
    case_test_suite = {'Manufacturer': '', 'Sandia ID': '', 'Material': '',
                       'IV Curves': []}
    for idx, (il, io, rs, rsh, n, ns) in enumerate(case_parameter_sets):
        voltages, currents = get_precise_i(il, io, rs, rsh, n, vth, ns, atol,
                                           num_pts)
        nstr16 = lambda x: mp.nstr(x, n=16)
        voltages_str_list = [nstr16(x) for x in voltages]
        currents_str_list = [nstr16(x) for x in currents]

        # find V_mp, I_mp, where P_mp = V_mp * I_mp is maximized
        V_mp, I_mp = voltages[0], currents[0]
        P_mp = V_mp * I_mp
        for v, c in zip(voltages, currents):
            if P_mp < v * c:
                V_mp, I_mp = v, c
                P_mp = v * c

        case_test_suite['IV Curves'].append({
            'Index': idx + 1, 'Voltages': voltages_str_list,
            'Currents': currents_str_list, 'V_oc': nstr16(voltages.max()),
            'I_sc': nstr16(currents.max()), 'V_mp': nstr16(V_mp),
            'I_mp': nstr16(I_mp), 'P_mp': nstr16(P_mp), 'Temperature': None,
            'Irradiance': None, 'Sweep direction': '',
            'Datetime': datetime.date.today().strftime('%m/%d/%Y')
        })

    with open(f'{case_filename}.json', 'w') as file:
        json.dump(case_test_suite, file, indent=2)


if __name__ == "__main__":
    mp.dps = 40 # set precision, 16*2 rounded up
    num_pts = 100
    atol = 1e-16

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = [1.380649e-23, 1.60217663e-19, 298.15]
    vth = (k * temp_cell) / q

    case_number = 1
    case_filename = f'case{case_number}'
    case_title = f'Case {case_number}'
    case_parameter_sets = read_case_parameters(case_filename)

    write_case_tests(case_filename, case_parameter_sets, vth, atol, num_pts)
    # plot_case_iv_curves(case_title, case_parameter_sets, vth, atol, num_pts)


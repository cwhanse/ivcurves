import argparse
import json
import pvlib
import numpy as np
import matplotlib.pyplot as plt

from ivcurves.utils import mp  # same instance of mpmath's mp imported in ivcurves/utils
import ivcurves.utils as utils


###################
# Max power point #
###################


def max_power_pt_finder(il, io, rs, rsh, n, vth, ns, atol):
    r"""
    Calculates power at maximum power point.

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

    atol : float
        The returned voltage will be at most ``atol`` from the voltage that
        yields the true maximum power of the given IV curve.

    Returns
    -------
    (max_voltage, max_current, max_power) : tuple of mpmath floats
        A good approximation for the voltage :math:`V`, current :math:`I`, and
        power :math:`P` of the maximum power point of the given IV curve.
    """
    params = il, io, rs, rsh, n, vth, ns
    v_oc = lambert_v_from_i(0, *params)
    if v_oc < atol:
        # means v_oc ~ 0 which causes problems when defining max_iters
        # this should only ever happen for very extreme parameters
        assert mp.chop(v_oc) == 0
        return 0, 0, 0

    max_iters = int(1 + mp.log(atol / v_oc) / mp.log((mp.sqrt(5) - 1) / 2))

    power_func = lambda v : v * lambert_i_from_v(v, *params)

    max_voltage = golden_search(0, v_oc, power_func, atol, max_iters)
    max_current = lambert_i_from_v(max_voltage, *params)

    # check max current and voltage is a precise solution to single diode eq
    diff = lambda v, i: diff_lhs_rhs(v, i, *params)
    if abs(diff(max_voltage, max_current)) >= atol:
        max_current = mp.findroot(lambda i: diff(max_voltage, i), max_current, tol=atol**2)
        # setting tol=atol**2 because findroot checks func(zero)**2 < tol

    max_power = max_current * max_voltage
    return max_voltage, max_current, max_power


########################
# Golden search method #
########################


def golden_search(a, b, func, atol, max_iters):
    r"""
    Finds a local maximizer of a function on an interval :math:`[a, b]` with at
    most ``atol`` error using golden-section serach.

    Parameters
    ----------
    a : float
        Left endpoint of interval.

    b : float
        Right endpoint of interval.

    func : function
        Single-variable function to maximize.

    atol : float
        The absolute tolerance between the true and calculated maximizer.

    max_iters : int
        Maximum number of iterations golden-section search before failing.

    Returns
    -------
    float
        The calculated maximizer.

    Notes
    -----
    For more information on the algorithm, see
    http://www.math.kent.edu/~reichel/courses/intr.num.comp.2/lecture16/lecture8.pdf.
    """
    rho = (1/2) * (3 - mp.sqrt(5))  # this value for rho is equivalent to using golden ratio
    for _ in range(max_iters):
        x_internal = lambda frac: a + frac * (b - a)
        xl = x_internal(rho)
        xr = x_internal(1 - rho)
        if func(xl) > func(xr):
            if abs(xr - a) < atol:
                return xl
            b = xr
        else:
            if abs(b - xl) < atol:
                return xr
            a = xl
    raise RuntimeError('Golden Search: maximum iteration count exceeded.')


#######################
# Auxiliary functions #
#######################


def lambert_i_from_v(v, il, io, rs, rsh, n, vth, ns):
    r"""
    Given a voltage, calculates the associated current using the Lambert W
    function.

    Parameters
    ----------
    v : numeric
        Voltage [V]

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
        Current associated to given voltage via the single diode equation.

    Notes
    -----
    See :external+pvlib:func:`pvlib.pvsystem.singlediode` for original function.
    This implemention differs only in that it does not accept Series as inputs
    and can take mpmath floats as inputs.
    """
    gsh = 1. / rsh
    if rs == 0:
        return il - io*mp.expm1(v / (n*vth*ns)) - gsh*v
    else:
        # handle overflow FIXME
        argW = rs * io / ((n*vth*ns) * (rs*gsh + 1)) * mp.exp((rs*(il + io) + v)/((n*vth*ns)*(rs*gsh + 1)))
        lambertw_term = mp.lambertw(argW).real
        return (il + io - v*gsh) / (rs*gsh + 1) - ((n*vth*ns) / rs)*lambertw_term


def lambert_v_from_i(i, il, io, rs, rsh, n, vth, ns):
    r"""
    Given a current, calculates the associated voltage using the Lambert W
    function.

    Parameters
    ----------
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
        Voltage associated to given current via the single diode equation.

    Notes
    -----
    See :external+pvlib:func:`pvlib.pvsystem.singlediode` for original function.
    This implemention differs only in that it does not accept Series as inputs
    and can take mpmath floats as inputs.
    """
    gsh = 1. / rsh
    if gsh == 0:
        return (n*vth*ns) * mp.log1p((il - i) / io) - i*rs
    else:
        # handle overflow FIXME
        argW = io / (gsh * (n*vth*ns)) * mp.exp((-i + il + io) / (gsh * (n*vth*ns)))
        lambertw_term = mp.lambertw(argW).real
        return (il + io - i) / gsh - i*rs - (n*vth*ns)*lambertw_term


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
    parameters, with an error of at most ``atol``.

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
        The error of each of the returned solution pairs is at most ``atol``.

        Each returned voltage, current pair is a solution to the single diode
        equation. Because we are working with inexact numbers, these pairs are
        rarely exact solutions to this equation. Here, an error of at most
        ``atol`` means that for a given :math:`(V, I)` pair,

        .. math::

           \left|
              I_L - I_0 \left[
                 \exp \left( \frac{V+I R_s}{n N_s V_{th}} \right) - 1
              \right] - \frac{V + I R_s}{R_{sh}} - I
           \right| < \text{atol}

    num_pts : int
        Number of points calculated on IV curve.

    Returns
    -------
    (vv, precise_i) : tuple of numpy arrays
        ``vv`` and ``precise_i`` are numpy arrays of mpmath floats of
        length ``num_pts``.

    Notes
    -----
    Uses :external+pvlib:func:`pvlib.pvsystem.singlediode` to generate solution
    pairs, then uses :func:`lambert_i_from_v` to sharpen the precision of the
    solutions if necessary.
    """
    # convert mpf to np.float64
    params_npfloat64 = map(lambda x: np.float64(x), [il, io, rs, rsh, n*vth*ns])
    res = pvlib.pvsystem.singlediode(*params_npfloat64, ivcurve_pnts=num_pts)

    # convert np.float64 to mpf
    vv = np.fromiter(map(mp.mpmathify, res['v']), dtype=mp.mpf)
    ii = np.fromiter(map(mp.mpmathify, res['i']), dtype=mp.mpf)

    # allocate array for new, more precise i's
    precise_i = np.zeros(ii.shape[0], dtype=mp.mpf)
    params = il, io, rs, rsh, n, vth, ns
    diff = lambda i: diff_lhs_rhs(v, i, *params)
    for idx, (v, i) in enumerate(zip(vv, ii)):
        # check if i val already precise enough
        if abs(diff(i)) < atol:
            new_i = i
        else:
            new_i = lambert_i_from_v(v, *params)
            assert abs(diff(new_i)) < atol, f'Index: {idx}'

        # updating array of i's
        precise_i[idx] = new_i

    # find precise v_oc and set as last coordinate
    precise_i[-1] = mp.mpmathify(0)
    vv[-1] = lambert_v_from_i(precise_i[-1], *params)

    assert vv[0] == 0, f'Must be zero: vv[0] = {vv[0]}'
    assert precise_i[-1] == 0, f'Must be zero: precise_i[-1] = {precise_i[-1]}'

    return vv, precise_i


def plot_iv_curves(test_set_filename, case_parameter_sets, vth, atol, num_pts,
                   show=True, savefig=False, stack_plots=True):
    """
    Plot IV curves.

    Parameters
    ----------
    test_set_filename : str
        The name of the test set CSV file that ``case_parameter_sets`` was
        loaded from. The filename must exclude its file extension.

    case_parameter_sets : dict
        A mapping of test case indices to a list of test case parameters.
        (See :func:`utils.read_iv_curve_parameter_sets`)

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as
        :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K),
        :math:`T_c` is the temperature of the p-n junction in Kelvin, and
        :math:`q` is the charge of an electron (coulombs).

    atol : numeric
        The absolute tolerance allowed when generating the IV curve data
        from the test case's parameters.

    num_pts : int
        Number of points calculated on IV curve.

    show : bool, default True
        Display the plot.

    savefig : bool, default False
        Save the plot as a PNG image.

    stack_plots : bool, default True
        Plot all of the IV curves generated from the parameters in
        ``case_parameter_sets`` on a single plot. If ``False``, only
        the last IV curve generated will be plotted.
    """
    plt.style.use('seaborn-darkgrid')
    for idx, (il, io, rs, rsh, n, ns) in case_parameter_sets.items():
        params = il, io, rs, rsh, n, vth, ns
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        v_vals, i_vals = get_precise_i(*params, atol=atol, num_pts=num_pts)
        plt.plot(v_vals, i_vals)
        if savefig:
            plt.savefig(utils.make_iv_curve_name(test_set_filename, idx),
                        bbox_inches='tight')
        if not stack_plots:
            plt.cla()
    if show:
        plt.show()


def build_test_set_json(case_parameter_sets, vth, temp_cell, atol, num_pts,
                        test_set_json=None):
    """
    Builds a dict of IV curve data compliant with the IV Curve JSON schema.

    Parameters
    ----------
    case_parameter_sets : dict
        A mapping of test case indices to a list of test case parameters.

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as
        :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K),
        :math:`T_c` is the temperature of the p-n junction in Kelvin, and
        :math:`q` is the charge of an electron (coulombs).

    temp_cell : numeric
        The test set's cell temperature.

    atol : numeric
        The absolute tolerance allowed when generating the IV curve data
        from the test case's parameters.

    num_pts : int
        Number of points calculated on IV curve.

    Returns
    -------
        cells_in_series, dict
    """
    test_set_json_template = {
        'Manufacturer': '', 'Model': '', 'Serial Number': '',
        'Module ID': '',  'Description': '', 'Material': '',
        'cells_in_series': None, 'IV Curves': []
    }
    if test_set_json is None:
        test_set_json = test_set_json_template
    else:
        # make sure the json only has the keys listed in the template
        for k in test_set_json_template:
            test_set_json_template[k] = test_set_json[k]
        test_set_json = test_set_json_template

    ivcurves = []

    # this is set in the loop below, and is the same for every iv curve
    cells_in_series = None

    for test_idx, (il, io, rs, rsh, n, ns) in case_parameter_sets.items():
        params = il, io, rs, rsh, n, vth, ns
        cells_in_series = int(params[-1])
        rs = params[2]

        vv, ii = get_precise_i(*params, atol=atol, num_pts=num_pts)
        v_oc = vv[-1]
        i_sc = ii[0]
        v_mp, i_mp, p_mp = max_power_pt_finder(*params, atol=atol)
        i_x = lambert_i_from_v(v_oc / 2, *params)
        i_xx = lambert_i_from_v((v_oc + v_mp) / 2, *params)

        nstr = utils.mp_nstr_precision_func
        vv_list = [nstr(x) for x in vv]
        ii_list = [nstr(x) for x in ii]
        diode_voltage_list = [nstr(dv) for dv in vv + rs * ii]
        ivcurves.append({
            'Index': test_idx, 'Voltages': vv_list,
            'Currents': ii_list, 'diode_voltage': diode_voltage_list,
            'v_oc': nstr(v_oc), 'i_sc': nstr(i_sc),
            'v_mp': nstr(v_mp), 'i_mp': nstr(i_mp), 'p_mp': nstr(p_mp),
            'i_x': nstr(i_x), 'i_xx': nstr(i_xx),
            'Temperature': mp.nstr(temp_cell, n=5), 'Irradiance': None,
            'Sweep direction': '', 'Datetime': '1970-01-01T00:00:00Z'
        })

    test_set_json['cells_in_series'] = cells_in_series
    test_set_json['IV Curves'] = ivcurves

    return test_set_json


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Generates precise IV curve data from the parameters of '
                    'the single diode equation.'
    )
    parser.add_argument('--test-set', dest='test_set_filename', type=str,
                        help='Test set filename (excluding file extension) to '
                             'generate curves for. If omitted, all test sets are used.')
    parser.add_argument('--save-json', dest='save_json_path', type=str,
                        help='Saves the test set JSON at the given path.')
    parser.add_argument('--save-images', dest='save_images_path', type=str,
                        help='Saves the test set IV curve plots at the given path.')
    parser.add_argument('--plot', action=argparse.BooleanOptionalAction,
                        help="Plot each test set's IV curves")
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()

    if args.test_set_filename:
        test_set_filenames = [args.test_set_filename]
    else:
        test_set_filenames = utils.get_filenames_in_directory(utils.TEST_SETS_DIR)

    constants = utils.constants()
    vth, temp_cell, atol, num_pts = (constants['vth'], constants['temp_cell'],
                                     constants['atol'], constants['num_pts'])
    for name in test_set_filenames:
        case_parameter_sets = utils.read_iv_curve_parameter_sets(utils.TEST_SETS_DIR / name)
        if args.save_json_path:
            test_set_json = None
            test_set_json_path = utils.TEST_SETS_DIR / f'{name}.json'
            if test_set_json_path.exists():
                test_set_json = utils.load_json(test_set_json_path)
                test_set_json = build_test_set_json(
                    name, case_parameter_sets, vth, temp_cell, atol, num_pts,
                    test_set_json=test_set_json
                )
                utils.save_json(test_set_json, args.save_json_path / f'{name}.json')
        if args.save_images_path:
            plot_iv_curves(
                f'{args.save_images_path}/{name}',
                case_parameter_sets, vth, atol, num_pts, show=False,
                savefig=True, stack_plots=False
            )
        if args.plot:
            plot_iv_curves(name, case_parameter_sets, vth, atol, num_pts)

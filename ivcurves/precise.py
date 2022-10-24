import pvlib
import numpy as np
import matplotlib.pyplot as plt
from utils import mp # same instance of mpmath's mp imported in ivcurves/utils
import utils
import argparse
import itertools
import json


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
    # make parameters mpf types for precision
    il, io, rs, rsh, n, vth = [mp.mpf(str(val)) for val in [il, io, rs, rsh, n, vth]]

    # function we want to maximize
    # x represents voltage
    power_func = lambda x : x * lambert_i_from_v(x, il, io, rs, rsh, n, vth, ns)

    # find initial interval endpts (using (0, I_sc) and (V_oc, 0) as starting
    # endpts, power=0 at both these points)
    # x-coord is voltage, y-coord is power (power = voltage * current)
    xl, yl = lambert_i_from_v(0, il, io, rs, rsh, n, vth, ns), 0
    xr, yr = lambert_v_from_i(0, il, io, rs, rsh, n, vth, ns), 0

    # run golden search on power function with above starting interval
    if (xr - xl) < atol: # this will cause problems when defining iterlimit
        # means xr == xl (roughly)
        # shouldn't happen unless both xr and xl are basically zero
        assert (mp.chop(xr) == 0 and mp.chop(xl) == 0)
        # this should only ever happen for very extreme parameters
        return 0, 0, 0
    else:
        iterlimit = 1 + mp.floor( mp.log(atol/(xr - xl)) / mp.log((mp.sqrt(5) - 1) / 2) )

    max_voltage, max_power = golden_search((xl, yl), (xr, yr), power_func, atol, iterlimit)

    # find current associated to max_voltage
    max_current = lambert_i_from_v(max_voltage, il, io, rs, rsh, n, vth, ns)

    # check that current, voltage pair is still a precise solution to single diode eq
    # if not precise enough, make precise using findroot
    dff = diff_lhs_rhs(max_voltage, max_current, il, io, rs, rsh, n, vth, ns)
    if abs(dff) > atol:
        max_current = mp.findroot(lambda x: diff_lhs_rhs(max_voltage, x, il, io, rs, rsh, n, vth, ns), max_current, tol=atol**2)
        # setting tol=atol**2 because findroot checks func(zero)**2 < tol

    return max_voltage, max_current, max_power



########################
# Golden search method #
########################


def golden_search(l_endpt, r_endpt, func, atol, iterlimit, int_pt=tuple(), is_right_int_pt=False, num_iter=0):
    r"""
    Uses golden-section search to recursively find maximum of the given
    function over the given interval, with at most ``atol`` error.

    Parameters
    ----------
    l_endpt : tuple of floats
        Left endpoint of interval in which a maximum occurs.

    r_endpt : tuple of floats
        Right endpoint of interval in which a maximum occurs.

    func : function
        Function we want to maximize (single-variable).

    atol : float
        The x-coordinate of the returned point will be at most ``atol`` from the
        x-coordinate that produces the true maximum in the given interval.

    iterlimit : int
        Maximum number of iterations for golden-section search. Should converge
        before we hit this.

    int_pt : tuple of floats, optional
        Coordinates of interior point in interval. The default value is an
        empty tuple.

    is_right_int_pt : bool, optional
        If ``int_pt`` is the right hand interior point, then True. If ``int_pt`` is
        given, this should also be passed in (default value is False).

    num_iter : int
        The number of iterations we've already done.

    Returns
    -------
    tuple of ints
        Coordinates of a point whose y-coordinate is a good approximation for
        the true maximum, and whose x-coordinate is within ``atol`` of the
        x-coordinate that yields the true maximum of ``func`` on the given
        interval.

    Notes
    -----
    This is a recursive function. When using, ``int_pt`` and ``is_right_int_pt``
    and ``num_iter`` should not be passed; they are only passed when the function
    recurses.

    For more information on the algorithm (and calculating the interior points
    in :func:`get_left_int_pt` and :func:`get_right_int_pt`), see
    http://www.math.kent.edu/~reichel/courses/intr.num.comp.2/lecture16/lecture8.pdf.
    """
    # overflow ? FIXME
    # take care of f(int_pt_1) == f(int_pt_2) (right now just pushing this case
    # into the f(int_pt_1) < f(int_pt_2) case) FIXME

    if num_iter >= iterlimit: raise Exception("Iterations exceeded maximum.")

    xl, yl = l_endpt
    xr, yr = r_endpt

    # find left and right interior points
    if int_pt == tuple(): # first iteration, no interior points yet
        l_int_pt = get_left_int_pt(xl, xr, func)
        r_int_pt = get_right_int_pt(xl, xr, func)
    else: # already have one interior point
        if is_right_int_pt:
            r_int_pt = int_pt
            l_int_pt = get_left_int_pt(xl, xr, func)
        else:
            l_int_pt = int_pt
            r_int_pt = get_right_int_pt(xl, xr, func)

    # recurse with new (smaller) interval and one interior point
    if l_int_pt[1] > r_int_pt[1]: # check which has higher y-coord (since we're searching for max)
        error = abs(r_int_pt[0] - l_endpt[0])
        if error < atol:
            return l_int_pt
        else:
            return golden_search(l_endpt, r_int_pt, func, atol, iterlimit, l_int_pt, is_right_int_pt=True, num_iter=num_iter+1)
            # `is_right_int_pt`=True because l_int_pt is now the right_int_pt of new interval
    else:
        error = abs(r_endpt[0] - l_int_pt[0])
        if error < atol:
            return r_int_pt
        else:
            return golden_search(l_int_pt, r_endpt, func, atol, iterlimit, r_int_pt, is_right_int_pt=False, num_iter=num_iter+1)
            # `is_right_int_pt`=False because r_int_pt is now the left_int_pt of new interval


def get_left_int_pt(left_x_endpt, right_x_endpt, func):
    r"""
    Calculates the left interior point of the interval.

    This is an auxiliary function for :func:`golden_search`.

    Parameters
    ----------
    left_x_endpt : float
        x-coordinate of the left endpoint.

    right_x_endpt : float
        x-coordinate of the right endpoint.

    func : function
        Function we want to maximize (single-variable).

    Returns
    -------
    tuple of mpmath floats
        Coordinate for left interior point.
    """
    xl, xr = left_x_endpt, right_x_endpt
    rho = (1/2) * (3 - mp.sqrt(5))
    # this value for rho is equivalent to using golden ratio

    l_int_x = xl + rho*(xr - xl) # voltage
    l_int_y = func(l_int_x)
    return (l_int_x, l_int_y)


def get_right_int_pt(left_x_endpt, right_x_endpt, func):
    r"""
    Calculates the right interior point of the interval.

    This is an auxiliary function for :func:`golden_search`.

    Parameters
    ----------
    left_x_endpt : float
        x-coordinate of the left endpoint.

    right_x_endpt : float
        x-coordinate of the right endpoint.

    func : function
        Function we want to maximize (single-variable).

    Returns
    -------
    tuple of mpmath floats
        Coordinate for right interior point.
    """
    xl, xr = left_x_endpt, right_x_endpt
    rho = (1/2) * (3 - mp.sqrt(5))
    # this value for rho is equivalent to using golden ratio

    r_int_x = xl + (1 - rho)*(xr - xl) # voltage
    r_int_y = func(r_int_x)
    return (r_int_x, r_int_y)



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
    parameters_npfloat64 = map(lambda x: np.float64(x), [il, io, rs, rsh, n*vth*ns])
    res = pvlib.pvsystem.singlediode(*parameters_npfloat64, ivcurve_pnts=num_pts)

    # convert np.float64 to mpf
    vv = np.fromiter(map(mp.mpmathify, res['v']), dtype=mp.mpf)
    ii = np.fromiter(map(mp.mpmathify, res['i']), dtype=mp.mpf)

    # allocate array for new, more precise i's
    precise_i = np.zeros(ii.shape[0], dtype=mp.mpf)
    i_precise_enough = lambda c: abs(diff_lhs_rhs(v, c, il, io, rs, rsh, n, vth, ns)) < atol
    for idx, (v, i) in enumerate(zip(vv, ii)):
        # check if i val already precise enough
        if i_precise_enough(i):
            new_i = i
        else:
            new_i = lambert_i_from_v(v, il, io, rs, rsh, n, vth, ns)
            assert i_precise_enough(new_i), f'Index: {idx}'

        # updating array of i's
        precise_i[idx] = new_i

    # find precise v_oc and set as last coordinate
    precise_i[-1] = mp.mpmathify(0)
    vv[-1] = lambert_v_from_i(precise_i[-1], il, io, rs, rsh, n, vth,
                                        ns)

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
        plt.xlabel('Voltage')
        plt.ylabel('Current')
        v_vals, i_vals = get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
        plt.plot(v_vals, i_vals)
        if savefig:
            plt.savefig(utils.make_iv_curve_name(test_set_filename, idx),
                        bbox_inches='tight')
        if not stack_plots:
            plt.cla()
    if show:
        plt.show()


def build_test_set_json(test_set_filename, case_parameter_sets, vth, temp_cell,
                        atol, num_pts):
    """
    Write JSON files of IV curve data.

    Parameters
    ----------
    test_set_filename : str
        The filename of the test set CSV file that ``case_parameter_sets`` was
        loaded from. The filename must exclude its file extension.

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
    """
    case_test_suite = {'Manufacturer': '', 'Sandia ID': '', 'Material': '',
                       'IV Curves': []}
    for test_idx, (il, io, rs, rsh, n, ns) in case_parameter_sets.items():
        vv, ii = get_precise_i(il, io, rs, rsh, n, vth, ns, atol,
                               num_pts)
        v_oc = vv.max()
        i_sc = ii.max()
        v_mp, i_mp, p_mp = max_power_pt_finder(il, io, rs, rsh, n,
                                               vth, ns, atol)

        nstr = utils.mp_nstr_precision_func
        vv_str_list = [nstr(x) for x in vv]
        ii_str_list = [nstr(x) for x in ii]
        case_test_suite['IV Curves'].append({
            'Index': test_idx, 'Voltages': vv_str_list,
            'Currents': ii_str_list, 'v_oc': nstr(v_oc),
            'i_sc': nstr(i_sc), 'v_mp': nstr(v_mp),
            'i_mp': nstr(i_mp), 'p_mp': nstr(p_mp),
            'cells_in_series': int(ns),
            'Temperature': mp.nstr(temp_cell, n=5), 'Irradiance': None,
            'Sweep direction': "", 'Datetime': ""
        })


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
        case_parameter_sets = utils.read_iv_curve_parameter_sets(f'{utils.TEST_SETS_DIR}/{name}')
        if args.save_json_path:
            with open(f'{test_set_filename}.json', 'w') as file:
                json.dump(build_test_set_json(f'{args.save_json_path}/{name}', case_parameter_sets, vth, temp_cell, atol, num_pts), file, indent=2)
        if args.save_images_path:
            plot_iv_curves(f'{args.save_images_path}/{name}',
                           case_parameter_sets, vth, atol, num_pts, show=False,
                           savefig=True, stack_plots=False)
        if args.plot:
            plot_iv_curves(name, case_parameter_sets, vth, atol, num_pts)


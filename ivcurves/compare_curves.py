import argparse
import csv
from pathlib import Path
import numpy as np

# from ivcurves repo
# same instance of mpmath's mp imported in ivcurves/utils
from ivcurves.utils import mp
import ivcurves.utils as utils
import ivcurves.precise as precise


#####################
# Find intersection #
#####################

def _find_x_intersection(interval, func, point, atol, maxsteps=100):
    r"""
    Finds x-coordinate of the intersection between the known IV curve and the
    line segment from the origin to the given point.

    This is an auxiliary function for :func:`total_score`.

    Parameters
    ----------
    interval : tuple(float, float)
        An interval on the x-axis that contains an intersection.

    func : function
        A two-variable scalar-valued function.

    point : tuple(float, float)
        The point defining the line from the origin.

    atol : float
        Absolute tolerance of how close the found x-coordiante must be to the
        true x-coordinate of the intersection point.

    maxsteps : int, default 100
        Maximum number of iterations for :func:`ivcurves.precise.golden_search`
        and :func:`mp.findroot` to find the x-coordinate of the intersection
        point. A ``ValueError`` by either iteration if ``atol`` is not met
        after ``maxsteps`` iterations.

    Returns
    -------
    float
        x-coordinate of an intersection.
    """
    xp, yp = point
    if xp == 0:
        # then line through origin and (`xp`, `yp`) sits on y-axis
        # so x coordinate at intersection must be zero
        return 0
    else:
        # line through origin and (`xp`, `yp`) has a defined slope
        line = lambda x: (yp / xp) * x

        # solve for intersection of line and single_diode
        solve_for_zero = lambda x: -abs(func(x, line(x)) - line(x))
        a, b = interval
        # golden_search finds a local maximum
        x_int = precise.golden_search(a, b, solve_for_zero, atol, maxsteps)
        try:
            # setting tol=atol**2 because findroot checks |func(zero)|**2 < tol
            x_int = mp.findroot(solve_for_zero, x_int, tol=atol**2,
                                maxsteps=maxsteps)
        except ValueError as e:
            raise ValueError("Can't find an intersection point. "
                             'Perhaps the curves are too far from each other? '
                             f'{e}')

        assert abs(solve_for_zero(x_int)) < atol

    return x_int


######################
# Calculate distance #
######################

def _find_distance(x, y, xp, yp):
    r"""
    Calculates the distance between two points :math:`(x, y)` and
    :math:`(x_p, y_p)` using a scaled Euclidean distance.

    The scaled Euclidean distance is defined as

    .. math::

       \sqrt{\left(\frac{x_p - x}{x}\right)^2 + \left(\frac{y_p - y}{y}\right)^2}

    If :math:`x` is zero, it is

    .. math::

       \sqrt{\left(\frac{y_p - y}{y}\right)^2}

    If :math:`y` is zero, it is

    .. math::

       \sqrt{\left(\frac{x_p - x}{x}\right)^2}

    It is assumed that :math:`(x, y) \neq 0`.

    Parameters
    ----------
    x : float
        x-coordinate of point on known curve.

    y : float
        y-coordinate of point on known curve.

    xp : float
        x-coordinate of point on fitted curve.

    yp : float
        y-coordinate of point on fitted curve.

    Returns
    -------
    mpmath float
        The calculated distance.
    """
    assert not (x == 0 and y == 0)  # this should never happen for these curves

    diff_x, diff_y = 0, 0
    if x != 0:
        diff_x = (xp - x) / x
    if y != 0:
        diff_y = (yp - y) / y

    return mp.sqrt(diff_x**2 + diff_y**2)


# Final score  per IV curve #

def score_curve(known_curve_params, fitted_curve_params, vth, num_pts, atol):
    r"""
    Calculates the total score for a given fitted curve.

    This score encodes how good an approximation the fitted curve is for the
    curve with the known parameters. If the score is small, then the fitted
    curve is close to the known curve.

    Parameters
    ----------
    known_curve_params : list
        A list of parameters representing a given IV curve. The list items
        should be in the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions.
            [ohm]

        rsh : numeric
            Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
            [ohm]

        n : numeric
            Diode ideality factor :math:`n`

        ns : numeric
            Number of cells in series :math:`N_s`

    fitted_curve_params : list
        A list of parameters representing a given IV curve. Should be passed in
        the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions.
            [ohm]

        rsh : numeric
            Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
            [ohm]

        n : numeric
            Diode ideality factor :math:`n`

        ns : numeric
            Number of cells in series :math:`N_s`

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as
        :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K),
        :math:`T_c` is the temperature of the p-n junction in Kelvin, and
        :math:`q` is the charge of an electron (coulombs).

    num_pts : int
        Number of points we want to compare between the two curves.

    atol : float
        The error of each of the solution pairs found is at most ``atol``. (See
        :func:`ivcurves.precise.get_precise_i`.)
        Each solution pair is a point on the curve.

    Returns
    -------
    score : mpmath float
        A measure for how close the two inputted curves are to each other.

    Notes
    -----
    The user inputs parameters that they’ve fitted to a particular known IV
    curve. The curve from these fitted parameters is compared to the known IV
    curve. A sampling of points on the fitted curve are chosen. To get an
    associated point on the known curve, we draw a line from the origin that
    passes through the point on the fitted curve. This line will intersect the
    known curve; this point of intersection is what we compare the fitted point
    to. We then find the distance between these two points, using the
    definition of distance given in find_distance. The sum of the distances for
    each pair of associated points is the score.
    """
    # get xs and ys for known and fitted curves
    known_xs, known_ys = calc_precise_curve(
        known_curve_params, vth, num_pts, atol)
    fit_xs, fit_ys = calc_precise_curve(
        fitted_curve_params, vth, num_pts, atol)

    il, io, rs, rsh, n, ns = known_curve_params
    single_diode = lambda v, i: il  - ((v + i*rs) / rsh) - \
        io * mp.expm1((v + i*rs) / (n * ns * vth))

    score = 0

    interval_current = min(known_xs), max(known_xs)
    for v, i in zip(fit_xs, fit_ys):
        # for each point (`v`, `i`) on the fitted curve, find the associated
        # point on known curve (`new_voltage`, `new_current`)
        new_voltage = _find_x_intersection(
            interval_current, single_diode, (v, i), atol)
         # find current associated to new_voltage
        new_current = precise.lambert_i_from_v(
            new_voltage, il, io, rs, rsh, n, vth, ns)

        # if voltage, current pair is not a precise enough solution to the
        # single diode equation, make more precise
        dff = precise.diff_lhs_rhs(
            new_voltage, new_current, il, io, rs, rsh, n, vth, ns)
        if abs(dff) > atol:
            new_current = mp.findroot(lambda y: precise.diff_lhs_rhs(
                new_voltage, y, il, io, rs, rsh, n, vth, ns), new_current,
                tol=atol**2)

        # calculate distance between these points, and add to score
        score += _find_distance(new_voltage, new_current, v, i)

    return score


def score_parameters(known_curve_params, fitted_curve_params):
    r'''
    Calculate score as a sum of absolute relative difference in each parameter

    Parameters
    ----------
    known_curve_params : list
        A list of parameters representing a given IV curve. The list items
        should be in the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions.
            [ohm]

        rsh : numeric
            Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
            [ohm]

        n : numeric
            Diode ideality factor :math:`n`

    fitted_curve_params : list
        A list of parameters representing a given IV curve. Should be passed in
        the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions.
            [ohm]

        rsh : numeric
            Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
            [ohm]

        n : numeric
            Diode ideality factor :math:`n`

    '''
    def _abs_rel_diff(x, y):
        return np.abs(x - y) / x

    score = 0.
    for x, y in zip(known_curve_params, fitted_curve_params):
        score += _abs_rel_diff(x, y)
    return score


def calc_precise_curve(curve_parameters, vth, num_pts, atol):
    r"""
    Gets precise voltage and current pairs for the given curve.

    Parameters
    ----------
    curve_parameters : list
        A list of parameters representing a given IV curve. Should be passed in
        the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions.
            [ohm]

        rsh : numeric
            Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
            [ohm]

        n : numeric
            Diode ideality factor :math:`n`

        ns : numeric
            Number of cells in series :math:`N_s`

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as
        :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K),
        :math:`T_c` is the temperature of the p-n junction in Kelvin, and
        :math:`q` is the charge of an electron (coulombs).

    num_pts : int
        Number of points to calculate on the given curve.

    atol : float
        The error of each of the solution pairs found is at most ``atol``.
        (See :func:`ivcurves.precise.get_precise_i`.) Each solution pair is a
        point on the curve.

    Returns
    -------
    (vv, ii) : tuple of numpy arrays
        ``vv`` is a numpy array of float64 and ``ii`` is a numpy array of
        mpmath floats. Each array has ``num_pts`` entries.
    """
    il, io, rs, rsh, n, ns = curve_parameters
    vv, ii = precise.get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
    return vv, ii


def _get_test_sets_to_score(tests='', directory=None):
    """
    Returns a list of valid test set names (excluding file extensions).

    Parameters
    ----------
    tests : str, default ''
        String composed of comma-separated test set names. If empty, all test
        set names from directory are returned.

    directory : pathlib.Path
        Directory that contains files whose filenames are test set filenames.

    Returns
    -------
    list
        A list of valid test set filenames (excluding file extensions).
    """
    test_set_names = utils.get_filenames_in_directory(utils.TEST_SETS_DIR)
    test_sets_to_score = []
    if tests != '':
        for t in tests.split(','):
            if t not in test_set_names:
                raise ValueError(f"'{t}' is not a test set")
            test_sets_to_score.append(t)
    elif directory is not None:
        # score all found in directory
        filenames = utils.get_filenames_in_directory(directory)
        test_sets_to_score = [f for f in filenames if f in test_set_names]
        test_sets_to_score.sort()
        if not test_sets_to_score:
            raise ValueError(f"no test sets found in '{directory}'")
    else:
        raise ValueError('Must supply either tests or directory')
    return test_sets_to_score


def write_test_set_score_per_curve_csvs(scores, csv_output_path):
    """
    Writes a CSV file containing a score for each test case in every test set.

    Parameters
    ----------
    scores : dict
        Dictionary of test set filenames (excluding file extensions) to a
        dictionary to test case indices to test case scores.

    csv_output_path : str
        Directory where the CSV files will be writen.
    """
    csv_columns = ['Index', 'score']
    nstr = utils.mp_nstr_precision_func
    for name, cases in scores.items():
        with open(csv_output_path / f'{name}_scores.csv', 'w') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(csv_columns)
            for idx, score in cases.items():
                writer.writerow([idx, nstr(score)])


def write_overall_scores_csv(scores, csv_output_path):
    """
    Writes a CSV file containing overall scores for each test set.
    An overall score is a sum of the test case scores.

    Parameters
    ----------
    scores : dict
        Dictionary of test set filenames (excluding file extensions) to a
        dictionary to test case indices to test case scores.

    csv_output_path : str
        Directory where the CSV files will be writen.
    """
    csv_columns = ['test_set', 'score']
    nstr = utils.mp_nstr_precision_func
    with open(csv_output_path / 'overall_scores.csv', 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(csv_columns)
        for name, cases in scores.items():
            test_set_score_sum = sum(s for s in cases.values())
            writer.writerow([name, nstr(test_set_score_sum)])


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Measure the distance between fitted and known parameter '
                    'sets.'
    )
    parser.add_argument(
        'directory', type=Path,
        help='Directory containing fitted parameter CSV files.')
    parser.add_argument(
        '--test-sets', dest='test_sets', type=str, default='',
        help='Comma-separated list of test sets to score.')
    parser.add_argument(
        '--csv-output-path', dest='csv_output_path', type=Path,
        default='.', help='Directory where to write output CSV files.')
    return parser


########
# MAIN #
########


if __name__ == '__main__':
    args = get_argparser().parse_args()
    test_sets_to_score = _get_test_sets_to_score(
        args.test_sets, args.directory)
    scores = {}
    num_compare_pts = 10
    num_total_pts = 200
    constants = utils.constants()
    vth, atol = constants['vth'], constants['atol']

    for name in sorted(utils.TEST_SETS_PRECISE.intersection(
            test_sets_to_score)):
        scores[name] = {}
        known_parameter_sets = utils.read_iv_curve_parameter_sets(
            utils.TEST_SETS_DIR / name)
        fitted_parameter_sets = utils.read_iv_curve_parameter_sets(
            args.directory / name)
        if known_parameter_sets.keys()==fitted_parameter_sets.keys():
            for idx, known_p in known_parameter_sets.items():
                fitted_p = fitted_parameter_sets[idx]
                scores[name][idx] = score_curve(known_p, fitted_p, vth,
                                                num_compare_pts, atol)
        else:
            raise ValueError(f'Fitted parameter Index in {name} must be the same as the Test Case')

    for name in sorted(utils.TEST_SETS_NOISY.intersection(
            test_sets_to_score)):
        scores[name] = {}
        known_parameter_sets = utils.read_iv_curve_parameter_sets(
            utils.TEST_SETS_DIR / name)
        fitted_parameter_sets = utils.read_iv_curve_parameter_sets(
            args.directory / name)
        if known_parameter_sets.keys()==fitted_parameter_sets.keys():
            for idx, known_p in known_parameter_sets.items():
                fitted_p = fitted_parameter_sets[idx]
                scores[name][idx] = score_parameters(known_p, fitted_p)
        else:
            raise ValueError(f'Fitted parameter Index in {name} must be the same as the Test Case')

    for name in utils.get_filenames_in_directory(
            utils.TEST_SETS_DIR).difference(test_sets_to_score):
        scores[name] = {0: float('NaN')}

    write_test_set_score_per_curve_csvs(scores, args.csv_output_path)
    write_overall_scores_csv(scores, args.csv_output_path)

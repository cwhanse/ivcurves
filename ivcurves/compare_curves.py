import argparse
import csv
import pvlib
import matplotlib.pyplot as plt

# from ivcurves repo
import utils
from utils import mp # same instance of mpmath's mp imported in ivcurves/utils
import precise


#####################
# Find intersection #
#####################

def find_x_intersection(single_diode, known_xs, known_ys, xp, yp, num_segments, atol, maxsteps=100):
    r"""
    Finds x-coordinate of the intersection between the known IV curve and the
    line segment from the origin to the given point :math:`(x_p, y_p)`.

    This is an auxiliary function for :func:`total_score`.

    Parameters
    ----------
    single_diode : function
        The single diode equation with two unknowns (the first being voltage,
        the second being current).

    known_xs : list of floats
        A list of x-coordinates.

    known_ys : list of floats
        A list of y-coordinates associated to ``known_xs`` that lie on the known
        curve.

    xp : float
        x-coordinate of point on fitted curve (voltage).

    yp : float
        y-coordinate of point on fitted curve (current).

    num_segments : int
        Number of segments to use when approximating known curve. (See
        :func:`get_guess_interval`.)

    atol : float
        Absolute tolerance of how close the found x-coordiante must be to the
        true x-coordinate of the intersection point.

    maxsteps : int
        Maximum number of iterations for ``mp.findroot`` to find the
        x-coordinate of the intersection point. A ``ValueError`` is thrown
        if ``atol`` is not met after ``maxsteps`` iterations.

    Returns
    -------
    float
        x-coordinate of the intersection of known IV curve and line through the
        origin and the given point.
    """
    if xp == 0:
        # then line through origin and (`xp`, `yp`) sits on y-axis
        # so x coordinate at intersection must be zero
        return 0

    else:
        # line through origin and (`xp`, `yp`) has a defined slope
        line = lambda x : (yp / xp) * x

        # solve for intersection of line and single_diode
        solve_for_zero = lambda x : single_diode(x, line(x)) - line(x)
        guess_int = get_guess_interval(known_xs, known_ys, (xp, yp), num_segments)
        try:
            # setting tol=atol**2 because findroot checks |func(zero)|**2 < tol
            x_int = mp.findroot(solve_for_zero, guess_int, tol=atol**2,
                                maxsteps=maxsteps)
        except ValueError as e:
            raise ValueError("Can't find an intersection point. "
                             'Perhaps the curves are too far from each other? '
                             f'{e}')

        assert abs(solve_for_zero(x_int)) < atol

    return x_int


def get_guess_interval(known_xs, known_ys, pt_on_line, num_segments):
    r"""
    Finds the interval in which the known curve intersects the given line.

    This is an auxiliary function for :func:`find_x_intersection`.

    Parameters
    ----------
    known_xs : list of floats
        A list of x-coordinates.

    known_ys : list of floats
        A list of y-coordinates associated to ``known_xs`` that lie on the known
        curve.

    pt_on_line : tuple of floats
        Point on fitted curve.

    num_segments : int
        Number of segments to use when approximating known curve.
        We will find ``num_segments`` points on the given curve, then consider
        consecutive pairs of points (which will define a line segment between
        them). If the line that passes through the origin and ``pt_on_line``
        crosses this particular segment, we return the x-coordinates of the
        endpoints of this segment.

    Returns
    -------
    tuple of mpmath floats
        The left and right x-coordinates of the interval that contains the
        intersection of the known curve with the line that passes through the
        origin and ``pt_on_line``.
    """
    # find slope and y-intercept of line, if finite
    if pt_on_line[0] != 0:
        line_slope, line_incpt = pt_on_line[1] / pt_on_line[0], 0

    pts = list(zip(known_xs, known_ys))
    # go through line segments
    for idx in range(len(pts)-1):
        if (pts[idx+1][0] == pts[idx][0]):
            # line segment is vertical, so x at intersection (`int_x`)
            # must be pts[idx][0] (==pts[idx+1][0])
            int_x = pts[idx][0]

            if pt_on_line[0] == 0: # line is on y-axis
                if int_x == 0: # segment is on y-axis
                    int_y = pts[idx][1] # segment is contained in line, just take an endpoint of segment as intersection pt
                else: # segment doesn't intersect y-axis (and so doesn't intersect line)
                    continue
            else: # slope of line is finite
                # so we can solve for y coordinate of intersection
                int_y = line_slope*int_x + line_incpt

        else: # segment is not vertical, and so has a finite slope
            # find slope and y-intercept of segment
            seg_slope = (pts[idx+1][1] - pts[idx][1]) / (pts[idx+1][0] - pts[idx][0])
            seg_incpt = pts[idx][1] - pts[idx][0]*seg_slope

            # intersection of line and segment
            if pt_on_line[0] != 0: # then line_slope and line_incpt are defined
                int_x = (seg_incpt - line_incpt) / (line_slope - seg_slope)
                int_y = line_slope*int_x + line_incpt

            else: # intersection is where segment crosses y-axis
                int_x = 0
                int_y = seg_incpt

        # check that found intersection point occurs within segment
        # if intersection point is within segment, return intersection point
        if min(pts[idx][0], pts[idx+1][0]) <= int_x and int_x <= max(pts[idx][0], pts[idx+1][0]):
            if min(pts[idx][1], pts[idx+1][1]) <= int_y and int_y <= max(pts[idx][1], pts[idx+1][1]):
                return pts[idx][0], pts[idx+1][0] # return x-coords of interval that contains intersection

    return pts[idx][0], pts[idx+1][0] # in case it misses the last interval because intersection occurs at (or near) last endpoint



######################
# Calculate distance #
######################

def find_distance(x, y, xp, yp):
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
    assert not (x == 0 and y == 0) # this should never happen for these curves

    if x == 0: # then xp == 0 too
        diff_x, diff_y = 0, (yp - y) / y
    elif y == 0: # then yp == 0 too
        diff_x, diff_y = (xp - x) / x, 0
    else: # x != 0 and y != 0
        diff_x, diff_y = (xp - x) / x, (yp - y) / y

    return mp.sqrt(diff_x**2 + diff_y**2)



###############
# Final score #
###############

def total_score(known_curve_params, fitted_curve_params, vth, num_pts, atol):
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
            Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

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
            Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

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
        :func:`precise.get_precise_i`.)
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
    known_xs, known_ys = get_curve(known_curve_params, vth, num_pts, atol)
    fit_xs, fit_ys = get_curve(fitted_curve_params, vth, num_pts, atol)

    il, io, rs, rsh, n, ns = known_curve_params
    single_diode = lambda v, i : il - io * mp.expm1((v + i*rs) / (n * ns * vth)) - ((v + i*rs) / rsh)

    score = 0

    for v, i in zip(fit_xs, fit_ys):
        # for each point (`v`, `i`) on the fitted curve, find the associated
        # point on known curve (`new_voltage`, `new_current`)
        new_voltage = find_x_intersection(single_diode, known_xs, known_ys, v, i, num_pts, atol)
        new_current = precise.lambert_i_from_v(new_voltage, il, io, rs, rsh, n, vth, ns) # find current associated to new_voltage

        # if voltage, current pair not a precise enough solution to single diode equation, make more precise
        dff = precise.diff_lhs_rhs(new_voltage, new_current, il, io, rs, rsh, n, vth, ns)
        if abs(dff) > atol:
            new_current = mp.findroot(lambda y : precise.diff_lhs_rhs(new_voltage, y, il, io, rs, rsh, n, vth, ns), new_current, tol=atol**2)

        # calculate distance between these points, and add to score
        score += find_distance(new_voltage, new_current, v, i)

    return score


########
# PLOT #
########

def get_curve(curve_parameters, vth, num_pts, atol):
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
            Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

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
        (See :func:`precise.get_precise_i`.) Each solution pair is a point
        on the curve.

    Returns
    -------
    (vv, ii) : tuple of numpy arrays
        ``vv`` is a numpy array of float64 and ``ii`` is a numpy array of
        mpmath floats. Each array has ``num_pts`` entries.
    """
    il, io, rs, rsh, n, ns = curve_parameters
    vv, ii = precise.get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
    return vv, ii


def iv_plotter(iv_known, iv_fitted, vth, num_pts, atol, pts=None, plot_lines=True):
    r"""
    Plots the fitted curve (green) and the known curve (cyan).

    Parameters
    ----------
    iv_known : list
        A list of parameters representing a given IV curve. Should be passed in
        the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

        rsh : numeric
            Shunt resistance :math:`R_{sh}` under desired IV curve conditions.
            [ohm]

        n : numeric
            Diode ideality factor :math:`n`

        ns : numeric
            Number of cells in series :math:`N_s`

    iv_fitted : list
        A list of parameters representing a given IV curve. Should be passed in
        the order [il, io, rs, rsh, n, ns].

        il : numeric
            Light-generated current :math:`I_L` (photocurrent) [A]

        io : numeric
            Diode saturation :math:`I_0` current under desired IV curve
            conditions. [A]

        rs : numeric
            Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

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
        Number of points to use when plotting curves.

    atol : float
        The error of each of the solution pairs found is at most ``atol``. (See
        :func:`precise.get_precise_i`.)
        Each solution pair is a point on the curve.

    pts : list, default []
        A list of points on the fitted curve that will be plotted, along with
        their associated points on the known curve.

    plot_lines : bool, default True
        If true, the lines connecting the points on the fitted curve and the
        associated points on the known curve will be plotted.

    Returns
    -------
    plot
        A matplotlib plot of the known curve (cyan) and the fitted curve
        (green). If ``pts`` is given, the image will also include points that
        were compared on the fitted curve (light green) and the known curve
        (magenta). If ``plot_lines`` is True, then lines are drawn that connect
        the associated points.
    """
    if not pts:
        pts = []

    plot = plt.plot()

    # plot known curve (known parameters)
    known_xs, known_ys = get_curve(iv_known, vth, num_pts, atol)
    plt.plot(known_xs, known_ys, 'cyan')

    # plot fitted curve (fitted parameters)
    fit_xs, fit_ys = get_curve(iv_fitted, vth, num_pts, atol)
    plt.plot(fit_xs, fit_ys, color='green')

    il, io, rs, rsh, n, ns = iv_known
    single_diode = lambda v, i : il - io * mp.expm1((v + i*rs) / (n * ns * vth)) - ((v + i*rs) / rsh)

    num_segments = len(pts)
    count = 0
    for vp, ip in pts:
        # plot point on fitted curve
        plt.plot(vp, ip, marker='o', color='lightgreen', markersize=3)

        # get intersection point on known curve
        try:
            new_voltage = find_x_intersection(single_diode, known_xs, known_ys, vp, ip, num_pts, atol)
            new_current = precise.lambert_i_from_v(new_voltage, il, io, rs, rsh, n, vth, ns)
        except:
            print("BAD PT @", count)
            count += 1
            continue
        else:
            count += 1

        # plot point on known curve
        plt.plot(new_voltage, new_current, marker='o', color='magenta', markersize=3)

        if plot_lines: # plot line intersecting curves
            xs = [0, min(vp, new_voltage), max(vp, new_voltage)]
            if xs[1] == vp:
                ys = [0, ip, new_current]
            else:
                ys = [0, new_current, ip]

            plt.plot(xs, ys, color='darkorchid', linewidth=0.4)

    return plot


def get_test_sets_to_score(fitted_files_directory, test_set=''):
    """
    Returns a list of valid test set filenames (excluding file extensions)
    based on the files found in ``fitted_files_directory``.

    Parameters
    ----------
    fitted_files_directory : str
        Directory that contains files whose filenames are test set filenames.

    test_set : str, default ''
        A singular test set filename to look for in ``fitted_files_directory``.
        The file extension must be excluded. This argument is ignored by
        default.

    Returns
    -------
    list
        A list of valid test set filenames (excluding file extensions).
    """
    test_set_names = utils.get_filenames_in_directory(utils.TEST_SETS_DIR)
    test_sets_to_score = []
    if test_set:
        if test_set not in test_set_names:
            raise ValueError(f'\'{test_set}\' is not a test set')
        test_sets_to_score = [test_set]
    else:
        filenames = utils.get_filenames_in_directory(fitted_files_directory)
        test_sets_to_score = [f for f in filenames if f in test_set_names]
        test_sets_to_score.sort()
        if not test_sets_to_score:
            raise ValueError(f'no test sets found in \'{fitted_files_directory}\'')
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
        with open(f'{csv_output_path}/{name}_scores.csv', 'w') as file:
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
    with open(f'{csv_output_path}/overall_scores.csv', 'w') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerow(csv_columns)
        for name, cases in scores.items():
            test_set_score_sum = sum(s for s in cases.values())
            writer.writerow([name, nstr(test_set_score_sum)])


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Measure the distance between IV curves generated from '
                    'the parameters of the single diode equation.'
    )
    parser.add_argument('fitted_files_directory', type=str,
                        help='Directory containing fitted parameter CSV files.')
    parser.add_argument('--test-set', dest='test_set', type=str, default='',
                        help='Name of test set to score.')
    parser.add_argument('--csv-output-path', dest='csv_output_path', type=str,
                        default='.', help='Directory where to write output CSV files.')
    parser.add_argument('--plot', action=argparse.BooleanOptionalAction,
                        help='Plot each IV curve fit.')
    return parser


########
# MAIN #
########


if __name__ == '__main__':
    args = get_argparser().parse_args()

    # # intersecting curves example
    # iv_known = list(map(mp.mpmathify, [6.0, 3.8500023e-06, 1.6816000000000002, 8832.800000000005, 1.4200000000000004, 72]))
    # iv_fitted = list(map(mp.mpmathify, [4.2, 6.500087e-07, 0.9453, 17881.40000000001, 1.6300000000000006, 72]))

    test_sets_to_score = get_test_sets_to_score(args.fitted_files_directory, args.test_set)
    scores = {}
    num_compare_pts = 10
    num_total_pts = 200
    constants = utils.constants()
    vth, atol = constants['vth'], constants['atol']

    for name in test_sets_to_score:
        scores[name] = {}
        known_parameter_sets = utils.read_iv_curve_parameter_sets(f'{utils.TEST_SETS_DIR}/{name}')
        fitted_parameter_sets = utils.read_iv_curve_parameter_sets(f'{args.fitted_files_directory}/{name}')
        for idx, known_p in known_parameter_sets.items():
            fitted_p = fitted_parameter_sets[idx]
            scores[name][idx] = total_score(known_p, fitted_p, vth, num_compare_pts, atol)

            if args.plot:
                fit_xs, fit_ys = get_curve(fitted_p, vth, num_compare_pts, atol)
                plot = iv_plotter(known_p, fitted_p, vth, num_total_pts, atol,
                                  pts=list(zip(fit_xs, fit_ys)), plot_lines=True)
                plt.show()
                plt.cla()

    write_test_set_score_per_curve_csvs(scores, args.csv_output_path)
    write_overall_scores_csv(scores, args.csv_output_path)

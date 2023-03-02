import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from ivcurves import compare_curves, precise, utils
from ivcurves.utils import mp  # same instance of mpmath's mp imported in ivcurves/utils


def plot_precise_iv_curves(test_set_filename, case_parameter_sets, vth, atol, num_pts):
    """
    Plot precise IV curves.

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
    """
    plt.style.use('seaborn-v0_8-darkgrid')
    for idx, (il, io, rs, rsh, n, ns) in case_parameter_sets.items():
        fig, ax = plt.subplots()
        params = il, io, rs, rsh, n, vth, ns
        ax.set_xlabel('Voltage')
        ax.set_ylabel('Current')
        v_vals, i_vals = precise.get_precise_i(*params, atol=atol, num_pts=num_pts)
        ax.plot(v_vals, i_vals)
        fig.savefig(utils.make_iv_curve_name(test_set_filename, idx),
                    bbox_inches='tight')
        plt.close(fig)


def scoring_visualization(test_set_filename, case_parameter_set_idx, iv_known,
                          iv_fitted, vth, num_pts, atol, pts=None,
                          plot_lines=True):
    r"""
    Plots the fitted curve (green) and the known curve (cyan) to visualize
    how compare_curves scores the fitted parameter sets.

    Parameters
    ----------
    test_set_filename : str
        The name of the test set CSV and JSON file. The filename must exclude
        its file extension.

    case_parameter_set_idx : int
        The index of the case parameter set loaded from a CSV file named
        ``test_set_filename``.

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
        :func:`ivcurves.precise.get_precise_i`.)
        Each solution pair is a point on the curve.

    pts : list, default []
        A list of points on the fitted curve that will be plotted, along with
        their associated points on the known curve.

    plot_lines : bool, default True
        If true, the lines connecting the points on the fitted curve and the
        associated points on the known curve will be plotted.
    """
    if not pts:
        pts = []

    fig, ax = plt.subplots()

    # plot known curve (known parameters)
    known_xs, known_ys = compare_curves.get_curve(iv_known, vth, num_pts, atol)
    ax.plot(known_xs, known_ys, 'cyan')

    # plot fitted curve (fitted parameters)
    fit_xs, fit_ys = compare_curves.get_curve(iv_fitted, vth, num_pts, atol)
    ax.plot(fit_xs, fit_ys, color='green')

    il, io, rs, rsh, n, ns = iv_known
    single_diode = lambda v, i: il - io * mp.expm1((v + i*rs) / (n * ns * vth)) - ((v + i*rs) / rsh)

    count = 0
    for vp, ip in pts:
        # plot point on fitted curve
        ax.plot(vp, ip, marker='o', color='lightgreen', markersize=3)

        # get intersection point on known curve
        try:
            interval_current = min(known_xs), max(known_xs)
            new_voltage = compare_curves.find_x_intersection(interval_current, single_diode, (vp, ip), atol)
            new_current = precise.lambert_i_from_v(new_voltage, il, io, rs, rsh, n, vth, ns)
        except ValueError:
            print("BAD PT @", count)
            count += 1
            continue
        else:
            count += 1

        # plot point on known curve
        ax.plot(new_voltage, new_current, marker='o', color='magenta', markersize=3)

        if plot_lines:  # plot line intersecting curves
            xs = [0, min(vp, new_voltage), max(vp, new_voltage)]
            if xs[1] == vp:
                ys = [0, ip, new_current]
            else:
                ys = [0, new_current, ip]

            ax.plot(xs, ys, color='darkorchid', linewidth=0.4)

        fig.savefig(f'{utils.make_iv_curve_name(test_set_filename, case_parameter_set_idx)}_compare_curves',
                    bbox_inches='tight')
        plt.close(fig)


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Plots precise IV curves, and visualizations for compare_curves scoring against them.'
    )
    parser.add_argument('save_images_path', type=Path,
                        help='Where to save the test set IV curve plots.')
    parser.add_argument('--fitted-files-path', dest='fitted_files_path', type=Path,
                        help='Where to find fitted parameter CSV files for compare_curves scoring visualization.')
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()

    constants = utils.constants()
    vth, temp_cell, atol, num_pts = (constants['vth'], constants['temp_cell'],
                                     constants['atol'], constants['num_pts'])

    # plot precise ivcurves
    for name in utils.get_filenames_in_directory(utils.TEST_SETS_DIR):
        case_parameter_sets = utils.read_iv_curve_parameter_sets(utils.TEST_SETS_DIR / name)
        plot_precise_iv_curves(
            args.save_images_path / name, case_parameter_sets, vth, atol,
            num_pts
        )

    # plot compare_curves scoring visualization
    if args.fitted_files_path:
        test_sets_to_score = compare_curves.get_test_sets_to_score(args.fitted_files_path)
        num_compare_pts = 10
        num_total_pts = 200
        for name in test_sets_to_score:
            known_parameter_sets = utils.read_iv_curve_parameter_sets(utils.TEST_SETS_DIR / name)
            fitted_parameter_sets = utils.read_iv_curve_parameter_sets(args.fitted_files_path / name)
            for idx, known_p in known_parameter_sets.items():
                fitted_p = fitted_parameter_sets[idx]
                fit_xs, fit_ys = compare_curves.get_curve(fitted_p, vth, num_compare_pts, atol)
                scoring_visualization(
                    args.save_images_path / name, idx,
                    known_p, fitted_p, vth, num_total_pts, atol,
                    pts=list(zip(fit_xs, fit_ys)), plot_lines=True
                )

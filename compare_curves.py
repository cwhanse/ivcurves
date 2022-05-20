import pvlib
from mpmath import mp
import matplotlib.pyplot as plt
from shapely.geometry import LineString

# from ivcurves repo
from max_power import lambert_i_from_v
from precise import diff_lhs_rhs, get_precise_i


#####################
# Find intersection #
#####################

def find_intersection(known_curve_params, vp, ip, vth, atol):
    # find intersection point on known curve
    il, io, rs, rsh, n, ns = known_curve_params
    single_diode = lambda v, i : il - io * mp.expm1((v + i*rs) / (n * ns * vth)) - ((v + i*rs) / rsh)  

    if vp == 0:
        new_volt = 0
        solve_for_zero = lambda i : single_diode(0, i) - i
        new_current = try_findroot(solve_for_zero, ip, atol)
        assert (single_diode(0, new_current) - new_current) < atol

    else:
        line = lambda v : (ip / vp) * v 
        solve_for_zero = lambda v : single_diode(v, line(v)) - line(v)

        guess_int = get_guess_interval(known_curve_params, (vp, ip), vth)
        new_volt = try_findroot(solve_for_zero, guess_int, atol)
        assert abs(solve_for_zero(new_volt)) < atol

        new_current = lambert_i_from_v(new_volt, il, io, rs, rsh, n, vth, ns)

    # if new_volt, new_current aren't a precise solution, make precise
    dff = diff_lhs_rhs(new_volt, new_current, il, io, rs, rsh, n, vth, ns)
    if abs(dff) > atol:
        new_current = mp.findroot(lambda i : diff_lhs_rhs(new_volt, i, il, io, rs, rsh, n, vth, ns), new_current, tol=atol**2)

    return new_volt, new_current


def try_findroot(func, guess, atol, max_steps=100):
    raisedException = True
    steps = 30 # initial default value for findroot secant method
    while raisedException:
        try:
            zero = mp.findroot(func, guess, tol=atol**2, maxsteps=steps)
        except:
            steps += 10
            if steps > max_steps:
                raise Exception("Can't find an intersection point. Perhaps the curves are too far from each other?")
        else:
            raisedException = False

    return zero


def get_guess_interval(known_curve_params, pt_on_line, vth, num_pts=10, atol=1e-16):
    # return interval in which curve intersects line (good guess for location of zero for findroot)
    line = LineString([(0,0), pt_on_line])
    vv, ii = get_curve(known_curve_params, vth, num_pts, atol)

    pts = list(zip(vv, ii))
    for idx in range(len(pts)-1): 
        segment = LineString([pts[idx], pts[idx+1]])
        if segment.intersects(line):
            return pts[idx][0], pts[idx+1][0] # return x-coords of interval that contains intersection

    return pts[idx][0], pts[idx+1][0] # in case it misses the last interval because intersection occurs at (or near) last endpoint



######################
# Calculate distance #
######################

def find_distance(v, i, vp, ip):
    # v, i : point on known curve & vp, ip : point on fitted curve
    # FIXME how to alter when v == 0 or i == 0
    diff_v, diff_i = (vp - v) / v, (ip - i) / i 
    return mp.sqrt(diff_v**2 + diff_i**2)



###############
# Final score #
###############

def total_score(known_curve_params, fitted_curve_params, vth, num_pts, atol=1e-16):
    fit_xs, fit_ys = get_curve(fitted_curve_params, vth, num_pts, atol)
    score = 0

    for x, y in list(zip(fit_xs, fit_ys)):
        new_volt, new_current = find_intersection(known_curve_params, x, y, vth, atol)
        if new_volt == 0 or new_current == 0: continue # FIXME
        score += find_distance(new_volt, new_current, x, y)

    return score 



########
# PLOT #
########

def get_curve(curve_parameters, vth, num_pts, atol):
    il, io, rs, rsh, n, ns = curve_parameters
    vv, ii = get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
    return vv, ii


def iv_plotter(iv_known, iv_fitted, vth, num_pts, pts=[], plot_lines=False, atol=1e-16):
    # iv_known, iv_fitted = [il, io, rs, rsh, n, ns]
    plot = plt.plot()

    # plot known curve (known parameters)
    known_xs, known_ys = get_curve(iv_known, vth, num_pts, atol)
    plt.plot(known_xs, known_ys, 'cyan')

    # plot fitted curve (fitted parameters)
    fit_xs, fit_ys = get_curve(iv_fitted, vth, num_pts, atol)
    plt.plot(fit_xs, fit_ys, color='green')

    for vp, ip in pts:
        # plot point on fitted curve
        plt.plot(vp, ip, marker='o', color='lightgreen', markersize=3)

        # get intersection point on known curve
        new_volt, new_current = find_intersection(iv_known, vp, ip, vth, atol)

        # plot point on known curve
        plt.plot(new_volt, new_current, marker='o', color='magenta', markersize=3)
        
        if plot_lines:
            # make sure line goes all the way to end of plot
            if known_xs[-1] > fit_xs[-1]: xs = known_xs
            else: xs = fit_xs

            # plot line intersecting curves
            if vp != 0:
                line_ys = [(ip / vp)*x for x in xs]
                plt.plot(xs, line_ys, color='darkorchid', linewidth=0.4)
            else:
                plt.plot([0]*len(xs), xs, color='darkorchid', linewidth=0.4)

    return plot



######## 
# MAIN #
######## 

if __name__ == "__main__":
    mp.dps = 40 # 16*2 rounded up (1e-16 is default tolerance)

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = [1.380649e-23, 1.60217663e-19, 298.15]
    vth = (k * temp_cell) / q
    num_pts = 10
    atol = 1e-16

    # intersecting curves example
    iv_known = [1.0, 5e-10, 0.1, 300, 1.01, 72]
    iv_fitted = [8.0, 3e-8, 0.1, 300, 1.01, 72] 

    print("Total score:", total_score(iv_known, iv_fitted, vth, num_pts, atol))

    fit_xs, fit_ys = get_curve(iv_fitted, vth, num_pts, atol)
    plot = iv_plotter(iv_known, iv_fitted, vth, 200, list(zip(fit_xs, fit_ys)), plot_lines=False)
    plt.show()


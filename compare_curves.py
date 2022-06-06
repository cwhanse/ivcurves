import pvlib
from mpmath import mp
import matplotlib.pyplot as plt

# from ivcurves repo
from max_power import lambert_i_from_v
from precise import diff_lhs_rhs, get_precise_i


#####################
# Find intersection #
#####################

def find_x_intersection(single_diode, known_xs, known_ys, xp, yp, num_segments, atol):
    if xp == 0: 
        return 0

    else:
        line = lambda x : (yp / xp) * x
        solve_for_zero = lambda x : single_diode(x, line(x)) - line(x)

        guess_int = get_guess_interval(known_xs, known_ys, (xp, yp), num_segments)
        x_int = try_findroot(solve_for_zero, guess_int, atol)
        assert abs(solve_for_zero(x_int)) < atol

    return x_int


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


def get_guess_interval(known_xs, known_ys, pt_on_line, num_segments):
    # return interval in which curve intersects line (good guess for location of zero for findroot)
    
    # find slope and y-intercept of line, if finite
    if pt_on_line[0] != 0:
        line_slope, line_incpt = pt_on_line[1] / pt_on_line[0], 0

    pts = list(zip(known_xs, known_ys))
    for idx in range(len(pts)-1): 
        if (pts[idx+1][0] == pts[idx][0]): # line segment is vertical
            int_x = pts[idx][0]

            if pt_on_line[0] == 0: # line is on y-axis
                if int_x == 0: # segment is on y-axis
                    int_y = pts[idx][1] # segment is contained in line, just take an endpoint of segment as intersection pt
                else: # segment doesn't intersect y-axis (and so doesn't intersect line)
                    continue
            else: # slope of line is finite
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

        # if intersection point is within segment, return intersection point
        if min(pts[idx][0], pts[idx+1][0]) <= int_x and int_x <= max(pts[idx][0], pts[idx+1][0]):
            if min(pts[idx][1], pts[idx+1][1]) <= int_y and int_y <= max(pts[idx][1], pts[idx+1][1]):
                return pts[idx][0], pts[idx+1][0] # return x-coords of interval that contains intersection

    return pts[idx][0], pts[idx+1][0] # in case it misses the last interval because intersection occurs at (or near) last endpoint



######################
# Calculate distance #
######################

def find_distance(x, y, xp, yp):
    # x, y is a point on known curve 
    # xp, yp is a point on fitted curve
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
    known_xs, known_ys = get_curve(known_curve_params, vth, num_pts, atol)
    fit_xs, fit_ys = get_curve(fitted_curve_params, vth, num_pts, atol)

    il, io, rs, rsh, n, ns = known_curve_params
    single_diode = lambda v, i : il - io * mp.expm1((v + i*rs) / (n * ns * vth)) - ((v + i*rs) / rsh)  

    score = 0

    for v, i in list(zip(fit_xs, fit_ys)):
        new_voltage = find_x_intersection(single_diode, known_xs, known_ys, v, i, num_pts, atol)
        new_current = lambert_i_from_v(new_voltage, il, io, rs, rsh, n, vth, ns) # find current associated to new_voltage

        # if voltage, current pair not a precise enough solution to single diode equation, make more precise
        dff = diff_lhs_rhs(new_voltage, new_current, il, io, rs, rsh, n, vth, ns)
        if abs(dff) > atol:
            new_current = mp.findroot(lambda y : diff_lhs_rhs(new_voltage, y, il, io, rs, rsh, n, vth, ns), new_current, tol=atol**2)

        score += find_distance(new_voltage, new_current, v, i)

    return score 



########
# PLOT #
########

def get_curve(curve_parameters, vth, num_pts, atol):
    il, io, rs, rsh, n, ns = curve_parameters
    vv, ii = get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
    return vv, ii


def iv_plotter(iv_known, iv_fitted, vth, num_pts, atol, pts=[], plot_lines=True):
    # iv_known, iv_fitted = [il, io, rs, rsh, n, ns]
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
            new_current = lambert_i_from_v(new_voltage, il, io, rs, rsh, n, vth, ns) 
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



######## 
# MAIN #
######## 

if __name__ == "__main__":
    mp.dps = 40 # 16*2 rounded up (1e-16 is default tolerance)

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = [1.380649e-23, 1.60217663e-19, 298.15]
    vth = (k * temp_cell) / q
    num_compare_pts = 10
    num_total_pts = 200
    atol = 1e-16

    # intersecting curves example
    iv_known = [6.0, 3.8500023e-06, 1.6816000000000002, 8832.800000000005, 1.4200000000000004, 72]
    iv_fitted = [4.2, 6.500087e-07, 0.9453, 17881.40000000001, 1.6300000000000006, 72]

    print("Total score:", total_score(iv_known, iv_fitted, vth, num_compare_pts, atol))

    fit_xs, fit_ys = get_curve(iv_fitted, vth, num_compare_pts, atol)
    plot = iv_plotter(iv_known, iv_fitted, vth, num_total_pts, atol, pts=list(zip(fit_xs, fit_ys)), plot_lines=True)
    plt.show()


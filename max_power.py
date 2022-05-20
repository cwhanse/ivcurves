from mpmath import mp
from itertools import product
import pvlib 


r"""
Parameters
--------------
il -> photocurrent 
io -> saturation_current
rs -> series resistance
rsh -> shunt resistance
n -> diode ideality factor 
ns -> number of cells in series
"""


###################
# Max power point # 
###################


def max_power_pt_finder(il, io, rs, rsh, n, vth, ns, atol):
    # calculate max power with at most atol error for IV curve with given parameters
    # returns max_voltage, max_current, max_power

    # make parameters mpf types for precision
    il, io, rs, rsh, n, vth = [mp.mpf(str(val)) for val in [il, io, rs, rsh, n, vth]]

    # function we want to maximize
    power_func = lambda x : x * lambert_i_from_v(x, il, io, rs, rsh, n, vth, ns) # x represents voltage

    # find initial interval endpts (using (0, I_sc) and (V_oc, 0) as starting endpts, power=0 at both these points)
    # x-coord is voltage, y-coord is power (power = voltage * current)
    xl, yl = lambert_i_from_v(0, il, io, rs, rsh, n, vth, ns), 0
    xr, yr = lambert_v_from_i(0, il, io, rs, rsh, n, vth, ns), 0

    # run golden search on power function with above starting interval
    max_voltage, max_power = golden_search((xl, yl), (xr, yr), power_func, atol)

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


def golden_search(l_endpt, r_endpt, func, atol, int_pt=None, is_right_int_pt=None):
    # func is function we want to maximize (single-variable)

    # add some sort of iterlimit FIXME
    # overflow ? FIXME
    # take care of f(int_pt_1) == f(int_pt_2) (right now just pushing this case into the f(int_pt_1) < f(int_pt_2) case) FIXME

    xl, yl = l_endpt
    xr, yr = r_endpt

    # find left and right interior points
    if int_pt == None: # first iteration, no interior points yet
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
            return golden_search(l_endpt, r_int_pt, func, atol, l_int_pt, True) # true because l_int_pt is now the right_int_pt of new interval
    else:
        error = abs(r_endpt[0] - l_int_pt[0])
        if error < atol:
            return r_int_pt
        else:
            return golden_search(l_int_pt, r_endpt, func, atol, r_int_pt, False) # false because r_int_pt is now the left_int_pt of new interval


def get_left_int_pt(left_x_endpt, right_x_endpt, func):
    xl, xr = left_x_endpt, right_x_endpt
    rho = (1/2) * (3 - mp.sqrt(5)) # from http://www.math.kent.edu/~reichel/courses/intr.num.comp.2/lecture16/lecture8.pdf
    # this value for rho is equivalent to using golden ratio

    l_int_x = xl + rho*(xr - xl) # voltage
    l_int_y = func(l_int_x)
    return (l_int_x, l_int_y)


def get_right_int_pt(left_x_endpt, right_x_endpt, func):
    xl, xr = left_x_endpt, right_x_endpt
    rho = (1/2) * (3 - mp.sqrt(5)) # from http://www.math.kent.edu/~reichel/courses/intr.num.comp.2/lecture16/lecture8.pdf
    # this value for rho is equivalent to using golden ratio

    r_int_x = xl + (1 - rho)*(xr - xl) # voltage
    r_int_y = func(r_int_x)
    return (r_int_x, r_int_y)



#######################
# Auxiliary functions #
#######################


def diff_lhs_rhs(v, i, il, io, rs, rsh, n, vth, ns):
    # returns the difference between the lhs and rhs of the single diode equation
    return (il - io*mp.expm1((v + i*rs)/(n*ns*vth)) - (v + i*rs)/(rsh) - i)


def lambert_i_from_v(v, il, io, rs, rsh, n, vth, ns):
    gsh = 1. / rsh
    if rs == 0:
        return il - io*mp.expm1(v / (n*vth*ns)) - gsh*v
    else:
        # handle overflow FIXME
        argW = rs * io / ((n*vth*ns) * (rs*gsh + 1)) * mp.exp((rs*(il + io) + v)/((n*vth*ns)*(rs*gsh + 1)))
        lambertw_term = mp.lambertw(argW).real
        return (il + io - v*gsh) / (rs*gsh + 1) - ((n*vth*ns) / rs)*lambertw_term


def lambert_v_from_i(i, il, io, rs, rsh, n, vth, ns):
    gsh = 1. / rsh
    if gsh == 0:
        return (n*vth*ns) * mp.log1p((il - i) / io) - i*rs
    else:
        # handle overflow FIXME
        argW = io / (gsh * (n*vth*ns)) * mp.exp((-i + il + io) / (gsh * (n*vth*ns)))
        lambertw_term = mp.lambertw(argW).real
        return (il + io - i) / gsh - i*rs - (n*vth*ns)*lambertw_term



########
# Main #
########


if __name__ == "__main__":
    mp.dps = 40 # 16*2 rounded up
    atol = 1e-16
    parameters = dict()
    max_vals = dict()

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = [1.380649e-23, 1.60217663e-19, 298.15]
    vth = (k * temp_cell) / q


    case1 = True

    if case1:
        IL = [1.0, 8.0]
        IO = [5e-10, 3e-8]
        RS = [0.1, 1.0]
        RSH = [300, 3000]
        N = [1.01, 1.3]
        ns = 72
    else: # we're in case 2
        IL = [0.5, 2.5] 
        IO = [1e-9, 1e-8]
        RS = [0.1, 1.0]
        RSH = [300, 3000]
        N = [1.3, 1.5]
        ns = 140

    count = 0
    for il, io, rs, rsh, n in product(IL, IO, RS, RSH, N):
        max_voltage, max_current, max_power = max_power_pt_finder(il, io, rs, rsh, n, vth, ns, atol)
        parameters[count] = [il, io, rs, rsh, n]
        max_vals[count] = [max_voltage, max_current, max_power]
        count += 1

    assert count == 32

    for idx in range(count):
        print(parameters[idx])
        print('Max voltage:', max_vals[idx][0])
        print('Max current:', max_vals[idx][1])
        print('Max power:', max_vals[idx][2])
        print('')


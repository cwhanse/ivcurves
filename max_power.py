from mpmath import mp
from itertools import product
import pvlib 

# from ivcurves repo
from precise import diff_lhs_rhs


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
        Diode saturation :math:`I_0` current under desired IV curve conditions. [A]

    rs : numeric
        Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

    rsh : numeric
        Shunt resistance :math:`R_{sh}` under desired IV curve conditions. [ohm]

    n : numeric
        Diode ideality factor :math:`n`

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K), :math:`T_c` is the temperature of the p-n junction in Kelvin, and :math:`q` is the charge of an electron (coulombs). 

    ns : numeric
        Number of cells in series :math:`N_s`

    atol : float
        The returned voltage will be at most `atol` from the voltage that yields the true maximum power of the given IV curve.

    Returns
    -------
    (max_voltage, max_current, max_power) : tuple of mpmath floats
        A good approximation for the voltage :math:`V`, current :math:`I`, and power :math:`P` of the maximum power point of the given IV curve.
    """
    # make parameters mpf types for precision
    il, io, rs, rsh, n, vth = [mp.mpf(str(val)) for val in [il, io, rs, rsh, n, vth]]

    # function we want to maximize
    power_func = lambda x : x * lambert_i_from_v(x, il, io, rs, rsh, n, vth, ns) # x represents voltage

    # find initial interval endpts (using (0, I_sc) and (V_oc, 0) as starting endpts, power=0 at both these points)
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
    Uses golden-section search to recursively find maximum of the given function over the given interval, with at most `atol` error.

    Parameters
    ----------
    l_endpt : tuple of floats
        Left endpoint of interval in which a maximum occurs.

    r_endpt : tuple of floats
        Right endpoint of interval in which a maximum occurs.

    func : function
        Function we want to maximize (single-variable).

    atol : float 
        The x-coordinate of the returned point will be at most `atol` from the x-coordinate that produces the true maximum in the given interval.

    iterlimit : int
        Maximum number of iterations for golden-section search. Should converge before we hit this.

    int_pt : tuple of floats, optional
        Coordinates of interior point in interval. The default value is an empty tuple.

    is_right_int_pt : bool, optional
        If `int_pt` is the right hand interior point, then True. If `int_pt` is given, this should also be passed in (default value is False).

    num_iter : int
        The number of iterations we've already done.

    Returns
    -------
    tuple of ints
        Coordinates of a point whose y-coordinate is a good approximation for the true maximum, and whose x-coordinate is within `atol` of the x-coordinate that yields the true maximum of `func` on the given interval.

    Notes
    -----
    This is a recursive function. When using, `int_pt` and `is_right_int_pt` and `num_iter` should not be passed; they are only passed when the function recurses.

    For more information on the algorithm (and calculating the interior points in get_left_int_pt and get_right_int_pt), see http://www.math.kent.edu/~reichel/courses/intr.num.comp.2/lecture16/lecture8.pdf.
    """
    # overflow ? FIXME
    # take care of f(int_pt_1) == f(int_pt_2) (right now just pushing this case into the f(int_pt_1) < f(int_pt_2) case) FIXME

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
            return golden_search(l_endpt, r_int_pt, func, atol, iterlimit, l_int_pt, True, num_iter+1) # true because l_int_pt is now the right_int_pt of new interval
    else:
        error = abs(r_endpt[0] - l_int_pt[0])
        if error < atol:
            return r_int_pt
        else:
            return golden_search(l_int_pt, r_endpt, func, atol, iterlimit, r_int_pt, False, num_iter+1) # false because r_int_pt is now the left_int_pt of new interval


def get_left_int_pt(left_x_endpt, right_x_endpt, func):
    r"""
    Calculates the left interior point of the interval.

    Auxiliary function for golden_search.

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

    Auxiliary function for golden_search.

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
    Given a voltage, calculates the associated current using the Lambert W function.

    Parameters
    ----------
    v : numeric
        Voltage [V]

    il : numeric
        Light-generated current :math:`I_L` (photocurrent) [A]

    io : numeric
        Diode saturation :math:`I_0` current under desired IV curve conditions. [A]

    rs : numeric
        Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

    rsh : numeric
        Shunt resistance :math:`R_{sh}` under desired IV curve conditions. [ohm]

    n : numeric
        Diode ideality factor :math:`n`

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K), :math:`T_c` is the temperature of the p-n junction in Kelvin, and :math:`q` is the charge of an electron (coulombs). 

    ns : numeric
        Number of cells in series :math:`N_s`

    Returns
    -------
    mpmath float
        Current associated to given voltage via the single diode equation.

    Notes
    -----
    See pvlib.singlediode for original function. This implemention differs only in that it does not accept Series as inputs and can take mpmath floats as inputs.
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
    Given a current, calculates the associated voltage using the Lambert W function.

    Parameters
    ----------
    i : numeric
        Current [A]

    il : numeric
        Light-generated current :math:`I_L` (photocurrent) [A]

    io : numeric
        Diode saturation :math:`I_0` current under desired IV curve conditions. [A]

    rs : numeric
        Series resistance :math:`R_s` under desired IV curve conditions. [ohm]

    rsh : numeric
        Shunt resistance :math:`R_{sh}` under desired IV curve conditions. [ohm]

    n : numeric
        Diode ideality factor :math:`n`

    vth : numeric
        Thermal voltage of the cell :math:`V_{th}` [V]
        The thermal voltage of the cell (in volts) may be calculated as :math:`k_B T_c / q`, where :math:`k_B` is Boltzmann's constant (J/K), :math:`T_c` is the temperature of the p-n junction in Kelvin, and :math:`q` is the charge of an electron (coulombs). 

    ns : numeric
        Number of cells in series :math:`N_s`

    Returns
    -------
    mpmath float
        Voltage associated to given current via the single diode equation.

    Notes
    -----
    See pvlib.singlediode for original function. This implemention differs only in that it does not accept Series as inputs and can take mpmath floats as inputs.
    """
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


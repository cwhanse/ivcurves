import pvlib
import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp
from itertools import product


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


def diff_lhs_rhs(v, i, il, io, rs, rsh, n, vth, ns):
    # returns the difference between the lhs and rhs of the single diode equation
    return (il - io*mp.expm1((v + i*rs)/(n*ns*vth)) - (v + i*rs)/(rsh) - i)


def get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts):
    res = pvlib.pvsystem.singlediode(il, io, rs, rsh, n*vth*ns, num_pts)
    vv = res['v']
    ii = res['i']
    prec_i = np.array([0 for _ in range(ii.shape[0])], dtype=mp.mpf) # initialize 'empty' array for new, more precise i's

    assert len(vv) == len(ii)
    for idx in range(len(vv)):
        i = mp.mpf(str(ii[idx]))
        # check if i val already precise enough
        if abs(diff_lhs_rhs(vv[idx], i, il, io, rs, rsh, n, vth, ns)) < atol: 
            new_i = i

        else:
            # findroot takes the function whose roots we want to find, a good guess for where the root is, and a tolerance
            # default solver uses secant method
            new_i = mp.findroot(lambda x: diff_lhs_rhs(vv[idx], x, il, io, rs, rsh, n, vth, ns), i, tol=atol)

        # check that mp.findroot did what we wanted 
        assert abs(diff_lhs_rhs(vv[idx], new_i, il, io, rs, rsh, n, vth, ns)) < atol

        # updating array of i's
        prec_i[idx] = new_i

    # vv is a np.array of float64 and prec_i is a np.array of mp.mpf 
    return vv, prec_i


def plotter(il, io, rs, rsh, n, vth, ns, atol, num_pts, case):
    # plot a single IV curve 
    plt.xlabel('Voltage')
    plt.ylabel('Current')

    if case == 1: plt.title('Case 1')
    else: plt.title('Case 2') # we're in case 2

    v_vals, i_vals = get_precise_i(il, io, rs, rsh, n, vth, ns, atol, num_pts)
    return plt.plot(v_vals, i_vals)



if __name__ == "__main__":
    mp.dps = 20 # set precision
    num_pts = 100 
    atol = 1e-16 

    # Boltzmann's const (J/K), electron charge (C), temp (K) 
    k, q, temp_cell = [1.380649e-23, 1.60217663e-19, 298.15]
    vth = (k * temp_cell) / q

    case1 = True
    case = [1 if case1 else 2][0]

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

    plt.style.use('seaborn-darkgrid')
    plot = plt.plot()
    for il, io, rs, rsh, n in product(IL, IO, RS, RSH, N):
        plot += plotter(il, io, rs, rsh, n, vth, ns, atol, num_pts, case)
    
    plt.show()


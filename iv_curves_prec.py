import pvlib
import numpy as np
from mpmath import *
from itertools import product
import matplotlib.pyplot as plt

plt.style.use('seaborn-darkgrid')


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


def diff_lhs_rhs(v, i, il, io, rs, rsh, n, ns, vth):
    # returns the difference between the lhs and rhs of the single diode equation
    return (il - io*(mp.exp((v + i*rs)/(n*ns*vth)) - 1) - (v + i*rs)/(rsh) - i)


# Boltzmann's const (J/K), electron charge (C), temp (K) 
k, q, temp_cell = [1.380649e-23, 1.60217663e-19, 298.15]
vth = (k * temp_cell) / q

num_pts = 100
mp.dps = 20
tol = 1e-16


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


# setting up plot
plot = plt.plot()
plt.xlabel('Voltage')
plt.ylabel('Current')
if case1: plt.title('Case 1')
else: plt.title('Case 2')


for il, io, rs, rsh, n in product(IL, IO, RS, RSH, N):
    res = pvlib.pvsystem.singlediode(il, io, rs, rsh, n*ns*vth, num_pts)
    vv = res['v']
    ii = res['i']
    prec_i = np.zeros_like(ii) # initialize 'empty' array for new, more precise i's

    assert len(vv) == len(ii)
    for idx in range(len(vv)):
        # check if i val already precise enough
        if abs(diff_lhs_rhs(vv[idx], ii[idx], il, io, rs, rsh, n, ns, vth)) < tol: 
            new_i = ii[idx]

        else:
            # findroot takes the function whose roots we want to find, a good guess for where the root is, and a tolerance
            # default solver uses secant method
            new_i = findroot(lambda i: diff_lhs_rhs(vv[idx], i, il, io, rs, rsh, n, ns, vth), ii[idx], tol=tol)

        # check that findroot did what we wanted 
        assert abs(diff_lhs_rhs(vv[idx], new_i, il, io, rs, rsh, n, ns, vth)) < tol

        # updating array of i's
        prec_i[idx] = new_i

    # graph
    plot += plt.plot(vv, prec_i)

plt.show()


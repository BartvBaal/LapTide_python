#!/bin/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

import Legendre

def shoot (m, lam, is_even):
    """Shoot ODE and return residual"""

    legendre_solver = Legendre.solver_t(m, is_even)

    return legendre_solver(lam)
        
## -------------------------------------------------------------------------- ##

def plot_coeffs (m, lam, is_even):
    """Plot the coefficients of ODE at set of even points"""

    leg = Legendre.ODE_t(m, lam)

    N = 21
    steps = np.linspace(Legendre.t0, Legendre.t1, N)[1:-1]
    coeffs = leg.coeffs(steps)

    f, axarr = plt.subplots(nrows=2, ncols=2, sharex=True)
    
    axarr[0, 0].plot(steps, coeffs[0], 'b.-', label='c00')
    axarr[0, 1].plot(steps, coeffs[1], 'b.-', label='c01')
    axarr[1, 0].plot(steps, coeffs[1], 'b.-', label='c10')
    axarr[1, 1].plot(steps, coeffs[1], 'b.-', label='c11')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_ODE (m, lam, is_even):
    """Plot a solution for ODE at integration points"""

    legendre_solver = Legendre.solver_t(m, is_even)

    steps, solun = legendre_solver.save(lam)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(steps, solun[:,0], 'b.-', label='Plm')
    axarr[1].plot(steps, solun[:,1], 'b.-', label='dPlm')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_interp (m, lam, is_even, N):
    """Plot a solution for ODE at set of even N points"""

    legendre_solver = Legendre.solver_t(m, is_even)

    steps = np.linspace(Legendre.t0, Legendre.t1, N)
    [P, Q] = legendre_solver.interp(lam, steps)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(steps, P, 'b.-', label='Plm interp')
    axarr[1].plot(steps, Q, 'b.-', label='dPlm interp')
    plt.show()

## -------------------------------------------------------------------------- ##

def main ():
    if len(sys.argv) != 3:
        raise RuntimeError("require m and l")

    m = int(sys.argv[1])
    l = int(sys.argv[2])

    if m > l:
        raise RuntimeError("Legendre polynomials undefined for m>l, input: ({}, {})".format(m, l))

    lam = Legendre.eigval(l)
    is_even = Legendre.check_is_even(m, l)

    print "Residual at {}: {}".format(Legendre.t1, shoot(m, lam, is_even))

    print "Plotting coefficients of ODE"
    plot_coeffs(m, lam, is_even)

    print "Plotting ODE solution at integration points"
    plot_ODE(m, lam, is_even)

    print "Plotting ODE solution at interpolation points"
    plot_interp(m, lam, is_even, 21)


if __name__ == "__main__":
    main()




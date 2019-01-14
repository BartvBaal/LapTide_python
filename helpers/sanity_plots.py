import matplotlib.pyplot as plt
import numpy as np

import classes.LaPlace as LaPlace
import classes.Legendre as Legendre

## -------------------------------------------------------------------------- ##

def plot_coeffs(ode_class, ode_args, is_even):
    """Plot the coefficients of ODE at set of even points"""

    ode = ode_class.ODE_t(*ode_args)

    N = 21
    steps = np.linspace(ode_class.t0, ode_class.t1, N)[1:-1]
    coeffs = ode.coeffs(steps)

    f, axarr = plt.subplots(nrows=2, ncols=2, sharex=True)
    
    axarr[0, 0].plot(steps, coeffs[0], 'b.-', label='c00')
    axarr[0, 1].plot(steps, coeffs[1], 'b.-', label='c01')
    axarr[1, 0].plot(steps, coeffs[1], 'b.-', label='c10')
    axarr[1, 1].plot(steps, coeffs[1], 'b.-', label='c11')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_ODE(ode_class, ode_args, lam, is_even):
    """Plot a solution for ODE at integration points"""
    ode_args.append(is_even)
    ode_solver = ode_class.solver_t(*ode_args)

    steps, solun = ode_solver.save(lam)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(steps, solun[:,0], 'b.-', label='Plm')
    axarr[1].plot(steps, solun[:,1], 'b.-', label='dPlm')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_interp(ode_class, ode_args, lam, is_even, N):
    """Plot a solution for ODE at set of even N points"""
    ode_args.append(is_even)
    ode_solver = ode_class.solver_t(*ode_args)

    steps = np.linspace(ode_class.t0, ode_class.t1, N)
    [P, Q] = ode_solver.interp(lam, steps)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(steps, P, 'b.-', label='Plm interp')
    axarr[1].plot(steps, Q, 'b.-', label='dPlm interp')
    plt.show()

## -------------------------------------------------------------------------- ##

#!/bin/bin/python

import seaborn as sns
import pandas as pd
import copy
import random
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy import optimize

import Legendre
import LaPlace
import plotting

def shoot (m, lam, is_even):
    """Shoot ODE and return residual"""

    legendre_solver = Legendre.solver_t(m, is_even)

    return legendre_solver(lam)

def alt_shoot(lam, m, is_even):
    """Returns the shoot method but with lam as first parameter
       so that one can use rootfinding algortims on lam (instead of m)

       alt_shoot(a,b,c) == shoot(b,a,c)"""
    return shoot(m, lam, is_even)

def rootfinder(m, lamlist, is_even):
    """Finds the root in the lamlow-lamhigh range. Uses the alt_shoot
       function as rootfinding works on the first specified parameter"""
    lamlow, lamhigh = lamlist
    return optimize.bisect(alt_shoot, lamlow, lamhigh, args=(m, is_even), full_output=True)
#TODO: use class with rootfinder

"""
Note about the rootfinder ^:
Need to rootfind on lambda, which thus needs to be the first parameter if we
want to use bisect - so the alt_shoot function is just shoot but with rotated
parameters to allow for the rootfinding on shoot.
Note that there was a small error in Legendre.check_if_even - (m+l)%2 was 
missing the () around the m+l term (so it did m+(l%2)) which is obviously not
going to return the expected True/False output
"""

def shoot_laplace(lam, m, q, is_even):
    laplace_solver = LaPlace.solver_t(m, q, is_even)
    return laplace_solver(lam)

def rootfinder_laplace(m, q, lamlist, is_even):
    lamlow, lamhigh = lamlist
    return optimize.bisect(shoot_laplace, lamlow, lamhigh, args=(m, q, is_even), full_output=True)

def multi_rootfind(m, qlist, lamlist, is_even):
    """
    does rootfinder_laplace for all values in q for a set m and k (and thus is_even)
    If the default new lamlist doesnt stradle the new point, algorithm will
    modify lamlist a bit before trying again. Limited options so far
    Now scales with stepsize of the qlist. Currently works for m=-2 setups ala Townsend,
    but that's because they dont have a decreasing slope when lambda is positive
        making root+diff -> root+(1-(10*shift))*diff for lower guess might/should
        fix that, but it might break the current setup from working (which takes
        ~1 hour to run so I am not testing that friday @1645)
    """
    found_lamlist = []
    oldroot = (lamlist[0]+lamlist[1])/2
    shift = np.abs(np.diff(qlist)[0])

    for q in qlist:

        root = rootfinder_laplace(m, q, lamlist, is_even)[0]
        diff = root - oldroot
        oldroot = root
        if diff < 1.5*shift:
            lamlist = [root-(1000*shift*diff), root+(1000*shift*diff)]
        else:
            lamlist = [root+diff, root+(1+(10*shift))*diff]
#        try:
#            root = rootfinder_laplace(m, q, lamlist, is_even)[0]
#        except:
#            print "excepting"
#            try:
#                lamlist = [root, root*(1 + 4*shift)]  # Try 4% if not working
#                root = rootfinder_laplace(m, q, lamlist, is_even)[0]
#            except:
#                try:
#                    print "bigger excepting"
#                    lamlist = [root, root*(1 + 7*shift)] # Try 7% next
#                    root = rootfinder_laplace(m, q, lamlist, is_even)[0]
#                except:
#                    print "negative excepting"
#                    lamlist = [root*(1 - 3*shift), root] # Try reduced 3%
#                    root = rootfinder_laplace(m, q, lamlist, is_even)[0]
#        lamlist = [root, root*(1 + 2*shift)]  # Take +.25-2% around last lambda as new guesses

        found_lamlist.append(root)
        print q, root, "checking", diff, lamlist

#    plt.plot(qlist, found_lamlist)
#    plt.yscale('log')
#    plt.show()
    found_lamlist = np.asarray(found_lamlist, dtype=float)
    return qlist, found_lamlist

def fullrange_multi_rootfind(m, qlists, lvals, townsendcomp=False):
    """
    Does multi_rootfind per qlist in qlists, for all values of l
    Will save the data as part of the process so the plotting can be done
    much quicker in later times

    TODO; make saving location function argument
    """
    if townsendcomp:
        asymptotic_plotting()
    for qlist in qlists:
        for l in lvals:
            k = l - np.abs(m)
            is_even = LaPlace.check_is_even(m, k)
            lamlist = [l*(l-.5025), l*(l+2.5)]

            print is_even

            qlist, found_lamlist = multi_rootfind(m, qlist, lamlist, is_even)
            plt.plot(qlist, found_lamlist, color="black", ls="--")
            savestring = "Numerics/Townsend2003/range_"+str(qlist[0])+"_"+str(qlist[-1])+\
                            "_steps_"+str(len(qlist))+"_kval_"+str(k)+".txt"
            np.savetxt(savestring, found_lamlist)
    plt.yscale('log')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_coeffs (m, q, lam, is_even):
    """Plot the coefficients of ODE at set of even points"""

#    leg = Legendre.ODE_t(m, lam)
    lap = LaPlace.ODE_t(m, q, lam)

    N = 21
    steps = np.linspace(Legendre.t0, Legendre.t1, N)[1:-1]
    coeffs = lap.coeffs(steps)

    f, axarr = plt.subplots(nrows=2, ncols=2, sharex=True)
    
    axarr[0, 0].plot(steps, coeffs[0], 'b.-', label='c00')
    axarr[0, 1].plot(steps, coeffs[1], 'b.-', label='c01')
    axarr[1, 0].plot(steps, coeffs[1], 'b.-', label='c10')
    axarr[1, 1].plot(steps, coeffs[1], 'b.-', label='c11')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_ODE (m, q, lam, is_even):
    """Plot a solution for ODE at integration points"""

#    legendre_solver = Legendre.solver_t(m, is_even)
    laplace_solver = LaPlace.solver_t(m, q, is_even)

    steps, solun = laplace_solver.save(lam)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(steps, solun[:,0], 'b.-', label='Plm')
    axarr[1].plot(steps, solun[:,1], 'b.-', label='dPlm')
    plt.show()

## -------------------------------------------------------------------------- ##

def plot_interp (m, q, lam, is_even, N):
    """Plot a solution for ODE at set of even N points"""

#    legendre_solver = Legendre.solver_t(m, is_even)
    laplace_solver = LaPlace.solver_t(m, q, is_even)

    steps = np.linspace(Legendre.t0, Legendre.t1, N)
    [P, Q] = laplace_solver.interp(lam, steps)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(steps, P, 'b.-', label='Plm interp')
    axarr[1].plot(steps, Q, 'b.-', label='dPlm interp')
    plt.show()

## -------------------------------------------------------------------------- ##

def asymptotic_plotting():
    """
    Plots the asymptotic functions to the Laplace Tidal Equations
    Solid lines for second order, dashdotted line for the first order
    """
#    l_eq_39 = lambda m, s, q : ( (m*q - m**2) / (q * (2.*s + 1.)) )**2  # < Need to update this to 2nd order term
#    l_eq_40_41 = lambda m, q : (q - m)**2

    eq38_fst = lambda m, s, q : (q**2) * ((2.*s + 1.)**2)
    eq38_snd = lambda m, s, q : (q**2) * ((2.*s + 1.)**2) - ( 2. * (m*q - m**2) )

    eq40 = lambda m, q : (q - m)**2

    eq55 = lambda m, q : (m**2) * 2.*m*q / (2*m*q + 1.)

    ## q > condition statement (eg q>1 for g_mode)
    g_mode_cond = np.abs(1)
    r_mode_cond = lambda m, s : ((2.*s + 1.)**2 + m**2) / np.abs(m)  # paper mistake in 2s+1 term
    k_mode_cond = lambda m : np.abs(3. / m)
    #y_mode_cond = lambda m, q : (m*q < m**2 and m*q > 0) ? r_mode_cond(m, 0) : g_mode_cond

    m = -2
    qneg = np.linspace(-10, 0, 10000)
    qpos = np.linspace(0, 10, 10000)

    plotting.plot_function(eq38_fst, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-.")
    plotting.plot_function(eq38_snd, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-")
    plotting.plot_function(eq38_fst, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-.")
    plotting.plot_function(eq38_snd, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-")
    plotting.plot_function(eq38_fst, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-.")
    plotting.plot_function(eq38_snd, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-")
    plotting.plot_function(eq38_fst, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-.")
    plotting.plot_function(eq38_snd, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-")
    plotting.plot_function(eq40, [m], qpos, [np.abs(qpos)>1.], "green", "-")  # its effecitvely a gmode constraint
    plotting.plot_function(eq55, [m], qpos, [np.abs(qpos)>np.sqrt(3/(-m*qpos))], "cyan", "-")
    plt.yscale('log')
##    plt.plot(left[left < -2.5], test_40_s0_l, color="red", label="40")
#    plt.plot(right, test_40_s0_r, color="green", label= "40")
#    plt.plot(right, test_55_s1_r, color="cyan", label="55")
#    plt.plot(left, test_38_s3_l, color="orange", label="38, s=3")
#    plt.plot(left, test_38_s2_l, color="yellow", label="38, s=2")
##    plt.plot(left[left < -6.], test_39_s1_l, color="black", label="39, s=1")
#    plt.yscale('log')
#    plt.legend(fontsize=16, frameon=True, fancybox=True, edgecolor="#000066")
#    plt.xlim([-10, 10])

def numerics_plotting():
    """
    Plots a dotted line for the loaded numerics values

    TODO; pass on a folder or conditions on how to load multiple plot_from_files
    """
    plotting.plot_from_file("Numerics/Townsend2003/range_0.0_10.0_steps_952_kval_0.txt", "black", "--")
    plotting.plot_from_file("Numerics/Townsend2003/range_0.0_10.0_steps_952_kval_1.txt", "black", "--")
    plotting.plot_from_file("Numerics/Townsend2003/range_0.0_10.0_steps_952_kval_2.txt", "black", "--")
    plotting.plot_from_file("Numerics/Townsend2003/range_0.0_-10.0_steps_951_kval_0.txt", "black", "--")
    plotting.plot_from_file("Numerics/Townsend2003/range_0.0_-10.0_steps_951_kval_1.txt", "black", "--")
    plotting.plot_from_file("Numerics/Townsend2003/range_0.0_-10.0_steps_951_kval_2.txt", "black", "--")
    plt.yscale('log')

## -------------------------------------------------------------------------- ##

def main ():
#    Relic from initial setups; need to look if this reworking setup requires this
    if len(sys.argv) != 3:
        raise RuntimeError("require m and l")

    m = int(sys.argv[1])
    l = int(sys.argv[2])

#    if m > l:
#        raise RuntimeError("Legendre polynomials undefined for m>l, input: ({}, {})".format(m, l))

#    lam = Legendre.eigval(l)
#    is_even = Legendre.check_is_even(m, l)
#    print is_even, m, l

#    lamlist = [l*(l-2), l*(l+3)]
#    root = rootfinder(m, lamlist, is_even)
#    print "Rootfinding solution: {}, original eigenvalue: {}".format(root, lam)
#    print "Using new root:"
#    print "Residual at {}: {}\n".format(Legendre.t1, shoot(m, root[0], is_even))

#    print "Using old root:"
#    print "Residual at {}: {}\n".format(Legendre.t1, shoot(m, lam, is_even))

    ## Works for m=1 setup - need to generalize this (ask Frank how is_even changes & how to make decent guesses for lambda)
    ## Could always start in q = 0 point since Legendre polynomials have known lambda
    ## Need to check what lambda does for negative m though - still l(l+1)?
    ## Also; k = l-|m| and Townsend uses k for different wave modes
    k = l - np.abs(m)
    is_even = LaPlace.check_is_even(m, k)
    lamlist = [30, 100]
    q = 2.  # If q=0 you have the legendre polynomials back

#    root = rootfinder_laplace(m, q, lamlist, is_even)
#    print "Rootfinding solution: {}, l-param: {}".format(root[0], l)  # root itself has much more info
#    print "Using rootfinder:"
#    print "Residual at {}: {}\n".format(Legendre.t1, shoot_laplace(root[0], m, q, is_even))

    qneg = np.linspace(0, -10, 9.5e2+1)  # Currently have split qpos and qneg since I know lambda for q=0 but not q=-10 or 10
    qpos = np.linspace(0, 10, 9.5e2+2)
    lamlist = [l*(l-1), l*(l+3)]  # Starting point since q=0 is legendre
    print is_even

#    multi_rootfind(m, qpos, lamlist, is_even)

    qlists = [qneg, qpos]
    lvals = [2, 3, 4]  # k=0,1,2

#    fullrange_multi_rootfind(m, qlists, lvals, townsendcomp=True)  # Mostly for plotting functionality
#    fullrange_multi_rootfind(m, [np.linspace(0, -10, 9.5e2+1)], [4])  # Testing just for k=2, negative part

    asymptotic_plotting()
    numerics_plotting()
    plt.xlim([-10, 10])
    plt.ylim([.1, 6800])  #6800 works for this setup
    plt.title("Asymptotic first and second order versus numerical solutions")
    plt.xlabel("q (spinparameter)")
    plt.ylabel(r"$\lambda$(q)")
    plt.show()

#    print "Plotting coefficients of ODE"
#    plot_coeffs(m, q, root[0], is_even)
#    plot_coeffs(m, lam, is_even)

#    print "Plotting ODE solution at integration points"
#    plot_ODE(m, q, root[0], is_even)
#    plot_ODE(m, lam, is_even)

#    print "Plotting ODE solution at interpolation points"
#    plot_interp(m, q, root[0], is_even, 21)
#    plot_interp(m, lam, is_even, 21)

if __name__ == "__main__":
    main()




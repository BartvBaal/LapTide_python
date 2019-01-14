#!/bin/bin/python

import copy
import random
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy import optimize

import classes.Legendre as Legendre
import classes.LaPlace as LaPlace

import helpers.Property as Property
import helpers.Straddle as Straddle
import helpers.rootfinder as roots
import helpers.plotting as plotting
import helpers.sanity_plots as sanplot


def fullrange_multi_rootfind(m, qlists, lvals, townsendcomp=False, saving=False):
    """
    Does multi_rootfind per qlist in qlists, for all values of l
    Will save the data as part of the process so the plotting can be done
    much quicker in later times, if the saving argument is used
    Will plot a comparison to Townsend2003's asymptotic approximations if that
    argument is used, otherwise it will only plot the found solutions
    """
    if townsendcomp:
        plotting.asymptotic_plotting()
    for qlist in qlists:
        for l in lvals:
            is_even = LaPlace.check_is_even(m, l)
            print is_even

            qlist, found_lamlist = roots.multi_rootfind(m, l, qlist, is_even)
            plt.plot(qlist, found_lamlist, color="black", ls="--")
            if saving:
                savestring = "Numerics/Townsend2003/range_"+str(qlist[0])+"_"+str(qlist[-1])+\
                                "_steps_"+str(len(qlist))+"_kval_"+str(k)+".txt"
                np.savetxt(savestring, found_lamlist)
    plt.yscale('log')
    plt.show()


def main():
    # Need to consider if I want the inputs as m, l or m, k - using l for now
    if len(sys.argv) != 3:
        raise RuntimeError("require m and l")

    m = int(sys.argv[1])
    l = int(sys.argv[2])

    mode_admin = Property.mode_admin(m, l)
    mode_admin.validate_values()
    is_even = mode_admin.is_even()
    k = mode_admin.get_k()
    wavemode = mode_admin.get_wavemode()

    # Currently have split qpos and qneg since I know lambda for q=0 but not q=-10 or 10
    qneg = np.linspace(0, -10, 9.5e2+1)  
    qpos = np.linspace(0, 10, 9.5e2+2)

    print is_even, k, wavemode

#    roots.multi_rootfind(m, l, qpos, is_even)  # Testing if the new function works

    qlists = [qneg, qpos]
    lvals = [2, 3, 4]  # k=0,1,2

    fullrange_multi_rootfind(m, qlists, lvals, townsendcomp=True)  # Mostly for plotting functionality
#    fullrange_multi_rootfind(m, [qneg], [4], townsendcomp=True)  # Testing just for k=2, negative part

#    plotting.asymptotic_plotting()
#    plotting.townsend_plotting()
#    plt.xlim([-10, 10])
#    plt.ylim([.1, 6800])  #6800 works for this setup
#    plt.title("Asymptotic first and second order versus numerical solutions")
#    plt.xlabel("q (spinparameter)")
#    plt.ylabel(r"$\lambda$(q)")
#    plt.show()

#    # Sanity plots; testing setup w/ q=0
#    q = 0
#    print "Plotting coefficients of ODE"
#    sanplot.plot_coeffs(LaPlace, [m, q, l*(l+1)], is_even)

#    print "Plotting ODE solution at integration points"
#    sanplot.plot_ODE(LaPlace, [m, q], l*(l+1), is_even)

#    print "Plotting ODE solution at interpolation points"
#    sanplot.plot_interp(LaPlace, [m, q], l*(l+1), is_even, 21)


if __name__ == "__main__":
    main()




#!/bin/bin/python
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import classes.Legendre as Legendre
import classes.LaPlace as LaPlace
import classes.Curvilinear as Curvilinear

import helpers.Property as Property
import helpers.Straddle as Straddle
import helpers.rootfinder as roots
import helpers.plotting as plotting
import helpers.sanity_plots as sanplot
import helpers.morsink_radius as oblate
import helpers.LaPlace_asymptotes as asym


def fullrange_multi_rootfind(m, kvals, qlists, aympcompare=False, saving=False):
    """
    Does multi_rootfind per qlist in qlists, for all values of k given in kvals
    Will save the data as part of the process so the plotting can be done
    much quicker in later times, if the saving argument is used
    Will plot a comparison to Townsend2003's asymptotic approximations if that
    argument is used, otherwise it will only plot the found solutions

    Should replace the is_even check with the mode_admin setup I created
    this will allow for more flexibility down the line with e.g. varying k *and* m
    """
    if aympcompare:
        plotting.asymptotic_plotting(m)
    for qlist in qlists:
        for k in kvals:
            is_even = LaPlace.check_is_even(m, k)
            print is_even

            qlist, found_lamlist = roots.multi_rootfind(m, k, qlist, is_even)
            plt.plot(qlist, found_lamlist, color="black", ls="--")
            if saving:
                savestring = "data/Townsend2003/range_{}_{}_steps_{}_kval_{}.txt"\
                                    .format(qlist[0],qlist[-1],len(qlist),str(k))
                np.savetxt(savestring, found_lamlist)
    plt.yscale('log')
    plt.show()


def fullrange_multi_rootfind_curvi(m, kvals, qlists, r_eq, mass, period, \
                                    aympcompare=False, saving=False):
    """
    Curvilinear version of fullrange_multi_rootfind

    Not yet in final version! and DO NOT USE the saving argument!
    """
    if aympcompare:
        plotting.asymptotic_plotting(m)
    for qlist in qlists:
        for k in kvals:
            is_even = Curvilinear.check_is_even(m, k)
            print is_even

            qlist, found_lamlist = roots.multi_rootfind_curvilinear(m, k, qlist, is_even, r_eq, mass, period)
            plt.plot(qlist, found_lamlist, color="black", ls="--")
            if saving:
                savestring = "data/Curvilinear/range_{}_{}_steps_{}_kval_{}.txt"\
                                    .format(qlist[0],qlist[-1],len(qlist),str(k))
                np.savetxt(savestring, found_lamlist)
    plt.yscale('log')
    plt.show()


def main():
    # Need to consider if I want the inputs as m, l or m, k - using k for now
    if len(sys.argv) != 3:
        raise RuntimeError("require m and l")

    m = int(sys.argv[1])
    k = int(sys.argv[2])

    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()  # Mode_admin now checks this upon init!
    is_even = mode_admin.is_even()
    l = mode_admin.get_l()
    wavemode = mode_admin.get_wavemode()

    # Currently have split qpos and qneg since I know lambda for q=0 but not q=-10 or 10
    qneg = np.linspace(0, -10, 9.5e2+4)
    qpos = np.linspace(0, 10, 9.5e2+4)

    print is_even, l, wavemode

#    roots.multi_rootfind(m, l, qpos, is_even)  # Testing if the new function works

    qlists = [qneg, qpos]
    kvals = [0, 1, 2]  # l=2,3,4

#    r_list = np.linspace(9e3, 16.5e3, 375)
#    m_list = np.linspace(1.2*1.9885e30, 2.75*1.9885e30, 325)
#    period = 1./581 # 4U 1636-536
##    period = 1./363 # 4U 1728-34
#    oblate.recover_radius_mass(r_list, m_list, period, 0.1)

    r_eq = 12000  # in meters
    mass = 1.8*1.9855e30  # in kg
    period = np.inf  # no period so no spin so **should match old equations** like this
#    fullrange_multi_rootfind_curvi(m, kvals, qlists, r_eq, mass, period, aympcompare=True)

    init_guess = asym.r_modes(m, k, -50.)
    print init_guess
    r_qlist = np.linspace(-50., -6.05, 250)  # r-modes are fine with far fewer steps really
    qlist, found_lamlist = roots.multi_rootfind_curvilinear_new(m, r_qlist, is_even, init_guess, r_eq, mass, period, verbose=False, inc=1.05)

    eq39 = lambda m, s, q : (m*q - m*m)**2 / (q*q * (2*s + 1)**2)
    eq39_2nd = lambda m, s, q : (m*q - m*m)**2 / (q*q * (2*s + 1)**2) + (2. * ((m*q - m*m)**3)) / (q**4 * ((2.*s + 1)**4))
    s = -k-1
    asy = eq39(m, s, r_qlist)
    asy_2nd = eq39_2nd(m, s, r_qlist)
    plt.plot(qlist, found_lamlist)
    plt.plot(qlist, asy, ls="--")
    plt.plot(qlist, asy_2nd, ls=":")
    plt.yscale('log')
    plt.show()

#    start_time_new = time.time()
#    roots.multi_rootfind_curvilinear_new(m, qneg, is_even, 5.5, r_eq, mass, period)
#    end_time_new = time.time()
#    start_time_old = time.time()
#    roots.multi_rootfind_curvilinear(m, 0, qneg, is_even, r_eq, mass, period)
#    end_time_old = time.time()
#    print "New function took: {} seconds, old function took: {} seconds".format(end_time_new - start_time_new, end_time_old - start_time_old)

    #TODO: need to fix the initial guess for higher spins
    # For 363/581 Hz the initial guess doesn't stradle a root!

#    fullrange_multi_rootfind(m, kvals, qlists, aympcompare=True)  # Mostly for plotting functionality
#    fullrange_multi_rootfind(m, [-2], [qneg], aympcompare=True)  # Testing just for k=-2, negative part

#    plotting.asymptotic_plotting(m)
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




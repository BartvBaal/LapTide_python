#!/bin/bin/python
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import classes.Legendre as Legendre
import classes.LaPlace as LaPlace

import helpers.Property as Property
import helpers.Straddle as Straddle
import helpers.rootfinder as roots
import helpers.plotting as plotting
import helpers.sanity_plots as sanplot
import helpers.morsink_radius as oblate


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


def main():
    # Need to consider if I want the inputs as m, l or m, k - using k for now
    if len(sys.argv) != 3:
        raise RuntimeError("require m and l")

    m = int(sys.argv[1])
    k = int(sys.argv[2])

    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()
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

    r_list = np.linspace(9e3, 16.5e3, 375)
    m_list = np.linspace(1.2*1.9885e30, 2.75*1.9885e30, 325)
    period = 1./581 # 4U 1636-536
#    period = 1./363 # 4U 1728-34
    oblate.recover_radius_mass(r_list, m_list, period, 0.1)

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




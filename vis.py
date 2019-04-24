# -*- coding: utf-8 -*-
#!/bin/bin/python
import sys
import numpy as np
import itertools
from functools import partial
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D

import classes.Legendre as Legendre
import classes.LaPlace as LaPlace
import classes.Curvilinear as Curvilinear

import helpers.Property as Property
import helpers.Straddle as Straddle
import helpers.rootfinder as roots
import helpers.plotting as plotting
import helpers.LaPlace_asymptotes as asym
import helpers.Curvilinear_asymptotes as curvasym
import helpers.sanity_plots as sanplot
import helpers.morsink_radius as oblate
import helpers.gravity_functions as grav

# Changing plotting parameters to bigger values for readable large plots
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 20, 16
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['axes.titlesize'] = 23
mpl.rcParams['axes.labelsize'] = 28
mpl.rcParams['axes.formatter.limits'] = -5, 5
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 1.6
mpl.rcParams['xtick.labelsize'] = 25
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.major.width'] = 1.6
mpl.rcParams['ytick.minor.size'] = 4.5
mpl.rcParams['ytick.minor.width'] = 1.2
mpl.rcParams['ytick.labelsize'] = 25


def fullrange_multi_rootfind(m, qlists, kvals, aympcompare=False, saving=False):
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


def multi_eccentricity_rootfind(m, k, qlist, is_even, r_eq, mass, periodlist):
    for period in periodlist:
        qlist, found_lamlist = roots.multi_rootfind_curvilinear(m, k, qlist, is_even, r_eq, mass, period)
        plt.plot(qlist, found_lamlist, label="{} Hz".format(1./period))
    plt.legend(fontsize=24, frameon=True, fancybox=True, edgecolor="#000000")
    plt.title(r"Wave m: {}, k: {}, with mass {} $M_\odot$ and radius {} km".format(\
                m, k, mass/1.9885e30, r_eq*1e-3))
    plt.yscale('log')
    plt.show()


def rootfind_dimless(m, k, qlist, ecc=0., chi=0., gravfunc=grav.chi_gravity_deriv, verbose=False, inc=1.0033):
    """
    For a given m, k and qlist, and potentially for different radius, mass and
    period as well, determines the wave mode and calculates the eigenvalues
    from the asymptotically calculated values.
    """
    dlngrav=partial(gravfunc, chi)
    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()
    is_even = mode_admin.is_even()
    l = mode_admin.get_l()

    if np.average(qlist)*m < 0:
        direction = "pro"
    else:
        direction = "retro"

    wavemode = mode_admin.get_wavemode(direction)
    print is_even, l, wavemode

    wavemode += "s"
    if wavemode[0] == "g":
        wavemode += "_list"
    wavemode = wavemode.replace(" ", "_")
    if wavemode[0] == "y" or wavemode[0] == "k":  # yanai and kelvin modes only have two arguments
        args = m, qlist, ecc, chi
    else:
        args = m, k, qlist, ecc, chi
    guesslist = getattr(curvasym, wavemode)(*args)

    qlist, found_lamlist = roots.multi_rootfind_fromguess_dimless(m, qlist, is_even, guesslist, ecc, dlngrav, verbose=False, inc=inc)
    plt.plot(qlist, found_lamlist)


def simple_numasyplot(m):
    custom_lines = [Line2D([0], [0], color="black", lw=2.5, linestyle="--", label="Numeric"),
                    Line2D([0], [0], color="blue", lw=2.5, label="2nd order"),
                    Line2D([0], [0], color="blue", lw=2.5, linestyle="-.", label="1st order")]
    plotting.asymptotic_plotting(m)
    plotting.townsend_plotting()
    plt.xlim([-5, 5])
    plt.ylim([3.8, 1400])  #5500 works for this setup
    plt.title("Asymptotic 1st and 2nd order compared with numerical solutions, m=-2")
    plt.xlabel(r"Spin parameter q = $2\Omega/\omega$")
    plt.ylabel(r"Eigenvalue $\lambda$")
#    plt.text(-4.25, 1000, "g-mode", fontsize=21)
#    plt.text(-4.75, 310, "g-mode", fontsize=21)
#    plt.text(-4.75, 95, "g-mode", fontsize=21)
#    plt.text(3.75, 220, "g-mode", fontsize=21)
#    plt.text(4., 28, "y-mode", fontsize=21)
#    plt.text(4, 4.5, "k-mode", fontsize=21)
#    plt.legend(handles=custom_lines, fontsize=24, frameon=True, fancybox=True, edgecolor="#000000", loc='upper right', bbox_to_anchor=(1, 1))
#    plt.savefig("TEST.png")
    plt.show()


def simple_numasyplot_rmodesINC(m, k, is_even):
    r_eq, mass, period = 10e3, 1.4*1.9885e30, np.inf
    r_qlist = np.linspace(-10., -6.025, 225)  # r-modes are fine with fewer steps really
    y_qlist = np.linspace(-10., -3.05, 325)
    r_guesslist = asym.r_modes(m, k, r_qlist)
    y_guesslist = asym.yanai_modes(m, y_qlist)
    r_qlist, r_found_lamlist = roots.multi_rootfind_fromguess(m, r_qlist, is_even, r_guesslist, r_eq, mass, period, verbose=False, inc=1.05)
    y_qlist, y_found_lamlist = roots.multi_rootfind_fromguess(m, y_qlist, not is_even, y_guesslist, r_eq, mass, period, verbose=False, inc=1.08)

    custom_lines = [Line2D([0], [0], color="black", lw=2.5, linestyle="--", label="Numeric"),
                    Line2D([0], [0], color="blue", lw=2.5, label="2nd order"),
                    Line2D([0], [0], color="blue", lw=2.5, linestyle="-.", label="1st order")]
    plotting.asymptotic_plotting(m)
    plotting.townsend_plotting()

    plt.plot(r_qlist, r_found_lamlist, ls="--", color='black')
    plt.plot(y_qlist, y_found_lamlist, ls="--", color='black')

    plt.xlim([-10, 10])
    plt.ylim([.1, 6000])
    plt.title("Asymptotic 1st and 2nd order compared with numerical solutions, m=-2")
    plt.xlabel(r"Spin parameter q = $2\Omega/\omega$")
    plt.ylabel(r"Eigenvalue $\lambda$")
    plt.show()


def compare_newfile(m):
    qneg = np.linspace(0, -10, 9.5e2+4)
    qpos = np.linspace(0, 10, 9.5e2+4)
    plotting.asymptotic_plotting(m)
    pro1 = asym.g_modes_list(m, 2, qpos)
    pro2 = asym.yanai_modes(m, qpos)
    pro3 = asym.kelvin_modes(m, qpos)
    ret1 = asym.g_modes_list(m, 2, qneg)
    ret2 = asym.g_modes_list(m, 1, qneg)
    ret3 = asym.g_modes_list(m, 0, qneg)
    plt.plot(qpos, pro1, label="new p1")
    plt.plot(qpos, pro2, label="new p2")
    plt.plot(qpos, pro3, label="new p3")
    plt.plot(qneg, ret1, label="new r1")
    plt.plot(qneg, ret2, label="new r2")
    plt.plot(qneg, ret3, label="new r3")
    plt.xlim([-5, 5])
    plt.ylim([3.8, 1400])
    plt.xlabel(r"Spin parameter q = $2\Omega/\omega$")
    plt.ylabel(r"Eigenvalue $\lambda$")
    plt.legend(fontsize=24, frameon=True, fancybox=True, edgecolor="#000000")
    plt.show()


def compare_curvasym(m):
    qneg = np.linspace(-1, -10, 9.5e2+4)
    qpos = np.linspace(1, 10, 9.5e2+4)
    for ecc in [0, .05, .1]:
        if ecc == .05:
            ls="dotted"
        elif ecc == .1:
            ls="dashed"
        else:
            ls="solid"
        chi = ecc*2
        pro1 = curvasym.g_modes_list(m, 2, qpos, ecc=ecc, chi=chi)
        pro2 = curvasym.yanai_modes(m, qpos, ecc=ecc, chi=chi)
        gpr2 = curvasym.g_modes_list(m, 1, qpos, ecc=ecc, chi=chi)
        pro3 = curvasym.kelvin_modes(m, qpos, ecc=ecc, chi=chi)
        ret1 = curvasym.g_modes_list(m, 2, qneg, ecc=ecc, chi=chi)
        ret2 = curvasym.g_modes_list(m, 1, qneg, ecc=ecc, chi=chi)
        ret3 = curvasym.g_modes_list(m, 0, qneg, ecc=ecc, chi=chi)
        ret4 = curvasym.yanai_modes(m, qneg, ecc=ecc, chi=chi)
        gre4 = curvasym.g_modes_list(m, -1, qneg, ecc=ecc, chi=chi)
        plt.plot(qpos, pro1, ls=ls)
        plt.plot(qpos, pro2, ls=ls)
        plt.plot(qpos, gpr2, ls=ls, color="black")
        plt.plot(qpos, pro3, ls=ls)
        plt.plot(qneg, ret1, ls=ls)
        plt.plot(qneg, ret2, ls=ls)
        plt.plot(qneg, ret3, ls=ls)
        plt.plot(qneg, ret4, ls=ls)
        plt.plot(qneg, gre4, ls=ls, color="black")
    plt.xlim([-10, 10])
    plt.ylim([3.8, 5500])
    plt.xlabel(r"Spin parameter q = $2\Omega/\omega$")
    plt.ylabel(r"Eigenvalue $\lambda$")
    plt.yscale('log')
    plt.legend(fontsize=16, frameon=True, fancybox=True, edgecolor="#000000")
    plt.show()


def compare_curvrmode(m):
    qneg = np.linspace(-50, -6.5, 9.5e2+4)
    rspheric = curvasym.r_modes(m, -2, qneg, ecc=0, chi=0)
    rcurv1 = curvasym.r_modes(m, -2, qneg, ecc=.05, chi=.1)
    rcurv2 = curvasym.r_modes(m, -2, qneg, ecc=.1, chi=.2)
    qnum = np.linspace(-50, -9.5, 250)
#    rootfind_dimless(m, -2, qnum, ecc=0., dlngrav=partial(grav.chi_gravity_deriv, 0.), verbose=False, inc=1.05)
#    rootfind_dimless(m, -2, qnum, ecc=0.05, dlngrav=partial(grav.chi_gravity_deriv, 0.1), verbose=False, inc=1.05)
#    rootfind_dimless(m, -2, qnum, ecc=0.1, dlngrav=partial(grav.chi_gravity_deriv, 0.2), verbose=False, inc=1.05)
    plt.plot(qneg, rspheric, ls='solid', color="black")
    plt.plot(qneg, rcurv1, ls='dotted', color="black")
    plt.plot(qneg, rcurv2, ls='dashed', color="black")
    plt.yscale('log')
    plt.show()


def error_numasyplot(m):
    plotting.difference_plotting(m)
    custom_lines = [Line2D([0], [0], color="black", lw=2.5, linestyle="-.", label="1st order"),
                    Line2D([0], [0], color="black", lw=2.5, label="2nd order")]
    plt.xlim([-5, 5])
    plt.ylim([.003, .6])
    plt.title("Asymptotic 1st and 2nd order relative errors with numerical solutions, m=-2")
    plt.xlabel(r"Spin parameter q = $2\Omega/\omega$")
    plt.ylabel(r"Relative error $\epsilon(\lambda)$")
#    plt.legend(handles=custom_lines, fontsize=18, frameon=True, fancybox=True, edgecolor="#000000")
    plt.show()


def oblate_plot(om_bar_sq):
    degrees = np.linspace(90, 0, 501)
    angles = np.radians(degrees)
    xlist = [.15, .18, .22, .37]
    clist = ["red", "orange", "blue", "black"]

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for x, c in zip(xlist, clist):
        rdata = oblate.calc_radius_14_dimless(x, om_bar_sq, angles)
        gdata = oblate.calc_grav_14_slow_dimless(x, om_bar_sq, angles)
        ax1.plot(degrees, rdata, ls="solid", c=c, label="x = {}".format(x))
        ax2.plot(degrees, gdata, ls="dotted", c=c)
    ax1.set_xlim([0, 90])
    ax1.set_xlabel(r"co-latitudinal angle $\cos(\theta)$")
    ax1.set_ylabel(r"radius $R(\theta) / R_eq$")
    ax2.set_ylabel(r"effective gravity $g(\theta) / g_0$")
    plt.title(r"Radial and gravitational changes as function of co-latitudinal angle $\cos(\theta)$ for $\bar\Omega^2 = {}$".format(om_bar_sq))
    ax1.legend(fontsize=24, frameon=True, fancybox=True, edgecolor="#000066", loc=8)
    plt.show()


def oblate_plot_alternative(om_bar_sq):
    degrees = np.linspace(90, 0, 501)
    angles = np.radians(degrees)
    xlist = [.15, .18, .22, .37]
    lslist = ["-", "--", "-.", ":"]

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for x, ls in zip(xlist, lslist):
        rdata = oblate.calc_radius_14_dimless(x, om_bar_sq, angles)
        gdata = oblate.calc_grav_14_slow_dimless(x, om_bar_sq, angles)
        ax1.plot(degrees, rdata, c="black", ls=ls, label="x = {}".format(x))
        ax2.plot(degrees, gdata, c="red", ls=ls)
    ax1.set_xlim([0, 90])
    ax1.set_xlabel(r"co-latitudinal angle $\cos(\theta)$")
    ax1.set_ylabel(r"radius $R(\theta) / R_eq$")
    ax2.set_ylabel(r"effective gravity $g(\theta) / g_0$", color="red")
    ax2.tick_params(axis='y', colors="red")
    ax2.spines['right'].set_color("red")
    plt.title(r"Radial and gravitational changes as function of co-latitudinal angle $\cos(\theta)$ for $\bar\Omega^2 = {}$".format(om_bar_sq))
    ax1.legend(fontsize=24, frameon=True, fancybox=True, edgecolor="#000066", loc=8)
    plt.show()


def plot_different_dimless(om_bar_sq_values, period=1./581, sn="4U 1636-536",\
                            cmap="Reds", rtol=5e-4, atol=0):
    r_list = np.linspace(9e3, 16.5e3, 3000)
    m_list = np.linspace(1.2*1.9885e30, 2.75*1.9885e30, 2500)

    # bit of code to create colors for the scatter plot, flexible with
    # the amount of values you wanted to use
    cmap = cm.get_cmap(cmap)
    if cmap.name in ["viridis", "inferno", "plasma", "magma"]:
        cmap = cmap.reversed()  # these four have the bright values at 1, so reverse them
    clist = []
    for i in range(len(om_bar_sq_values)):
        clist.append(cmap((i+1.)/len(om_bar_sq_values))) 

    # Run the recover function to find out which of the radius-mass combinations
    # gets within the wanted accuracy of the wanted omega_bar_squared value
    for om_bar_sq, c in zip(om_bar_sq_values, clist):
        data = oblate.recover_radius_mass(r_list, m_list, period, om_bar_sq, rtol=rtol, atol=atol)
        if data.size:
            radius, mass = data[:,0], data[:,1]
            plt.scatter(radius*1e-3, mass/1.9885e30, c=c, label=r"$\bar\Omega^2$: {}".format(om_bar_sq))

    # Some extra layout for the plot, including a rescaling of the legend dotsizes
    plt.xlabel("Equatorial radius (km)")
    plt.ylabel(r"NS mass ($M_\odot$)")
    lgnd = plt.legend(fontsize=22, frameon=True, fancybox=True, edgecolor="#660000")
    for handle in lgnd.legendHandles:
        handle._sizes = [200]
    plt.title(r"Radii and masses for given $\bar\Omega^2$ at {:.0f} Hz ({})".format(1/period, sn))
    plt.show()


def compare_dimless(period=1./581, sn="4U 1636-536", r_cmap="inferno", m_cmap="viridis"):
    r_list = np.linspace(9e3, 16.5e3, 375)
    m_list = np.linspace(1.2*1.9885e30, 2.75*1.9885e30, 325)

    const_R = [9e3, 10.5e3, 12e3, 13.5e3, 15e3, 16.5e3]
    const_M = [1.2*1.9885e30, 1.51*1.9885e30, 1.82*1.9885e30, 2.13*1.9885e30, 2.44*1.9885e30, 2.75*1.9885e30]

    # Some code for the color maps, to keep the brightest colors away
    r_cmap = cm.get_cmap(r_cmap)
    if r_cmap.name in ["viridis", "inferno", "plasma", "magma"]:
        r_cmap = r_cmap.reversed()
    m_cmap = cm.get_cmap(m_cmap)
    if m_cmap.name in ["viridis", "inferno", "plasma", "magma"]:
        m_cmap = m_cmap.reversed()

    for i, M in enumerate(const_M):
        x, om_bar_sq = oblate.find_x_ombarsq(r_list, M, period)
        plt.plot(om_bar_sq, x, ls="-", lw=2.5, label="M={} $M_\odot$".format(M/1.9885e30), c=m_cmap((i+1.)/len(const_M)))
    for i, R in enumerate(const_R):
        x, om_bar_sq = oblate.find_x_ombarsq(R, m_list, period)
        plt.plot(om_bar_sq, x, ls="--", lw=2.5, label="R={} km".format(R*1e-3), c=r_cmap((i+1.)/len(const_R)))
    plt.xlabel(r"$\bar\Omega^2 = \Omega^2 R_{eq}^3 / GM$")
    plt.ylabel(r"$x = GM/R_{eq}c^2$")
    plt.legend(fontsize=22, frameon=True, fancybox=True, edgecolor="#660000")
    plt.title("Dimensionless parameters for constant values of R and M at {:.0f} Hz ({})".format(1/period, sn))
    plt.show()


def rootfind_dimless(m, k, qlist, ecc=0., dlngrav=partial(grav.chi_gravity_deriv, 0.), verbose=False, inc=1.0033):
    """
    For a given m, k and qlist, and potentially for different radius, mass and
    period as well, determines the wave mode and calculates the eigenvalues
    from the asymptotically calculated values.
    """
    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()
    is_even = mode_admin.is_even()
    l = mode_admin.get_l()

    if np.average(qlist) < 0:
        direction = "retro"
    else:
        direction = "pro"

    wavemode = mode_admin.get_wavemode(direction)
    print is_even, l, wavemode

    wavemode += "s"
    if wavemode[0] == "g":
        wavemode += "_list"
    wavemode = wavemode.replace(" ", "_")
    if wavemode[0] == "y" or wavemode[0] == "k":  # yanai and kelvin modes only have two arguments
        args = m, qlist
    else:
        args = m, k, qlist
    guesslist = getattr(asym, wavemode)(*args)

    qlist, found_lamlist = roots.multi_rootfind_fromguess_dimless(m, qlist, is_even, guesslist, ecc, dlngrav, verbose=False, inc=inc)
    if ecc == .05:
        ls="dotted"
    elif ecc == .1:
        ls="dashed"
    else:
        ls="solid"
    plt.plot(qlist, found_lamlist, ls=ls)


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

#    fullrange_multi_rootfind(m, qlists, kvals, aympcompare=True)  # Mostly for plotting functionality
#    fullrange_multi_rootfind(m, [qneg], [2], aympcompare=True)  # Testing just for k=-1, negative part

#    for ecc in [.0, .05, .1]:
#        dlngrav=partial(grav.chi_gravity_deriv, ecc*2)
#        rootfind_dimless(m, 2, np.linspace(-10., -1.25, 170), ecc=ecc, dlngrav=dlngrav, inc=1.0075)
#        rootfind_dimless(m, 1, np.linspace(-10., -1.25, 170), ecc=ecc, dlngrav=dlngrav, inc=1.05)
#        rootfind_dimless(m, 0, np.linspace(-10., -1.25, 170), ecc=ecc, dlngrav=dlngrav, inc=1.05)
#        rootfind_dimless(m, 2, np.linspace(10., 1.25, 170), ecc=ecc, dlngrav=dlngrav, inc=1.0075)
#        rootfind_dimless(m, 1, np.linspace(10., 1.25, 170), ecc=ecc, dlngrav=dlngrav, inc=1.05)
#        rootfind_dimless(m, 0, np.linspace(10., 1.25, 170), ecc=ecc, dlngrav=dlngrav, inc=1.05)
#        rootfind_dimless(m, -1, np.linspace(-10., -3.25, 100), ecc=ecc, dlngrav=dlngrav, inc=1.05)
##        rootfind_dimless(m, -2, np.linspace(-10., -6.1, 50), ecc=ecc, dlngrav=dlngrav, inc=1.05)
#    plt.yscale('log')
#    plt.show()
    compare_curvasym(m)
    compare_curvrmode(m)


    r_eq, mass, period = 1.2e4, 1.6*1.9885e30, 1./581
    x, om_bar_sq = oblate.find_x_ombarsq(r_eq, mass, period)
    r_polar = r_eq*oblate.calc_radius_14_dimless(x, om_bar_sq, 0.)
    eps = 1 - (r_polar / r_eq)
    ecc = 1 - (r_polar / r_eq)**2
    sigma = (r_polar / r_eq)**2
    Gamma = (2*om_bar_sq + 4*eps) / (1-om_bar_sq)
    print "sigma: {}, Gamma: {}, ecc: {}".format(sigma, Gamma, ecc)

#    plotting.asymptotic_plotting(m)
#    plotting.curvi_asymptotic_plotting(m, sigma, Gamma, ecc)
#    plt.show()
#    simple_numasyplot_rmodesINC(m, k, is_even)
#    r_eq, mass = 1e4, 1.4*1.9885e30
#    for period in [np.inf, 1./100, 1./363, 1./581]:
#        r_qlist = np.linspace(-100., -7.85, 450)  # r-modes are fine with fewer steps really
#        r_guesslist = asym.r_modes(m, k, r_qlist)
#        r_qlist, r_found_lamlist = roots.multi_rootfind_fromguess(m, r_qlist, is_even, r_guesslist, r_eq, mass, period, verbose=False, inc=1.15)
#        plt.plot(r_qlist, r_found_lamlist, ls="--")
#    plt.show()

#    error_numasyplot(m)
#    compare_newfile(m)

#    r_eq = 10000
#    mass = 1.4*1.9885e30
#    periodlist = [np.inf, 1./100, 1./363, 1./581]
#    multi_eccentricity_rootfind(m, k, qpos, is_even, r_eq, mass, periodlist)

#    rtol = 5e-4
#    period = 1./581
#    sn = "4U 1636-536"
##    period = 1./363
##    sn = "4U 1728-34"
##    period = 1./620
##    sn = "4U 1608-522"
##    period = 1./294
##    sn = "IGR J17191-2821"
#    om_bar_sq_values = [.03, .065, .1, .15, .20, .27]
#    plot_different_dimless(om_bar_sq_values, period=period, sn=sn, cmap="inferno", rtol=rtol)
#    compare_dimless(period=period, sn=sn)

##    oblate_plot(.1)
#    oblate_plot_alternative(.1)

#    for om_bar_sq in [.04, .10, .16, .22]:
#        oblate_plot(om_bar_sq)
#    for om_bar_sq in [.04, .10, .16, .22]:
#        oblate_plot_alternative(om_bar_sq)


if __name__ == "__main__":
    main()




# -*- coding: utf-8 -*-
#!/bin/bin/python
import sys
import numpy as np
import itertools
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D

import classes.Legendre as Legendre
import classes.LaPlace as LaPlace

import helpers.Property as Property
import helpers.Straddle as Straddle
import helpers.rootfinder as roots
import helpers.plotting as plotting
import helpers.LaPlace_asymptotes as asym
import helpers.sanity_plots as sanplot
import helpers.morsink_radius as oblate

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

#    simple_numasyplot(m)
#    error_numasyplot(m)
#    compare_newfile(m)


    rtol = 5e-4
    period = 1./581
    sn = "4U 1636-536"
#    period = 1./363
#    sn = "4U 1728-34"
#    period = 1./620
#    sn = "4U 1608-522"
#    period = 1./294
#    sn = "IGR J17191-2821"
    om_bar_sq_values = [.03, .065, .1, .15, .20, .27]
    plot_different_dimless(om_bar_sq_values, period=period, sn=sn, cmap="inferno", rtol=rtol)
    compare_dimless(period=period, sn=sn)

#    oblate_plot(.1)
    oblate_plot_alternative(.1)

#    for om_bar_sq in [.04, .10, .16, .22]:
#        oblate_plot(om_bar_sq)
#    for om_bar_sq in [.04, .10, .16, .22]:
#        oblate_plot_alternative(om_bar_sq)


if __name__ == "__main__":
    main()




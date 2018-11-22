import numpy as np
import matplotlib.pyplot as plt

def plot_function(func, params, values, condition, color, ls):
    plotrange = values[condition]
    params.append(plotrange)
    plotline = func(*params)
    plt.plot(plotrange, plotline, color=color, ls=ls)

def plot_from_file(file_loc, color, ls):
    splitstring = file_loc.split("_")
    data = np.loadtxt(file_loc)
    recon_qvals = np.linspace(float(splitstring[1]), float(splitstring[2]), float(splitstring[4]))
    plt.plot(recon_qvals, data, color=color, ls=ls)

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

    plot_function(eq38_fst, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-.")
    plot_function(eq38_snd, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-")
    plot_function(eq38_fst, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-.")
    plot_function(eq38_snd, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-")
    plot_function(eq38_fst, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-.")
    plot_function(eq38_snd, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-")
    plot_function(eq38_fst, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-.")
    plot_function(eq38_snd, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-")
    plot_function(eq40, [m], qpos, [np.abs(qpos)>1.], "green", "-")  # its effecitvely a gmode constraint
    plot_function(eq55, [m], qpos, [np.abs(qpos)>np.sqrt(3/(-m*qpos))], "cyan", "-")
    plt.yscale('log')


def numerics_plotting():
    """
    Plots a dotted line for the loaded numerics values

    TODO; pass on a folder or conditions on how to load multiple plot_from_files
    """
    plot_from_file("Numerics/Townsend2003/range_0.0_10.0_steps_952_kval_0.txt", "black", "--")
    plot_from_file("Numerics/Townsend2003/range_0.0_10.0_steps_952_kval_1.txt", "black", "--")
    plot_from_file("Numerics/Townsend2003/range_0.0_10.0_steps_952_kval_2.txt", "black", "--")
    plot_from_file("Numerics/Townsend2003/range_0.0_-10.0_steps_951_kval_0.txt", "black", "--")
    plot_from_file("Numerics/Townsend2003/range_0.0_-10.0_steps_951_kval_1.txt", "black", "--")
    plot_from_file("Numerics/Townsend2003/range_0.0_-10.0_steps_951_kval_2.txt", "black", "--")
    plt.yscale('log')

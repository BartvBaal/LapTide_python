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

def plot_func_file_diff(file_loc, func, params, values, condition, color, ls):
    splitstring = file_loc.split("_")
    data = np.loadtxt(file_loc)
    recon_qvals = np.linspace(float(splitstring[1]), float(splitstring[2]), float(splitstring[4]))
    plotrange = recon_qvals[np.abs(recon_qvals)>1.]
    data = data[np.abs(recon_qvals)>1.]
    params.append(plotrange)
    asyline = func(*params)
    plt.plot(plotrange, np.abs((data-asyline)/data), color=color, ls=ls)

def difference_plotting(m):
    eq38_fst = lambda m, s, q : (q**2) * ((2.*s + 1.)**2)
    eq38_snd = lambda m, s, q : (q**2) * ((2.*s + 1.)**2) - ( 2. * (m*q - m**2) )

    if m < 0:
        qneg = np.linspace(-10, 0, 10000)
        qpos = np.linspace(0, 10, 10000)
    else:
        qneg = np.linspace(0, 10, 10000)
        qpos = np.linspace(-10, 0, 10000)

    plot_func_file_diff("data/Townsend2003/range_0.0_-10.0_steps_951_kval_0.txt", eq38_fst, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-.")
    plot_func_file_diff("data/Townsend2003/range_0.0_-10.0_steps_951_kval_0.txt", eq38_snd, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-")
    plot_func_file_diff("data/Townsend2003/range_0.0_10.0_steps_952_kval_2.txt", eq38_fst, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-.")
    plot_func_file_diff("data/Townsend2003/range_0.0_10.0_steps_952_kval_2.txt", eq38_snd, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-")
    plot_func_file_diff("data/Townsend2003/range_0.0_-10.0_steps_951_kval_1.txt", eq38_fst, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-.")
    plot_func_file_diff("data/Townsend2003/range_0.0_-10.0_steps_951_kval_1.txt", eq38_snd, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-")
    plot_func_file_diff("data/Townsend2003/range_0.0_-10.0_steps_951_kval_2.txt", eq38_fst, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-.")
    plot_func_file_diff("data/Townsend2003/range_0.0_-10.0_steps_951_kval_2.txt", eq38_snd, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-")
    plt.yscale('log')


def asymptotic_plotting(m):
    """
    Plots the asymptotic functions to the Laplace Tidal Equations
    Solid lines for second order, dashdotted line for the first order
    """
    eq39_fst = lambda m, s, q : (m*q - m*m)**2 / (q*q * (2*s + 1)**2)
    eq39_snd = lambda m, s, q : (m*q - m*m)**2 / (q*q * (2*s + 1)**2) + (2. * ((m*q - m*m)**3)) / (q**4 * ((2.*s + 1)**4))
#    l_eq_40_41 = lambda m, q : (q - m)**2

    eq38_fst = lambda m, s, q : (q**2) * ((2.*s + 1.)**2)
    eq38_snd = lambda m, s, q : (q**2) * ((2.*s + 1.)**2) - ( 2. * (m*q - m**2) )

    eq40 = lambda m, q : (q - m)**2

    eq55 = lambda m, q : (m**2) * 2.*m*q / (2*m*q + 1.)

    ## q > condition statement (eg q>1 for g_mode)
    g_mode_cond = np.abs(1)
    r_mode_cond = lambda m, s : ((2.*s + 1.)**2 + m**2) / np.abs(m)  # paper mistake in 2s+1 term
    k_mode_cond = lambda m : np.abs(3. / m)
    y_mode_cond = lambda m : np.abs(m) + 1./np.abs(m)
    #y_mode_cond = lambda m, q : (m*q < m**2 and m*q > 0) ? r_mode_cond(m, 0) : g_mode_cond

    if m < 0:
        qneg = np.linspace(-10, 0, 10000)
        qpos = np.linspace(0, 10, 10000)
    else:
        qneg = np.linspace(0, 10, 10000)
        qpos = np.linspace(-10, 0, 10000)

    plot_function(eq38_fst, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-.")
    plot_function(eq38_snd, [m, 1], qneg, [np.abs(qneg)>1.], "blue", "-")
    plot_function(eq39_fst, [m, 1], qneg, [np.abs(qneg)>r_mode_cond(m, 1)], "brown", "-.")
    plot_function(eq39_snd, [m, 1], qneg, [np.abs(qneg)>r_mode_cond(m, 1)], "brown", "-")
    plot_function(eq38_fst, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-.")
    plot_function(eq38_snd, [m, 1], qpos, [np.abs(qpos)>1.], "purple", "-")
    plot_function(eq38_fst, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-.")
    plot_function(eq38_snd, [m, 2], qneg, [np.abs(qneg)>1.], "red", "-")
    plot_function(eq38_fst, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-.")
    plot_function(eq38_snd, [m, 3], qneg, [np.abs(qneg)>1.], "orange", "-")
    plot_function(eq40, [m], qpos, [np.abs(qpos)>1.], "green", "-")  # its effecitvely a gmode constraint
    plot_function(eq40, [m], qneg, [np.abs(qneg)>y_mode_cond(m)], "green", "-")
    plot_function(eq55, [m], qpos, [np.abs(qpos)>k_mode_cond(m)], "cyan", "-")
    plt.yscale('log')


def curvi_asymptotic_plotting(m, sigma, Gamma, ecc):
    eq39_fst = lambda m, s, sigma, Gamma, q : (m*q - (sigma**2) * m*m - Gamma)**2 / (q*q * (2*s + 1)**2)
    eq39_snd = lambda m, s, sigma, Gamma, ecc, q : (sigma*sigma/q/q) * ( (m*q/sigma/sigma - m*m + Gamma/sigma/sigma)**2 / ((2*s + 1)**2) - (ecc*ecc/(sigma**4)*m*q - m*q*Gamma/sigma/sigma - ecc*ecc*Gamma/(sigma**4)) )

#    l_eq_40_41 = lambda m, q : (q - m)**2

    eq38_fst = lambda m, s, sigma, q : (q**2) * ((2.*s + 1.)**2) / (sigma**4)
    eq38_snd = lambda m, s, sigma, Gamma, ecc, q : ((2*s + 1)**2)*q*q/(sigma**4) - 2* ( m*q/sigma/sigma - m*m + Gamma/sigma/sigma ) + ecc*ecc/sigma/sigma*m/q - m*Gamma/q - ecc*ecc/sigma/sigma*Gamma/q/q

    eq40 = lambda m, q : (q - m)**2

    eq55 = lambda m, q : (m**2) * 2.*m*q / (2*m*q + 1.)

    ## q > condition statement (eg q>1 for g_mode)
    g_mode_cond = np.abs(1)
    r_mode_cond = lambda m, s : ((2.*s + 1.)**2 + m**2) / np.abs(m)  # paper mistake in 2s+1 term
    k_mode_cond = lambda m : np.abs(3. / m)
    y_mode_cond = lambda m : np.abs(m) + 1./np.abs(m)
    #y_mode_cond = lambda m, q : (m*q < m**2 and m*q > 0) ? r_mode_cond(m, 0) : g_mode_cond

    if m < 0:
        qneg = np.linspace(-10, 0, 10000)
        qpos = np.linspace(0, 10, 10000)
    else:
        qneg = np.linspace(0, 10, 10000)
        qpos = np.linspace(-10, 0, 10000)

    plot_function(eq38_fst, [m, 1, sigma], qneg, [np.abs(qneg)>1.], "blue", "-.")
    plot_function(eq38_snd, [m, 1, sigma, Gamma, ecc], qneg, [np.abs(qneg)>1.], "blue", "-")
    plot_function(eq39_fst, [m, 1, sigma, Gamma], qneg, [np.abs(qneg)>r_mode_cond(m, 1)], "brown", "-.")
    plot_function(eq39_snd, [m, 1, sigma, Gamma, ecc], qneg, [np.abs(qneg)>r_mode_cond(m, 1)], "brown", "-")
    plot_function(eq38_fst, [m, 1, sigma], qpos, [np.abs(qpos)>1.], "purple", "-.")
    plot_function(eq38_snd, [m, 1, sigma, Gamma, ecc], qpos, [np.abs(qpos)>1.], "purple", "-")
    plot_function(eq38_fst, [m, 2, sigma], qneg, [np.abs(qneg)>1.], "red", "-.")
    plot_function(eq38_snd, [m, 2, sigma, Gamma, ecc], qneg, [np.abs(qneg)>1.], "red", "-")
    plot_function(eq38_fst, [m, 3, sigma], qneg, [np.abs(qneg)>1.], "orange", "-.")
    plot_function(eq38_snd, [m, 3, sigma, Gamma, ecc], qneg, [np.abs(qneg)>1.], "orange", "-")
    plot_function(eq40, [m], qpos, [np.abs(qpos)>1.], "green", "-")  # its effecitvely a gmode constraint
    plot_function(eq40, [m], qneg, [np.abs(qneg)>y_mode_cond(m)], "green", "-")
    plot_function(eq55, [m], qpos, [np.abs(qpos)>k_mode_cond(m)], "cyan", "-")
    plt.yscale('log')



def numerics_plotting(plotlist):
    """
    Plots a dotted line for the numerics values files in plotlist
    """
    for location in plotlist:
        plot_from_file(location, "black", "--")
    plt.yscale('log')


def townsend_plotting():
    plotlist = ["data/Townsend2003/range_0.0_10.0_steps_952_kval_0.txt", 
    "data/Townsend2003/range_0.0_10.0_steps_952_kval_1.txt", 
    "data/Townsend2003/range_0.0_10.0_steps_952_kval_2.txt", 
    "data/Townsend2003/range_0.0_-10.0_steps_951_kval_0.txt", 
    "data/Townsend2003/range_0.0_-10.0_steps_951_kval_1.txt", 
    "data/Townsend2003/range_0.0_-10.0_steps_951_kval_2.txt"]
    numerics_plotting(plotlist)




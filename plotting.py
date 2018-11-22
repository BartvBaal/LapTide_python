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


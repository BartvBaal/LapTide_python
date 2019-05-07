import numpy as np
from functools import partial
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors


import classes.Curvilinear as Curvilinear

import helpers.Curvilinear_asymptotes as curvasym
import helpers.Property as Property
import helpers.rootfinder as roots
import helpers.gravity_functions as grav


def create_grid(m, k, size, qmin, qmax, ecc, chi, dlngrav, saving=False):
    qarr = np.linspace(qmin, qmax, size)

    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()
    is_even = mode_admin.is_even()
    l = mode_admin.get_l()

    if np.average(qarr) < 0:
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
        args = m, qarr, ecc, chi
    else:
        args = m, k, qarr, ecc, chi
    guesslist = getattr(curvasym, wavemode)(*args)
    lamarr = np.linspace(min(guesslist)*.9, max(guesslist)*1.1, size)

    qdim = np.repeat(qarr, size)
    lamdim = np.tile(lamarr, size)

    values = []
    for q, lam in zip(qdim, lamdim):
        value = roots.shoot_curvi_dimless(lam, m, q, is_even, ecc, dlngrav)
        values.append(value)
        if lam in [lamdim[0]]:
            print q, lam, value
    heatmap = np.asarray(values).reshape(size, size).T

    all_data = np.asarray(zip(qdim, lamdim, values))
    if saving:
        folderloc = "data/gridplot/"
        if not os.path.exists(folderloc):
            os.makedirs(folderloc)
        savestring = folderloc+"qrange_{}_{}_steps_{}_kval_{}_ecc_{}_chi_{}.txt"\
                            .format(qarr[0],qarr[-1],size,str(k), str(ecc), str(chi))
        np.savetxt(savestring, all_data)
    return all_data


def lambda_grid(m, k, size, qmin, qmax, lammin, lammax, ecc, chi, dlngrav, saving=False):
    qarr = np.linspace(qmin, qmax, size)
    lamarr = np.linspace(lammin, lammax, size)
    print lamarr

    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()
    is_even = mode_admin.is_even()
    if is_even:
        evenstr = "Even"
    else:
        evenstr = "Odd"

    qdim = np.repeat(qarr, size)
    lamdim = np.tile(lamarr, size)

    values = []
    for q, lam in zip(qdim, lamdim):
        value = roots.shoot_curvi_dimless(lam, m, q, is_even, ecc, dlngrav)
        values.append(value)
        if lam in [lamdim[0]]:
            print q, lam, value
    heatmap = np.asarray(values).reshape(size, size).T

    all_data = np.asarray(zip(qdim, lamdim, values))
    if saving:
        folderloc = "data/lamgrid/"
        if not os.path.exists(folderloc):
            os.makedirs(folderloc)
        savestring = folderloc + "qrange_{}_{}_lamrange_{}_{}_steps_{}_ecc_{}_chi_{}_{}.txt"\
                    .format(qarr[0],qarr[-1],lamarr[0],lamarr[-1],size,str(ecc),str(chi),evenstr)
        np.savetxt(savestring, all_data)
    return all_data


def load_grid(filename):
    if os.path.isfile(filename):
        all_data = np.loadtxt(filename)
        return all_data
    else:
        raise ValueError("Given filename does not exist. Create grid first!")


def plot_grid(all_data, asymcompare=None, title=None, interpolation=None, logscale=False):
    size = int(np.sqrt(len(all_data)))
    qvals = all_data[:,0]
    lamvals = all_data[:,1]
    solns = all_data[:,2].reshape(size, size).T
    print solns

    fig, ax = plt.subplots()
    im = ax.imshow(solns, cmap="RdGy", extent=[qvals[0], qvals[-1], lamvals[0], lamvals[-1]], aspect="auto", origin="lower", norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5), interpolation=interpolation)
    ax.set_xlabel(r"Spin parameter q = 2$\Omega / \omega$")
    ax.set_ylabel(r"Eigenvalue $\lambda$")
    if logscale:
        ax.set_yscale('log')
    if title:
        ax.set_title(title)

    if asymcompare:
        ax.plot(*asymcompare, color="black", lw=4.5)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Value", rotation=-45, va="bottom")
    plt.show()


def plot_pcolormesh(all_data, asymcompare=None, title=None, logscale=True):
    """
    Alternative version from above, for if you have non-linear boxed data
    Cannot use interpolation for this, unfortunately :(
    """
    size = int(np.sqrt(len(all_data)))
    qvals = all_data[:,0][::size]
    lamvals = all_data[:,1][:size]
    solns = all_data[:,2].reshape(size, size).T
    solns = solns[:-1, :-1]
    print solns

    fig, ax = plt.subplots()
    im = ax.pcolormesh(qvals, lamvals, solns, cmap="RdGy", norm=colors.SymLogNorm(linthresh=0.01, linscale=0.5))
    ax.set_xlabel(r"Spin parameter q = 2$\Omega / \omega$")
    ax.set_ylabel(r"Eigenvalue $\lambda$")
    if logscale:
        ax.set_yscale('log')
    if title:
        ax.set_title(title)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Value", rotation=-45, va="bottom")
    plt.show()


#def main():
#    m, k = -2, 2
#    qmin, qmax, size = -3., -2., 8
#    ecc, chi = 0.0, 0.0
#    dlngrav = partial(grav.chi_gravity_deriv, chi)

#    qneg = np.linspace(qmin, qmax, 200)
#    asymcompare = [qneg, curvasym.g_modes_list(m, 2, qneg, ecc=ecc, chi=chi)]

#    all_data = load_grid("data/gridplot/qrange_-3.0_-2.0_steps_8_kval_2.txt")
##    all_data = create_grid(m, k, size, qmin, qmax, ecc, dlngrav, True)
#    plot_grid(all_data)  # , asymcompare=asymcompare


if __name__ == "__main__":
    main()

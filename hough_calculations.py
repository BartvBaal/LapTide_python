# -*- coding: utf-8 -*-
import numpy as np
import scipy.special as special
from functools import partial
import matplotlib.pyplot as plt

import classes.Curvilinear as Curvilinear

import helpers.rootfinder as roots
import helpers.Property as Property
import helpers.LaPlace_asymptotes as asym
import helpers.gravity_functions as grav

def rootfind_dimless(m, k, q, ecc=0, dlngrav=partial(grav.chi_gravity_deriv, 0.), verbose=False, inc=1.0033):
    """
    For a given m, k and qlist, and potentially for different radius, mass and
    period as well, determines the wave mode and calculates the eigenvalues
    from the asymptotically calculated values.
    """
    mode_admin = Property.Mode_admin(m, k)
    mode_admin.validate_values()
    is_even = mode_admin.is_even()
    l = mode_admin.get_l()

    if q*m < 0:
        direction = "pro"
    else:
        direction = "retro"

    wavemode = mode_admin.get_wavemode(direction)
    print is_even, l, wavemode, direction

    wavemode += "s"
    if wavemode[0] == "g":
        wavemode += "_list"
    wavemode = wavemode.replace(" ", "_")
    if wavemode[0] == "y" or wavemode[0] == "k":  # yanai and kelvin modes only have two arguments
        args = m, q
    else:
        args = m, k, q
    guess = getattr(asym, wavemode)(*args)

    qlist, found_lamlist = roots.multi_rootfind_fromguess_dimless(m, q, is_even, guess, ecc, dlngrav, verbose=False, inc=inc)
    return guess, found_lamlist, wavemode.split("_")[0], direction


def houghHat(s, sig):
    hermite = special.eval_hermite(s, sig)
    exp = np.exp(-(sig**2)/2.)
    return hermite*exp


def hough(s, sig, m, lam, q):
    L = lam**.5
    Lnu = L*q
    prefactor = np.sqrt(Lnu) / (lam - m**2)
    minhermite = s*((m/L) + 1) * special.eval_hermite(s-1, sig)
    maxhermite = .5*((m/L) - 1) * special.eval_hermite(s+1, sig)
    exp = np.exp(-(sig**2)/2.)
    hermites = minhermite+maxhermite
    return prefactor * hermites * exp


def houghTilde(s, sig, m, lam, q):
    L = lam**.5
    Lnu = L*q
    prefactor = m*np.sqrt(Lnu) / (m**2 - L**2)
    minhermite = s*((L/m) + 1) * special.eval_hermite(s-1, sig)
    maxhermite = .5*((L/m) - 1) * special.eval_hermite(s+1, sig)
    exp = np.exp(-(sig**2)/2.)
    hermites = minhermite+maxhermite
    return prefactor * hermites * exp


def kelvinhough(mu, m, q):
    min_mq_root = (-m*q)**.5
    tau = min_mq_root * mu
    return np.exp(-(tau**2)/2.)
    

def kelvinhoughHat(mu, m, q):
    min_mq_root = (-m*q)**.5
    tau = min_mq_root * mu
    fst = 1./min_mq_root
    snd = -(m*m) / (2*m*q + 1.)
    thr = tau * np.exp(-(tau**2)/2.)
    return fst * snd * thr


def kelvinhoughTilde(mu, m, q):
    min_mq_root = (-m*q)**.5
    tau = min_mq_root * mu
    fst = -m
    brckt = tau*tau / (2*m*q + 1) + 1
    exp = np.exp(-(tau**2)/2.)
    return fst * brckt * exp


def numerics(ode_class, ode_args, lam, is_even, N):
    ode_args.insert(2, is_even)
    ode_solver = ode_class.solver_t_dimless(*ode_args)

    steps = np.linspace(ode_class.t0, ode_class.t1, N)
    [P, Q] = ode_solver.interp(lam, steps)
    return P, Q


def normalize(array):
    return array/max(np.abs(array))


def townsendify(hough, houghhat, houghtilde, k, m):
    if k%2 == 0:
        divisor = hough[-1]
    else:
        divisor = -houghhat[-1]
    if m < 0:
        houghtilde *= -1
    return hough/divisor, houghhat/divisor, houghtilde/divisor


def make_houghs(m, k, s, q, ecc, chi):
    N = 125
    dlngrav = partial(grav.chi_gravity_deriv, chi)
    guess, lam, wavename, direc = rootfind_dimless(m, k, np.asarray([q]), ecc, dlngrav, verbose=False, inc=1.075)
    mu = np.linspace(1., 0., N)
    L = lam**.5
    Lnu = L * q  # it should just be q but that breaks if q is negative?
    sig = np.sqrt(Lnu) * mu
    guesssig = np.sqrt(guess**.5 * q) * mu

    print "Lnu: {}, guessLnu: {}".format(Lnu, (guess**.5) * q)
    print "k: {}, s: {}, q: {}".format(k, s, q)

    # Numeric values
    is_even = Curvilinear.check_is_even(m, k)
    num_hough, num_houghHat = numerics(Curvilinear, [m, q, ecc, dlngrav], lam, is_even, N)
    num_houghTilde = -m * num_hough - q*mu * num_houghHat

    ## Analytic values
    ana_houghHat = houghHat(s, guesssig)
    ana_hough = hough(s, guesssig, m, guess, q)
    ana_houghTilde = houghTilde(s, guesssig, m, guess, q)
    if wavename == "kelvin":
        ana_houghHat = kelvinhoughHat(mu, m, q)
        ana_hough = kelvinhough(mu, m, q)
        ana_houghTilde = kelvinhoughTilde(mu, m, q)

    num_hough, num_houghHat, num_houghTilde = townsendify(num_hough, num_houghHat, num_houghTilde, k, m)
    ana_hough, ana_houghHat, ana_houghTilde = townsendify(ana_hough, ana_houghHat, ana_houghTilde, k, m)

    fig = plt.figure()
    ax1 = fig.add_subplot(3, 1, 1)
    ax1.set_title("m: {}, k: {}, s: {}, q: {}  ({}grade {}); ecc: {}, $\chi$: {}".format(m, k, s, q, direc, wavename, ecc, chi))
    ax1.plot(mu, num_hough, ls="--", label="Numeric")
    ax1.plot(mu, ana_hough, label="Analytic")
    ax1.set_ylabel(r"$\Theta(\sigma)$")
    ax1.set_xlim([0, 1])
    ax2 = fig.add_subplot(3, 1, 2)
    ax2.plot(mu, num_houghHat, ls="--")
    ax2.plot(mu, ana_houghHat)
    ax2.set_ylabel(r"$\hat\Theta(\sigma)$")
    ax2.set_xlim([0, 1])
    ax3 = fig.add_subplot(3, 1, 3)
    ax3.plot(mu, num_houghTilde, ls="--")
    ax3.plot(mu, ana_houghTilde)
    ax3.set_ylabel(r"$\tilde\Theta(\sigma)$")
    ax3.set_xlabel(r"$\mu \equiv \cos(\theta)$")
    ax3.set_xlim([0, 1])

    ax1.tick_params(axis='y', which='both', left='on', right='on')
    ax2.tick_params(axis='y', which='both', left='on', right='on')
    ax3.tick_params(axis='y', which='both', left='on', right='on')
    ax1.legend()
    plt.show()


if __name__ == "__main__":
#    m, k, s, q = -2, 2, 1, 3  # Pro g mode
    m, k, s, q = 2, 2, 3, 2.3  # Retro g mode
#    m, k, s, q = 2, 1, 2, 3  # Retro g mode
#    m, k, s, q = 2, 0, 1, 3  # Retro g mode
#    m, k, s, q = -2, 1, 0, 3  # Pro Yanai
#    m, k, s, q = 2, -1, 0, 6  # Retro Yanai
#    m, k, s, q = 2, -2, 1, 12.5  # Retro r mode
#    m, k, s, q = -2, 0, -1, 3  # Kelvin check - this should use different functions!
#    m, k, s, q = -2, 0, -1, 10  # LeeSaio1997 check - these do not require new functions ?

    ecc = 0.
    chi = 0.

    make_houghs(m, k, s, q, ecc, chi)

    mlist = [-2, 2, -2, 2, -2, 2, 2, 2]
    klist = [2, 2, 1, 1, 0, 0, -1, -2]
    slist = [1, 3, 0, 2, -1, 1, 0, 1]
    qlist = [3, 3, 3, 3, 3, 3, 6, 15]
    ecc, chi = 0.15, 0.3

    for m, k, s, q in zip(mlist, klist, slist, qlist):
        make_houghs(m, k, s, q, ecc, chi)




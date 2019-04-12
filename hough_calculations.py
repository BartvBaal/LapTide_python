import numpy as np
import scipy.special as special

import matplotlib.pyplot as plt

import classes.Curvilinear as Curvilinear

import helpers.rootfinder as roots
import helpers.Property as Property
import helpers.LaPlace_asymptotes as asym

def rootfind_any(m, k, q, r_eq=1e4, mass=1.4*1.9885e30, period=np.inf, verbose=False, inc=1.0033):
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

    qlist, found_lamlist = roots.multi_rootfind_fromguess(m, q, is_even, guess, r_eq, mass, period, verbose=verbose, inc=inc)
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


def numerics(ode_class, ode_args, lam, is_even, N):
    ode_args.insert(2, is_even)
    ode_solver = ode_class.solver_t(*ode_args)

    steps = np.linspace(ode_class.t0, ode_class.t1, N)
    [P, Q] = ode_solver.interp(lam, steps)
    return P, Q


def normalize(array):
    return array/max(np.abs(array))


if __name__ == "__main__":
#    m, k, s, q = -2, 2, 1, 3  # Pro g mode
#    m, k, s, q = 2, 2, 3, 3  # Retro g mode
#    m, k, s, q = -2, 1, 0, 3  # Pro Yanai
#    m, k, s, q = 2, -1, 0, 6  # Retro Yanai
    m, k, s, q = 2, -2, 1, 6.005  # Retro r mode
#    m, k, s, q = -2, 0, -1, 3  # Kelvin check - this should use different functions!
#    m, k, s, q = -2, 0, -1, 10  # LeeSaio1997 check - these do not require new functions ?

    N = 125
    guess, lam, wavename, direc = rootfind_any(m, k, np.asarray([q]), verbose=False, inc=1.075)
    mu = np.linspace(1., 0., N)
    L = lam**.5
    Lnu = L * q  # it should just be q but that breaks if q is negative?
    sig = np.sqrt(Lnu) * mu
    guesssig = np.sqrt(guess**.5 * q) * mu

    print "Lnu: {}, guessLnu: {}".format(Lnu, (guess**.5) * q)
    print "k: {}, s: {}, q: {}".format(k, s, q)

    # Numeric values
    is_even = Curvilinear.check_is_even(m, k)
    num_hough, num_houghHat = numerics(Curvilinear, [m, q, 1e4, 1.4*1.9885e30, np.inf], lam, is_even, N)
    num_houghTilde = -m * num_hough - q*mu * num_houghHat

    ## Analytic values
    ana_houghHat = houghHat(s, guesssig)
    ana_hough = hough(s, guesssig, m, guess, q)
    ana_houghTilde = houghTilde(s, guesssig, m, guess, q)

    norm = True
    if norm:
        num_hough, num_houghHat, num_houghTilde = normalize(num_hough), normalize(num_houghHat), normalize(num_houghTilde)
        ana_hough, ana_houghHat, ana_houghTilde = normalize(ana_hough), normalize(ana_houghHat), normalize(ana_houghTilde)

    # This doesn't work for Kelvin modes I think
    if np.abs(m) == m or not is_even:
        ana_hough *= -1
        ana_houghHat *= -1
        ana_houghTilde *= -1
    
    plt.subplot(3, 1, 1)
    plt.title("m: {}, k: {}, s: {}, q: {}  ({}grade {})".format(m, k, s, q, direc, wavename))
    plt.plot(mu, num_hough, ls="--")
    plt.plot(mu, ana_hough)
    plt.ylabel(r"$\Theta(\sigma)$")
    plt.xlim([0, 1])
    plt.subplot(3, 1, 2)
    plt.plot(mu, num_houghHat, ls="--")
    plt.plot(mu, ana_houghHat)
    plt.ylabel(r"$\hat\Theta(\sigma)$")
    plt.xlim([0, 1])
    plt.subplot(3, 1, 3)
    plt.plot(mu, num_houghTilde, ls="--")
    plt.plot(mu, ana_houghTilde)
    plt.ylabel(r"$\tilde\Theta(\sigma)$")
    plt.xlabel(r"$\mu \equiv \cos(\theta)$")
    plt.xlim([0, 1])
    plt.show()


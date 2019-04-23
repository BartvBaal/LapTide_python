import numpy as np

import helpers.morsink_radius as oblate

def AM14_gravity_deriv(r_eq, mass, period, mu):
    """
    Returns the eccentricity and gamma as set by algendy&morsink2014
    With those two you can describe the Curvilinear ODE (and all other LaPlace terms)
    """
    x, om_bar_sq = oblate.find_x_ombarsq(r_eq, mass, period)
    factor = oblate.calc_radius_14_dimless(x, om_bar_sq, 0)
    ecc = 1-(factor**2)
    gamma = calc_gamma(om_bar_sq, ecc)

    top = 2.*gamma*mu
    bot = 2. + mu*mu*gamma
    return top/bot


def chi_gravity_deriv(chi, mu):
    """
    Using a more general form of gravity, g(mu) = g_E(1 + chi*mu^2)
    so dlngdmu = 2chimu / (1+chi*mu^2)
    """
    top = 2.*chi*mu
    bot = 1. + chi*mu*mu
    return top/bot


def calc_gamma(om_bar_sq, ecc):
    return (2*om_bar_sq + 4*ecc) / (1-om_bar_sq)

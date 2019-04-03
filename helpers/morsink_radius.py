# -*- coding: utf-8 -*-
import numpy as np

def calc_radius_07(r_eq, mass, period, angle):
    """
    Calculates the radius at the angle, for the given parameters
    r_eq should be given in meters, and is the equatorial radius
    mass should be given in kilograms, and is the stellar mass
    period should be given in seconds (or 1/freq)
    angle should be given in radians, is the co-latitudinal angle!
    """
    GM = 6.67408e-11 * mass
    epsi = 4*np.pi*np.pi * r_eq**3 / (period*period * GM)
    zeta = GM / (r_eq * 299792458*299792458)
    p0 = 1
    p2 = .5*(3*np.cos(angle)**2 - 1)
    p4 = .125*(35*np.cos(angle)**4 - 30*np.cos(angle)**2 + 3)
    a0 = -.18*epsi + .23*zeta*epsi - .05*epsi**2
    a2 = -.39*epsi + .29*zeta*epsi - .13*epsi**2
    a4 = 0.04*epsi - .15*zeta*epsi + .07*epsi**2
    factor = 1 + a0*p0 + a2*p2 + a4*p4
    return factor


def calc_radius_14(r_eq, mass, period, angle):
    """
    Uses new functions to calculate ellipticity.
    Parameters the same as for calc_radius_07
    """
    GM = 6.67408e-11 * mass
    x = GM / (r_eq * 299792458*299792458) 
    om_bar_sq = 4*np.pi*np.pi * r_eq**3 / (period*period * GM)
    o20 = -.788
    o21 = 1.03 * x
    factor = 1 + om_bar_sq * (o20 + o21) * np.cos(angle)**2
    return factor


def calc_grav_14_slow(r_eq, mass, period, angle):
    """
    Uses equation 49 from AlGendy_Morsink(2014) to calculate the effective
    gravity at different angles for the "slow rotation" limit
    """
    GM = 6.67408e-11 * mass
    x = GM / (r_eq * 299792458*299792458) 
    om_bar_sq = 4*np.pi*np.pi * r_eq**3 / (period*period * GM)
    ce0 = -.791
    ce1 = .776 * x
    cp0 = 1.138
    cp1 = -1.431 * x
    factor = 1 + om_bar_sq * (ce0 + ce1) * np.sin(angle)**2 \
                + om_bar_sq * (cp0 + cp1) * np.cos(angle)**2
    return factor


def find_x_ombarsq(r_eq, mass, period):
    """
    Returns x and omega_bar_squared for the given radius, mass and period
    x = (GM/(Rc^2)) and om_bar_sq = (om_sqR^3/(GM))
    """
    GM = 6.67408e-11 * mass
    x = GM / (r_eq * 299792458*299792458) 
    om_bar_sq = 4*np.pi*np.pi * r_eq**3 / (period*period * GM)
    return x, om_bar_sq
    

def calc_radius_14_dimless(x, om_bar_sq, angle):
    """
    Uses the dimensionless parameters x (GM/(Rc^2)) and om_bar_sq
    (om_sqR^3/(GM)) to calculate the radius at the given angle
    """
    o20 =  -.788
    o21 = 1.03 * x
    factor = 1 + om_bar_sq * (o20 + o21) * np.cos(angle)**2
    return factor


def calc_grav_14_slow_dimless(x, om_bar_sq, angle):
    """
    Uses the dimensionless parameters x (GM/(Rc^2)) and om_bar_sq
    (om_sqR^3/(GM)) to calculate the effective gravity at the given angle
    """
    ce0 = -.791
    ce1 = .776 * x
    cp0 = 1.138
    cp1 = -1.431 * x
    factor = 1 + om_bar_sq * (ce0 + ce1) * np.sin(angle)**2 \
                + om_bar_sq * (cp0 + cp1) * np.cos(angle)**2
    return factor


def recover_radius_mass(r_list, m_list, period, om_bar_sq_target, rtol=5e-4, atol=0):
    """
    Much faster version of recover_radius_mass_slow(*params)

    For all combinations of values in r_list and m_list, this function will
    calculate x and om_bar_sq at the given period. It will return all of the
    combinations that are close to om_bar_sq_target, with rtol deciding how
    accurate the search has to be.
    
    Returns the result as a numpy array of tuples (radius, mass)
    """
    rlen = len(r_list)
    mlen = len(m_list)

    newr = np.repeat(r_list, mlen)
    newm = np.tile(m_list, rlen)
    x, om_bar_sq = find_x_ombarsq(newr, newm, period)

    result = []
    loc = np.nonzero(np.isclose(om_bar_sq, om_bar_sq_target, rtol=rtol, atol=atol))[0]

    for data in zip(newr[loc], newm[loc], x[loc], om_bar_sq[loc]):
        r, m, x, om_bar_sq = data
        result.append((r, m))
        print r"x: {:.4f}, Omega_bar²: {:.4f} from r: {:.2f} km, mass: {:.5f} Msol".format(x, om_bar_sq, r*1e-3, m/1.9885e30)
    return np.asarray(result)


if __name__ == "__main__":
    r_eq = 12000  # in meters
    mass = 1.8*1.9885e30  # in kg
    period = 1./581
    degrees = [90, 75, 62.5, 52.5, 45, 37.5, 27.5, 15, 0]
    print "Equatorial radius: {} km, mass: {} Msol, spin frequency: {:.1f} Hz\n".format(r_eq*1e-3, mass/1.9885e30, 1/period)
    angles = np.radians(degrees)
    polar = calc_radius_07(r_eq, mass, period, angles)
    data = zip(r_eq*polar*1e-3, 1-polar, degrees)
    print "Morsink et al 2007 calculations:"
    for elem in data:
        print u"Radius of the star: {:.2f} km,\tdiff: {:.5f},\tco-latitudinal angle: {}\u00b0".format(*elem)
    print "Ellipciticy: {}, e²: {}".format(1 - polar[-1]/polar[0], 1 - (polar[-1]/polar[0])**2)

    print "\n ======================================== \n"
    newpolar = calc_radius_14(r_eq, mass, period, angles)
    grav = calc_grav_14_slow(r_eq, mass, period, angles)
    newdata = zip(r_eq*newpolar*1e-3, 1-newpolar, grav, degrees)
    print "AlGendy & Morsink 2014 calculations (including slow rotation GR-gravity):"
    for elem in newdata:
        print u"Radius of the star: {:.2f} km,\tdiff: {:.5f},\teff_gravity: {:.3f},\tco-latitudinal angle: {}\u00b0".format(*elem)
    print "Ellipciticy: {}, e²: {}, relative gravity change: {}".format(1 - newpolar[-1]/newpolar[0], 1 - (newpolar[-1]/newpolar[0])**2, grav[-1]/grav[0])


#    om_bar_sq = find_x_ombarsq(r_eq, mass, period)[1]
#    eps = 1 - newpolar[-1]/newpolar[0]
#    gamma = (2*om_bar_sq + 4*eps) / (1 - om_bar_sq)
#    c = newpolar[-1]
#    a = newpolar[0]
#    lhs = gamma*c/a
#    rhs = (a**2 - c**2)/c**2
#    print lhs, rhs, lhs/rhs 




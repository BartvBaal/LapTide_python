import numpy as np

def radius(r_eq, mass, period, angle):
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

polar = radius(18000, 2.5*1.989e30, 1./600, np.pi/2)
print polar, 1-polar

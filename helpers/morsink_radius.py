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


def test_calc(r_eq, mass, freq, angle):
    """
    with G=c=1 so units don't match => values are garbage
    """
    x = mass / (r_eq)
    angvel = (2*np.pi*r_eq * freq)**2
    om_bar_sq = angvel * r_eq**3 / (mass)
    o20 = -.788
    o21 = 1.03 * x
    factor = 1 + om_bar_sq * (o20 + o21) * np.cos(angle)**2
    print "values", angvel, om_bar_sq, x
    return factor

if __name__ == "__main__":
    r_eq = 12000  # in meters
    mass = 2.2*1.989e30  # in kg
    period = 1./600
    degrees = [90, 75, 62.5, 52.5, 45, 37.5, 27.5, 15, 0]
    angles = np.radians(degrees)
    polar = calc_radius_07(r_eq, mass, period, angles)
    data = zip(r_eq*polar*1e-3, 1-polar, degrees)
    print "Morsink et al 2007 calculations:"
    for elem in data:
        print u"Radius of the star: {:.2f} km,\tdiff: {:.5f},\tco-latitudinal angle: {}\u00b0".format(*elem)
    print "Ellipciticy: {}, e: {}".format(1 - polar[-1]/polar[0], 1 - (polar[-1]/polar[0])**2)

    print "\n ======================================== \n"
    newpolar = calc_radius_14(r_eq, mass, period, angles)
    newdata = zip(r_eq*newpolar*1e-3, 1-newpolar, degrees)
    print "AlGendy & Morsink 2014 calculations:"
    for elem in newdata:
        print u"Radius of the star: {:.2f} km,\tdiff: {:.5f},\tco-latitudinal angle: {}\u00b0".format(*elem)
    print "Ellipciticy: {}, e: {}".format(1 - newpolar[-1]/newpolar[0], 1 - (newpolar[-1]/newpolar[0])**2)

#    print "\n ======================================== \n"
#    newr = 10.
#    newpolar = test_calc(newr, 1.4, 600., angles)
#    newdata = zip(newr*newpolar, 1-newpolar, degrees)
#    for elem in newdata:
#        print u"Radius of the star: {:.2f} km,\tdiff: {:.5f},\tco-latitudinal angle: {}\u00b0".format(*elem)
#    print "Ellipciticy: {}".format(1 - newpolar[-1]/newpolar[0])




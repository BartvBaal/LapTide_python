def r_modes(m, k, q):
    """
    First validates if the given parameters can create an r-mode wave.
    Returns the second order approximation.
    """
#    if type(m) != int or type(k) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")
    if m*q <= m*m:
        raise ValueError("Invalid parameters; m*q must be greater than m*m")
    if q == 0:
        raise ValueError("Invalid parameter; q has to be non-zero")
    if k > m:  # Not sure on this constraint yet
        raise ValueError("Invalid parameter; k has to be lower than m")
    n1 = (m*q - m*m)**2.
    d1 = q*q * ((2.*k - 1)**2)
    n2 = 2. * ((m*q - m*m)**3)
    d2 = q**4 * ((2.*k - 1)**4)
    approx = (n1 / d1) + (n2 / d2)
    return approx


def g_modes(m, k, q):
    """
    Figures out which k should be used, since it differs for prograde and
    retrograde modes. Not an optimal way of doing this yet.
    Returns the second order approximation.
    """
#    if type(m) != int or type(k) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")
    if q[-1] < 0:
        s = k + 1.
        print "retrograde"
    else:
        s = k - 1.
        print "prograde"
    fo = q*q * ((2.*s + 1)**2)
    so = -2. * (m*q - m*m)
    approx = fo + so
    return approx


def yanai_modes(m, q):
    """
    Returns the approximation for the Yanai modes.
    """
#    if type(m) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")
    approx = (q - m)**2
    return approx


def kelvin_modes(m, q):
    """
    Returns the approximation for the Kelvin modes.
    """
#    if type(m) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")
    approx = (m**2) * 2.*m*q / (2*m*q + 1.)
    return approx

import numpy as np

def r_modes(m, k, q, ecc=0, chi=0):
    """
    First validates if the given parameters can create an r-mode wave.
    Returns the second order approximation.

    TODO; figure out a better pattern to switch between single entry q and
    numpy arrays of q values
    """
    if type(q) == int:
        if m*q <= m*m:
            raise ValueError("Invalid parameters; m*q must be greater than m*m")
        if q == 0:
            raise ValueError("Invalid parameter; q has to be non-zero")
    elif type(q) == np.ndarray:
        if m*q[m*q <= m*m]:
            raise ValueError("Invalid parameters; all m*q must be greater than m*m")
        if q[q == 0.]:
            raise ValueError("Invalid parameter; all q have to be non-zero")
    if type(m) != int or type(k) != int:
        raise TypeError("Input parameters m and k should be given as integers")
    if k > m:  # Not sure on this constraint yet
        raise ValueError("Invalid parameter; k has to be lower than m")
    sigma = np.sqrt(1.-ecc)
    sigsq = sigma**2
    s = -k-1
    eccdivsigsq = ecc / (sigsq)
    Gamma = 2.*chi
    Upsilon = eccdivsigsq*m*q - m*q*Gamma - eccdivsigsq*Gamma

    n1 = (m*q - m*m*sigsq - Gamma)**2.
    d1 = q*q * ((2.*s + 1)**2)
    n2 = Upsilon / (q*q)
    n3 = 2*n2 * sigsq * (m*q - m*m*sigsq + Gamma)
    approx = (n1 / d1) - (n2) + (n3 / d1)
    return approx


def g_modes(m, k, q, ecc=0, chi=0):
    """
    Figures out which k should be used, since it differs for prograde and
    retrograde modes.
    Returns the second order approximation.
    """
#    if type(m) != int or type(k) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")

    # Figure out pro/retro; q=0 case does not matter since the s terms drop
    if q*m < 0:  # prograde
        s = k - 1.
    else:  # retrograde
        s = k + 1.
    sigma = np.sqrt(1.-ecc)
    sigsq = sigma**2
    eccdivsigsq = ecc / (sigsq)
    Gamma = 2.*chi

    fo = q*q * ((2.*s + 1)**2) / (sigsq**2)
    so = -2. * (m*q/sigsq - m*m - Gamma/sigsq)
    to = eccdivsigsq*m/q - Gamma*m/q - Gamma*eccdivsigsq/(q*q)
    approx = fo + so + to
    return approx


def g_modes_list(m, k, qlist, ecc=0, chi=0):
    """
    Will call g_modes(m, k, q) for all values of q in qlist.
    """
    approx = []
    for q in qlist:
        approx.append(g_modes(m, k, q, ecc=ecc, chi=chi))
    return np.asarray(approx)


def yanai_modes(m, q, ecc=0, chi=0):
    """
    Returns the approximation for the Yanai modes.
    """
#    if type(m) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")
    approx = (q - m)**2
    return approx


def kelvin_modes(m, q, ecc=0, chi=0):
    """
    Returns the approximation for the Kelvin modes.
    """
#    if type(m) != int or type(q) != int:
#        raise TypeError("Input parameters should be given as integers")
    approx = (m**2) * 2.*m*q / (2*m*q + 1.)
    return approx

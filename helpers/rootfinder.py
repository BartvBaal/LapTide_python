import numpy as np
from scipy import optimize

import classes.LaPlace as LaPlace
import helpers.Straddle as Straddle


def shoot_laplace(lam, m, q, is_even):
    laplace_solver = LaPlace.solver_t(m, q, is_even)
    return laplace_solver(lam)


def rootfinder_laplace(m, q, lamlist, is_even):
    lamlow, lamhigh = lamlist
    return optimize.bisect(shoot_laplace, lamlow, lamhigh, args=(m, q, is_even), full_output=True)


def straddlefinder(fn, x0, neg_allowed=False):
    """
    For a given function fn and starting guess x0, it will look for a straddling
    point through linear steps, always in the positive direction, and also in
    the negative direction if the neg_allowed argument is True.

    Considering to change to logarithmic stepping to increase speed
    """
    inc = 1.005  # qneg, k=2 will break if inc is too large!
    N_steps = 100
    straddle = Straddle.straddle_t(fn, inc, N_steps)
    bisec = straddle.search_log(x0, neg_allowed=neg_allowed)

    return bisec


def multi_rootfind(m, l, qlist, is_even):
    """
    does rootfinder_laplace for all values in q for a set m and l (and thus is_even)
    Takes a starting value of "root" just below the Legendre solution
    Then for every q in qlist, the function will look for the first straddle it
    can find, and pass the straddling list on to the rootfinder_laplace() function
    the straddlefinder will take the previous root as the next guess
    
    Currently testing if this works for all previously known cases!
    Had to work in the neg_allowed searching for the k=2, qneg cases
    """
    found_lamlist = []
    root = l*(l+.99)  # Initial guessing point, called "root" since that's what its called later
    neg_allowed = False
    if qlist[-1] > 0:
        neg_allowed = True

    for q in qlist:
        fn = LaPlace.solver_t(m, q, is_even)
        lamlist = straddlefinder(fn, root, neg_allowed)

        root = rootfinder_laplace(m, q, lamlist, is_even)[0]

        found_lamlist.append(root)
        print q, root, lamlist

    found_lamlist = np.asarray(found_lamlist, dtype=float)
    return qlist, found_lamlist

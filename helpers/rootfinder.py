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


def straddlefinder(fn, x0, verbose=False, neg_allowed=False):
    """
    For a given function fn and starting guess x0, it will look for a straddling
    point with log steps, always in the positive direction, and also in
    the negative direction if the neg_allowed argument is True.

    Returns the points which straddle the root
    """
    ## qneg, k=3 still jumps at this inc rate (so need to make it smaller still? That slows down the code a lot :( )
    inc = 1.0033  # qneg, k=2 will break if inc is too large!
    N_steps = 100
    straddle = Straddle.straddle_t(fn, inc, N_steps)
    bisec = straddle.search_log(x0, verbose=verbose, neg_allowed=neg_allowed)

    return bisec


def multi_rootfind(m, k, qlist, is_even):
    """
    Does rootfinder_laplace for all values in q for a set m and k (and thus is_even)
    Currently calculates l from k and m; should pass in Property.Mode_admin()
    since that will already come with a stored m, k and l.
    Takes a starting value of "root" just below the Legendre solution
    Then for every q in qlist, the function will look for the first straddle it
    can find, and pass the straddling list on to the rootfinder_laplace() function
    the straddlefinder will take the previous root as the next guess
    
    Currently testing if this works for all previously known cases!
    Had to work in the neg_allowed searching for the k=2, qneg cases

    Need to work the neg_allowed into the admin_mode settings because it's really not
    that easy to determine when I should allow the code to look for decreasing solutions
    (eg for m=2, qnegs are allowed to go negative?
    NOTE; I **THINK** that m*q positive is only allowed postive searched)

    Do note; the k=2 cases (m=-2, negative q & m=2, positive q) are still bumpy around |2.2|
    for the current qlists - 9.5e2+1 but it breaks for 9.5e+2 => its still magic numbery
    """
    l = k + np.abs(m)
    found_lamlist = []
    root = l*(l+.99)  # Initial guessing point, called "root" since that's what its called later
    neg_allowed = False
    verbose=False
    if qlist[-1]*m < 0:
        print "\nAllowing negatives\n"
        neg_allowed = True

    for q in qlist:
        fn = LaPlace.solver_t(m, q, is_even)
        lamlist = straddlefinder(fn, root, verbose, neg_allowed)

        root = rootfinder_laplace(m, q, lamlist, is_even)[0]

        found_lamlist.append(root)
        print q, root, lamlist

    found_lamlist = np.asarray(found_lamlist, dtype=float)
    return qlist, found_lamlist

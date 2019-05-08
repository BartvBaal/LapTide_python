import numpy as np
from scipy import optimize

import classes.LaPlace as LaPlace
import classes.Curvilinear as Curvilinear
import helpers.Straddle as Straddle
import helpers.Property as Property
import helpers.Curvilinear_asymptotes as curvasym


def shoot_laplace(lam, m, q, is_even):
    laplace_solver = LaPlace.solver_t(m, q, is_even)
    return laplace_solver(lam)


def rootfinder_laplace(m, q, lamlist, is_even):
    lamlow, lamhigh = lamlist
    return optimize.bisect(shoot_laplace, lamlow, lamhigh, args=(m, q, is_even), full_output=True)


def shoot_curvilinear(lam, m, q, is_even, r_eq, mass, period):
    curvi_solver = Curvilinear.solver_t(m, q, is_even, r_eq, mass, period)
    return curvi_solver(lam)


def rootfinder_curvilinear(m, q, lamlist, is_even, r_eq, mass, period):
    return optimize.brentq(shoot_curvilinear, lamlist[0], lamlist[1], args=(m, q, is_even, r_eq, mass, period), full_output=True)


def shoot_curvi_dimless(lam, mode_admin, q):
    curvi_solver = Curvilinear.solver_t_dimless(mode_admin, q)  # mode_admin.mode should be the strings of the wave, not the strings of the function calls
    return curvi_solver(lam)


def rootfinder_curvi_dimless(mode_admin, q, lamlist):
    return optimize.brentq(shoot_curvi_dimless, lamlist[0], lamlist[1], args=(mode_admin, q), full_output=True)



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
    verbose = False
    if qlist[-1]*m < 0:
        print "\nAllowing negatives!\n"
        neg_allowed = True

    for q in qlist:
        fn = LaPlace.solver_t(m, q, is_even)
        lamlist = straddlefinder(fn, root, verbose, neg_allowed)

        root = rootfinder_laplace(m, q, lamlist, is_even)[0]

        found_lamlist.append(root)
        print q, root, lamlist

    found_lamlist = np.asarray(found_lamlist, dtype=float)
    return qlist, found_lamlist


def multi_rootfind_curvilinear_new(m, qlist, is_even, init_guess, r_eq, mass, period, verbose=False, inc=1.0033):
    """
    New way of calculating all the eigenvalues for q's in qlist
    Not yet faster than the old method, unfortunately. But at least it's not slower
    """
    root = init_guess
    neg_allowed = True  # init_guess might be too high or low so should allow pos&neg searching
    inc = inc  # qneg, k=2 will break if inc is too large!
    N_steps = 100
#    if qlist[-1]*m < 0:
#        print "\nAllowing negatives!\n"
#        neg_allowed = True

    found_lamlist = np.arange(float(len(qlist)))
    fn = Curvilinear.solver_t(m, qlist[0], is_even, r_eq, mass, period)
    straddle = Straddle.straddle_t(fn, inc, N_steps)
    qlam = zip(qlist, found_lamlist)

    for q, lam in qlam:
        fn.set_q(q)
        bisec = straddle.search_log(root, verbose=verbose, neg_allowed=neg_allowed)
        root = rootfinder_curvilinear(m, q, bisec, is_even, r_eq, mass, period)[0]

#        np.put(found_lamlist, lam, root)
        found_lamlist[lam] = root
        print q, root, bisec

    return qlist, found_lamlist


def multi_rootfind_fromguess(m, qlist, is_even, guesslist, r_eq, mass, period, verbose=False, inc=1.0033):
    neg_allowed = True
    inc = inc
    N_steps = 100
    found_lamlist = np.zeros(len(qlist))
    fn = Curvilinear.solver_t(m, qlist[0], is_even, r_eq, mass, period)
    straddle = Straddle.straddle_t(fn, inc, N_steps)

    for i in range(len(qlist)):
        q = qlist[i]
        fn.set_q(q)
        bisec = straddle.search_log(guesslist[i], verbose=verbose, neg_allowed=neg_allowed)
        root = rootfinder_curvilinear(m, q, bisec, is_even, r_eq, mass, period)[0]

        found_lamlist[i] = root
        print q, root, bisec

    return qlist, found_lamlist


def multi_rootfind_curvilinear(m, k, qlist, is_even, r_eq, mass, period):
    """
    multi_rootfind but for the curvilinear function instead. Should reduce to
    the multi_rootfind (laplace version) in the no-spin limit (so period=np.inf)

    should check this over with Frank, since the function 'fn' never changes
    and there is lots of room for optimization probably, which would be very
    helpful
    """
    l = k + np.abs(m)
    found_lamlist = []
    root = l*(l+.99)
    neg_allowed = False
    verbose = False
    if qlist[-1]*m < 0:
        print "\nAllowing negatives!\n"
        neg_allowed = True

    for q in qlist:
        fn = Curvilinear.solver_t(m, q, is_even, r_eq, mass, period)
        lamlist = straddlefinder(fn, root, verbose, neg_allowed)

        root = rootfinder_curvilinear(m, q, lamlist, is_even, r_eq, mass, period)[0]

        found_lamlist.append(root)
        print q, root, lamlist

    found_lamlist = np.asarray(found_lamlist, dtype=float)
    return qlist, found_lamlist


#def multi_rootfind_fromguess_dimless(m, qlist, is_even, guesslist, ecc, dlngrav, verbose=False, inc=1.0033):
def multi_rootfind_fromguess_dimless(mode_admin, verbose=False, inc=1.0033):
#    m = mode_admin.m
#    is_even = mode_admin.is_even
#    qlist = mode_admin.qlist
#    ecc = mode_admin.ecc
#    dlngrav = mode_admin.dlngrav
#    mode = mode_admin.mode

    direction = mode_admin.get_direction()
    wavemode = mode_admin.get_wavemode(direction)
    print mode_admin.is_even, mode_admin.l, wavemode

    wavemode += "s"
    if wavemode[0] == "g":
        wavemode += "_list"
    wavemode = wavemode.replace(" ", "_")
    if wavemode[0] == "y" or wavemode[0] == "k":  # yanai and kelvin modes only have two arguments
        args = mode_admin.m, mode_admin.qlist, mode_admin.ecc, mode_admin.chi
    else:
        args = mode_admin.m, mode_admin.k, mode_admin.qlist, mode_admin.ecc, mode_admin.chi
    guesslist = getattr(curvasym, wavemode)(*args)
    mode_admin.guess = guesslist

    neg_allowed = True
    inc = inc
    N_steps = 250
    found_lamlist = np.zeros(len(mode_admin.qlist))
    fn = Curvilinear.solver_t_dimless(mode_admin, mode_admin.qlist[0])
    straddle = Straddle.straddle_t(fn, inc, N_steps)

    for i in range(len(mode_admin.qlist)):
        q = mode_admin.qlist[i]
        fn.set_q(q)
        bisec = straddle.search_log(guesslist[i], verbose=verbose, neg_allowed=neg_allowed)
        root = rootfinder_curvi_dimless(mode_admin, q, bisec)[0]

        found_lamlist[i] = root
        print q, root, bisec

    return mode_admin.qlist, found_lamlist


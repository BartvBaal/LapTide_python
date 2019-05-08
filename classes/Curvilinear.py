"""Solve ODE for curvilinear, oblate variant of LaPlace Tidal Equations:

  (1 - x^2) \diffrac{f}{x} = (2*alpha*x + m*q*x)*f + sigma*((q^2 * x^2)/sigma^2 - 1)*g
  (1 - x^2) \diffrac{g}{x} = (2*alpha*x - m*q*x)*g + sigma*(lambda*(1-x^2) - m^2)*f + (1-x^2)*\diffrac{\ln{gravity}}*f

with

  P = (1 - x^2) f
  Q = (1 - x^2) g


to avoid singularity and implement boundary condition at \pm 1"""

import numpy as np

from scipy.integrate import ode
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from functools import partial

import Observers
import helpers.morsink_radius as oblate

## -------------------------------------------------------------------------- ##

eps = 1.0e-6
t0 = 1. - eps
t1 = 0.

def norm(arr):
    fmax = np.amax(arr)
    fmin = np.amin(arr)
    if abs(fmax) > abs(fmin):
        return fmax
    else:
        return fmin

def check_is_even(m, l):
    if (m+l) % 2 == 0:
        return True
    else:
        return False

def is_even_str(is_even):
    if is_even is True:
        return "even"
    else:
        return "odd"

def eigval(l):
    return l*(l+1.0)


## -------------------------------------------------------------------------- ##

class ODE_t:
    """Return RHS of ODE, and helpful functions for transformation"""
    # TODO: add in the parameters required for the new ODE, which are:
        # sigma (-> eccentricity; e^2 - x is called "t" and thus present already)
        # dxln(grav) (= 2\Gamma x / (2 + x^2*\Gamma))
            # \Gamma \equiv (2\bar\Omega^2 + 4\epsilon) / (1 - \bar\Omega^2)
        # With eccentricity, epsilon and om_bar_sq the parameters above are set
            # We need R_{eq} and R_{polar} for ecc/eps
            # We need R_{eq}, \Omega *and* M for om_bar_sq&R_{polar}
            # Fortunately we can choose these at the start of the simulation
            # and they will remain constants, but the initial conditions of the
            # star will now matter for the outcome
    # TODO: decide if I will input R_eq, R_polar, M, Omega *or* sigma, Gamma *or* x, om_bar_sq
    ## NOTE: mathematical approach with ecc, dlngrav as input is in comment at EoF!
    def __init__(self, m, q, lam, r_eq, mass, period):
        self.m = 1. * m
        self.q = 1. * q
        self.msq = m*m
        self.qsq = q*q

        self.r_eq = r_eq
        self.mass = mass
        self.period = period
        self.r_polar = oblate.calc_radius_14(r_eq, mass, period, 0) * r_eq

        self.x, self.om_bar_sq = oblate.find_x_ombarsq(r_eq, mass, period)
        self.ecc = 1 - (self.r_polar / r_eq)**2  # Eccentricity, e**2
        self.ell = 1 - (self.r_polar / r_eq)  # Ellipticity, epsilon

        self.alpha = .5 * abs(self.m)
        self.lam = 1. * lam

    def init_y(self):
        y0 = 1.0e-4
        sig = np.sqrt(1 - self.ecc * (1 - y0**2))
        y1 = (2 * self.alpha + self.m * self.q)*t0*y0 / (sig*(1 - t0*t0*self.qsq/sig/sig))
        return [y0, y1]  # Starting point, variation of RHS of eq (10)

    def coeffs(self, t):
        ## TODO: update as this is now Curvilinear!, not Legendre!
        sinsq = 1. - t*t
        twoax = 2. * self.alpha * t
        return np.array([twoax / sinsq, -1. / sinsq, self.lam - self.msq / sinsq, twoax / sinsq])

    def __call__(self, t, y):
    # Comments show what variables are called in legendre-ode_derivation.pdf
        sinsq = 1. - t*t  # 1-(x**2) (\equiv 1-\mu^2 \equiv sin^2 => name)
        twoax = 2. * self.alpha * t  # 2*alpha*x
        mqx = self.m * self.q * t  # m*q*x
        sig = np.sqrt(1 - self.ecc * sinsq)
        gam = (2*self.om_bar_sq + 4*self.ell) / (1-self.om_bar_sq)
        qxsqmo_sigcor = ((self.qsq * t*t / (sig*sig)) - 1) *sig  # [(x^2*q^2/sig^2)-1]*sig
        dlng = 2*gam*t / (2 + t*t*gam)
        dy0dt = ( (twoax + mqx)*y[0] + qxsqmo_sigcor*y[1] ) / sinsq  # eq (8), rewritten
#        dy1dt = ( (self.lam*sinsq - self.msq)*y[0] + (twoax - mqx)*y[1] ) / sinsq
        dy1dt = (dlng + sig*self.lam)*y[0] - \
                (sig*self.msq*y[0] - (twoax - mqx)*y[1]) / sinsq  # eq(9), rewritten
        return [dy0dt, dy1dt]

    def transform(self, steps, solun):
        for t, y in zip(steps, solun):
            one_m_xsq_a = np.power(1. - t*t, self.alpha)
            y *= one_m_xsq_a

def run_ode(cur, obs):
    """Shoot the ode for given m and observer. Could implement 
    this as an object"""
    # stepper = 'dopri5'
    stepper = 'dop853'
    atol = 0.
    rtol = 2.**-32  # 1./(2**30), 1./(2**48)
    nsteps = 9600  # 2000, 9600
    y0 = cur.init_y()
#    first_step = 2.**-30
    solver = ode(cur)
    solver.set_integrator(stepper, atol=atol, rtol=rtol, nsteps=nsteps)
    solver.set_solout(obs)
    solver.set_initial_value(y0, t0)
    return solver.integrate(t1)

## -------------------------------------------------------------------------- ##

class score_t:
    """Provide f or g depending on odd or even"""
    def __init__(self, is_even):
        self.idx = self.set_idx(is_even)
    def set_idx (self, is_even):
        if is_even is True:
            return 1
        return 0
    def __call__(self, y1):
        return y1[self.idx]

## -------------------------------------------------------------------------- ##

class solver_t:
    """Shoot for x=0 from x=1-eps, can also save eigenvalues of solution"""
    def __init__(self, m, q, is_even, r_eq, mass, period):
        self.m = m
        self.q = q
        self.score = score_t(is_even)
        self.r_eq = r_eq
        self.mass = mass
        self.period = period

    def set_m(self, m):
        self.m = m

    def set_q(self, q):
        self.q = q

    def set_is_even(self, is_even):
        self.score.set_idx(is_even)

    def shoot(self, lam):
        cur = ODE_t(self.m, self.q, lam, self.r_eq, self.mass, self.period)
        obs = Observers.max_t()
        y1 = run_ode(cur, obs)
        return self.score(y1 / obs.max_f)

    def __call__(self, lam):
        return self.shoot(lam)

    def save(self, lam):
        cur = ODE_t(self.m, self.q, lam, self.r_eq, self.mass, self.period)
        obs = Observers.save_t()
        run_ode(cur, obs)
        steps = obs.steps()
        solun = obs.solun()
        cur.transform(steps, solun)
        N = norm(solun[:,0])
        return steps, solun / N

    def interp(self, lam, steps):
        """note: interpolation is performed using cubic spline.
        Should probably check how scipy.integrate.solve_ivp interpolates;
        possibly uses dense output stepper for interpolation"""
        t_interp, y_interp = self.save(lam)
        P = interp1d(t_interp, y_interp[:,0], kind='cubic')
        Q = interp1d(t_interp, y_interp[:,1], kind='cubic')
        return np.array([P(steps), Q(steps)])

## -------------------------------------------------------------------------- ##


## Version with ecc and dlngrav as input parameters, rather than physical properties
class ODE_t_dimless:
    """Return RHS of ODE, and helpful functions for transformation"""
    # TODO: add in the parameters required for the new ODE, which are:
        # sigma (-> eccentricity; e^2 - x is called "t" and thus present already)
        # dxln(grav) (= 2\Gamma x / (2 + x^2*\Gamma))
            # \Gamma \equiv (2\bar\Omega^2 + 4\epsilon) / (1 - \bar\Omega^2)
        # With eccentricity, epsilon and om_bar_sq the parameters above are set
            # We need R_{eq} and R_{polar} for ecc/eps
            # We need R_{eq}, \Omega *and* M for om_bar_sq&R_{polar}
            # Fortunately we can choose these at the start of the simulation
            # and they will remain constants, but the initial conditions of the
            # star will now matter for the outcome
    # TODO: decide if I will input R_eq, R_polar, M, Omega *or* sigma, Gamma *or* x, om_bar_sq
    def __init__(self, m, q, lam, ecc, dlngrav):  #r_eq, mass, period):
        self.m = 1. * m
        self.q = 1. * q
        self.msq = m*m
        self.qsq = q*q

        self.ecc = ecc
        self.dlngrav = dlngrav  # this stores the gravity function with chi set

        self.alpha = .5 * abs(self.m)
        self.lam = 1. * lam

    def init_y(self):
        y0 = 1.0e-4
        sig = np.sqrt(1 - self.ecc * (1 - y0**2))
        y1 = (2 * self.alpha + self.m * self.q)*t0*y0 / (sig*(1 - t0*t0*self.qsq/sig/sig))
        return [y0, y1]  # Starting point, variation of RHS of eq (10)

    def coeffs(self, t):
        ## TODO: update as this is now Curvilinear!, not Legendre!
        sinsq = 1. - t*t
        twoax = 2. * self.alpha * t
        return np.array([twoax / sinsq, -1. / sinsq, self.lam - self.msq / sinsq, twoax / sinsq])

    def __call__(self, t, y):
    # Comments show what variables are called in legendre-ode_derivation.pdf
        sinsq = 1. - t*t  # 1-(x**2) (\equiv 1-\mu^2 \equiv sin^2 => name)
        twoax = 2. * self.alpha * t  # 2*alpha*x
        mqx = self.m * self.q * t  # m*q*x
        sig = np.sqrt(1 - self.ecc * sinsq)
        qxsqmo_sigcor = ((self.qsq * t*t / (sig*sig)) - 1) *sig  # [(x^2*q^2/sig^2)-1]*sig
        dlng = self.dlngrav(t)
        dy0dt = ( (twoax + mqx)*y[0] + qxsqmo_sigcor*y[1] ) / sinsq  # eq (8), rewritten
#        dy1dt = ( (self.lam*sinsq - self.msq)*y[0] + (twoax - mqx)*y[1] ) / sinsq
        dy1dt = (dlng + sig*self.lam)*y[0] - \
                (sig*self.msq*y[0] - (twoax - mqx)*y[1]) / sinsq  # eq(9), rewritten
        return [dy0dt, dy1dt]

    def transform(self, steps, solun):
        for t, y in zip(steps, solun):
            one_m_xsq_a = np.power(1. - t*t, self.alpha)
            y *= one_m_xsq_a


class solver_t_dimless:
    """Shoot for x=0 from x=1-eps, can also save eigenvalues of solution"""
    def __init__(self, mode_admin, q):
        self.m = mode_admin.m
        self.q = q
        self.score = score_t(mode_admin.is_even)
        self.ecc = mode_admin.ecc
        self.dlngrav = mode_admin.dlngrav
        self.wavemode = mode_admin.mode
        self.idx = 1
        if self.wavemode == "kelvin mode":
            self.idx = 0
#        if self.m*self.q > 0:
#            self.is_prograde = False
#        else:
#            self.is_prograde = True
#        self.is_even_int = 1  # Should only be 0 for Kelvin modes?
##        if self.m*self.q < 0 and ## Need more info to recognize wavemode

    def set_m(self, m):
        self.m = m

    def set_q(self, q):
        self.q = q

    def set_is_even(self, is_even):
        self.score.set_idx(is_even)

    def shoot(self, lam):
        cur = ODE_t_dimless(self.m, self.q, lam, self.ecc, self.dlngrav)
        obs = Observers.max_t(self.idx)
        y1 = run_ode(cur, obs)
#        print y1, obs.max_f
        return self.score(y1 / obs.max_f)  # TODO: find alternative for normalization!

    def __call__(self, lam):
#        print self.m*self.m*.75 < lam < self.m*self.m*1.25
#        if self.m*self.m*.75 < lam < self.m*self.m*1.25 and self.is_prograde:
#            self.is_even_int = 0
        return self.shoot(lam)

    def save(self, lam):
        cur = ODE_t_dimless(self.m, self.q, lam, self.ecc, self.dlngrav)
        obs = Observers.save_t()
        run_ode(cur, obs)
        steps = obs.steps()
        solun = obs.solun()
        cur.transform(steps, solun)
        N = norm(solun[:,0])
        return steps, solun / N

    def interp(self, lam, steps):
        """note: interpolation is performed using cubic spline.
        Should probably check how scipy.integrate.solve_ivp interpolates;
        possibly uses dense output stepper for interpolation"""
        t_interp, y_interp = self.save(lam)
        P = interp1d(t_interp, y_interp[:,0], kind='cubic')
        Q = interp1d(t_interp, y_interp[:,1], kind='cubic')
        return np.array([P(steps), Q(steps)])


#!/bin/bin/python

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt


def logistic_fn (self, t, y, r):
    return r * y * (1.0 - y)

class logistic_t:
    def __init__ (self, r):
        self.r = r
    def __call__ (self, t, y):
        return self.r * y * (1.0 - y)

class observer_t:
    """Save the current step"""
    def __init__ (self):
        self.t = []
        self.y = []
    def __call__ (self, t, y):
        self.t.append(t)
        self.y.append(list(y))  ## to use copy by value, not reference
    def steps (self):
        return np.array(self.t)
    def solun (self):
        return np.array(self.y)

def logistic_ode ():
    r = .01
    t0 = 0
    y0 = [1e-5]
    t1 = 5000.0

    ## if using the class (has r saved locally)
    logistic = logistic_t(r)
    solver = ode(logistic)

    # ## if using that function (requires r on each call)
    # solver = ode(logistic_fn)
    # solver.set_f_params(r)

    ## which stepper algorithm
    stepper = 'dop853'  # 'dopri5'
    solver.set_integrator(stepper)

    ## save the steps and solution here
    obs = observer_t()
    solver.set_solout(obs)

    ## set initial condition
    solver.set_initial_value(y0, t0)

    ## integrate:
    solver.integrate(t1)

    ## output the solution to nice np arrays
    t = obs.steps()
    solun = obs.solun()

    plt.plot(t, solun[:,0], 'b.-')
    plt.show()


def main ():
    logistic_ode()
    

if __name__ == "__main__":
    main()

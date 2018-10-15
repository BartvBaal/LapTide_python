#!/bin/bin/python

import numpy as np
from scipy.integrate import ode
from scipy.integrate import odeint
from scipy import interpolate
import matplotlib.pyplot as plt

class projectile_t:
    """Solving the vector ode: 
    \ddot{\vec{r}} = -\gamma \dot{\vec{r}} + \vec{g}
    use y = [x, y, v_x, v_y]
    and \vec{g} = - g \hat{\mathrm{e}}_y
    """
    def __init__ (self, g, gamma):
        self.g = g
        self.gamma = gamma
    def __call__ (self, t, y):
        dxdt = y[2]
        dydt = y[3]
        dvxdt = - (self.gamma * y[2])
        dvydt = - (self.gamma * y[3]) - self.g
        return [dxdt, dydt, dvxdt, dvydt]

def v0 (theta):
    v_init = 12.
    v = np.zeros(4)
    v[2] = v_init * np.cos(theta)
    v[3] = v_init * np.sin(theta)
    return v


def test_theta_range (theta_lower, theta_upper):
    g = 9.81    # gravity acceleration
    gamma = .5  # drag / mass
    projectile = projectile_t(g, gamma)

    t_steps = np.linspace(0., 3., 30)

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))

    for theta in np.linspace(theta_lower, theta_upper, 5):
        v = odeint(projectile, v0(theta), t_steps, tfirst=True)
        ax.plot(v[:, 0], v[:, 1], 'o-', label='theta={:.1e}'.format(theta))

    ax.legend()
    ax.set_ylim(0, 6)
    plt.show()


class shoot:

    def __init__ (self, gamma, g, target=[]):
        ## ode
        self.projectile = projectile_t(g, gamma)
        self.target = target
        ## 
        self.t0 = 0.
        self.tmax = 3.

    def __call__ (self, theta):
        print theta

def main ():

    tc = np.pi / 4.
    td = np.pi / 10.
    test_theta_range(tc-td, tc+td)


if __name__ == "__main__":
    main()




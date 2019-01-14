#!/bin/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from helpers.Straddle import straddle_t
from helpers.Straddle import obs_print_t


class polynomial_t:
    def __init__(self, coeffs):
        self.coeffs = coeffs

    def __call__(self, x):
        f  = 0.0
        xp = 1.0
        for ci in self.coeffs:
            f  += ci * xp
            xp *= x
        return f


def main():

    fn = polynomial_t([-4.0, 2.0, -1.5, 1.0, .5])

    x_grid = np.linspace(-4, 4, 70)
    f_grid = fn(x_grid)

    inc = .1
    N_steps = 100
    straddle = straddle_t(fn, inc, N_steps)

    x0 = -1.0
    verbose = True
    print "Straddling to find bisection:"
    bisec = straddle.search_lin(x0, verbose)

    print "root in range: {}".format(bisec)

    plt.plot(x_grid, f_grid)
    plt.plot(bisec, [0, 0], '-o')
    plt.show()


if __name__ == "__main__":
    main()




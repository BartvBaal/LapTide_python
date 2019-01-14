# I think a useful tool might be one that can work out some mode properties given a (m, l) or (m, k), for instance: whether it is odd/even, type of mode (g, r, Kelvin, Yanai), the asymptotic approximation with q, location of the equatorial waveguide with q. I suggest writing a header that can work this out and implements the approximations - you could call this something like mode-admin.
import numpy as np

class mode_admin:
    def __init__(self, m, l):
        self.m = m
        self.l = l

    def get_k(self):
        """
        Sets self.k and returns it according to the following:
        k = l - np.abs(m)
        """
        self.k = self.l - np.abs(self.m)
        return self.k

    def is_even(self):
        """
        Check if m+l is even (True) or odd (False)
        """
        if (self.m + self.l) % 2 == 0:
            return True
        return False

    def get_wavemode(self):
        """
        Returns the wavemode for the given parameters. #TODO; expand this to include all wavemodes
        """
        if not self.k:
            self.get_k()

        if self.k == -2:
            self.mode = "r mode"
        else:
            self.mode = "Cannot determine the other modes *just* from k & l ?"
        return self.mode

    def validate_values(self):
        if self.m > self.l:
            raise RuntimeError("Legendre polynomials undefined for m>l, input: ({}, {})".format(self.m, self.l)) 

import numpy as np

class print_t:
    """Print the current step"""
    def __call__ (self, t, y):
        print "t = {:.10e}, y = {}".format(t, y)

class save_t:
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

class max_t:
    """Hold the max value of the ODE for normalisation"""
    def __init__ (self, idx=1):
        self.max_f = 0.
        self.idx = idx
#    def __call__(self, t, y):
#        checkmax = max(y.max(), y.min(), key=abs)
#        if abs(checkmax) > abs(self.max_f):
#            self.max_f = checkmax
    def __call__ (self, t, y):
        if abs(y[self.idx]) > abs(self.max_f):
            self.max_f = y[self.idx]
    
class compose_t:
    """Collect a set of independent observers"""
    def __init__ (self, list_obs):
        self.list_obs = list_obs
    def __call__ (self, t, y):
        for obs in self.list_obs:
            obs(t, y)

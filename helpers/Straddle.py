
class lin:
    @staticmethod
    def inc(x, fact):
        return x + fact
    @staticmethod
    def red(x, fact):
        return x - fact
class log:
    @staticmethod
    def inc(x, fact):
        return x * fact
    @staticmethod
    def red(x, fact):
        return x / fact

class obs_empty_t:
    def __call__(self, point, itt):
        pass
class obs_print_t:
    # in the case of outputing to file, could use somethings like this:
    # def __init__ (self, output):
    #     self.output = output
    def __call__(self, point, itt):
        print "{} {} {}".format(itt, point[0], point[1])

class straddle_t:

    def __init__(self, fn, fact, N_steps, root=None):
        self.fn = fn
        self.fact = fact
        self.N_steps = N_steps
        if root is None:
            self.root = 0.0
        else:
            self.root = root

    def guess(self, x):
        return [x, self.fn(x)]

    def search_general(self, init_x, op, obs, neg_allowed=False):
        itt = 0
        pt_inc = self.guess(init_x)
        if neg_allowed:
            pt_red = self.guess(init_x)
        f0 = self.guess(init_x)[1]
        while itt < self.N_steps:
            pt_inc = self.guess(op.inc(pt_inc[0], self.fact))
            obs(pt_inc, itt)
            if pt_inc[1]*f0 < 0: 
                return [init_x, pt_inc[0]]

            if neg_allowed:
                pt_red = self.guess(op.red(pt_red[0], self.fact))
                obs(pt_red, itt)
                if pt_red[1]*f0 < 0:
                    return [pt_red[0], init_x]
            itt += 1

    def search_lin(self, init_x, verbose=False, neg_allowed=False):
        if verbose is True:
            return self.search_general(init_x, lin, obs_print_t(), neg_allowed)
        else:
            return self.search_general(init_x, lin, obs_empty_t(), neg_allowed)

    def search_log(self, init_x, verbose=False, neg_allowed=False):
        if verbose is True:
            return self.search_general(init_x, log, obs_print_t(), neg_allowed)
        else:
            return self.search_general(init_x, log, obs_empty_t(), neg_allowed)

def search_lin(fn, init_x, fact, N_steps, verbose=False):
    straddle = straddle_t(fn, fact, N_steps)
    return straddle.search_lin(init_x, verbose)

def search_log(fn, init_x, fact, N_steps, verbose=False):
    straddle = straddle_t(fn, fact, N_steps)
    return straddle.search_log(init_x, verbose)



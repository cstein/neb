import numpy

class SteepestDescent(object):
    """ The Steepest Descent method takes a step along
        the direction of the force

        R_i+1 = R_i + k * F_i

        where k is the stepsize.
    """
    def __init__(self, stepsize=1.0e-3, eps=1.0e-2, verbose=False):
        self._stepsize = stepsize
        self._eps = eps
        self._verbose = verbose

    def step(self, energy, force):
        return self._stepsize * force

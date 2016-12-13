import numpy

import util

class Atom(object):
    """ An atom

        The minimum amount of information required is the
        nuclear charge Z of the atom.

        Arguments:
        Z -- nuclear charge of atom. This argument is mandatory.

        Keyword Arguments:
        mass -- the mass of the atom in atomic units. Default is specified by using the nuclear charge.
        coords -- the Cartesian coordinate of the atom in Angstrom. Default is origin.
    """
    def __init__(self, Z, **kwargs):
        assert Z > 0, "Nuclear charge of atom must be greater than zero."
        self._z = Z
        self._c = numpy.array(kwargs.get('xyz', [0, 0, 0]))
        self._mass = kwargs.get('mass', util.MASSES[Z])
        self._vdw_radius = kwargs.get('vwdradius', util.VDWRADII[Z])
        self._cov_radius = kwargs.get('covradius', util.COVALENTRADII[Z])
        self._coordination = kwargs.get('coordination', util.COORDINATION[Z])
        self._label = util.Z2LABEL[Z]

    def getMass(self):
        return self._mass

    def getNuclearCharge(self):
        return self._z

    def getLabel(self):
        return self._label

    def getCoordinate(self):
        return self._c

    def setCoordinate(self, value):
        (n, ) = numpy.shape(value)
        assert n == 3, "Dimensions of data do not match. Expected 3 but got {}".format(n)
        self._c = value

    def getVDWRadius(self):
        return self._vdw_radius

    def getCovalentRadius(self):
        return self._cov_radius

    def setCoordination(self, value):
        if not isinstance(int, value):
            raise TypeError

        max_coordination = util.COORDINATION[self._z]
        if value > max_coordination:
            raise ValueError("Coordination number too large.")

        self._coordination = value

    def getCoordination(self):
        return self._coordination

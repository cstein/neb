from ..molecule import Molecule

import path

class Linear(path.Path):
    """ A linear interpolator that generates n-2 new molecules
    """
    def __init__(self, initial, final, nsteps=10):
        path.Path.__init__(self)

        assert isinstance(nsteps, int)

        self._molecules = [initial]

        ci = initial.getCoordinates()
        cf = final.getCoordinates()
        delta = (cf - ci) / (nsteps - 1)

        # only generate the inner range
        for k in range(1, nsteps-1):
            m2 = Molecule.fromMolecule(initial)
            m2.setCoordinates(ci + k*delta)
            self._molecules.append(m2)

        self._molecules.append(final)
        assert self.getNumBeads() == nsteps

import copy
import numpy

import atom
import bond
import angle

class Molecule(object):
    """ A molecule.

        A molecule is at the minimum a collection of atoms.

        The molecule class can also be asked to identify all bonds. This can be
        quite costly since we use a brute force approach.
    """
    _bond_threshold = 0.45 # Added threshold for bonds. Replicates openbabel

    def __init__(self):
        self._charge = 0
        self._multiplicity = 1
        self._atoms = []
        self._bonds = []
        self._name = ""

    # class methods
    @classmethod
    def fromMolecule(cls, m):
        M = cls()
        M.setCharge(m.getCharge())
        M.setMultiplicity(m.getMultiplicity())
        M.setName(m.getName())
        if m.getNumAtoms() > 0:
            M.addAtoms(*m.getAtoms())

        # currently we do not transfer bond information
        return M

    # getters and setters for various properties
    def addAtom(self, _atom):
        #assert isinstance(_atom, atom.Atom), "You attempted to add something that was not an atom."
        self._atoms.append(copy.deepcopy(_atom))

    def addAtoms(self, *args):
        for _atom in args:
            self.addAtom(_atom)

    def getNumAtoms(self):
        """ Returns the number of atoms in the molecule """
        return len(self._atoms)

    def getAtoms(self):
        for _atom in self._atoms:
            yield _atom

    def getBonds(self):
        """ Returns all bonds (as an iterator) in the molecule

            If the bond list has not been calculated before, the bonds are
            percieved through the percieveBonds method
        """
        if len(self._bonds) == 0:
            self._bonds = list(self.percieveBonds())

        for _bond in self._bonds:
            yield _bond

    def getName(self):
        return self._name

    def setName(self, value):
        assert isinstance(value, str)
        self._name = value

    def getCharge(self):
        return self._charge

    def setCharge(self, value):
        assert isinstance(value, int)
        self._charge = value

    def getMultiplicity(self):
        return self._multiplicity

    def setMultiplicity(self, value):
        assert isinstance(value, int)
        self._multiplicity = value

    # properties that are lazily evaluated such as bonds and angles
    def percieveBonds(self):
        """ This method attempts to percieve bonds

            It works by comparing atom distances to covalent radii of the atoms.
            It is not optimized in any way.
        """

        for iat, atom1 in enumerate(self.getAtoms()):
            for jat, atom2 in enumerate(self.getAtoms()):
                if iat <= jat: continue
                dr = atom2.getCoordinate() - atom1.getCoordinate()
                R2 = dr.dot(dr)

                dr_cov = atom1.getCovalentRadius() + atom2.getCovalentRadius() + self._bond_threshold
                R2_cov = dr_cov**2
                if R2 < R2_cov:
                    yield bond.Bond(id1=iat, id2=jat)

    def percieveAngles(self):
        """ This method attemps to percieve angles

            It works by iterating through all bonds in the molecule
        """
        for ibd, bond1 in enumerate(self.getBonds()):
            for jbd, bond2 in enumerate(self.getBonds()):
                if ibd <= jbd: continue
                jatm = bond1.sharesAtom(bond2)
                if jatm >= 0:
                    iatm = bond1.getNbrAtomIdx(jatm)
                    katm = bond2.getNbrAtomIdx(jatm)
                    yield angle.Angle(iatm, jatm, katm)

    # specialized options to extract information stored in
    # other classes related to molecule
    def getCoordinates(self):
        """ Returns a numpy array with all the coordinates
            of all the atoms in the molecule
        """
        c = numpy.zeros((self.getNumAtoms(), 3))
        for iat, _atom in enumerate(self.getAtoms()):
            c[iat] = _atom.getCoordinate()

        return c

    def setCoordinates(self, c):
        """ Sets the coordinates of all atoms in the molecule from
            the numpy array
        """
        assert isinstance(c, numpy.ndarray)
        (n,k) = numpy.shape(c)
        assert n == self.getNumAtoms()
        for iat, _atom in enumerate(self.getAtoms()):
            _atom.setCoordinate(c[iat])

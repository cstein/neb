import numpy


class Bond(object):
    """ A bond between two Atom objects.

        Currently, there is no reference to the actual Atom objects but instead
        they are indexed with integers in from the parent molecule class.
    """
    def __init__(self, id1, id2):
        self._id1 = id1
        self._id2 = id2
        assert self._id1 != -1
        assert self._id2 != -1
        assert self._id1 != self._id2, "indices cannot refer to same atom."

    def sharesAtom(self, other):
        """ Returns an atom index if two bonds shares an atom.
            If an atom is not found or the bonds are the same a -1 is returned.

            Arguments:
            other -- the other bond

            Returns:
            integer of atom the bonds share. -1 of None.
        """
        if self == other:
            return -1

        if self._id1 == other._id1 and self._id2 != other._id2:
            return self._id1

        if self._id1 == other._id2 and self._id2 != other._id1:
            return self._id1

        if self._id2 == other._id1 and self._id1 != other._id2:
            return self._id2

        if self._id2 == other._id2 and self._id1 != other._id1:
            return self._id2

        return -1

    def getNbrAtomIdx(self, value):
        """ Returns the neighboring atom index in the bond """
        if self._id1 == value: return self._id2
        if self._id2 == value: return self._id1
        raise ValueError("The atom index {0:d} is not in the bond.".format(value))

    def __eq__(self, other):
        """ Bonds are equal if they refer to the same atoms """
        c1 = self._id1 == other._id1 and self._id2 == other._id2
        c2 = self._id1 == other._id2 and self._id2 == other._id1
        return c1 or c2

    def __repr__(self):
        return("Bond({0:d},{1:d})".format(self._id1, self._id2))

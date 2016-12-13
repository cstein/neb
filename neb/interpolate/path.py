from .. import molecule

class Path(object):
    """ Represents a path of molecules along some coordinate """
    def __init__(self):
        self._molecules = []

    def __iter__(self):
        for c in self._molecules:
            yield c

    def __getitem__(self, index):
        return self._molecules[index]

    def getNumBeads(self):
        return len(self._molecules)


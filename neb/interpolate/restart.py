import path

class Restart(path.Path):
    """ A path which uses molecules from a previous run """
    def __init__(self, *args):
        path.Path.__init__(self)
        for _molecule in args:
            self._molecules.append(_molecule)


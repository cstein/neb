import math

import numpy

class NEB(object):
    """ A Nudged Elastic Band implementation

        This NEB implementation is based on http://dx.doi.org/10.1063/1.1323224
        by Henkelman et al.
    """
    def __init__(self, path, k):
        """ Initialize the NEB with a predefined path and force
            constants between images.

            Typical use-case might look like:

            >>> m1 = molecule_from_xyz('m1.xyz')
            >>> m2 = molecule_from_xyz('m2.xyz')
            >>> apath = neb.interpolate.Linear(m1, m2, 10)
            >>> neb = neb.Neb(apath, 5.0)
            >>> eandg = somefunction
            >>> minimizer = neb.minimizers.SteepestDescent
            >>> neb.minimize(100, 0.01, eandg, minimizer)

            Arguments:
            path -- Path between two endpoints to be optimized
            k -- force constant in units of eV / A^2 between each bead in the path
        """
        self._path = path
        self._k = k

        # set bead energies, tangents, forces and spring forces to zero initially
        self._tangents = []
        self._beadgradients = []
        self._springforces = []
        self._forces = []
        self._energies = []

        # accounting variables
        self._grms = []

        for bead in path:
            (n, k) = numpy.shape(bead.getCoordinates())
            self._tangents.append(numpy.zeros((n,k)))
            self._springforces.append(numpy.zeros((n,k)))
            self._beadgradients.append(numpy.zeros((n,k)))
            self._forces.append(numpy.zeros((n,k)))
            self._energies.append(0.0)
            self._grms.append(-1.0)

        # now we calculate the tangents and springforces
        # for the initial beads
        self._beadTangents()
        self._springForces()

    def innerBeads(self):
        """ an iterator over the inner beads """
        n = self._path.getNumBeads()
        for i, bead in enumerate(self._path):
            if i > 0 and i < n-1:
                yield bead

    def innerBeadForces(self):
        """ iterator over the forces of the inner beads """
        for i, bead in enumerate(self.innerBeads(), start=1):
            yield self._forces[i]

    def _beadTangents(self):
        """ Evaluates all tangents for all the inner beads """
        for ibead, bead in enumerate(self.innerBeads(), start=1):
            self._tangents[ibead] = self._beadTangent(bead, self._path[ibead-1], self._path[ibead+1])

    def _beadTangent(self, ibead, mbead, pbead):
        """ Calculates the tangent for ibead given the bead
            indexed by i-1 (mbead) and i+1 (pbead).

            Calculated according to eq 2 in http://dx.doi.org/10.1063/1.1323224

            Arguments:
            ibead -- the current (i'th) bead
            mbead -- the (i-1)'th bead to use in the calculation of the tanget
            pbead -- the (i+1)'th bead to use in the calculation of the tanget

            Returns:
            tanget of the bead
        """
        Ri = ibead.getCoordinates()
        Rm = mbead.getCoordinates()
        Rp = pbead.getCoordinates()

        vm = Ri - Rm
        vp = Rp - Ri
        ti = vm / numpy.linalg.norm(numpy.ravel(vm)) + vp / numpy.linalg.norm(numpy.ravel(vp));
        return ti / numpy.linalg.norm(ti)

    def _springForces(self):
        """ Evaluates all spring forces between the beads """
        for ibead, bead in enumerate(self.innerBeads(), start=1):
            self._springforces[ibead] = self._springForce(bead, self._path[ibead-1], self._path[ibead+1], self._tangents[ibead])

    def _springForce(self, ibead, mbead, pbead, tangent):
        """ Calculates the spring force for ibead given the bead
            indexed by i-1 (mbead) and i+1 (pbead).


        """
        Ri = numpy.ravel(ibead.getCoordinates())
        Rm = numpy.ravel(mbead.getCoordinates())
        Rp = numpy.ravel(pbead.getCoordinates())

        # old spring force calculated according
        # to eq 5 in http://dx.doi.org/10.1063/1.1323224
        r = numpy.dot(numpy.ravel(Rp + Rm - 2*Ri), numpy.ravel(tangent))

        return self._k * r * tangent

    def _beadGradients(self, func):
        """ Calculates the forces on each bead using the func supplied

            Calculated according to eq 4 in http://dx.doi.org/10.1063/1.1323224

            Arguments:
            bead -- the bead whose internal force is to be evaluated
            func -- function that returns energy and forces for a bead

            Returns:
            e, g -- internal energy and force with component projected out
        """
        if func is None:
            return

        for ibead, bead in enumerate(self.innerBeads(), start=1):
            energy, gradient = func(bead)
            tangent = self._tangents[ibead]

            grad_perp = numpy.dot(numpy.ravel(gradient), numpy.ravel(tangent))

            # calculate regular NEB bead gradient
            self._beadgradients[ibead] = gradient - grad_perp * tangent

            self._energies[ibead] = energy

    def beadForces(self, func):
        """ Calculates the forces of all 'inner' beads

            Arguments:
            func -- function that returns energy and forces for a bead
        """
        self._beadTangents()
        self._springForces()
        self._beadGradients(func)

        for ibead, bead in enumerate(self.innerBeads(), start=1):
            bead_force = - self._beadgradients[ibead]

            bead_force += self._springforces[ibead]

            self._forces[ibead] = bead_force[:]

            # Accounting and statistics
            f = numpy.ravel(bead_force)
            self._grms[ibead] = math.sqrt(f.dot(f)/len(f))

    def minimize(self, nsteps, opttol, func, minimizer):
        """ Minimizes the NEB path

            The minimization is carried out for nsteps to a tolerance
            of opttol with the energy and gradients calculated
            for each bead by func. The minimizer used is suppplied
            via the minimizers argument.

            When the method ends, one can iterate over all the beads
            in this class to get the states and continue from there.

            NOTE: The opttol argument is not active

            Arguments:
            nstesp -- perform a maximum of nsteps steps
            opttol -- the maximum rms gradient shall be below this value
            func -- energy and gradient function
            minimizer -- a minimizer
        """
        for i in range(1, nsteps):
            self.beadForces(func)

            s  = "-"*89 + "\nI={0:3d} ENERGY={1:12.6f} G RMS={2:13.9f}"
            s2 = " E     ="
            s3 = " F RMS ="
            s4 = " F SPR ="

            maxerg = max(self._energies[1:-1])
            grms = 0.0
            grmsnrm = 0
            for ibead, bead in enumerate(self.innerBeads(), start=1):
                c = bead.getCoordinates()
                (n, k) = numpy.shape(c)
                bead.setCoordinates(c + minimizer.step(self._energies[ibead], self._forces[ibead]))

                f = numpy.ravel(self._forces[ibead])
                grms += numpy.linalg.norm(f)
                grmsnrm += len(f)

                s2 += "{0:9.4f}".format(self._energies[ibead])
                s3 += "{0:9.4f}".format(self._grms[ibead])
                s4 += "{0:9.4f}".format(numpy.max(self._springforces[ibead]))

            print s.format(i, maxerg, math.sqrt(grms/grmsnrm))
            print s2
            print s3
            print s4

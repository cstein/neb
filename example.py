""" Implementation of the 3-atom LEPS potential

    URL: http://theory.cm.utexas.edu/henkelman/pubs/jonsson98_385.pdf
"""

import numpy
import matplotlib
import matplotlib.pyplot as plt

import neb
from neb.util import idamax
from neb.minimizers import SteepestDescent
from neb.methods import LEPSEnergyAndGradient
from neb.interpolate import Linear

def minimize(_mol, nsteps, opttol, func, minimizer):
    """ Minimizes a single molecule

        Arguments:
        nstesp -- perform a maximum of nsteps steps
        opttol -- the maximum rms gradient shall be below this value
        func -- energy and gradient function
        minimizer -- a minimizer
    """

    for k in range(nsteps):
        cc = _mol.getCoordinates()
        rab, rbc, v, g = func(_mol)
        _mol.setCoordinates(cc + minimizer.step(v, -g))

        gval = numpy.ravel(g)
        i = idamax(gval)
        gmax = gval[i]
        grms = numpy.sqrt(gval.dot(gval)/9)

        s = "Step = {0:04d} rAB = {1:6.2f} rBC = {2:6.2f} E = {3:9.4f} Gmax = {4:9.4f} Grms = {5:9.4f}".format(k, rab, rbc, v, gmax, grms)
        if k % 100==0:
            print s

        if grms < opttol:
            break

        if gmax/3 >= grms:
            break

    print "------ OPTIMIZATION CONVERGED ------"
    print s
    print
    print
    return _mol

if __name__ == '__main__':
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

    VLeps = neb.methods.leps.VLeps
    # The following example code
    # 1) plots the potential energy surface
    # 2) creates two molecules on the LEPS potential energy surface (crosses) and minimizes them (squares)
    # 3) creates a linear interpolated path between them (squares connected by black line)
    # 4) minimizes that path using NEB
    #
    # NOTE: The code is a mess because we do plotting
    #       at the same time as simulation which
    #       should be considered an abhorrent thing
    #       to do. Please do not take this code as
    #       any form of embrace of good coding
    #       practice.
    #

    # ---------------------------------
    # plot the potential energy surface
    nx = 60
    ny = 60
    x = numpy.linspace(0.4,4.0,nx)
    y = numpy.linspace(0.4,4.0,ny)
    z = numpy.zeros((nx,ny))

    for ia, xa in enumerate(x):
        for ic, yc in enumerate(y):
            m = neb.Molecule()
            m.addAtoms(neb.Atom(1, xyz=[xa, 0.0, 0.0]), neb.Atom(1, xyz=[0.0, 0.0, 0.0]), neb.Atom(1, xyz=[0.0, yc, 0.0]))
            rab, rbc, v, g = VLeps(m)
            z[ia,ic] = v

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.set_xlabel(r'$r_\mathrm{AB}$', fontsize=16)
    ax.set_ylabel(r'$r_\mathrm{BC}$', fontsize=16)
    ax.contourf(x,y,z, levels=numpy.linspace(-5.0, 0.0, 50), cmap='Blues_r')
    c = ax.contour(x,y,z, levels=numpy.linspace(-5.0, 0.0, 6), colors='white', alpha=0.3, linewidths=2)
    # ---------------------------------

    sd = SteepestDescent(stepsize=0.01)
    # ---------------------------------
    # Setup molecule 1, minimize it
    # and plot it's initial and minimized
    # coordinates
    m1 = neb.Molecule()
    xa = 0.6
    yc = 2.0
    m1.addAtoms(
        neb.Atom(1, xyz=[xa, 0.0, 0.0]),
        neb.Atom(1, xyz=[0.0, 0.0, 0.0]),
        neb.Atom(1, xyz=[0.0, yc, 0.0])
    )
    rab, rbc, v, g = VLeps(m1)
    ax.scatter([rab], [rbc], s=30, marker='x', linewidth=2, c='k')

    m1opt = minimize(m1, 1000, 0.02, VLeps, sd)
    rab, rbc, v, g = VLeps(m1opt)
    ax.scatter([rab], [rbc], s=20, marker='s', linewidth=0, c='k')

    # ---------------------------------
    # Setup molecule 2, minimize it
    # and plot it's initial and minimized
    # coordinates
    m2 = neb.Molecule()
    xa = 2.0
    yc = 0.8
    m2.addAtoms(
        neb.Atom(1, xyz=[xa, 0.0, 0.0]),
        neb.Atom(1, xyz=[0.0, 0.0, 0.0]),
        neb.Atom(1, xyz=[0.0, yc, 0.0])
    )
    rab, rbc, v, g = VLeps(m2)
    ax.scatter([rab], [rbc], s=30, marker='x', linewidth=2, c='k')

    m2opt = minimize(m2, 1000, 0.02, VLeps, sd)
    rab, rbc, v, g = VLeps(m2opt)
    ax.scatter([rab], [rbc], s=20, marker='s', linewidth=0, c='k')

    # ---------------------------------
    # linear interpolation between m1 and m2
    # optimized molecular geometries
    l = Linear(m1opt, m2opt, 20)
    n = neb.NEB(l, 1.0)
    RAB = []
    RBC = []
    for b in n.innerBeads():
        rab, rbc, v, g = VLeps(b)
        RAB.append(rab)
        RBC.append(rbc)

    ax.plot(RAB, RBC, 'k--', marker='o')

    # ---------------------------------
    # minimize path using NEB
    n.minimize(200, 0.2, LEPSEnergyAndGradient, sd)

    RAB = []
    RBC = []
    for b in n.innerBeads():
        rab, rbc, v, g = VLeps(b)
        RAB.append(rab)
        RBC.append(rbc)

    ax.plot(RAB, RBC, 'k-', marker='x')
    ax.set_xlim(0.5,4.0)
    ax.set_ylim(0.5,4.0)

    plt.show()

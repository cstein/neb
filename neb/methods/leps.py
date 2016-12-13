import numpy

def VLeps(molecule):
    """ Energy and gradient of the LEPS potential """

    def Q(r,d,dr):
        alpha = 1.942
        value = alpha * (r - 0.742)
        energy_factor = 1.5 * numpy.exp(-2.0*value) - numpy.exp(-value)
        gradient_factor = alpha*dr/r*(-3.0*numpy.exp(-2.0*value) + numpy.exp(-value))
        return 0.5*d*energy_factor, 0.5*d*gradient_factor

    def J(r,d,dr):
        alpha = 1.942
        value = alpha * (r - 0.742)
        energy_factor = numpy.exp(-2.0*value) - 6.0*numpy.exp(-value)
        gradient_factor = alpha*dr/r*(6.0*numpy.exp(-value) - 2.0*numpy.exp(-2.0*value))
        return 0.25*d*energy_factor, 0.25*d*gradient_factor

    c = molecule.getCoordinates()
    drab = c[1]-c[0]
    drbc = c[2]-c[1]
    drac = c[2]-c[0]

    rab = numpy.linalg.norm(drab)
    rbc = numpy.linalg.norm(drbc)
    rac = numpy.linalg.norm(drac)

    opai = 1.0 / (1.0 + 0.05)
    opbi = 1.0 / (1.0 + 0.30)
    opci = 1.0 / (1.0 + 0.05)
    opai2 = opai*opai
    opbi2 = opbi*opbi
    opci2 = opci*opci

    dab = 4.476
    dbc = 4.476
    dac = 3.445

    qab, gqab = Q(rab, dab, drab)
    qbc, gqbc = Q(rbc, dbc, drbc)
    qac, gqac = Q(rac, dac, drac)

    qvalue = qab*opai + qbc*opbi + qac*opci
    qgrad = gqab*opai + gqbc*opbi + gqac*opci

    jab, gjab = J(rab, dab, drab)
    jbc, gjbc = J(rbc, dbc, drbc)
    jac, gjac = J(rac, dac, drac)

    jvalue  = jab*jab*opai2
    jvalue += jbc*jbc*opbi2
    jvalue += jac*jac*opci2
    jvalue -= jab*jbc*opai*opbi
    jvalue -= jbc*jac*opbi*opci
    jvalue -= jab*jac*opai*opci

    qgrad = numpy.zeros((3,3))
    jgrad = numpy.zeros((3,3))

    qgrad[0,:] = gqab*opai + gqac*opci
    qgrad[1,:] = gqbc*opbi - gqab*opai
    qgrad[2,:] =-gqbc*opbi - gqac*opci

    djab = gjab*opai*(2*jab*opai - jac*opci - jbc*opbi)
    djbc = gjbc*opbi*(2*jbc*opbi - jab*opai - jac*opci)
    djac = gjac*opci*(2*jac*opci - jbc*opbi - jab*opai)
    jgrad[0,:] =  djab + djac
    jgrad[1,:] =  djbc - djab
    jgrad[2,:] = -djbc - djac

    return rab, rbc, qvalue - numpy.sqrt(jvalue), -(qgrad - 0.5/numpy.sqrt(jvalue)*jgrad)

def LEPSEnergyAndGradient(molecule):
    """ Wrapper for the LEPS potential

    """
    rab, rbc, e, g = VLeps(molecule)
    return e, g

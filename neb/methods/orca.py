def OrcaEnergyAndGradient(bead):
    """ Calculates the energy and gradient of a bead using ORCA

        NOTE:
        This method requires a folder called orca_scratch
        in the directory for scratch files.

        Arguments:
        bead -- the current bead / molecule to calculate
    """

    def parse_gradient(file, n):
        """ Gradient parser """
        g = numpy.zeros((n,3))
        for k in range(n):
            tokens = (file.readline()).split()
            try:
                g[k] = numpy.array(map(float, tokens[3:]))
            except ValueError:
                print tokens
                raise

        return g * util.aa2au  # Convert from Eh/bohr to Eh/AA


    """ Returns the gradient in Eh/angstrom """
    (n,k) = numpy.shape(bead.getCoordinates())

    if not os.path.exists('orca_scratch'):
        raise ValueError("ORCA scratch directory 'orca_scratch' does not exist.")

    e = 0.0
    g = numpy.zeros((n,k))

    s = "! PM3 ENGRAD\n* xyz 0 1\n"
    for _atom in bead.getAtoms():
        s += "{0:6>s}{1[0]:16.9f}{1[1]:16.9f}{1[2]:16.9f}\n".format(_atom.getLabel(), _atom.getCoordinate())
    s += "*"
    with open('orca_scratch/bead.inp', 'w') as orcafile:
        orcafile.write(s)

    # go to scratch directory
    os.chdir('orca_scratch')
    orcacmd = shlex.split("orca bead.inp") # | grep -A16 \"The cartesian gradient:\"")
    orcajob = subprocess.Popen(orcacmd, stdout=open('bead.out', 'w')) # subprocess.PIPE)

    # to avoid the python code continuing before
    time.sleep(0.5)

    with open('bead.out', 'r') as orcafile:
        line = orcafile.readline()
        while True:
            line = orcafile.readline()
            if "TOTAL RUN TIME:" in line:
                break

            tokens = line.split()

            # Semi-Empirical gradient
            if "The cartesian gradient:" in line:
                g = parse_gradient(orcafile, n)

            # HF or DFT gradient
            if "CARTESIAN GRADIENT" in line:
                line = orcafile.readline()
                line = orcafile.readline()
                g = parse_gradient(orcafile, n)

            # energy
            if "Total Energy       :" in line:
                e = float(tokens[3])

    cleancmd = shlex.split("rm -f bead.*")
    subprocess.Popen(cleancmd)
    os.chdir('..')

    return e, g

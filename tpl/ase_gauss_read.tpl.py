"""
Functions to read Gaussian output.
"""

def read_geom(outfile, dummy):
    """
    Read the final geometry from a Gaussian file.
    """

    with open(outfile) as f:
        lines = f.readlines()
    geom = np.zeros((len(mol),3))
    for index, line in enumerate(reversed(lines)):
        if re.search('Input orientation:', line) != None:
            for n in range(len(mol)):
                geom[n][0:3] = np.array(lines[-index+4+n].split()[3:6]).astype(float)
            break

    for i,d in enumerate(dummy):
        geom[-(i+1)][0:3] = d[0:3]

    return geom



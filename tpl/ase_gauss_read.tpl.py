"""
Functions to read Gaussian output.
"""

def read_geom(outfile, mol, dummy):
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


def read_zpe(outfile):
    """
    Read the zpe
    """

    with open(outfile) as f:
        lines = f.readlines()

    for line in reversed(lines):  
        if re.search('Zero-point correction=', line) != None:
            zpe = float(line.split()[2])

    return zpe


def read_freq(outfile, atom):
    """
    Read the frequencies
    """

    with open(outfile) as f:
        lines = f.readlines()

    natom = len([at for at in atom if at !='X']) #filter out the dummy atoms
    if natom == 1:
        freq = []
    else:
        freq = []
        for line in lines:
            if re.search('Frequencies', line) != None:
                if natom == 2:
                    freq.append(np.array(line.split()[2]).astype(float))
                    break
                else:
                    f = np.array(line.split()[2:5]).astype(float)
                    freq.extend(f)

    return freq

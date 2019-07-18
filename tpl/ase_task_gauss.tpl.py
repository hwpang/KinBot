# Gaussian ASE calculator

from ase.calculators.gaussian import Gaussian
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)

mol.set_calculator(calc)

# Code to convert KinBot's constraints into Gaussian constraints

bonds = []
angles = []
dihedrals = []
for fi in fix:
    if len(fi) == 2:
        #careful: atom indices in the fix lists start at 1
        bondlength = mol.get_distance(fi[0]-1, fi[1]-1)
        bonds.append([bondlength,[fi[0]-1, fi[1]-1]])
    if len(fi) == 3:
        #careful: atom indices in the fix lists start at 1
        angle = mol.get_angle(fi[0]-1, fi[1]-1, fi[2]-1) * np.pi / 180.
        angles.append([angle,[fi[0]-1, fi[1]-1, fi[2]-1]])
    if len(fi) == 4:
        #careful: atom indices in the fix lists start at 1
        dihed = mol.get_dihedral(fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1) * np.pi / 180.
        dihedrals.append([dihed,[fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1]])
for ci in change:
    if len(ci) == 3:
        #careful: atom indices in the fix lists start at 1
        bondlength = ci[2]
        bonds.append([bondlength,[ci[0]-1, ci[1]-1]])
    if len(ci) == 4:
        #careful: atom indices in the fix lists start at 1
        angle = ci[3] * np.pi / 180.
        angles.append([angle, [ci[0]-1, ci[1]-1, ci[2]-1]])
    if len(ci) == 5:
        #careful: atom indices in the fix lists start at 1
        dihed = ci[4] * np.pi / 180.
        dihedrals.append([dihed,[ci[0]-1, ci[1]-1, ci[2]-1, ci[3]-1]])


# Define the task
# ASE does not have access to geometries, ZPE and frequencies
# KinBot's reader is needed

for attempt in range(maxattempt): 
    success = 0
    try:
        if attempt > 0:  # read the last geom
            outfile = '{label}.log'
            mol.positions = reader_{qc}.read_geom(outfile, mol, dummy)

        e = mol.get_potential_energy()

        #read the geometry from the output file
        outfile = '{label}.log'
        mol.positions = reader_{qc}.read_geom(outfile, mol, dummy)

        # read frequency and zpe if available:
        if freq:
            zpe = reader_{qc}.read_zpe(outfile)
            freq = reader_{qc}.read_freq(outfile, atom)

        for d in dummy:  #remove the dummy atom before writing to the database
            mol.pop()
 
        db = connect('kinbot.db')

        if freq:
            db.write(mol,name = label,data={{'energy': e,'frequencies': np.asarray(freq), 'zpe':zpe, 'status' : 'normal'}})
        else:
            db.write(mol, name=label, data={{'energy': e,'status' : 'normal'}})

        success = 1
        break

    except:
        pass


if success == 0:
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'error'}})


with open(label + '.log', 'a') as f:
    f.write('done\n')

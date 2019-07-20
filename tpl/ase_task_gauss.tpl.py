# Gaussian ASE calculator

from ase.calculators.gaussian import Gaussian
Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
calc = Gaussian(**kwargs)

mol.set_calculator(calc)

# only used if change is invoked, otherwise modredundant is used
bonds, angles, dihedrals = reader_{qc}.constraint(mol, fix, change)


# Define the task
# ASE does not have access to geometries, ZPE and frequencies
# KinBot's reader is needed

e = 0.

for attempt in range(maxattempt): 
    success = 0
    try:
        if attempt > 0:  # read the last geom
            outfile = '{label}.log'
            mol.positions = reader_{qc}.read_geom(outfile, mol, dummy)

        e = mol.get_potential_energy()
        success = 1
        break
    except:
        pass

# for IRC we don't care if the run has an error
if success or task == 'ircf' or task == 'ircr':
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
        db.write(mol,name = label,data={{'energy': e, 
                                         'geometry': mol.positions, 
                                         'frequencies': np.asarray(freq), 
                                         'zpe':zpe, 
                                         'status' : 'normal'}})
    else:
        db.write(mol, name=label, data={{'energy': e, 
                                         'geometry': mol.positions, 
                                         'status' : 'normal'}})
else:
    db = connect('kinbot.db')
    db.write(mol, name = label, data = {{'status' : 'error'}})


with open(label + '.log', 'a') as f:
    f.write('done\n')

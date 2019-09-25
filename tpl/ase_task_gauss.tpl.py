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


"""
Template to define the task
"""

for attempt in range({maxattempt}): 
    success = 0
    try:
        if attempt > 0:  # read the last geom
            outfile = '{label}.log'
            mol.positions = read_geom(outfile, mol, dummy)

        e = mol.get_potential_energy()

        #read the geometry from the output file
        outfile = '{label}.log'
        mol.positions = read_geom(outfile, mol, dummy)

        # read frequency and zpe if available:
        if {freq}:
            zpe = read_zpe(outfile)
            freq = read_freq(outfile, atom)

        for d in dummy:  #remove the dummy atom before writing to the database
            mol.pop()
 
        db = connect('kinbot.db')

        if {freq}:
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


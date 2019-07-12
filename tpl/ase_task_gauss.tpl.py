"""
Template to define the task
"""


for attempt in range({maxattempt}): 
    success = 0
    try:
        if attempt > 0:  # read the last geom
            outfile = '{label}.log'
            mol.positions = read_geom(outfile, dummy)

        e = mol.get_potential_energy()

        #read the geometry from the output file
        outfile = '{label}.log'
        mol.positions = read_geom(outfile, dummy)

        for d in dummy:  #remove the dummy atom before writing to the database
            mol.pop()
 
        db = connect('kinbot.db')
        db.write(mol, name=label, data={{'energy': e,'status' : 'normal'}})

        success = 1
        break

    except:
        pass


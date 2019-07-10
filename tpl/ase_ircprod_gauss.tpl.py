"""
Optimize the products when using Gaussian with Gaussian's optimizer.
"""

if success:
    # start the product optimization
    pr_kwargs = {prod_kwargs}
    label = '{label}_prod'
    
    calc = Gaussian(**pr_kwargs)
    mol = Atoms(symbols=atom, positions=geom)
    mol.set_calculator(calc)

    for attempt in {maxattempt}:
        try:
            e = mol.get_potential_energy() 
            outfile = '{label}_prod.log'
            mol.positions = read_gauss(outfile) 

            db = connect('kinbot.db')
            db.write(mol, name=label, data={{'energy': e,'status': 'normal'}})

    except RuntimeError: 

        #read the geometry from the output file
        outfile = '{label}_prod.log'
        mol.positions = read_gauss(outfile) 
                
        if new_geom:
            mol.positions = geom
            db = connect('kinbot.db')
            db.write(mol, name = label, data = {{'status' : 'normal'}}) #although there is an error, continue from the final geometry
        else:
            db = connect('kinbot.db')
            db.write(mol, name = label, data = {{'status' : 'error'}})



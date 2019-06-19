"""
Template to define the task
"""

import sella  # TODO

for attempt in {maxattempt}: 
    constraints = {x}
    sella(mol, calc, {task}, constraints)

    if sella_success:
        if {freq}:
            sella(mol, calc, 'freq', None)

        for d in dummy:  #remove the dummy atom before writing to the database
            mol.pop()

        db = connect('kinbot.db')
        db.write(mol, name = label, data = {{'energy': e,'status' : 'normal'}})

        if {prodopt}:
            sella(mol, calc, 'opt', None)

        break
    
    if attempt == {maxattempt}:
         db.write(mol, name = label, data = {{'status' : 'error'}}) 




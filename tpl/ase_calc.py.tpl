"""
Setting up the approproate calculator
"""
if {qc} == 'gauss':
    from ase.calculators.gaussian import Gaussian
    Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
    calc = Gaussian(**kwargs)
elif {qc} == 'nwchem':
    from ase.calculators.nwchem import NWChem
    NWChem.command = 'mpirun -np {ppn} -path /usr/local/bin nwchem PREFIX.nw > PREFIX.out'
    calc = NWChem(**kwargs)

mol.set_calculator(calc)



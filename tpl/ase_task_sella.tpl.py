#Provide the constraints in Sella format.
#Sella is zero-indexed, i.e., first atom is 0.

constraints = {{'fix': [], 'bonds': [], 'angles': [], 'dihedrals': []}}  # fix is different for sella than for kinbot

for f in fix:
    if len(f) == 2:  # bond
        constraints['bonds'].append((f[0], f[1]))
    if len(f) == 3:  # angle
        constraints['angles'].append((f[0], f[1], f[2]))
    if len(f) == 4:  # dihedral
        constraints['dihedrals'].append((f[0], f[1], f[2], f[3]))

for f in change:
    if len(f) == 3:  # bond
        constraints['bonds'].append((f[0], f[1]), f[2])
    if len(f) == 4:  # angle
        constraints['angles'].append((f[0], f[1], f[2]), f[3])
    if len(f) == 5:  # dihedral
        constraints['dihedrals'].append((f[0], f[1], f[2], f[3]), f[4])

#Template to define the task

#!/usr/bin/env python3


from sella import Sella

#myatoms.calc = calc

e = mol.get_potential_energy()  # do a first calculation to create a chk file

if qc == 'gauss':
    kwargs['guess'] = 'Read'
    calc = Gaussian(**kwargs)
elif qc == 'nwchem':
    calc = NWChem(**kwargs)

mol.set_calculator(calc)
outfile = '{label}.log'

# delta0: the initial step size, default is 1e-1 for minimization 
#         and 1.3e-3 for saddle point
#         likely need to change the minimization default to smaller
# gamma: the convergence criterion for the Hessian, default is 0.4
#        For molecular systems it needs to be smaller. Max is about 100, min is 1e-15
if order: # saddle point
    dyn = Sella(mol, trajectory='{label}.traj', order=order, eig=False, gamma=1e-2)
else:
    dyn = Sella(mol, trajectory='{label}.traj', order=order, eig=False, delta0=1e-3)

# fmax: optimization ends if the forces are below the threshold
#       If set to zero, then it is not likely to cut the loop
# steps: optimization ends if the number of steps reaches it
for converged in dyn.irun(fmax=0., steps=3000):
    # x = mol.get_positions()
    # f = mol.get_forces()

    # the Hessian matrix
    # H = dyn.mm.H

    # its eigenvalues excluding translation/rotation
    # lams = dyn.mm.lams

    # the corresponding eigenvectors
    # vecs = dyn.mm.vecs

    if reader_{qc}.read_convergence(outfile) == 1:
        success = 1
        break



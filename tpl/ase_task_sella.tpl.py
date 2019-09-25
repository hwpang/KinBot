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

dyn = Sella(mol, trajectory='{label}.traj', order=order, eig=False)

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
        break

#from sella import MinModeAtoms, optimize
#
## Create a Sella MinMode object
#myminmode = MinModeAtoms(mol,  # Your Atoms object
#                         calc,  # Your calculator
#                         constraints=constraints,  # Your constraints
#                         trajectory='{label}.traj',  # Optional trajectory TODO remove option
#                         )
#
## These need to be keywords in KinBot TODO
#x1 = optimize(myminmode,    # Your MinMode object
#              maxiter=500,  # Maximum number of force evaluations
#              ftol=1e-3,    # Norm of the force vector, convergence threshold
#              r_trust=5e-4, # Initial trust radius (Angstrom per d.o.f.)
#              order=order,  # Order of saddle point to find (set to 0 for minimization)
#              dxL=1e-4,     # Finite difference displacement magnitude (Angstrom)
#              maxres=0.1,   # Maximum residual for eigensolver convergence (should be <= 1)
#              eig=(order!=0), # Switch on eigensolve for saddle points only
#              )
#

# if freq:
    # TODO what to do?



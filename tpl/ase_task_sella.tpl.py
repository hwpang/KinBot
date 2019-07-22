#Provide the constraints in Sella format.
#Sella is zero-indexed, i.e., first atom is 0.

fix = {fix}
change = {change}

constraints = {'fix': [], 'bonds': [], 'angles': [], 'dihedrals': []}  # fix is different for sella than for kinbot

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

from sella import MinModeAtoms, optimize

# Create a Sella MinMode object
myminmode = MinModeAtoms(mol,  # Your Atoms object
                         calc,  # Your calculator
                         constraints=dict(fix=fix),  # Your constraints
                         trajectory='{label}.traj',  # Optional trajectory
                         )

# These need to be keywords in KinBot TODO
x1 = optimize(myminmode,    # Your MinMode object
              maxiter=500,  # Maximum number of force evaluations
              ftol=1e-3,    # Norm of the force vector, convergence threshold
              r_trust=5e-4, # Initial trust radius (Angstrom per d.o.f.)
              order=order,  # Order of saddle point to find (set to 0 for minimization)
              dxL=1e-4,     # Finite difference displacement magnitude (Angstrom)
              maxres=0.1,   # Maximum residual for eigensolver convergence (should be <= 1)
              )


if freq:
    # TODO what to do?



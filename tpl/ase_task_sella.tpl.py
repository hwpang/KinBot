#Provide the constraints in Sella format.
#Sella is zero-indexed, i.e., first atom is 0.

constraints = {{'fix': [], 'bonds': [], 'angles': [], 'dihedrals': []}}  # fix is different for sella than for kinbot

for f in fix:
    if len(f) == 2:  # bond
        constraints['bonds'].append((f[0]-1 , f[1]-1))
    if len(f) == 3:  # angle
        constraints['angles'].append((f[0]-1, f[1]-1, f[2]-1))
    if len(f) == 4:  # dihedral
        constraints['dihedrals'].append((f[0]-1, f[1]-1, f[2]-1, f[3]-1))

for f in change:
    if len(f) == 3:  # bond
        constraints['bonds'].append(((f[0]-1, f[1]-1), f[2]))
    if len(f) == 4:  # angle
        constraints['angles'].append(((f[0]-1, f[1]-1, f[2]-1), f[3]))
    if len(f) == 5:  # dihedral
        constraints['dihedrals'].append(((f[0]-1, f[1]-1, f[2]-1, f[3]-1), f[4]))

#Template to define the task

#!/usr/bin/env python3

from sella import IRC
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

# IRC
if task[:3] == 'irc':
    job = '{label}'[:-6] 
    if qc == 'gauss':
        hess = reader_{qc}.read_hess(job, natom)
        for i, row in enumerate(hess):
            for j, _ in enumerate(row):
                hess[i][j] *= Hartree / Bohr / Bohr

if task[:4] == 'ircf':
    dyn = IRC(mol, trajectory='{label}.traj', H0=hess, dx=0.1, eta=1e-4, gamma=1e-2)
    for converged in dyn.irun(fmax=0.01, steps=100, direction='forward'):
        if converged:
            success = 1
            break

elif task[:4] == 'ircr':
    dyn = IRC(mol, trajectory='{label}.traj', H0=hess, dx=0.1, eta=1e-4, gamma=1e-2)
    for converged in dyn.irun(fmax=0.01, steps=100, direction='reverse'):
        if converged:
            success = 1
            break

# optimization
else:
    try:
        hess = np.load('{label}.npy')
        if hess == None:
            hess = None
    except:
       hess = None 
    
    if len(constraints['fix']) + len(constraints['bonds']) + len(constraints['angles']) + len(constraints['dihedrals']):
        constrained = True
    else:
        constrained = False

    if constrained:
        if tight:
            fmax = 0.005
            constraints_tol = 1e-5
        else:
            fmax = 0.05
            constraints_tol = 1e-3
    else:
        fmax = 0. # use Gaussian criterion TODO implement for other codes in a general sense
        constraints_tol = 1e-5

    if order: # saddle point
        dyn = Sella(mol, constraints=constraints, trajectory='{label}.traj', H0=hess, 
                order=order, eig=True, gamma=1e-2, append_trajectory={app_traj}, constraints_tol=constraints_tol)
    else:
        dyn = Sella(mol, constraints=constraints, trajectory='{label}.traj', H0=hess, 
                order=order, eig=False, delta0=1e-3, append_trajectory={app_traj}, constraints_tol=constraints_tol)

    for converged in dyn.irun(fmax=fmax, steps=100):
        if constrained:
            if converged:
                success = 1
                break
        else:
            if reader_{qc}.read_convergence(outfile) > 0:
                success = 1
                break
    if task == 'preopt' or task == 'preopt0':
        success = 1

    np.save('{label}', dyn.pes.H)

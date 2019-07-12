"""
Provide the constraints in Sella format.
Sella is zero-indexed, i.e., first atom is 0.
"""

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



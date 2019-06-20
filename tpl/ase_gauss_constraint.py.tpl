"""
Code to convert KinBot's constraints into Gaussian constraints
"""

fix = {fix}
change = {change}

bonds = []
angles = []
dihedrals = []
for fi in fix:
    if len(fi) == 2:
        #careful: atom indices in the fix lists start at 1
        bondlength = mol.get_distance(fi[0]-1, fi[1]-1)
        bonds.append([bondlength,[fi[0]-1, fi[1]-1]])
    if len(fi) == 3:
        #careful: atom indices in the fix lists start at 1
        angle = mol.get_angle(fi[0]-1, fi[1]-1, fi[2]-1) * np.pi / 180.
        angles.append([angle,[fi[0]-1, fi[1]-1, fi[2]-1]])
    if len(fi) == 4:
        #careful: atom indices in the fix lists start at 1
        dihed = mol.get_dihedral(fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1) * np.pi / 180.
        dihedrals.append([dihed,[fi[0]-1, fi[1]-1, fi[2]-1, fi[3]-1]])
for ci in change:
    if len(ci) == 3:
        #careful: atom indices in the fix lists start at 1
        bondlength = ci[2]
        bonds.append([bondlength,[ci[0]-1, ci[1]-1]])
    if len(ci) == 4:
        #careful: atom indices in the fix lists start at 1
        angle = ci[3] * np.pi / 180.
        angles.append([angle, [ci[0]-1, ci[1]-1, ci[2]-1]])
    if len(ci) == 5:
        #careful: atom indices in the fix lists start at 1
        dihed = ci[4] * np.pi / 180.
        dihedrals.append([dihed,[ci[0]-1, ci[1]-1, ci[2]-1, ci[3]-1]])


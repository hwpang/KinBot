"""
Translation template to use ASE for Gaussian when using Sella.
No keywords that affect forces or geometry should appear here.
"""

kwargs = {
'nprocshared' : {ppn},
'mem' : '1000MW',
'chk' : {label},
'label': {label}, 
'NoSymm' : 'NoSymm',
'multiplicity': {mult},
'charge': {charge},
'scf' : 'xqc',
}

# level of theory keywords
if {level} == 0:
    kwargs['method']: 'am1' 
    kwargs['basis']: ''
if {level} == 1:
    kwargs['method']: {method}
    kwargs['basis']: {basis}
if {level} == 2:
    kwargs['method'] = {high_level_method}
    kwargs['basis'] = {high_level_basis}
    if len({integral}) > 0:
        kwargs['integral'] = {integral}
if {level} == 'm': # if {scan} or 'R_Addition_MultipleBond' in {job}:
    kwargs['method'] = 'mp2'
    kwargs['basis']: {basis}

dummy = {dummy}


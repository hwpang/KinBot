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

if {delchk}:
   del kwargs['chk'] 

if {ts}:
    if {step} < {max_step}:
        kwargs['method'] = 'am1'
        kwargs['basis'] = ''
    else:
        kwargs['method'] = {method}
        kwargs['basis'] = {basis}

if {scan} or 'R_Addition_MultipleBond' in {job}:
        kwargs['method'] = 'mp2'
        kwargs['basis'] = {basis}

if {irc} is not None: 
    kwargs['method']:i {method}
    kwargs['basis']: {basis}

if {high_level}:
    kwargs['method'] = {high_level_method}
    kwargs['basis'] = {high_level_basis}
    if len({integral}) > 0:
        kwargs['integral'] = {integral}

if {hir}:
    kwargs['method'] = {high_level_method}
    kwargs['basis'] = {high_level_basis}
    if len({integral}) > 0:
        kwargs['integral'] = {integral}

dummy = {dummy}


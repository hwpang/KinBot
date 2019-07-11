"""
Translation template to use ASE for Gaussian when using Sella.
No keywords that effect forces or geometry should appear here.
"""

kwargs = {{
'nprocshared' : {ppn},
'mem' : '1000MW',
'label': '{label}', 
'method': '{method}',
'basis': '{basis}',
'NoSymm' : 'NoSymm',
'multiplicity': {mult},
'charge': {charge},
'scf' : 'xqc',
}}

if {chk}:
    kwargs['chk'] = '{label}'
if {guess}:
    kwargs['guess'] = 'Read'
if len('{integral}') > 0:
    kwargs['integral'] = '{integral}'

dummy = {dummy}


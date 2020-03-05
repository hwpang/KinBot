# Translation template to use ASE for Gaussian when using Sella.
# No keywords that affect forces or geometry should appear here.

kwargs = {{
'nprocshared' : ppn,
'mem' : mem + memu,
'label': label, 
'method': method,
'basis': basis,
'NoSymm' : 'NoSymm',
'multiplicity': mult,
'charge': charge,
'scf' : 'xqc',
'population' : 'None',
}}

if chk:
    try:
        kwargs['chk'] = label[label.index('/')+1:] # remove the directory if it was there
    except:
        kwargs['chk'] = label
if guess:
    kwargs['guess'] = 'Read'
if len(integral) > 0:
    kwargs['integral'] = integral
if guessmix and not guess:  # only for fresh start jobs
    kwargs['guess'] = 'Mix' 



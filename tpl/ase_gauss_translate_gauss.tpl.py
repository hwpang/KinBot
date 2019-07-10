"""
Translation template to use ASE for Gaussian when using Gaussian's optimizer.
This only contains items that are in addition to 
ase_gauss_translate_sella
"""

if not {sella}

    if {ts}:
        if {step} < {max_step}: 
            kwargs['opt'] = 'ModRedun,Tight,CalcFC,MaxCycle=999'
        else:
            kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'

    if {freq}:
        kwargs['freq'] = 'freq'

    if {irc}: 
        kwargs['geom'] = 'AllCheck,NoKeepConstants'
        kwargs['irc'] = 'RCFC,{},MaxPoints={},StepSize={}'.format({irc}, {irc_maxpoints}, {irc_stepsize})

    if {hir}:
        kwargs['opt'] = 'ModRedun,CalcFC'
        if {ts}:
            kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest,MaxCycle=999'




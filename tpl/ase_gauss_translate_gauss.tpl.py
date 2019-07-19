# Translation template to use ASE for Gaussian when using Gaussian's optimizer.
# This only contains items that are in addition to 
# ase_gauss_translate_sella

if not sella:

    if order:
        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'
    else:
        kwargs['opt'] = 'ModRedun,Tight,CalcFC,MaxCycle=999'

    if freq:
        kwargs['freq'] = 'freq'

    if task == 'ircf':
        kwargs['geom'] = 'AllCheck,NoKeepConstants'
        kwargs['irc'] = 'RCFC,Forward,MaxPoints={{}},StepSize={{}}'.format(irc_maxpoints, irc_stepsize)
    if task == 'ircr':
        kwargs['geom'] = 'AllCheck,NoKeepConstants'
        kwargs['irc'] = 'RCFC,Reverse,MaxPoints={{}},StepSize={{}}'.format(irc_maxpoints, irc_stepsize)

    if task == 'hir':
        if order:
            kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest,MaxCycle=999'
        else:
            kwargs['opt'] = 'ModRedun,CalcFC'

    if fix:
        kwargs['fix'] = fix

    if release:
        kwargs['release'] = release

    if change:
        kwargs['change'] = change



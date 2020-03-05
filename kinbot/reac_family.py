###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
"""
Generic methods for the reaction families
"""

import os
import numpy as np
import copy
import time
import pkg_resources

from kinbot import modify_geom

def carry_out_reaction(rxn, step, command, sella):
    """
    Verify what has been done and what needs to be done
    
    skip: boolean which tells to skip the first 12 steps in case of an instance shorter than 4
    
    scan: boolean which tells if this is part of an energy scan along a bond length coordinate
    """
    if step > 0:
        status = rxn.qc.check_qc(rxn.instance_name)
        if status != 'normal' and status != 'normal freq' and status != 'error': return step

    skipped = 0

    if step == 0:
        if rxn.qc.is_in_database(rxn.instance_name):
            if rxn.qc.check_qc(rxn.instance_name) == 'normal' or rxn.qc.check_qc(rxn.instance_name) == 'normal freq': 
                err, freq = rxn.qc.get_qc_freq(rxn.instance_name, rxn.species.natom)
                if err == 0 and len(freq) > 0.:
                    err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom)
                    step = rxn.max_step + 1
                    return step
        if rxn.skip and len(rxn.instance) < 4: 
            step = 12
            skipped = 1
        geom = rxn.species.geom
    else:
        err, geom = rxn.qc.get_qc_geom(rxn.instance_name, rxn.species.natom, allow_error = 1)

    #the the constraints for this step
    step, fix, change, release = rxn.get_constraints(step, geom)

    if step > rxn.max_step:
        return step
    
    pcobfgs = 0
    if pcobfgs == 0:
        #apply the geometry changes internally and fix the coordinates that changed
        #if not rxn.qc.sella or step < rxn.max_step - 1:
        if 1:
            change_starting_zero = []
            for c in change:
                c_new = [ci - 1 for ci in c[:-1]]
                c_new.append(c[-1])
                change_starting_zero.append(c_new)
            if len(change_starting_zero) > 0:
                success, geom = modify_geom.modify_coordinates(rxn.species, rxn.instance_name, geom, change_starting_zero, rxn.species.bond)
                if not rxn.qc.sella or step < rxn.max_step -1:
                    for c in change:
                        fix.append(c[:-1])
                    change = []
        if step > 0 and rxn.qc.sella:
            app_traj = True
        else:
            app_traj = None
        if step == rxn.max_step - 1 or skipped:
            tight = True
        else:
            tight = False
        if rxn.scan or 'R_Addition_MultipleBond' in rxn.instance_name:
            step += rxn.qc.assemble_ase_template(rxn.instance_name, 'optmp2', rxn.species, geom, 1, rxn.qc.sella, 
                    fix=fix, change=change, release=release, app_traj=app_traj, tight=True, singlejob=False)
        elif step == 0 or skipped:
            step += rxn.qc.assemble_ase_template(rxn.instance_name, 'preopt0', rxn.species, geom, 0, rxn.qc.sella, 
                    fix=fix, change=change, release=release, app_traj=app_traj, tight=tight)
        elif step < rxn.max_step:
            step += rxn.qc.assemble_ase_template(rxn.instance_name, 'preopt', rxn.species, geom, 0, rxn.qc.sella, 
                    fix=fix, change=change, release=release, app_traj=app_traj, tight=tight, singlejob=False)
        else:
            step += rxn.qc.assemble_ase_template(rxn.instance_name, 'opt', rxn.species, geom, 1, rxn.qc.sella, 
                    fix=fix, change=change, release=release, app_traj=app_traj, singlejob=False)
        
    else:
        # use the pcobfgs algorithm for the geometry update, currently disabled and abandoned
        if step < rxn.max_step:
            del kwargs['opt']
            conv_crit = 0.01  # force convergence criterion 
            template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_search_pcobfgs.py.tpl'.format(qc = rxn.qc.qc))
            template = open(template_file,'r').read()
            template = template.format(label=rxn.instance_name, kwargs=kwargs, atom=list(rxn.species.atom), 
                                       geom=list([list(gi) for gi in geom]), ppn=rxn.qc.ppn, fix=fix,
                                       change=change, conv_crit=conv_crit)
        else:
            template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_end.py.tpl'.format(qc = rxn.qc.qc))
            template = open(template_file,'r').read()
            template = template.format(label = rxn.instance_name, kwargs = kwargs, atom = list(rxn.species.atom), 
                                       geom = list([list(gi) for gi in geom]), ppn = rxn.qc.ppn, qc_command=command)
    
    return step
    

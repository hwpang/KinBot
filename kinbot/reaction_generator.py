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
from __future__ import print_function
import numpy as np
import sys
import os
import copy
import time
import logging

from kinbot import constants
from kinbot import geometry
from kinbot import pes
from kinbot import postprocess
from kinbot import reac_family
from kinbot.irc import IRC
from kinbot.optimize import Optimize
from kinbot.stationary_pt import StationaryPoint


class ReactionGenerator:
    """
    This class generates the reactions using the qc codes
    It builds an initial guess for the ts, runs that ts towards a saddle point
    and does IRC calculations 
    """
    
    def __init__(self, species, par, qc):
        self.species = species
        self.par = par
        self.qc = qc
    
    def generate(self):
        """ 
        Creates the input for each reaction, runs them, and tests for success.
        If successful, it creates the barrier and product objects.
        It also then does the conformational search, and finally, the hindered rotor scans.
        To make the code the most efficient, all of these happen in parallel, in a sense that
        the jobs are not waiting for each other. E.g., one reaction can still be in the stage
        of TS search, while the other can be already at the hindered rotor scan. This way, 
        all cores are occupied efficiently.

        The switching between the various stages are done via the reac_ts_done variable.
        TODO: write the list, the old one was very outdated
        
        If at any times the calculation fails, reac_ts_done is set to -999.
        If all steps are successful, reac_ts_done is set to -1.
        """

        deleted = []
        if len(self.species.reac_inst) > 0:
            alldone = 1
        else: 
            alldone = 0

        while alldone:
            for index, instance in enumerate(self.species.reac_inst):
                obj = self.species.reac_obj[index]
                instance_name = obj.instance_name

                # START REATION SEARCH
                if self.species.reac_ts_done[index] == 0 and self.species.reac_step[index] == 0:
                    #verify after restart if search has failed in previous kinbot run
                    status = self.qc.check_qc(instance_name)
                    if status == 'error' or status == 'killed':  # TODO remove killed everywhere??
                        logging.info('\tRxn search failed (error or killed) for {}'.format(instance_name))
                        self.species.reac_ts_done[index] = -999
                
                if self.species.reac_ts_done[index] == 0: # ts search is ongoing
                    
                    if obj.scan == 0: #don't do a scan of a bond
                        if self.species.reac_step[index] == obj.max_step + 1:  # reached last step, no freq yet though
                            status = self.qc.get_qc_energy(instance_name, self.species.natom)[0]
                            if status == 0 and self.species.reac_ts_freq_done[index] == 0:  # meaning success and freq was not started
                                # start freq calculation
                                opt_geom = self.qc.get_qc_geom(instance_name, self.species.natom, wait=1, allow_error=0)[1]
                                self.qc.qc_freq(self.species, instance_name, opt_geom, 1)
                                self.species.reac_ts_freq_done[index] = 1
                            elif status == 0 and self.species.reac_ts_freq_done[index] == 1:  # meaning success and freq is running
                                status, freq = self.qc.get_qc_freq(instance_name, self.species.natom)
                                if status == 0:
                                    self.species.reac_ts_freq_done[index] = 2  # won't return here any more
                                    if freq[0] < 0. and freq[1] > 0.:
                                        self.species.reac_ts_done[index] = 1
                                    else:
                                        logging.info('\tNot the right number of imeginary frequencies for {}'.format(instance_name))
                                        self.species.reac_ts_done[index] = -999
                                else:
                                    self.species.reac_ts_done[index] = -999
                                    logging.info('\tRxn search failed for {}'.format(instance_name))
                            elif status == -1:
                                logging.info('\tRxn search failed for {}'.format(instance_name))
                                self.species.reac_ts_done[index] = -999
                        else: 
                            self.species.reac_step[index] = reac_family.carry_out_reaction(obj, self.species.reac_step[index], self.par.par['qc_command'], self.par.par['sella'])
                    
                    else: # do a bond scan
                        if self.species.reac_step[index] == self.par.par['scan_step'] + 1:
                            status = self.qc.get_qc_energy(instance_name, self.species.natom)[0]
                            if status == 0 and self.species.reac_ts_freq_done[index] == 0:  # meaning success and freq was not started
                                # start freq calculation
                                opt_geom = self.qc.get_qc_geom(instance_name, self.species.natom, wait=1, allow_error=0)[1]
                                self.qc.qc_freq(self.species, instance_name, opt_geom, 1)
                                self.species.reac_ts_freq_done[index] = 1
                            elif status == 0 and self.species.reac_ts_freq_done[index] == 1:  # meaning success and freq is running
                                status, freq = self.qc.get_qc_freq(instance_name, self.species.natom)
                                if status == 0:
                                    self.species.reac_ts_freq_done[index] = 2  # won't return here any more
                                    if freq[0] < 0. and freq[1] > 0.:
                                        self.species.reac_ts_done[index] = 1
                                    else:
                                        logging.info('\tNot the right number of imeginary frequencies for {}'.format(instance_name))
                                        self.species.reac_ts_done[index] = -999
                                else:
                                    self.species.reac_ts_done[index] = -999
                                    logging.info('\tRxn search failed for {}'.format(instance_name))
                            elif status == -1:
                                logging.info('\tRxn search failed for {}'.format(instance_name))
                                self.species.reac_ts_done[index] = -999
                        else:        
                            if self.species.reac_step[index] == 0:
                                self.species.reac_step[index] = reac_family.carry_out_reaction(obj, self.species.reac_step[index], self.par.par['qc_command'], self.par.par['sella'])
                            elif self.species.reac_step[index] > 0:
                                status = self.qc.check_qc(instance_name)
                                if status == 'error' or status == 'killed':
                                    logging.info('\tRxn search failed for {}'.format(instance_name))
                                    self.species.reac_ts_done[index] = -999
                                else:
                                    err, energy = self.qc.get_qc_energy(instance_name)
                                    if err == 0:
                                        self.species.reac_scan_energy[index].append(energy)
                                        if len(self.species.reac_scan_energy[index]) > 1:
                                            if self.species.reac_scan_energy[index][-1] < self.species.reac_scan_energy[index][-2]:
                                                self.species.reac_step[index] = self.par.par['scan_step'] 
                                        self.species.reac_step[index] = reac_family.carry_out_reaction(obj, self.species.reac_step[index], self.par.par['qc_command'], self.par.par['sella'])

                elif self.species.reac_ts_done[index] == 1:
                    status = self.qc.check_qc(instance_name)
                    if status == 'running': continue
                    elif status == 'error': 
                        logging.info('\tRxn search failed (Gaussian error) for {}'.format(instance_name))
                        self.species.reac_ts_done[index] = -999
                    else: 
                        #check the barrier height:
                        if self.species.reac_type[index] == 'R_Addition_MultipleBond':
                            sp_energy = self.qc.get_qc_energy(str(self.species.chemid) + '_well_mp2')[1]
                            barrier = (self.qc.get_qc_energy(instance_name)[1] - sp_energy) * constants.AUtoKCAL
                        else:
                            sp_energy = self.qc.get_qc_energy(str(self.species.chemid) + '_well')[1]
                            barrier = (self.qc.get_qc_energy(instance_name)[1] - sp_energy) * constants.AUtoKCAL
                        if barrier > self.par.par['barrier_threshold']:
                            logging.info('\tRxn barrier too high ({val}) for {name}'.format(val=barrier,name=instance_name))
                            self.species.reac_ts_done[index] = -999
                        else:
                            obj.irc = IRC(obj, self.par) #TODO: this doesn't seem like a good design
                            irc_status = obj.irc.check_irc()
                            if 0 in irc_status:
                                # No IRC started yet, start the IRC now
                                logging.info('\tStarting IRC calculations for {}'.format(instance_name))
                                obj.irc.do_irc_calculations()
                            elif irc_status[0] == 'running' or irc_status[1] == 'running':
                                continue
                            else: 
                                #The main IRC calculation is done, now onto IRC products
                                #It was either 'normal' or 'error'
                                #TODO not allowing unfinished jobs!!!
                                irc_prod_status = obj.irc.check_irc_prod()
                                if 0 in irc_prod_status:
                                    logging.info('\t\tStarting IRC product calculations for {}'.format(instance_name))
                                    obj.irc.do_irc_prod_calculations()
                                elif irc_prod_status[0] == 'running' or irc_status[1] == 'running':
                                    continue
                                else:
                                    #read the geometries and try to make products out of them
                                    #verify which of the ircs leads back to the reactant, if any
                                    prod = obj.irc.irc2stationary_pt()
                                    if prod == 0:
                                        logging.info('\t\tNo product found for {}'.format(instance_name))
                                        self.species.reac_ts_done[index] = -999
                                    else:
                                        #IRC's are done
                                        obj.products = prod
                                        obj.product_bonds = prod.bond
                                        self.species.reac_ts_done[index] = 2

                elif self.species.reac_ts_done[index] == 2:
                    #identify bimolecular products and wells
                    fragments, maps = obj.products.start_multi_molecular()
                    obj.products = []
                    for frag in fragments:
                        logging.info('\t\tFor reaction {} product {} is identified, now running optimization.'.format(instance_name, frag.chemid))
                        obj.products.append(frag)
                        self.qc.qc_opt(frag, frag.geom)
                    self.species.reac_ts_done[index] = 3

                elif self.species.reac_ts_done[index] == 3:
                    #do freq calculation if successful
                    for frag in obj.products:
                        fragname = str(frag.chemid) + '_well'
                        stat = self.qc.check_qc(fragname)
                        if stat == 'normal' or 'normal freq':
                            #start freq calculation if not done
                            logging.info('\t\tFor reaction {} product {} is identified, now running frequency calculation.'.format(instance_name, fragname))
                            self.qc.qc_freq(frag, fragname, frag.geom, 0) 
                        elif stat == 'error':
                            self.species.reac_ts_done[index] = -999
                            logging.info('\tProduct optimization failed for {}, product {}'.format(instance_name, fragname))

                    prodstat = 0 # status of all product optimization + freq calculation in this channel
                    for frag in obj.products:
                        fragname = str(frag.chemid) + '_well'
                        if self.qc.check_qc(fragname) != 'normal freq':
                            prodstat = 1
                    if prodstat == 0:
                        # now we can move on to the next step
                        self.species.reac_ts_done[index] = 4

                elif self.species.reac_ts_done[index] == 4:
                    #wait for the freq calculation to finish 
                    err = 0
                    for st_pt in obj.products:
                        chemid = st_pt.chemid
                        orig_geom = copy.deepcopy(st_pt.geom)
                        e, st_pt.geom = self.qc.get_qc_geom(str(st_pt.chemid) + '_well', st_pt.natom)
                        if e < 0:
                            logging.info('\tProduct optimization failed for {}, product {}'.format(instance_name, st_pt.chemid))
                            self.species.reac_ts_done[index] = -999
                            err = -1
                        elif e != 0:
                            err = -1
                        else:
                            e2, st_pt.energy = self.qc.get_qc_energy(str(st_pt.chemid) + '_well')
                            e2, st_pt.zpe = self.qc.get_qc_zpe(str(st_pt.chemid) + '_well')
                            st_pt.bond_mx()
                            st_pt.characterize(0)  # not allowed to use the dimer option here
                            st_pt.calc_chemid()
                            if chemid != st_pt.chemid:
                                # product was optimized to another structure, give warning and remove this reaction
                                logging.info('\tProduct optimizatied to other structure for {}, product {} to {}'.format(instance_name,chemid,st_pt.chemid))
                                self.species.reac_ts_done[index] = -999
                                err = -1
                    if err == 0:
                        self.species.reac_ts_done[index] = 5

                elif self.species.reac_ts_done[index] == 5:
                    # Do the TS and product optimization
                    
                    #make a stationary point object of the ts
                    bond_mx = np.zeros((self.species.natom, self.species.natom), dtype=int)
                    for i in range(self.species.natom):
                        for j in range(self.species.natom):
                            bond_mx[i][j] = max(self.species.bond[i][j],obj.product_bonds[i][j])
                    err, geom = self.qc.get_qc_geom(instance_name, self.species.natom)
                    ts = StationaryPoint(   instance_name, self.species.charge, self.species.mult,
                                            atom = self.species.atom, geom = geom, wellorts = 1)
                    err, ts.energy = self.qc.get_qc_energy(instance_name)
                    err, ts.zpe = self.qc.get_qc_zpe(instance_name)
                    ts.bond = bond_mx
                    ts.find_cycle()
                    ts.find_conf_dihedral()
                    obj.ts = ts
                    #do the ts optimization
                    obj.ts_opt = Optimize(obj.ts, self.par, self.qc)
                    obj.ts_opt.do_optimization()
                    #do the products optimizations
                    for st_pt in obj.products:
                        #check for products of other reactions that are the same as this product
                        #in the case such products are found, use the same Optimize object for both
                        new = 1
                        for i, inst_i in enumerate(self.species.reac_inst):
                            if not i == index:
                                obj_i = self.species.reac_obj[i]
                                if self.species.reac_ts_done[i] > 4:
                                    for j,st_pt_i in enumerate(obj_i.products):
                                        if st_pt_i.chemid == st_pt.chemid:
                                            if len(obj_i.prod_opt) > j:
                                                prod_opt = obj_i.prod_opt[j]
                                                new = 0
                                                break
                        if new:
                            prod_opt = Optimize(st_pt,self.par,self.qc)
                            prod_opt.do_optimization()
                        obj.prod_opt.append(prod_opt)
                    self.species.reac_ts_done[index] = 6

                elif self.species.reac_ts_done[index] == 6:
                    #check up on the TS and product optimizations 
                    opts_done = 1
                    fails = 0
                    #check if ts is done
                    if not obj.ts_opt.shir == 1:
                        opts_done = 0
                        obj.ts_opt.do_optimization()
                    if obj.ts_opt.shigh == -999:
                        fails = 1
                    for pr_opt in obj.prod_opt:
                        if not pr_opt.shir == 1:
                            opts_done = 0
                            pr_opt.do_optimization()
                        if pr_opt.shigh == -999:
                            fails = 1
                    if fails:
                        self.species.reac_ts_done[index] = -999
                    elif opts_done:
                        self.species.reac_ts_done[index] = 7

                elif self.species.reac_ts_done[index] == 7:
                    #Finilize the calculations
                    
                    #continue to PES search in case a new well was found
                    if self.par.par['pes']:
                        #verify if product is monomolecular, and if it is new
                        if len(obj.products) ==1:
                            st_pt = obj.prod_opt[0].species
                            chemid = st_pt.chemid
                            energy = st_pt.energy
                            well_energy = self.species.energy
                            new_barrier_threshold = self.par.par['barrier_threshold'] - (energy-well_energy)*constants.AUtoKCAL
                            dir = os.path.dirname(os.getcwd()) 
                            jobs = open(dir+'/chemids','r').read().split('\n')
                            jobs = [ji for ji in jobs]
                            if not str(chemid) in jobs:
                                #this well is new, add it to the jobs
                                while 1:
                                    try:
                                        #try to open the file and write to it
                                        pes.write_input(self.par,obj.products[0],new_barrier_threshold,dir)
                                        f = open(dir+'/chemids','a')
                                        f.write('{}\n'.format(chemid))
                                        f.close()
                                        break
                                    except IOError:
                                        #wait a second and try again
                                        time.sleep(1)
                                        pass
                                        
                    #check for wrong number of negative frequencies
                    neg_freq = 0
                    for st_pt in obj.products:
                        if any([fi < 0. for fi in st_pt.reduced_freqs]):
                            neg_freq = 1
                    if any([fi < 0. for fi in obj.ts.reduced_freqs[1:]]): 
                        neg_freq = 1
                    
                    if neg_freq:
                        logging.info('\tFound negative frequency for ' + instance_name)
                        self.species.reac_ts_done[index] = -999
                    else:
                        #the reaction search is finished
                        self.species.reac_ts_done[index] = -1 # this is the success code
                        
                        # write a temporary pes input file
                        # remove old xval and im_extent files
                        if os.path.exists('{}_xval.txt'.format(self.species.chemid)):
                            os.remove('{}_xval.txt'.format(self.species.chemid))
                        if os.path.exists('{}_im_extent.txt'.format(self.species.chemid)):
                            os.remove('{}_im_extent.txt'.format(self.species.chemid))
                        postprocess.createPESViewerInput(self.species, self.qc, self.par)
                elif self.species.reac_ts_done[index] == -999:
                    if not self.species.reac_obj[index].instance_name in deleted:
                        if self.par.par['delete_intermediate_files'] == 1:
                            self.delete_files(self.species.reac_obj[index].instance_name)
                            deleted.append(self.species.reac_obj[index].instance_name)
                        
            alldone = 1
            for index, instance in enumerate(self.species.reac_inst):
                if any(self.species.reac_ts_done[i] >= 0 for i in range(len(self.species.reac_inst))):
                    alldone = 1
                    break 
                else: 
                    alldone = 0
            
            # write a small summary while running
            wr = 1
            if wr:
                f_out = open('kinbot_monitor.out','w')
                for index, instance in enumerate(self.species.reac_inst):
                    f_out.write('{}\t{}\t{}\n'.format(self.species.reac_ts_done[index],self.species.reac_step[index],self.species.reac_obj[index].instance_name))
                f_out.close()
            time.sleep(1)
        
        s = []
        for index, instance in enumerate(self.species.reac_inst):
            obj = self.species.reac_obj[index]
            instance_name = obj.instance_name
            # Write a summary on the combinatorial exploration
            if 'combinatorial' in instance_name:
                s.append('NAME\t' + instance_name)
                
                # Write the bonds that were broken and formed
                s.append('BROKEN_BONDS\t' + '\t'.join('[{}, {}]'.format(re[0], re[1]) for re in obj.reac))
                s.append('FORMED_BONDS\t' + '\t'.join('[{}, {}]'.format(pr[0], pr[1]) for pr in obj.prod))
                
                # Populate the ts_bond_lengths dict with the values
                # of this reaction
                if self.species.reac_ts_done[index] == -1:
                    for i in range(self.species.natom - 1):
                        for j in range(i + 1, self.species.natom):
                            if self.species.bond[i][j] != obj.product_bonds[i][j]:
                                if (self.species.bond[i][j] == 0 or
                                        obj.product_bonds[i][j] == 0):
                                    syms = []
                                    syms.append(self.species.atom[i])
                                    syms.append(self.species.atom[j])
                                    syms = ''.join(sorted(syms))
                                    dist = np.linalg.norm(obj.ts.geom[i] - obj.ts.geom[j])
                                    s.append('TS_BOND_LENGTHS\t{}\t{}'.format(syms, dist))
                # write the expected inchis
                s.append('EXPECTED_INCHIS\t' + '\t'.join(inchi for inchi in obj.prod_inchi))
                # get the inchis the reaction found
                if self.species.reac_ts_done[index] == -1:
                    inchis = obj.get_final_inchis()
                    s.append('FOUND_INCHIS\t' + '\t'.join(inchis))
                s.append('\n')
            with open('combinatorial.txt', 'w') as f:
                f.write('\n'.join(s) + '\n')

        logging.info("Reaction generation done!")


    def delete_files(self, name):
        # job names
        names = []
        zf = self.par.par['zf']

        names.append(name)
        names.append(name + '_high')
        names.append(name + '_IRC_F')
        names.append(name + '_IRC_R')
        names.append(name + '_IRC_F_prod')
        names.append(name + '_IRC_R_prod')
        extensions = ['chk', 'py', 'sbatch']
        
        for name in names:
            for ext in extensions:
                # delete file
                file = '.'.join([name, ext])
                try:
                    os.remove(file)
                except FileNotFoundError:
                    pass

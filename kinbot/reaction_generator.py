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
from kinbot import filecopying
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
    
    def __init__(self,species,par,qc):
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
        0: initiate the TS search
        1: check barrier height and errors in TS, and initiates normal mode displacement test, start the irc calculations 
        2: submit product optimization
        3: submit the frequency calculation 
        4: do the optimization of the ts and the products
        5: follow up on the optimizations
        6: finalize calculations, check for wrong number of negative frequencies
        
        If at any times the calculation fails, reac_ts_done is set to -999.
        If all steps are successful, reac_ts_done is set to -1.
        """

        deleted = []
        if len(self.species.reac_inst) > 0:
            alldone = 1
        else: 
            alldone = 0
        
        # status to see of kinbot needs to wait for the product optimizations
        # from another kinbot run, to avoid duplication of calculations
        products_waiting_status = [[] for i in self.species.reac_inst]
        count=0
        for i in self.species.reac_inst:
             count=count+1
        all_unique_prod=[]
        frag_unique=[]
        while alldone:
            for index, instance in enumerate(self.species.reac_inst):
                obj = self.species.reac_obj[index]
                instance_name = obj.instance_name

                # START REACTION SEARCH
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
                            if status == 0 and self.species.reac_ts_freq[index] == 0:  # meaning success and freq was not started
                                # start freq calculation
                                self.qc.qc_freq(instance_name, self.geom)
                                self.species_ts_freq[index] = 1
                            elif status == 0 and self.species.reac_ts_freq[index] == 1:  # meaning success and freq is running
                                status, freq = self.qc.get_qc_freq(instance_name, self.species.natom)
                                if status == 0:
                                    self.species.reac_ts_freq[index] = 2  # won't return here any more
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
                            self.species.reac_step[index] = reac_family.carry_out_reaction(obj, self.species.reac_step[index], self.par.par['qc_command'])
                    
                    else: # do a bond scan
                        if self.species.reac_step[index] == self.par.par['scan_step'] + 1:
                            status = self.qc.get_qc_energy(instance_name, self.species.natom)[0]
                            if status == 0 and self.species.reac_ts_freq[index] == 0:  # meaning success and freq was not started
                                # start freq calculation
                                self.qc.qc_freq(instance_name, self.geom)
                                self.species_ts_freq[index] = 1
                            elif status == 0 and self.species.reac_ts_freq[index] == 1:  # meaning success and freq is running
                                status, freq = self.qc.get_qc_freq(instance_name, self.species.natom)
                                if status == 0:
                                    self.species.reac_ts_freq[index] = 2  # won't return here any more
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
                                self.species.reac_step[index] = reac_family.carry_out_reaction(obj, self.species.reac_step[index], self.par.par['qc_command'])
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
                                        self.species.reac_step[index] = reac_family.carry_out_reaction(obj, self.species.reac_step[index], self.par.par['qc_command'])

                elif self.species.reac_ts_done[index] == 1:
                    status = self.qc.check_qc(instance_name)
                    if status == 'running': continue
                    elif status == 'error': 
                        logging.info('\tRxn search failed (gaussian error) for {}'.format(instance_name))
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
                                    #IRC's are done
                                    obj.products = prod
                                    obj.product_bonds = prod.bond
                                    self.species.reac_ts_done[index] = 2

                elif self.species.reac_ts_done[index] == 2:
                    if len(products_waiting_status[index]) == 0:
                        #identify bimolecular products and wells
                        fragments, maps = obj.products.start_multi_molecular()
                        obj.products = []

                        a=[]
                        for frag in fragments:
                            a.append(frag)
                            if len(frag_unique) == 0:
                                frag_unique.append(frag)
                            elif len(frag_unique) > 0:
                                i=1
                                for fragb in frag_unique:
                                    if frag.chemid == fragb.chemid:
                                        a.pop()
                                        frag=fragb
                                        a.append(frag)
                                        break
                                    if i == len(frag_unique):
                                        frag_unique.append(frag)
                                        break
                                    i=i+1

                        n=1
                        for frag in a:
                            obj.products.append(frag)
                            self.qc.qc_opt(frag, frag.geom)
                            n=n+1

                        #check products make sure they are the same
                        for i, st_pt_i in enumerate(obj.products):
                            for j, st_pt_j in enumerate(obj.products):
                                if st_pt_i.chemid == st_pt_j.chemid and i < j:
                                    obj.products[j] = obj.products[i]



                    #print products generated by IRC
                    products=[] 
                    for st_pt in obj.products: 
                        products.append(st_pt.chemid) 
                    products.append(' ') 
                    products.append(' ') 
                    products.append(' ') 
                    barrier = (self.qc.get_qc_energy(instance_name)[1] - sp_energy) * constants.AUtoKCAL
                    logging.info('\tReaction {} has a barrier of {} and lead to products {} {} {}'.format(instance_name,barrier,products[0],products[1],products[2]))

                    #check geom post optimization
                    obj.products_final = []

                    for st_pt in obj.products:
                        obj.products_final.append(st_pt)
                        chemid = st_pt.chemid
                        orig_geom = copy.deepcopy(st_pt.geom)
                        e, st_pt.geom = self.qc.get_qc_geom(str(st_pt.chemid) + '_well', st_pt.natom)
                        if e < 0:
                            logging.info('\tProduct optimization failed for {}, product {}'.format(instance_name,st_pt.chemid))
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
                                obj.products_final.pop()
                                newfrags, newmaps = st_pt.start_multi_molecular()
                                products_waiting_status[index] = [0 for frag in newfrags]
                                fragChemid=[]
                                for a in newfrags:
                                    for prod in frag_unique:
                                        if a.chemid == prod.chemid:
                                            newfrags.pop()
                                            a=prod
                                            newfrags.append(a)
                                            break
                                    obj.products_final.append(a)
                                    self.qc.qc_opt(a, a.geom)
                                    fragChemid.append(a.chemid)
                                if len(fragChemid) == 1:
                                    fragChemid.append(" ")
                                for i, frag in enumerate(newfrags):
                                    products_waiting_status[index][i] = 1
                                
                                logging.info('\ta) Product optimized to other structure for {}, product {} to {} {}'.format(instance_name,chemid, fragChemid[0], fragChemid[1]))                
                    
                    obj.products=[]

                    for prod in obj.products_final:
                        obj.products.append(prod)
                    
                    obj.products_final=[] 


                    if all([pi == 1 for pi in products_waiting_status[index]]):
                        self.species.reac_ts_done[index] = 3
 
>>>>>>> master
                elif self.species.reac_ts_done[index] == 3:
                    # wait for the optimization to finish 
                    # if two st_pt are the same in the products, we make them exactly identical otherwise
                    # the different ordering of the atoms causes the chemid of the second to be seemingly wrong
                    for i, st_pt_i in enumerate(obj.products):
                        for j, st_pt_j in enumerate(obj.products):
                            if st_pt_i.chemid == st_pt_j.chemid and i < j:
                                obj.products[j] = obj.products[i]
                    err = 0
                    for st_pt in obj.products:
                        chemid = st_pt.chemid
                        orig_geom = copy.deepcopy(st_pt.geom)
                        e, st_pt.geom = self.qc.get_qc_geom(str(st_pt.chemid) + '_well', st_pt.natom)
                        if e < 0:
                            logging.info('\tProduct optimization failed for {}, product {}'.format(instance_name,st_pt.chemid))
                            self.species.reac_ts_done[index] = -999
                            rr = -1
                        elif e != 0:
                            err = -1
                        else:
                            e2, st_pt.energy = self.qc.get_qc_energy(str(st_pt.chemid) + '_well')
                            e2, st_pt.zpe = self.qc.get_qc_zpe(str(st_pt.chemid) + '_well')
                            st_pt.bond_mx()
                            st_pt.characterize(0)  # not allowed to use the dimer option here
                            st_pt.calc_chemid()
                            if chemid != st_pt.chemid:
                                # product was optimized to another structure, give warning but don't remove reaction
                                logging.info('\t b) Product optimized to other structure for {}, product {} to {}'.format(instance_name,chemid,st_pt.chemid))
                                e, st_pt.geom = self.qc.get_qc_geom(str(st_pt.chemid) + '_well', st_pt.natom)
                                if e < 0:
                                    err = -1
                    if err == 0:
                        self.species.reac_ts_done[index] = 4

                elif self.species.reac_ts_done[index] == 4:
                    # Do the TS and product optimization
                    # make a stationary point object of the ts
                    bond_mx = np.zeros((self.species.natom, self.species.natom))
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
                    obj.ts_opt = Optimize(obj.ts,self.par,self.qc)
                    obj.ts_opt.do_optimization()
                    #do the products optimizations
                    for st_pt in obj.products:
					#do the products optimizations
                        #check for products of other reactions that are the same as this product
                        #in the case such products are found, use the same Optimize object for both
                        for i, inst_i in enumerate(self.species.reac_inst):
                            new=1
                            if not i == index:
                                obj_i = self.species.reac_obj[i]
                                if self.species.reac_ts_done[i] > 3:
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

                    for st_pt in obj.products:                        
                    #section where comparing products in same reaction occurs
                        if len(obj.prod_opt) > 0:
                            for j, st_pt_opt in enumerate(obj.prod_opt):
                                if st_pt.chemid == st_pt_opt.species.chemid:
                                    if len(obj.prod_opt) > j:
                                        prod_opt = obj.prod_opt[j]
                                        break

                    elog=open("energy.log",'a')
                    for prod_opt in obj.prod_opt:
                        elog.write("prod_opt: {} |\tenergy: {}\n".format(prod_opt.species.chemid, prod_opt.species.energy))
                    elog.close()

                    self.species.reac_ts_done[index] = 5

                elif self.species.reac_ts_done[index] == 5:
                    #check up on the TS and product optimizations 
                    opts_done = 1
                    fails = 0
                    #check if ts is done
                    if not obj.ts_opt.shir == 1:
                        opts_done = 0
                        obj.ts_opt.do_optimization()
                    if obj.ts_opt.shigh == -999:
                        logging.info("Reaction {} ts_opt_shigh failure".format(instance_name))
                        fails = 1
                    for pr_opt in obj.prod_opt:
                        if not pr_opt.shir == 1:
                            opts_done = 0
                            pr_opt.do_optimization()
                        if pr_opt.shigh == -999:
                            logging.info("Reaction {} pr_opt_shigh failure".format(instance_name))
                            fails = 1
                    if fails:
                        self.species.reac_ts_done[index] = -999
                    elif opts_done:
                        self.species.reac_ts_done[index] = 6

                elif self.species.reac_ts_done[index] == 6:
                    #Finilize the calculations
                    
                    #continue to PES search in case a new well was found
                    if self.par.par['pes']:
                        #verify if product is monomolecular, and if it is new
                        if len(obj.products)==1:
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
           
                        # copy the files of the species to an upper directory
                        frags = obj.products
                        for frag in frags:
                            filecopying.copy_to_database_folder(self.species.chemid, frag.chemid, self.qc)

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
                    if self.par.par['delete_intermediate_files'] == 1:
                        if not self.species.reac_obj[index].instance_name in deleted:
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
                except OSError:
                    pass

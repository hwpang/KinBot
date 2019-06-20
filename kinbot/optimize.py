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
import sys
import os
import copy
import logging
import time

from kinbot import frequencies
from kinbot import geometry
from kinbot import symmetry
from kinbot.conformers import Conformers
from kinbot.hindered_rotors import HIR
from kinbot.molpro import Molpro


class Optimize:
    """
    This class does the following:

    1. Conformational search of the species
    2. High level optimization and freq calc of the species
    3. Hindered rotor scans
    4. Repeat steps 2-3 as long as lower energy structures are found

    TODO: find better name for this module and class
    """

    def __init__(self, species, par, qc, wait=0):
        self.species = species
        self.par = par
        self.qc = qc

        # wait for all calcualtions to finish before returning
        self.wait = wait

        # high level job name
        if self.species.wellorts:
            self.job_high = self.species.name + '_high'
            self.job_hir = 'hir/' + self.species.name + '_hir_'
        else:
            self.job_high = str(self.species.chemid) + '_well_high'
            self.job_hir = 'hir/' + str(self.species.chemid) + '_hir_'

        # status of the various parts
        # -1: not yet started
        #  0: running
        #  1: finished
        # -999:failed
        self.scycconf = -1
        self.sconf = -1
        self.shigh = -1
        self.shir = -1

        # restart counter: number of times the high-level and hir calculations
        # has been restarted in case a lower energy structure has been found
        self.restart = 0
        # maximum restart count
        self.max_restart = par.par['rotation_restart']

    def do_optimization(self):
        while 1:
            # do the conformational search
            if self.par.par['conformer_search'] == 1:
                if self.scycconf == -1 and self.sconf == -1:
                    # conformational analysis has to be started
                    logging.info('\tStarting conformational search of {}'.format(self.species.name))
                    self.species.confs = Conformers(self.species, self.par, self.qc)

                # first do the cyclic part of the molecule
                if self.scycconf == -1:
                    # start the ring conf search
                    if len(self.species.cycle_chain) > 0:
                        # there are rings in the molecule, do a search
                        self.species.confs.generate_ring_conformers(copy.deepcopy(self.species.geom))
                        # set the cyclic conf status to running
                        self.scycconf = 0
                    else:
                        # there are no rings in the molecule, continue from the current geometry
                        self.species.confs.cyc_conf_geoms = [copy.deepcopy(self.species.geom)]
                        # no ring conf search has to be done, set status to finished
                        self.scycconf = 1
                if self.scycconf == 0:
                    # ring conf search is running, check if finished
                    status, self.species.confs.cyc_conf_geoms = self.species.confs.check_ring_conformers()
                    if status:
                        # ring conf search is finished
                        self.scycconf = 1
                # do the open chain par of the molecule
                if self.scycconf == 1:
                    # do open chain part if cyclic part is done
                    if self.sconf == -1:
                        # open chain part has not started yet
                        for geom in self.species.confs.cyc_conf_geoms:
                            # take all the geometries from the cyclic part
                            # generate the conformers for the current geometry
                            self.species.confs.generate_conformers(0, geom)
                        # set conf status to running
                        self.sconf = 0
                    if self.sconf == 0:
                        # conformational search is running
                        # check if the conformational search is done
                        status, geom = self.species.confs.check_conformers(wait=self.wait)
                        if status == 1:
                            # conf search is done
                            # save lowest energy conformer as species geometry
                            self.species.geom = geom
                            # set conf status to finished
                            self.sconf = 1
            else:
                # no conf search necessary, set status to finished
                self.sconf = 1
            if self.sconf == 1:  # conf search is finished
                while self.restart < self.max_restart:
                    # do the high level calculations
                    if self.par.par['high_level'] == 1:
                        if self.shigh == -1:
                            # high level calculation did not start yet
                            logging.info('\tStarting high level optimization of {}'.format(self.species.name))
                            if self.species.wellorts:
                                # do the high level optimization of a ts
                                self.qc.qc_opt_ts(self.species, self.species.geom, high_level=1)
                            else:
                                # do the high level optimization of a well
                                self.qc.qc_opt(self.species, self.species.geom, high_level=1)
                            self.shigh = 0  # set the high status to running
                        if self.shigh == 0:
                            # high level calculation is running
                            # check if it is finished
                            status = self.qc.check_qc(self.job_high)
                            if status == 'error':
                                # found an error
                                logging.info('\tHigh level optimization failed for {}'.format(self.species.name))
                                self.shigh = -999
                            if status == 'normal':
                                # finished successfully
                                err, new_geom = self.qc.get_qc_geom(self.job_high, self.species.natom, wait=self.wait)
                                if geometry.equal_geom(self.species.bond, self.species.geom, new_geom, 0.1):
                                    # geometry is as expected
                                    err, self.species.geom = self.qc.get_qc_geom(self.job_high, self.species.natom)
                                    err, self.species.energy = self.qc.get_qc_energy(self.job_high)
                                    err, self.species.freq = self.qc.get_qc_freq(self.job_high, self.species.natom)
                                    err, self.species.zpe = self.qc.get_qc_zpe(self.job_high)
                                    self.shigh = 1
                                else:
                                    # geometry diverged to other structure
                                    logging.info('\tHigh level ts optimization converged to different structure for {}'.format(self.species.name))
                                    self.shigh = -999
                    else:
                        # no high-level calculations necessary, set status to finished
                        self.shigh = 1
                    if self.shigh == 1:
                        # do the HIR calculation
                        if self.par.par['rotor_scan'] == 1:
                            if self.shir == -1:
                                # hir not stated yet
                                logging.info('\tStarting hindered rotor calculations of {}'.format(self.species.name))
                                self.species.hir = HIR(self.species, self.qc, self.par)
                                self.species.hir.generate_hir_geoms(copy.deepcopy(self.species.geom))
                                self.shir = 0
                            if self.shir == 0:
                                # hir is running
                                # check if it is done:
                                status = self.species.hir.check_hir(wait=self.wait)
                                if status:
                                    if len(self.species.hir.hir_energies) > 0:
                                        # check if along the hir potential a structure was found with a lower energy
                                        min = self.species.hir.hir_energies[0][0]
                                        min_rotor = -1
                                        min_ai = -1
                                        for rotor in range(len(self.species.dihed)):
                                            for ai in range(self.species.hir.nrotation):
                                                # use a 0.1kcal/mol cutoff for numerical noise
                                                if self.species.hir.hir_energies[rotor][ai] < min - 1.6E-4:
                                                    min = self.species.hir.hir_energies[rotor][ai]
                                                    min_rotor = rotor
                                                    min_ai = ai
                                        if min_rotor > -1:
                                            self.restart += 1
                                            if self.restart < self.max_restart:
                                                # lower energy structure found
                                                logging.info('\t\tLower energy found during hindered rotor scan for {}'.format(self.species.name))
                                                logging.info('\t\tRestart number: ' + str(self.restart))
                                                logging.info('\t\tRotor: ' + str(min_rotor))
                                                logging.info('\t\tScan point: ' + str(min_ai))
                                                job = self.job_hir + str(min_rotor) + '_' + str(min_ai).zfill(2)

                                                err, self.species.geom = self.qc.get_qc_geom(job, self.species.natom)
                                                # delete the high_level log file and the hir log files
                                                if os.path.exists(self.job_high + '.log'):
                                                    # logging.info("\t\t\tRemoving file " + self.job_high + '.log')
                                                    os.remove(self.job_high + '.log')
                                                for rotor in range(len(self.species.dihed)):
                                                    for ai in range(self.species.hir.nrotation):
                                                        if os.path.exists(self.job_hir + str(rotor) + '_' + str(ai).zfill(2) + '.log'):
                                                            # logging.info("\t\t\tRemoving file " + self.job_hir + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                                                            os.remove(self.job_hir + str(rotor) + '_' + str(ai).zfill(2) + '.log')
                                                # set the status of high and hir back to not started
                                                self.shigh = -1
                                                self.shir = -1
                                            else:
                                                logging.info('\t\tLower energy found, but readched max restart for {}'.format(self.species.name))
                                                self.shir = 1
                                        else:
                                            self.shir = 1
                                    else:
                                        self.shir = 1
                        else:
                            # no hir calculations necessary, set status to finished
                            self.shir = 1
                    if not self.wait or self.shir == 1 or self.shigh == -999:
                        # break the loop if no waiting is required or
                        # if the hir calcs are done or
                        # if the high level calc failed
                        break
                    else:
                        time.sleep(1)
            if self.shir == 1:
                # finilize if no waiting is required or if the hir calcs are done
                # calculate the symmetry numbers
                symmetry.calculate_symmetry(self.species)

                # calculate the new frequencies with the internal rotations projected out
                fr_file = self.species.name
                if not self.species.wellorts:
                    fr_file += '_well'
                if self.par.par['high_level']:
                        fr_file += '_high'
                hess = self.qc.read_qc_hess(fr_file, self.species.natom)
                self.species.kinbot_freqs, self.species.reduced_freqs = frequencies.get_frequencies(self.species, hess, self.species.geom)

                # write the molpro input and read the molpro energy, if available
                if self.par.par['single_point_qc'] == 'molpro':
                    molp = Molpro(self.species, self.par)
                    molp.create_molpro_input()
                    molp.create_molpro_submit()
                    status, molpro_energy = molp.get_molpro_energy()
                    if status:
                        self.species.energy = molpro_energy
                
                # delete unnecessary files
                delete = 1
                if delete:
                    self.delete_files()
            if self.wait:
                if self.shir == 1 or self.shigh == -999:
                    return 0
                time.sleep(1)
            else:
                return 0
    
    def delete_files(self):
        # job names
        names = []
        zf = self.par.par['zf']
        if self.species.wellorts:
            names.append(self.species.name)
            names.append(self.species.name + '_high')
            names.append(self.species.name + '_IRC_F')
            names.append(self.species.name + '_IRC_R')
            names.append(self.species.name + '_IRC_F_prod')
            names.append(self.species.name + '_IRC_R_prod')
            
            if self.par.par['high_level'] == 1:
                for count in range(self.species.hir.nrotation):
                    for rot_num in range(self.par.par['nrotation']):
                        names.append('hir/' + self.species.name + '_hir_' + str(count) + '_' + str(rot_num).zfill(2))
            if self.par.par['conformer_search'] == 1:
                for count in range(self.species.confs.conf):
                    names.append('conf/' + self.species.name + '_' + str(count).zfill(zf))
                for count in range(self.species.confs.cyc_conf):
                    for num in range(self.species.conf.cyc_conf_index[count]):
                        names.append('conf/' + self.species.name + '_r' + str(count).zfill(zf) + '_' + str(num).zfill(zf))
        else:
            names.append(str(self.species.chemid) + '_well')
            names.append(str(self.species.chemid) + '_well_high')
            if self.par.par['high_level'] == 1:
                for count in range(self.species.hir.nrotation):
                    for rot_num in range(self.par.par['nrotation']):
                        names.append('hir/' + str(self.species.chemid) + '_hir_' + str(count) + '_' + str(rot_num).zfill(2))
            if self.par.par['conformer_search'] == 1:
                for count in range(self.species.confs.conf):
                    names.append('conf/' + str(self.species.chemid) + '_' + str(count).zfill(zf))
                for count in range(self.species.confs.cyc_conf):
                    for num in range(self.species.conf.cyc_conf_index[count]):
                        names.append('conf/' + self.species.name + '_r' + str(count).zfill(zf) + '_' + str(num).zfill(zf))
        extensions = ['chk', 'py', 'sbatch']
        
        for name in names:
            for ext in extensions:
                # delete file
                file = '.'.join([name, ext])
                # print(file)
                try:
                    os.remove(file)
                except FileNotFoundError:
                    pass



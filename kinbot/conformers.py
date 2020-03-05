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
from __future__ import division
import random
import time
import copy
import logging
import numpy as np

from shutil import copyfile
from kinbot import geometry
from kinbot import zmatrix


class Conformers:
    """
    Class that does all the steps for the conformers of one species
    """
    def __init__(self, species, par, qc):
        """
        species: instance of StationaryPoint
        qc: instance of QuantumChemistry
        """

        self.species = species
        self.qc = qc

        # status of the conformational analysis
        # -1: not yet started
        # 0: running
        # 1: finished
        # -999:failed
        self.scycconf = -1
        self.sconf = -1

        # final geometries of the cyclic search
        self.cyc_conf_geoms = []

        # for each conformer, check the index of its progress
        self.cyc_conf_index = []
        # dihedral atoms of the cyclic conformers
        self.cyc_dih_atoms = []
        # dihedral values of the cyclic conformers
        self.cyc_dih_values = []
        # number of cyclic conformers generated
        self.cyc_conf = 0

        # number of open chain conformers generated
        self.conf = 0

        # -1 (not finished), 0 (successful) or
        # 1 (failed) for each cyclic conformer
        self.cyc_conf_status = []
        # -1 (not finished), 0 (successful) or
        # 1 (failed) for each open chain conformer
        self.conf_status = []
        self.zf = par.par['zf']

        # Maximum number of diherals for which exhaustive
        # conformation searches are done
        self.max_dihed = par.par['max_dihed']
        # Number of random conformers in case no
        # exhaustive search is done
        self.nconfs = par.par['random_conf']


    def generate_ring_conformers(self, cart):
        """
        Generate the conformers of a cyclic structure
        by randomly sampling the dihedrals of the ring
        """

        # iterate the different rings in the species
        for cyc in self.species.cycle_chain:
            if len(cyc) > 3:  # three membered rings don't have conformers
                dihs = []  # list of the ring dihedrals
                for i, at in enumerate(cyc):
                    dihs.append([cyc[i-3], cyc[i-2], cyc[i-1], cyc[i]])

                # define the flatness of the ring by the sum of the
                # absolute values of the dihedrals along the ring
                # divided by the number of atoms in the ring
                cycdih = 0.
                for dih in dihs:
                    val = geometry.calc_dihedral(cart[dih[0]], cart[dih[1]],
                                                 cart[dih[2]], cart[dih[3]])[0]
                    cycdih += np.abs(val)
                cycdih /= float(len(cyc))

                # randomly select N-3 dihedrals,
                # with N the number of dihedrals in the ring
                random_dihs = random.sample(dihs, len(dihs) - 3)
                # number of independent dihedrals
                nd = len(dihs) - 3
                # number of conformers for this ring:
                nc = np.power(3, nd)
                for i in range(nc):
                    self.cyc_dih_atoms.append(random_dihs)
                    # values the dihedrals will be modified to
                    values = []
                    for j in range(nd):
                        values.append(cycdih*(np.mod(i // np.power(3, j), 3) - 1))
                    self.cyc_dih_values.append(values)
                    self.cyc_conf_index.append(-1)
                    self.cyc_conf += 1
            else:
                self.cyc_conf_geoms.append(copy.deepcopy(cart))

        for ci in range(self.cyc_conf):
            self.start_ring_conformer_search(ci, copy.deepcopy(self.species.geom))


    def start_ring_conformer_search(self, index, cart):
        """
        index: number of the conformer
        """

        if self.cyc_conf_index[index] == len(self.cyc_dih_atoms[index]) - 1:
            # this conformer has finished
            return 0
        else:
            self.cyc_conf_index[index] += 1
            fix = []
            change = []
            for j, da in enumerate(self.cyc_dih_atoms[index]):
                if j == self.cyc_conf_index[index]:
                    new_dih = self.cyc_dih_values[index][j]
                    change.append([da[0] + 1, da[1] + 1, da[2] + 1, da[3] + 1, new_dih])
                    break
                else:
                    fix.append([da[0] + 1, da[1] + 1, da[2] + 1, da[3] + 1])
            for i in range(self.species.natom - 1):
                for j in range(i+1, self.species.natom):
                    if self.species.bond[i][j] > 0:
                        fix.append([i + 1, j + 1])
            self.qc.qc_ring_conf(self.species, cart, fix, change, index, self.cyc_conf_index[index])
        return 1


    def test_ring_conformer(self, index):
        """
        Test whether a conformer has the same bond matrix as the original structure.
        Returns the conformer object and -1 if not yet finished, 0 if same, and 1 if not.
        """

        if self.species.wellorts:
            job = 'conf/' + self.species.name + '_r' + str(index).zfill(self.zf) + '_' + str(self.cyc_conf_index[index]).zfill(self.zf)
        else:
            job = 'conf/' + str(self.species.chemid) + '_r' + str(index).zfill(self.zf) + '_' + str(self.cyc_conf_index[index]).zfill(self.zf)

        status, geom = self.qc.get_qc_geom(job, self.species.natom)
        if status == 1:  # still running
            return np.zeros((self.species.natom, 3)), -1
        elif status == -1:  # conformer search failed
            logging.debug('Conformer search failed for scan point {}'.format(job))
            return np.zeros((self.species.natom, 3)), 1
        else:
            if self.start_ring_conformer_search(index, geom):
                logging.debug('Running the next dihedral for conformer {}'.format(job))
                return geom, -1
            else:
                # check if all the bond lenghts are withing 10% or the original bond lengths
                if geometry.equal_geom(self.species.bond, self.species.geom, geom, 0.10):
                    logging.debug('Successfullly finished conformer {}'.format(job))
                    return geom, 0
                else:
                    logging.debug('Conformer too far from original structure {}'.format(job))
                    return np.zeros((self.species.natom, 3)), 1


    def check_ring_conformers(self, wait=0):
        """
        Check if the conformer optimizations finished.
        Test them, and submit frequency calculations.
        Then select the lowest energy one.

        returns:
        *status: 0 if still running, 1 if finished
        *geometries of all the conformers

        wait: wait for all the conformer calculations to finish before returning anything
        """

        if len(self.cyc_conf_status) < self.cyc_conf:
            for i in range(self.cyc_conf):
                self.cyc_conf_status.append(-1)
        while 1:
            # check if conformational search is finished
            for i, si in enumerate(self.cyc_conf_status):
                if si == -1:
                    self.cyc_conf_status[i] = self.test_ring_conformer(i)[1]
            if all([si >= 0 for si in self.cyc_conf_status]):
                geoms = [self.species.geom]  # list used for intermediate ring conformer generation
                final_geoms = []
                for ci in range(self.cyc_conf):
                    si = self.cyc_conf_status[ci]
                    if si == 0:  # this is a valid confomer
                        if self.species.wellorts:
                            job = 'conf/' + self.species.name + '_r' + str(ci).zfill(self.zf) + '_' + str(self.cyc_conf_index[ci]).zfill(self.zf)
                        else:
                            job = 'conf/' + str(self.species.chemid) + '_r' + str(ci).zfill(self.zf) + '_' + str(self.cyc_conf_index[ci]).zfill(self.zf)
                        err, geom = self.qc.get_qc_geom(job, self.species.natom)
                        geoms.append(geom)
                        final_geoms.append(geom)
                    else:
                        final_geoms.append(np.zeros((self.species.natom, 3)))
                self.write_profile(self.cyc_conf_status, final_geoms, [0 for gi in final_geoms], ring=1)
                return 1, geoms
            else:
                if wait:
                    time.sleep(1)
                else:
                    return 0, np.zeros((self.species.natom, 3))


    def generate_conformers(self, rotor, cart):
        """
        Generate guesses for all of the canonical conformers.
        This is a recursive routine to generate them.
        rotor: the rotor number in the order it was discovered
        """

        if len(self.species.conf_dihed) > self.max_dihed:
            self.generate_conformers_random_sampling(cart)
            return 0

        if rotor == len(self.species.conf_dihed):  # the end point of the recursion
            self.qc.qc_conf(self.species, cart, self.conf, freq=0)  # TODO add MP2 flags!!!
            self.conf += 1
            return 0

        # TODO: it would be nice to write this in a loop, and then the user could
        # control how many divisions to take in the 360 degree space
        # low priority, since these 3 steps are adequate for most cases
        cart = np.asarray(cart)
        zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(self.species, rotor, cart, 1)

        rotor += 1
        cart0 = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
        self.generate_conformers(rotor, cart0)

        zmat[3][2] += 120.
        for i in range(4, self.species.natom):
            if zmat_ref[i][2] == 4:
                zmat[i][2] += 120.
            if zmat_ref[i][2] == 1:
                zmat[i][2] += 120.
        cart1 = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
        self.generate_conformers(rotor, cart1)

        zmat[3][2] += 120.
        for i in range(4, self.species.natom):
            if zmat_ref[i][2] == 4:
                zmat[i][2] += 120.
            if zmat_ref[i][2] == 1:
                zmat[i][2] += 120.
        cart2 = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
        self.generate_conformers(rotor, cart2)

        return 0


    def generate_conformers_random_sampling(self, ini_cart):
        """
        Generate a random sampling of each dihedral for a number nconfs of conformers
        """

        for ni in range(self.nconfs):
            cart = copy.deepcopy(ini_cart)
            if ni == 0:
                sample = [0. for di in self.species.conf_dihed]
            else:
                sample = [random.choice([0., 120., 240.]) for di in self.species.conf_dihed]
            for rotor in range(len(self.species.conf_dihed)):
                zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(self.species, rotor, cart, 1)
                zmat[3][2] += sample[rotor]
                for i in range(4, self.species.natom):
                    if zmat_ref[i][2] == 4:
                        zmat[i][2] += sample[rotor]
                    if zmat_ref[i][2] == 1:
                        zmat[i][2] += sample[rotor]
                cart = zmatrix.make_cart_from_zmat(zmat, zmat_atom, zmat_ref, self.species.natom, self.species.atom, zmatorder)
            self.qc.qc_conf(self.species, cart, self.conf, freq=0)  # TODO also add MP2 flag
            self.conf += 1

        return 0


    def test_conformer(self, conf, freq=0):
        """
        Test whether a conformer has the same bond matrix as the original structure.
        Returns 
        -1 if not yet finished, 
        0 if same,
        10 if not
        1 if freq running
        11 if freq success
        12 if freq fail
        and the geometry.
    
        if freq=1 then it is looking at the frequencies

        """

        if self.species.wellorts:
            job = 'conf/' + self.species.name + '_' + str(conf).zfill(self.zf)
        else:
            job = 'conf/' + str(self.species.chemid) + '_' + str(conf).zfill(self.zf)

        status, geom = self.qc.get_qc_geom(job, self.species.natom)
        if status == 1:  # still running
            if freq:
                return np.zeros((self.species.natom, 3)), 1
            else:
                return np.zeros((self.species.natom, 3)), -1
        elif status == -1:  # conformer search failed
            return np.zeros((self.species.natom, 3)), 10
        else:  # normal termination
            if freq:
                status, zpe = self.qc.get_qc_zpe(job)
                status, frequencies = self.qc.get_qc_freq(job)
                if self.wellorts:
                    if frequencies[0] > 0. or frequencies[1] < 0.:
                        return np.zeros((self.species.natom, 3)), 10  # wrong number for imag freq
                elif frequencies[0] < 0.:
                    return np.zeros((self.species.natom, 3)), 10  # wrong number for imag freq
                else:
                    return geom, 0
            else:
                # check if all the bond lenghts are withing 10% or the original bond lengths
                if geometry.equal_geom(self.species.bond, self.species.geom, geom, 0.10):
                    return geom, 0
                else:
                    return np.zeros((self.species.natom, 3)), 10


    def check_conformers(self, wait=0):
        """
        Check if the conformer optimizations finished.
        Test them, and submit frequency calculations.
        Then select the lowest energy one.

        returns:
        *status: 0 if still running, 1 if finished
        *geometry of lowest energy conformer

        wait: wait for all the conformer calculations to finish before returning anything
        """

        # TODO why is this like this?
        if len(self.conf_status) < self.conf:
            for i in range(len(self.conf_status), self.conf):
                self.conf_status.append(-1)
        status = self.conf_status

<<<<<<< HEAD
        geom = [np.zeros((self.species.natom, 3)) for st in status] # create placeholder for geoms
=======
        lowest_conf = str(0).zfill(self.zf) # the index of the lowest conf, to be updated as we go
>>>>>>> master

        while 1:
            # check if conformational search is finished
            for i, si in enumerate(status):
                if si == -1:
                    geom[i], status[i] = self.test_conformer(i)
            for i, si in enumerate(status):
                if si == 0:  # was a success but there is not frequency yet
                    # submit a frequency calculation
                    self.qc.qc_conf(self.species, geom[i], i, freq=1)  # TODO add MP2 flags!!!
                    status[i] = 1  # meaning that frequency calculation has started
            for i, si in enumerate(status):
                if si == 1:
                    status[i] = self.test_conformer(i, freq=1)[1]

            if all([si >= 10 for si in status]):  # all are done
                lowest_energy = self.species.energy + self.species.zpe
                lowest_e_geom = self.species.geom
                lowest_job = ''  # name of the job that is the lowest energy conformer
                final_geoms = []  # list of all final conformer geometries
                energies = []
                for ci in range(self.conf):
                    si = status[ci]
                    if si == 0:  # this is a valid confomer
                        if self.species.wellorts:
                            job = 'conf/' + self.species.name + '_' + str(ci).zfill(self.zf)
                        else:
                            job = 'conf/' + str(self.species.chemid) + '_' + str(ci).zfill(self.zf)
                        err, energy = self.qc.get_qc_energy(job)
                        err, zpe = self.qc.get_qc_zpe(job)
                        err, geom = self.qc.get_qc_geom(job, self.species.natom)
                        final_geoms.append(geom)
                        energies.append(energy + zpe)
                        if energy < lowest_energy:
                            lowest_conf = str(ci).zfill(self.zf) 
                            lowest_energy = energy
                            lowest_e_geom = geom
                            lowest_job = job
                    else:
                        energies.append(0.)
                        final_geoms.append(np.zeros((self.species.natom, 3)))
               
                self.write_profile(status, final_geoms, energies)
<<<<<<< HEAD
                # copy the lowest energy geometry to the main directory
                copyfile('{}.log'.format(job), '../{}.log'.format(job))
                # the original file is not overwritten, but the geometry is updated in KinBot's memory
                return 1, lowest_e_geom
=======
                return 1, lowest_conf, lowest_e_geom, lowest_energy

>>>>>>> master
            else:
                if wait:
                    time.sleep(1)
                else:
                    return 0, lowest_conf, np.zeros((self.species.natom, 3)), self.species.energy


    def write_profile(self, status, final_geoms, energies, ring=0):
        """
        Write a molden-readable file with the CONF analysis (geometries and energies)
        """

        r = ''
        if ring:
            r = 'r'
        if self.species.wellorts:
            file = open('conf/' + self.species.name + r + '.xyz', 'w')
        else:
            file = open('conf/' + str(self.species.chemid) + r + '.xyz', 'w')
        for i, st in enumerate(status):
            s = str(self.species.natom) + '\n'
            s += 'energy = ' + str(energies[i]) + '\n'
            for j, at in enumerate(self.species.atom):
                x, y, z = final_geoms[i][j]
                s += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            if st == 0:  # valid conformer:
                file.write(s)
        file.close()

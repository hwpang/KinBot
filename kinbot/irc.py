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
import numpy as np
import os
import copy
import time
import logging
from shutil import copyfile
import pkg_resources

from kinbot.stationary_pt import StationaryPoint


class IRC:
    """
    Class to run the IRC's for one specific reaction
    """
    def __init__(self, rxn, par):
        # instance of the reac_family object
        # the family this reaction belongs to
        self.rxn = rxn
        self.par = par

    def irc2stationary_pt(self):
        """
        Read the irc files
        There are three possible scenarios:
        1. One of the ircs leads the initial well and
           the other to another well or bimolecular product
        2. Neither of the ircs lead to the inital well,
           transition state structure is not the one
           kinbot was looking for
        3. Both the ircs lead to the initial well,
           KinBot found either an identical reaction
           or the ts is not correct
        """
        instance_name = self.rxn.instance_name

        directions = ['Forward', 'Reverse']

        ini_well_hits = 0
        prod_hit = -1
        st_pts = [-1, -1]
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}_prod'.format(instance_name, direction[0])
            err, geom = self.rxn.qc.get_qc_geom(irc_name,
                                                self.rxn.species.natom,
                                                allow_error=1)
            if err == -1:
                return 0
            if self.problem_in_geom(geom):
                # this happens seldomly that all the atoms are
                # very close to one another (problem in Gaussian)
                logging.info('\tProblem with product geometry for {}'.format(instance_name))
                return 0

            temp = StationaryPoint(irc_name,
                                   self.rxn.species.charge,
                                   self.rxn.species.mult,
                                   atom=self.rxn.species.atom,
                                   geom=geom)
            temp.calc_chemid()

            st_pts[i] = temp
            if temp.chemid == self.rxn.species.chemid:
                ini_well_hits += 1
            else:
                prod_hit = i

        if ini_well_hits == 0:
            logging.info('\tNeither IRC leads to the well for {}'.format(instance_name))
            return 0
        elif ini_well_hits == 2:
            logging.info('\tBoth IRCs lead to the well, identical reaction found: {}'.format(instance_name))
            return 0
        else:
            # ircs OK: well and product found
            logging.info('\tIRCs successful for {}'.format(instance_name))
            return st_pts[prod_hit]

    def problem_in_geom(self, geom):
        # check if interatomic distances are closer than 0.3 Angstrom
        for i in range(len(geom)):
            for j in range(i+1, len(geom)):
                dist = np.linalg.norm(geom[i] - geom[j])
                if dist < 0.3:
                    return 1
        return 0

    def check_irc(self):
        instance_name = self.rxn.instance_name
        directions = ['Forward', 'Reverse']
        status = [-1, -1]
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])
            status[i] = self.rxn.qc.check_qc(irc_name)
        return status

    def do_irc_calculations(self):
        """
        Carry out the IRC calculation.
        """
        instance_name = self.rxn.instance_name
        err, geom = self.rxn.qc.get_qc_geom(instance_name,
                                            self.rxn.species.natom)
        directions = ['Forward', 'Reverse']
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])
            if self.rxn.qc.qc == 'gauss':
                # copy the chk file
                if os.path.exists(instance_name + '.chk'):
                    copyfile(instance_name + '.chk', irc_name + '.chk')

            if self.rxn.qc.qc == 'nwchem' and direction == 'Reverse':
                direction = 'Backward'

            odft = self.rxn.species.mult > 1

            if direction == 'Forward':
                if rxn.scan or 'R_Addition_MultipleBond' in rxn.instance_name:
                    self.rxn.qc.assemble_ase_template(irc_name, 'ircfmp2', self.rxn.species, geom, 0, self.par.par['sella'])
                else:
                    self.rxn.qc.assemble_ase_template(irc_name, 'ircf', self.rxn.species, geom, 0, self.par.par['sella'])
            else:
                if rxn.scan or 'R_Addition_MultipleBond' in rxn.instance_name:
                    self.rxn.qc.assemble_ase_template(irc_name, 'ircrmp2', self.rxn.species, geom, 0, self.par.par['sella'])
                else:
                    self.rxn.qc.assemble_ase_template(irc_name, 'ircr', self.rxn.species, geom, 0, self.par.par['sella'])

        return 0


    def check_irc_prod(self):
        instance_name = self.rxn.instance_name
        directions = ['Forward', 'Reverse']
        status = [-1, -1]
        for i, direction in enumerate(directions):
            irc_prod_name = '{}_IRC_{}_prod'.format(instance_name, direction[0])
            status[i] = self.rxn.qc.check_qc(irc_prod_name)
        return status


    def do_irc_prod_calculations(self):
        """
        Carry out the product optimization for the IRCs in both directions
        """
        instance_name = self.rxn.instance_name
        directions = ['Forward', 'Reverse']
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])
            err, geom = self.rxn.qc.get_qc_geom(irc_name,
                                                self.rxn.species.natom, wait=0, allow_error=1)
            irc_prod_name = '{}_IRC_{}_prod'.format(instance_name, direction[0])
            if self.rxn.qc.qc == 'gauss':
                # copy the chk file
                try: 
                    copyfile(irc_name + '.chk', irc_prod_name + '.chk')
                except:
                    logging.info('\tIRC checkpoint file is not present for {}.'.format(irc_name))

            if self.rxn.qc.qc == 'nwchem' and direction == 'Reverse':
                direction = 'Backward'

            odft = self.rxn.species.mult > 1

            if rxn.scan or 'R_Addition_MultipleBond' in rxn.instance_name:
                self.rxn.qc.assemble_ase_template(irc_prod_name, 'prodircmp2', self.rxn.species, geom, 0, self.par.par['sella'])
            else:
                self.rxn.qc.assemble_ase_template(irc_prod_name, 'prodirc', self.rxn.species, geom, 0, self.par.par['sella'])

        return 0

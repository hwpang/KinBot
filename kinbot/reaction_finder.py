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
# -*- coding: utf-8 -*-
import numpy as np
import sys
import os
import copy
import time
import logging

from kinbot import bond_combinations
from kinbot import find_motif
from kinbot.reac_Cyclic_Ether_Formation import CyclicEtherFormation
from kinbot.reac_Diels_alder_addition import DielsAlder
from kinbot.reac_Intra_Diels_alder_R import IntraDielsAlder
from kinbot.reac_12_shift_S_F import S12ShiftF
from kinbot.reac_12_shift_S_R import S12ShiftR
from kinbot.reac_cpd_H_migration import CpdHMigration
from kinbot.reac_intra_H_migration import IntraHMigration
from kinbot.reac_intra_H_migration_suprafacial import IntraHMigrationSuprafacial
from kinbot.reac_intra_OH_migration import IntraOHMigration
from kinbot.reac_Intra_R_Add_Endocyclic_F import IntraRAddEndocyclicF
from kinbot.reac_Intra_R_Add_Exocyclic_F import IntraRAddExocyclicF
from kinbot.reac_Intra_R_Add_ExoTetCyclic_F import IntraRAddExoTetCyclicF
from kinbot.reac_intra_R_migration import IntraRMigration
from kinbot.reac_Retro_Ene import RetroEne
from kinbot.reac_r22_cycloaddition import R22Cycloaddition
from kinbot.reac_r12_insertion_R import R12Insertion
from kinbot.reac_r13_insertion_RSR import R13InsertionRSR
from kinbot.reac_r13_insertion_ROR import R13InsertionROR
from kinbot.reac_r13_insertion_CO2 import R13InsertionCO2
from kinbot.reac_r12_cycloaddition import R12Cycloaddition
from kinbot.reac_R_Addition_MultipleBond import RAdditionMultipleBond
from kinbot.reac_R_Addition_CSm_R import RAdditionCS
from kinbot.reac_R_Addition_COm3_R import RAdditionCO
from kinbot.reac_Korcek_step2 import KorcekStep2
from kinbot.reac_ketoenol import KetoEnol
from kinbot.reac_Intra_RH_Add_Exocyclic_R import IntraRHAddExoR
from kinbot.reac_Intra_RH_Add_Exocyclic_F import IntraRHAddExoF
from kinbot.reac_Intra_RH_Add_Endocyclic_R import IntraRHAddEndoR
from kinbot.reac_Intra_RH_Add_Endocyclic_F import IntraRHAddEndoF
from kinbot.reac_HO2_Elimination_from_PeroxyRadical import HO2Elimination

from kinbot.reac_combinatorial import Combinatorial

from kinbot.reac_birad_recombination_F import BiradRecombinationF
from kinbot.reac_birad_recombination_R import BiradRecombinationR
from kinbot.reac_Intra_disproportionation_R import IntraDisproportionationR
from kinbot.reac_Intra_disproportionation_F import IntraDisproportionationF
from kinbot.reac_r14_birad_scission import R14BiradScission
from kinbot.reac_r14_cyclic_birad_scission_R import R14CyclicBiradScission

class ReactionFinder:
    """
    Class to find all the potential reactions starting from a well
    """
    
    def __init__(self,species,par,qc):
        self.species = species
        self.qc = qc
        self.par = par
        self.families = par.par['families']
        self.specific_reaction = par.par['specific_reaction']
        self.break_bond = par.par['break_bonds']
        self.form_bond = par.par['form_bonds']

        self.one_reaction_comb = par.par['one_reaction_comb']
        self.one_reaction_fam = par.par['one_reaction_fam']
        # make a set of frozen sets from the breaking and forming bond lists
        self.reac_bonds = set()
        for i, bond in enumerate(par.par['break_bonds']):
            self.reac_bonds.add(frozenset(par.par['break_bonds'][i]))
        self.prod_bonds = set()
        for i, bond in enumerate(par.par['form_bonds']):
            self.prod_bonds.add(frozenset(par.par['form_bonds'][i]))

        #keys: names of the families
        #values: list of instances
        #this dict is used to keep track of the unique reactions found,
        #and to verify whether a new reaction is indeed unique 
        self.reactions = {}

    def find_reactions(self):
        """
        List all reaction types available, and find the key atoms for them 
        for the current structure.
        """

        atom = self.species.atom
        natom = self.species.natom
        
        for i, bond in enumerate(self.species.bonds):
            rad = self.species.rads[i]
            
            if self.one_reaction_comb:
                # search for just one reaction, given by the list of bonds to be 
                # broken or formed
                
                # based on the combinatorial reaction family, because they are also
                # defined by the list of bonds to be broken or formed
                name = 'combinatorial'
                self.reactions[name] = []
                
                self.reac_bonds = self.par.par['break_bonds']
                self.prod_bonds = self.par.par['form_bonds']
                ts = bond_combinations.generate_ts(self.reac_bonds, self.prod_bonds, self.species.bond)
                self.reactions[name].append([self.reac_bonds, self.prod_bonds, ts, 1])
                
            else:
                if 'intra_H_migration' in self.families or 'all' in self.families:
                    self.search_intra_H_migration(natom,atom,bond,rad)
                    
                if 'intra_H_migration_suprafacial' in self.families or 'all' in self.families:
                    self.search_intra_H_migration_suprafacial(natom,atom,bond,rad)
                    
                if 'intra_R_migration' in self.families or 'all' in self.families:
                    self.search_intra_R_migration(natom,atom,bond,rad)
                    
                if 'intra_OH_migration' in self.families or 'all' in self.families:
                    self.search_intra_OH_migration(natom,atom,bond,rad)
                    
                if 'cpd_H_migration' in self.families or 'all' in self.families:
                    self.search_cpd_H_migration(natom,atom,bond,rad)
                    
                if 'Intra_RH_Add_Endocyclic_F' in self.families or 'all' in self.families:
                    self.search_Intra_RH_Add_Endocyclic_F(natom,atom,bond,rad)
                    
                if 'Intra_RH_Add_Endocyclic_R' in self.families or 'all' in self.families:
                    self.search_Intra_RH_Add_Endocyclic_R(natom,atom,bond,rad)
                    
                if 'Cyclic_Ether_Formation' in self.families or 'all' in self.families:
                    self.search_Cyclic_Ether_Formation(natom,atom,bond,rad)
                    
                if 'Intra_RH_Add_Exocyclic_F' in self.families or 'all' in self.families:
                    self.search_Intra_RH_Add_Exocyclic_F(natom,atom,bond,rad)
                    
                if 'Intra_RH_Add_Exocyclic_R' in self.families or 'all' in self.families:
                    self.search_Intra_RH_Add_Exocyclic_R(natom,atom,bond,rad)
                    
                if 'Retro_Ene' in self.families or 'all' in self.families:
                    self.search_Retro_Ene(natom,atom,bond,rad)
                    
                if 'Intra_R_Add_Endocyclic_F' in self.families or 'all' in self.families:
                    self.search_Intra_R_Add_Endocyclic_F(natom,atom,bond,rad)
                    
                if 'Intra_R_Add_ExoTetCyclic_F' in self.families or 'all' in self.families:
                    self.search_Intra_R_Add_ExoTetCyclic_F(natom,atom,bond,rad)
                    
                if 'Intra_R_Add_Exocyclic_F' in self.families or 'all' in self.families:
                    self.search_Intra_R_Add_Exocyclic_F(natom,atom,bond,rad)
                    
                if 'Korcek_step2' in self.families or 'all' in self.families:
                    self.search_Korcek_step2(natom,atom,bond,rad)
                    
                if 'r22_cycloaddition' in self.families or 'all' in self.families:
                    self.search_r22_cycloaddition(natom,atom,bond,rad)
                    
                if 'r12_cycloaddition' in self.families or 'all' in self.families:
                    self.search_r12_cycloaddition(natom,atom,bond,rad)
                    
                if 'r12_insertion_R' in self.families or 'all' in self.families:
                    self.search_r12_insertion_R(natom,atom,bond,rad)
                    
                if 'r13_insertion_CO2' in self.families or 'all' in self.families:
                    self.search_r13_insertion_CO2(natom,atom,bond,rad)
                    
                if 'r13_insertion_ROR' in self.families or 'all' in self.families:
                    self.search_r13_insertion_ROR(natom,atom,bond,rad)
                    
                if 'Diels_alder_addition' in self.families or 'all' in self.families:
                    self.search_Diels_alder_addition(natom,atom,bond,rad)
                    
                if 'Intra_Diels_alder_R' in self.families or 'all' in self.families:
                    self.search_Intra_Diels_alder_R(natom,atom,bond,rad)
                    
                if 'ketoenol' in self.families or 'all' in self.families:
                    self.search_ketoenol(natom,atom,bond,rad)
                    
                if 'HO2_Elimination_from_PeroxyRadical' in self.families or 'all' in self.families:
                    self.search_HO2_Elimination_from_PeroxyRadical(natom,atom,bond,rad)
                    
                if 'R_Addition_COm3_R' in self.families or 'all' in self.families:
                    self.search_R_Addition_COm3_R(natom,atom,bond,rad)
                    
                if 'R_Addition_MultipleBond' in self.families or 'all' in self.families:
                    self.search_R_Addition_MultipleBond(natom,atom,bond,rad)
                    
                if '12_shift_S_F' in self.families or 'all' in self.families:
                    self.search_12_shift_S_F(natom,atom,bond,rad)
                    
                if '12_shift_S_R' in self.families or 'all' in self.families:
                    self.search_12_shift_S_R(natom,atom,bond,rad)
                    
                if 'R_Addition_CSm_R' in self.families or 'all' in self.families:
                    self.search_R_Addition_CSm_R(natom,atom,bond,rad)
                    
                if 'r13_insertion_RSR' in self.families or 'all' in self.families:
                    self.search_r13_insertion_RSR(natom,atom,bond,rad)

                if 'combinatorial' in self.families:
                    self.search_combinatorial(natom,atom,bond,rad)
            
            
            #if 'birad_recombination_F' in self.families or 'all' in self.families:
            #    self.search_birad_recombination_F(natom,atom,bond,rad)
            #if 'birad_recombination_R' in self.families or 'all' in self.families:
            #    self.search_birad_recombination_R(natom,atom,bond,rad)
            #if 'Intra_disproportionation_F' in self.families or 'all' in self.families:
            #    self.search_Intra_disproportionation_F(natom,atom,bond,rad)
            #if 'Intra_disproportionation_R' in self.families or 'all' in self.families:
            #    self.search_Intra_disproportionation_R(natom,atom,bond,rad)
            #if 'r14_birad_scission' in self.families or 'all' in self.families:
            #    self.search_r14_birad_scission(natom,atom,bond,rad)
            #if 'r14_cyclic_birad_scission_R' in self.families or 'all' in self.families:
            #    self.search_r14_cyclic_birad_scission_R(natom,atom,bond,rad)

        
        for name in self.reactions:
            self.reaction_matrix(self.reactions[name], name) 
        
        #verify if every name is unique
        for index in range(len(self.species.reac_name)-1):
            if self.species.reac_name[index] in self.species.reac_name[index+1:]:
                logging.error('Found reaction name "%s" more than once'%self.species.reac_name[index])
                logging.error('Exiting')
                sys.exit()

        #write the reactions that were found to the log
        logging.info('\tFound the following reactions:')
        for rxn in self.species.reac_name:
            logging.info('\t\t%s'%rxn)
        
        return 0  
    
    def search_combinatorial(self, natom, atom, bond, rad):
        """ 
        This is a method to create all possible combinations of maximum 3 bond breakings 
        and maximum 3 bond formations.
        TODO: allow bond breaking without the atoms forming new bond (only for radicals)
        """
        
        name = 'combinatorial'

        if not name in self.reactions:
            self.reactions[name] = []

        instances = bond_combinations.generate_all_product_bond_matrices(self.species, self.par)
        for inst in instances:
            self.reactions[name].append(inst)
        #~ self.reactions[name] = []
        #~ reac = [[0, 5], [1, 2], [3, 4]]
        #~ prod = [[0, 1], [2, 3], [4, 5]]
        #~ ts = self.species.bond
        #~ self.reactions[name].append([reac, prod, ts])

        return 0


    def search_intra_H_migration(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        H-R~~~~~~~R* <==> R*~~~~~~~R-H

        Find all unique cases for ring sizes between 3 and 9. Works in both directions.
        H is moved to
        * radical site
        * multiple bond
        * lone pair
        """
        
        name = 'intra_H_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        if np.sum(rad) == 0: 
            #find H-migrations over double bonds and to lone pairs
            
            for ringsize in range(3, 9):
                # double bonds 
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
                for instance in instances:
                    if any([bi > 1 for bi in bond[instance[0]]]):
                        rxns += [instance]
                # lone pairs
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
                for instance in instances:
                    if (self.species.atom[instance[0]] == 'O' or  
                       self.species.atom[instance[0]] == 'S' or 
                       self.species.atom[instance[0]] == 'N'):
                        rxns += [instance]

        else:
            instances = []
            for ringsize in range(3, 9):
                motif = ['X' for i in range(ringsize)]
                motif[-1] = 'H'
                for rad_site in np.nonzero(rad)[0]:
                    instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            for instance in instances: 
                rxns.append(instance)
        
        #filter for the same reactions
        for inst in rxns:
            new = 1
            # first filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    def search_intra_H_migration_suprafacial(self, natom, atom, bond, rad):
        """ 
        This is a special case of H migration reactions over a double bond 
        (keto-enol type) that proceeds through a suprafacial instead of the
        common antrafacial TS.
        """
        
        name = 'intra_H_migration_suprafacial'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        # search for keto-enol type reactions
        motif = ['X', 'X', 'X', 'H']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        
        # filter for the double bond
        for instance in instances:
            if bond[instance[0]][instance[1]] == 2:
                rxns += [instance]

        for inst in rxns:
            new = 1
            # first filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    
    def search_intra_R_migration(self, natom, atom, bond, rad):
        """ 
        This is an class that covers several RMG classes.
        
        R cannot be an H, this is already taken care of in the intra_H_migration
        
        TODO: merge this with intra H migration families?
        yes, because it is the same rule
        no, because then it's hard to search for just one of the types
        TODO: this should also include migration to lone pair electrons?
        currently it moves atoms to radical sites only
        """
        
        if np.sum(rad) != 1: return 

        name = 'intra_R_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        instances = []
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize)]
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        for instance in instances: 
            if not atom[instance[-1]] == 'H':
                rxns.append(instance)
        
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    def search_cpd_H_migration(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        H-C1-C=C-C=C-1 <==> C1=C-C=C-C(-H)-1

        """
        
        if not any([len(ci) == 5 for ci in self.species.cycle_chain]) : return
        
        name = 'cpd_H_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        bondsum = 0
        
        for cycle in self.species.cycle_chain:
            if len(cycle) == 5:
                for index, atomi in enumerate(cycle):
                    if index < 4:
                        atomj = cycle[index + 1]
                    else:
                        atomj = cycle[0]
                    if index == 0:
                        atomk = cycle[-1]
                    else:
                        atomk = cycle[index - 1]
                    bondsum += bond[atomi][atomj]
                    if bond[atomi][atomj] == 1 and bond[atomi][atomk] == 1:
                        start = atomi
                        startindex = index
                if bondsum != 7: return # exactly two double bonds
                ring_forw = np.ndarray.tolist(np.roll(cycle, 5 - startindex))
                ring_rev = ring_forw[::-1] # look at the ring in the reverse direction for an H-shift to the other side
                ring_rev = np.ndarray.tolist(np.roll(ring_rev, 1))
                rings = [ring_forw,ring_rev]
                
                Hatomi = -1
                for atomi in range(natom):
                    if atom[atomi] == 'H':
                        if bond[atomi][start] == 1:
                            Hatomi = atomi
                if Hatomi > -1:
                    for ring in rings:
                        instance = ring[:]
                        instance.append(Hatomi)
                        rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                # TODO need to check if this is correct
                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0
        
    

    def search_intra_OH_migration(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R*~~~~~~~O-OH <==> HOR~~~~~~~O*

        Find all unique cases for ring sizes between 3 and 9. The H atom is not counted in the cycle size but has to be there.
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'intra_OH_migration'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            instances = []
            # forward direction
            motif = ['X' for i in range(ringsize+1)]
            motif[-1] = 'H'
            motif[-2] = 'O'
            motif[-3] = 'O'
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            # reverse direction
            motif = ['X' for i in range(ringsize+1)]
            motif[-1] = 'H'
            motif[-2] = 'O'
            motif[0] = 'O'
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            for ins in instances:
                rxns.append(ins)

        for case in range(len(rxns)):
            rxns[case] = rxns[case][:-1] #cut off H
            
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-3], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0



    def search_Intra_RH_Add_Endocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                                  H
                                  | 
        H-R~~~~~~~R=R ==> R~~~~~~~R-R
                          |         |
                           ---------

        Find all unique cases for ring sizes between 3 and 9. This is for the forward direction.
        """
        
        if np.sum(rad) != 0: return
        if len(self.species.cycle_chain) > 0: return
        
        name = 'Intra_RH_Add_Endocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize + 1)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            bondpattern = ['X' for i in range(ringsize)]
            bondpattern[0] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 
            

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-2] == instance[-2] and len(inst) == len(instance):
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]}), frozenset({inst[-1], inst[1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
                
        return 0
        


    def search_Intra_RH_Add_Endocyclic_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                H
                | 
        R~~~~~~~R-R ==> H-R~~~~~~~R=R
        |         |
         ---------

        Find all unique cases for ring sizes between 3 and 9. This is for the reverse direction.
        """
        
        if len(self.species.cycle_chain) == 0: return
        if np.sum(rad) != 0: return
        
        name = 'Intra_RH_Add_Endocyclic_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ci in self.species.cycle_chain:
            motif = ['X' for i in range(len(ci) + 1)]
            motif[-1] = 'H'
            
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            # check if there is a bond between the first and second to last atom
            for instance in instances:
                if bond[instance[0]][instance[-2]] > 0:
                    rxns += [instance[-4:]]
        
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-3], inst[-4]}), frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[-1], inst[-4]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        self.reactions[name]

        return 0



    def search_Cyclic_Ether_Formation(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R*~~~~~~~O-OR ==> R~~~~~~~O + OR
                          |_______|

        Find all unique cases for ring sizes between 3 and 9. The OR groups are not counted in the cycle size but have to be there.
        Only the forward direction is included.
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'Cyclic_Ether_Formation'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(4, 10):
            motif = ['X' for i in range(ringsize)]
            motif[-2] = 'O'
            motif[-3] = 'O'
            motif[0] = 'C'
            for rad_site in np.nonzero(rad)[0]:
                rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        for instance in range(len(rxns)):
            rxns[instance] = rxns[instance][:-2] #cut off OR
            
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-2], inst[-3]})} or self.prod_bonds != {frozenset({inst[0], inst[-3]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_Intra_R_Add_Endocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~~~~R=R ==> R~~~~~~~~R*-R
                          |___________|

        """
        
        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_Endocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize)]
            instances = []
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            bondpattern = ['X' for i in range(ringsize-1)]
            bondpattern[-1] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]
            
            bondpattern[-1] = 3
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
                
        return 0


    def search_Intra_R_Add_ExoTetCyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~~~~R-R ==> R~~~~~~~~R + R*
                          |________|

        """

        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_ExoTetCyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize + 1)]
            for rad_site in np.nonzero(rad)[0]:
                rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst) 

        return 0


    def search_Intra_R_Add_Exocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~~~~R=R ==> R~~~~~~~~R-R*
                          |________|

        """
        
        if np.sum(rad) == 0: return

        name = 'Intra_R_Add_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize + 1)]
            instances = []
            for rad_site in np.nonzero(rad)[0]:
                instances += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
            bondpattern = ['X' for i in range(ringsize)]
            bondpattern[-1] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]
                    
            bondpattern[-1] = 3
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]


        #filter for the same reactions
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-2]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst) 

        return 0


    def search_Intra_RH_Add_Exocyclic_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        The general scheme is:

        H-R~~~~~~R=R ==> R~~~~~~R-R-H
                         |      |
                          ------

        The special case of this reaction is Korcel_step1:
              
                R        R   OH 
        R      /          \ /
         \ / \C=O          C
          |       ==>     / \
          O   H          |   O
           \ /          / \ /
            O          R   O

        Implemented as:

                                   --O--O--
                                  |        |
        O=C~~~~~~~~C-O-O-H ==> HO-C~~~~~~~~C
          |                       |
          R                       R

        Find all unique cases for final ring sizes between 3 and 9. The carbonyl dangling R and the
        tail H are included, but are not counted as the ring size, but these two atoms are kept
        because they are needed in the geometry manipulation step.
        """
        
        if len(self.species.cycle_chain) > 0: return
        
        name = 'Intra_RH_Add_Exocyclic_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize+2)]
            motif[-1] = 'H'

            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            bondpattern = ['X' for i in range(ringsize+1)]
            bondpattern[0] = 2
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[1], inst[-2]}), frozenset({inst[-1], inst[0]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_Intra_RH_Add_Exocyclic_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                                    H
                                    | 
        H-R~~~~~~~R=R <== R~~~~~~~R-R
                          |_______|

        Find all unique cases for ring sizes between 3 and 9. This is for the reverse direction.
        """
        
        name = 'Intra_RH_Add_Exocyclic_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer


        if len(self.species.cycle_chain) == 0: return
        if np.sum(rad) != 0: return

        for ci in self.species.cycle_chain:
            motif = ['X' for i in range(len(ci) + 2)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            # check if there is a bond between the first and second to last atom
            for instance in instances:
                if bond[instance[0]][instance[-3]] > 0:
                    rxns += [instance[-4:]]


        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]}), frozenset({inst[-3], inst[0]})} or self.prod_bonds != {frozenset({inst[-1], inst[0]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0


    def search_Retro_Ene(self, natom, atom, bond, rad):
        """ 
        This is not an RMG class.

        R-R-R-R=R ==> R=R + R=R-R

        """
        
        
        if np.sum(rad) != 0: return
        
        name = 'Retro_Ene'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X' for i in range(6)]
        motif[-1] = 'H'
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        bondpattern = ['X' for i in range(5)]
        bondpattern[0] = 2
        for instance in instances:
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance] 
        
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0

    def search_Korcek_step2_judit(self, natom, atom, bond, rad):
        """ 
        This is an RMG class, but that only has the 5-mem ring case (?).
              
        Implemented for 4, 5, and 6-membered rings (the R groups are not necessarily identical):

        4-membered ring:

            --O--O--
           |        |
        HO-C--------C-R  ==> RCOOH + R2CO
           |        |
           R        R

        5-membered ring:

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

        6-membered ring:

            ----O---O----
           |             |
        HO-C---C- ---C---C-R  ==> RCOOH + C2R2 + R2CO 
           |  / \   / \  |
           R R   R R   R R

        FIXME: need to generalize to larger structures
        Only the forward direction is included.
        
        This family is not currently used, superseded by the next one

        """
        
        
        name = 'Korcek_step2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #  reactions found with the current resonance isomer
        ring_var = [] #  a helper variable to temporarily mark the ring size within this function
        
        for ringsize in range(6, 8):
            motif = ['C' for i in range(ringsize)]
            motif[0] = 'O'
            motif[1] = 'O'
            motif[-1] = 'H'
            motif[-2] = 'O'

            for atomi in range(natom):
                if atom[atomi] == 'O':
                    for atomj in range(natom):
                        if atom[atomj] == 'O':
                            if bond[atomi][atomj] == 1:
                                korcek_chain =  find_motif.start_motif(motif, natom, bond, atom, atomi, self.species.atom_eqv)
                                for case in range(len(korcek_chain)):
                                    if bond[korcek_chain[case][0]][korcek_chain[case][-3]] == 1:
                                        for ringbond in range(len(korcek_chain[0]) - 2 - 3): # FIXME, assuming just one Korcek hit
                                            marked_chain = korcek_chain[case] + [ringbond + 2] # mark the first atom of the second bond fission
                                            rxns += [marked_chain]
                                            ring_var.append(ringsize)

        for n, inst in enumerate(rxns):
            new = 1
            #filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if ring_var[n] == 6: 
                    if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[-3], inst[-4]})} or self.prod_bonds != {frozenset()}:
                        new = 0
                if ring_var[n] == 7: 
                    #  TODO, incomplete
                    if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[-3], inst[-4]})} or self.prod_bonds != {frozenset()}:
                        new = 0
                if ring_var[n] == 8: 
                    if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]}), frozenset({inst[-3], inst[-4]})} or self.prod_bonds != {frozenset()}:
                        new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_Korcek_step2(self, natom, atom, bond, rad):
        """ 
        Generalized Korcek step 
        
        The 4 membered ring equals a 2,2 cycloaddition and is not considered here (no H shift involved)
        
        The 5 membered ring proceeds through a 6 membered transition state (including a 1,2 H migration):

            --O--O--
           |        |
        HO-C---C----C-R  ==> RCOOH + R3CC(R)O
           |  / \   |
           R R   R  R

        6-membered ring: TODO

        Only the forward direction is included.

        """
        
        
        name = 'Korcek_step2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        ring_var = [] #  a helper variable to temporarily mark the ring size within this function

        for ringsize in range(5, 6):
            motif = ['X' for i in range(ringsize + 1)]
            #motif[-1] = 'H'  #  deleted because atom types are no longer checked
            korcek_chain =  find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            for ins in korcek_chain:
                if bond[ins[0]][ins[-2]] == 1:
                    rxns += [ins]
                    ring_var.append(ringsize)

        for n, inst in enumerate(rxns):
            new = 1
            #filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if ring_var[n] == 7: 
                    if (not {frozenset({inst[-2], inst[-3]}), frozenset({inst[0], inst[1]})}.issubset(self.reac_bonds)) or self.prod_bonds != {frozenset()}:
                        new = 0
                if ring_var[n] == 8: 
                    #  TODO this is an incomplete check
                    if self.reac_bonds != {frozenset({inst[-2], inst[-3]}), frozenset({inst[-4], inst[-5]}), frozenset({inst[0], inst[1]})}:
                        new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_r22_cycloaddition(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R      R         R---R
        ||  +  ||  <==   |   |
        R      R         R---R

        N.B.: only the reverse direction is available. Also, the 3 related RMG classes are treated as one.

        """
        
        
        name = 'r22_cycloaddition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        if not any([len(ci) == 4 for ci in self.species.cycle_chain]): return 
        
        for ci in self.species.cycle_chain:
            if len(ci) == 4:
                # there are two ways to slice a 4-mem ring
                ring1 = ci
                ring2 = np.ndarray.tolist(np.roll(ring1, 1))

                # FIXME only works for 1 cycle
                rxns += [ring1]
                rxns += [ring2]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                # TODO need to make sure that these are the bonds that are broken, see the reaction details
                if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0



    def search_r12_cycloaddition(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                       R--R
        R=R + R: <==   \  /
                        R 

        N.B.: only the reverse direction is available. 

        """
        
        
        name = 'r12_cycloaddition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        if not any([len(ci) == 3 for ci in self.species.cycle_chain]): return 
        
        for ci in self.species.cycle_chain:
            if len(ci) == 3:
                # there are three ways to slice a 3-mem ring
                ring1 = self.species.cycle_chain
                ring2 = np.ndarray.tolist(np.roll(ring1, 1))
                ring3 = np.ndarray.tolist(np.roll(ring1, 2))

                # FIXME only works for 1 cycle
                rxns += ring1
                rxns += ring2
                rxns += ring3 

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                # TODO need to make sure that these are the bonds that are broken, see the reaction details
                if self.reac_bonds != {frozenset({inst[0], inst[2]}), frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0



    def search_r12_insertion_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                          X
                          |
        X-P + R-R <==   R-P-R

        """
        
        #if np.sum(rad) != 0: return
        
        name = 'r12_insertion_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X','X','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        
        for instance in instances:
            #if all([atom[atomi] != 'H' for atomi in instance]):
            rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[1], inst[2]})}) or self.prod_bonds != {frozenset({inst[0], inst[2]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_r13_insertion_CO2(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

                          O
                          ||
        O=C=O + R-R <== R-C-O-R


        """
        
        #if np.sum(rad) != 0: return
        
        name = 'r13_insertion_CO2'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        motif = ['X','C','O','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        for instance in instances:
            for atomi in range(natom):
                if not atomi in instance:
                    if atom[atomi] == 'O':
                        if bond[atomi][instance[1]] == 2:
                            rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_r13_insertion_ROR(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R1-O-R2 + R=R <== R1-R-R-O-R2


        """
        
        #if np.sum(rad) != 0: return
        name = 'r13_insertion_ROR'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X','X','O']
        rxns = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})}) or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_Diels_alder_addition(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

          R                  R
        //                 /   \
        R       R         R     R
        |  +    ||  <==   ||    |
        R       R         R     R
         \\                \   /
           R                 R

        N.B.: only the reverse direction is available. 

        """
        
        
        name = 'Diels_alder_addition'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        if not any([len(ci) == 6 for ci in self.species.cycle_chain]): return 

        for ci in self.species.cycle_chain:
            if len(ci) == 6:
                bondsum = 0
                for index, atomi in enumerate(ci):
                    if index < 5:
                        atomj = ci[index + 1]
                    else:
                        atomj = ci[0]
                    bondsum += bond[atomi][atomj]
                    if bond[atomi][atomj] == 2:
                        start = atomi
                        startindex = index
                if bondsum != 7: return # exactly one double bond
                ring = np.ndarray.tolist(np.roll(ci, 6 - startindex))

                rxns += [ring] # FIXME only works for 1 cycle

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != set({frozenset({inst[2], inst[3]}), frozenset({inst[4], inst[5]})}) or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0



    def search_Intra_Diels_alder_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.
        TODO it seems like this is the forward reaction, but the naming is confusing.

                             C
                            / \\
                           C   C
        C=C-C=C~~~C=C <==  |   |
                           C   C
                            \ //
                              C

        """
        
        name = 'Intra_Diels_alder_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(3, 9):  # TODO what is the meaning of these larger rings?
            motif = ['X' for i in range(ringsize + 4)]
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

            bondpattern = ['X' for i in range(ringsize + 3)]
            bondpattern[0] = 2
            bondpattern[2] = 2
            bondpattern[-1] = 2

            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    #inst = instance[:4] + instance[-2:]
                    rxns += [instance] 
     

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
                
        return 0



    def search_ketoenol(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R=R-O-R <==> R-R-R=O
        """
        
        name = 'ketoenol'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        # enol to keto
        motif = ['C', 'C', 'O', 'X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        # keto to enol
        motif = ['O', 'C', 'C', 'X']
        instances += find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        bondpattern = [2, 'X', 'X', 'X']
        for instance in instances:
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                rxns += [instance]
            
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if (inst[0] == instance[0] and inst[1] == instance[1]
                    and inst[2] == instance[2] and inst[3] == instance[3]):
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[2], inst[3]})} or self.prod_bonds != {frozenset({inst[0], inst[1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0
 


    def search_HO2_Elimination_from_PeroxyRadical(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        H-R-R-O-O* ==> R=R + HO2

        N.B.: only the forward direction is available.
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'HO2_Elimination_from_PeroxyRadical'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['H', 'X', 'X', 'O', 'O']
        rxns += find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[2], inst[3]})}) or self.prod_bonds != {frozenset({inst[0], inst[4]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0
        

    def search_R_Addition_COm3_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        C#O + R* <== R-C*=O

        N.B.: only the reverse direction is available. 
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_COm3_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'C', 'O']

        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        for instance in instances:
            bondpattern = [1, 2]
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                if rad[instance[1]] == 1:
                    rxns += [instance]
        
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0


        
    def search_R_Addition_MultipleBond(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R=R + R* <== R*-R-R

        N.B.: only the reverse direction is available. 
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_MultipleBond'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'X', 'X']
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
    
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_12_shift_S_F(self, natom, atom, bond, rad):
        """
        This is an RMG class.
        """

        if np.sum(rad) != 1: return
        
        name = '12_shift_S_F'
        
        if not name in self.reactions:
            self.reactions[name] = []
        
        motif = ['X','S','X']
        rxns = []
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)

        #filter for identical reactions
        for inst in rxns:
            new = 1
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        return 0


    def search_12_shift_S_R(self, natom, atom, bond, rad):
        """
        This is an RMG class.

        C-S-R* <== *S-R-C

        TODO: why not forward??

        """
        
        if np.sum(rad) != 1: return
        
        name = '12_shift_S_R'
        
        if not name in self.reactions:
            self.reactions[name] = []
        
        motif = ['S','X','X']
        rxns = []
        for rad_site in np.nonzero(rad)[0]:
            rxns += find_motif.start_motif(motif, natom, bond, atom, rad_site, self.species.atom_eqv)
        
        for inst in rxns:
            new = 1
            # filter for identical reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset({inst[0], inst[2]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        return 0


    def search_r13_insertion_RSR(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R-S-R + R1=R2 <== R-R1-R2-S-R


        """
        
        #if np.sum(rad) != 0: return
        name = 'r13_insertion_RSR'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X','X','S']
        rxns = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != set({frozenset({inst[0], inst[1]}), frozenset({inst[0], inst[1]})}) or self.prod_bonds != {frozenset({inst[0], inst[3]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst) 
        
        return 0


    def search_R_Addition_CSm_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        C#S + R* <== R-C*=S

        N.B.: only the reverse direction is available. 
        """
        
        
        if np.sum(rad) == 0: return
        
        name = 'R_Addition_CSm_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X', 'C', 'S']

        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)

        for instance in instances:
            bondpattern = [1, 2]
            if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                if rad[instance[1]] == 1:
                    rxns += [instance]
        
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)

        return 0


    def search_r14_birad_scission(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        It is now renamed to 1,4_Linear_birad_scission on the RMG website,

        *R-R-R-R* ==> R=R + R=R

        Problematic reaction because of the biradical character.

        """

        if np.sum(rad) != 2: return
        
        
        name = 'r14_birad_scission'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        
        motif = ['X','X','X','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        for instance in instances: 
            if rad[instance[0]] == 1 and rad[instance[-1]] == 1:
                rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[1] == instance[1] and inst[2] == instance[2]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[1], inst[2]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_r14_cyclic_birad_scission_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        R1-R*~~~~~~R*-R2   <==  R1=R~~~~~~R=R2
        |______________|
        (this one bond)

        TODO forward?

        """

        if np.sum(rad) != 0: return
        
        name = 'r14_cyclic_birad_scission_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            bondpattern[-1] = 2
            
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 
                    
        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_birad_recombination_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~~~~R* ==> R~~~~~~~~R
                         |________|

        """

        if np.sum(rad) != 2: return
        
        name = 'birad_recombination_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(3, 9):
            motif = ['X' for i in range(ringsize)]
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
            for instance in instances: 
                if rad[instance[0]] == 1 and rad[instance[-1]] == 1:
                    rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset()} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_birad_recombination_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~~~~R* <== R~~~~~~~~R
                         |________|

        """

        if np.sum(rad) != 0: return
        if len(self.species.cycle_chain) == 0: return
        
        name = 'birad_recombination_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        motif = ['X','X']
        instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
        for instance in instances: 
            if instance[0] in self.cycle and instance[1] in self.cycle :
                rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[1] == instance[1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[0], inst[1]})} or self.prod_bonds != {frozenset()}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_Intra_disproportionation_F(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~R*-R-H ==> H-R~~~~~R=R

        """

        if np.sum(rad) != 2: return
        
        name = 'Intra_disproportionation_F'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer

        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
           
            for instance in instances: 
                if rad[instance[0]] == 1 and rad[instance[-3]] == 1:
                    rxns += [instance]

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def search_Intra_disproportionation_R(self, natom, atom, bond, rad):
        """ 
        This is an RMG class.

        *R~~~~~R*-R-H <== H-R~~~~~R=R

        """

        if np.sum(rad) != 0: return
        
        name = 'Intra_disproportionation_R'
        
        if not name in self.reactions:
            self.reactions[name] = []

        rxns = [] #reactions found with the current resonance isomer
        
        for ringsize in range(5, 9):
            motif = ['X' for i in range(ringsize)]
            motif[-1] = 'H'
            
            instances = find_motif.start_motif(motif, natom, bond, atom, -1, self.species.atom_eqv)
            
            bondpattern = ['X' for i in range(ringsize - 1)]
            bondpattern[0] = 2
            
            for instance in instances:
                if find_motif.bondfilter(instance, bond, bondpattern) == 0:
                    rxns += [instance] 

        for inst in rxns:
            new = 1
            # filter for the same reactions
            for instance in self.reactions[name]:
                if inst[0] == instance[0] and inst[-1] == instance[-1]:
                    new = 0
            # filter for specific reaction after this
            if self.one_reaction_fam and new:
                if self.reac_bonds != {frozenset({inst[-1], inst[-2]})} or self.prod_bonds != {frozenset({inst[0], inst[-1]})}:
                    new = 0
            if new:
                self.reactions[name].append(inst)
        
        return 0


    def reaction_matrix(self, reac_list, reac_id):
        """ 
        Create arrays to store all reactions for species.
        input: 
        reac_list: atom motifs from individual searches
        reac_id: reaction name (e.g., HO2_Elimination_from_PeroxyRadical) from individual searc functions
        Every reaction type just makes the below arrays longer, generated as reactions are found.

        generated:
        reac_type: reaction class identifier
        reac_inst: reaction instance defined by the important atoms
        reac_step: the step at which the search is at
        reac_scan_energy: for each reaction the energy as a function of steps, only used for scanning type searches, e.g. R_Addition_MultipleBond
        rec_ts_done: the last calculations is submitted in the sequence
        reac_ts_geom: the geometry of the TS
        reac_ts_freq: the freqencies of the TS
        reac_name: the base name of the file to run - created for each reaction later
        """
        
        
        self.species.reac_type += [reac_id for i in range(len(reac_list))]
        self.species.reac_inst += reac_list
        self.species.reac_step += [0 for i in range(len(reac_list))]
        self.species.reac_scan_energy += [[] for i in range(len(reac_list))]
        self.species.reac_ts_done += [0 for i in range(len(reac_list))] 
        self.species.reac_ts_geom += [0 for i in range(len(reac_list))]
        self.species.reac_ts_freq += [0 for i in range(len(reac_list))]
        
        for i in range(len(reac_list)):
            if reac_id == 'intra_H_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraHMigration(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'intra_H_migration_suprafacial':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraHMigrationSuprafacial(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'intra_R_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRMigration(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'intra_OH_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraOHMigration(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'cpd_H_migration':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1) + '_' + str(reac_list[i][-2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(CpdHMigration(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_RH_Add_Endocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(len(reac_list[i])) + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddEndoF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_RH_Add_Endocyclic_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddEndoR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Cyclic_Ether_Formation':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(CyclicEtherFormation(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_RH_Add_Exocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddExoF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_RH_Add_Exocyclic_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRHAddExoR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Retro_Ene':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RetroEne(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_R_Add_Endocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRAddEndocyclicF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_R_Add_ExoTetCyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRAddExoTetCyclicF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_R_Add_Exocyclic_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraRAddExocyclicF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Korcek_step2':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(KorcekStep2(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r22_cycloaddition':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R22Cycloaddition(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r12_cycloaddition':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R12Cycloaddition(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r12_insertion_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R12Insertion(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r13_insertion_CO2':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R13InsertionCO2(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r13_insertion_ROR':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R13InsertionROR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r14_birad_scission':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R14BiradScission(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r14_cyclic_birad_scission_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R14CyclicBiradScission(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'birad_recombination_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BiradRecombinationF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'birad_recombination_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(BiradRecombinationR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_disproportionation_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraDisproportionationF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_disproportionation_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraDisproportionationR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Diels_alder_addition':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(DielsAlder(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'Intra_Diels_alder_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(IntraDielsAlder(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'ketoenol':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1)  + '_' + str(reac_list[i][2] + 1)  + '_' + str(reac_list[i][3] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(KetoEnol(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'HO2_Elimination_from_PeroxyRadical':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][-1] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(HO2Elimination(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'R_Addition_COm3_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RAdditionCO(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'R_Addition_MultipleBond':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RAdditionMultipleBond(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == '12_shift_S_F':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(S12ShiftF(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == '12_shift_S_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(S12ShiftR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'R_Addition_CSm_R':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(RAdditionCS(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'r13_insertion_RSR':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(reac_list[i][0] + 1) + '_' + str(reac_list[i][1] + 1) + '_' + str(reac_list[i][2] + 1) + '_' + str(reac_list[i][3] + 1)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(R13InsertionRSR(self.species,self.qc,self.par,reac_list[i],name))
            elif reac_id == 'combinatorial':
                name = str(self.species.chemid) + '_' + reac_id + '_' + str(i)
                self.species.reac_name.append(name)
                self.species.reac_obj.append(Combinatorial(self.species,self.qc,self.par,reac_list[i],name))
            else:
                self.species.reac_name.append(0)
        return 0


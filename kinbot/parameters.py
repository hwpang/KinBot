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
File that contains all the default and user-defined parameters for the
reaction search.
The defaults are listed in this file and the user defined parameters
are defined in a json file that needs to be given as an argument to the
initializer.
"""
from __future__ import with_statement
import os
import sys
import json
import logging;

import numpy as np


class Parameters:
    """
    This class initiates all parameters to their defaults and reads in the
    user-defined variables, which overwrite the defaults
    """
    def __init__(self, file=None):
        """
        Initialize all the variable and read the file which is the user input
        file
        """
        # user input file
        self.input_file = file

        self.par = {
            # GENERAL INFO
            # title of the current calculations
            'title': '',
            # verbose log file
            'verbose': 0,

            # INPUT SPECIES INFOR
            # SMILES of the species
            'smiles': '',
            # geometry of the species
            'structure': [],
            # Charge of the species
            'charge': 0,
            # Multiplicity of the species
            'mult': 0,
            # Whether it is a dimer
            'dimer': 0,

            # WHICH STEPS TO TAKE
            # Do a reaction search
            'reaction_search': 1,
            # Which reaction families to include in the search
            'families': ['all'],
            # break all single bonds to find the barriers
            # of potential homolytic scissions
            'homolytic_scissions': 0,
            # if requested with specific_reaction = 1
            # then only these bonds are broken and formed
            'specific_reaction': 0,
            'break_bonds': [],
            'form_bonds': [],
            # Threshold above which barriers are deemed unimportant
            'barrier_threshold': 100.,
            # Number of 0.1 Angstrom steps in bond scans
            'scan_step': 30,
            # Do a full PES scan instead of one well
            'pes': 0,
            # Maximum number of simultaneous kinbot runs in a pes search
            'simultaneous_kinbot': 5,
            # Perform high level optimization and freq calculation (L2)
            'high_level': 0,
            # Do a conformational search
            'conformer_search': 0,
            # Do a hindered rotor scan
            'rotor_scan': 0,
            # Number of points along the rotor scan
            'nrotation': 12,
            # Make figures of the HIR profiles
            'plot_hir_profiles': 0,
            # Do master equation calculations
            'me': 0,
            # Number of HIR restarts in case a lower energy point gets found
            'rotation_restart': 3,
            # Maximum number of diherals for which exhaustive
            # comformation searches are done
            'max_dihed': 5,
            # Number of random conformers in case no exhaustive search is done
            'random_conf': 500,
            # For the combinatorial search, minimum number of bonds to break
            # this value is decreased by 1 for radical reactions
            'min_bond_break': 2,
            # For the combinatorial search, maximum number of bonds to break
            # this value is decreased by 1 for radical reactions
            'max_bond_break': 3,
            # For the combinatorial search, include molecular pathways
            'comb_molec': 1,
            # For the combinatorial search, include reactions with radicals
            'comb_rad': 1,
            # For the combinatorial search, include reactions with lone electron pairs
            'comb_lone': 1,
            # For the combinatorial search, include reactions with pi electrons
            'comb_pi': 1,
            # For the combinatorial search, allow the breaking of valence
            'break_valence': 1,
            # Search for one specific reaction using combinatorial approach
            'one_reaction_comb' : 0,
            # Search for one specific reaction using family approach
            'one_reaction_fam' : 0,
            # break the following bonds
            'break_bonds' : [],
            # form the following bonds
            'form_bonds' : [],

            # QUANTUM CHEMISTRY INFO
            # Which quantum chemistry code to use
            'qc': 'gauss',  # or nwchem
            # nwchem-specific parameter
            'methodclass': 'dft',  # or scf or mp2
            # Command for the quantum chemistry code
            'qc_command': 'g09',
            # Quantum chemistry method to use as L1
            'method': 'b3lyp',
            # Basis set to use
            'basis': '6-31G',
            # Quantum chemistry method to use for high-level L2
            'high_level_method': 'M062X',
            # Basis set to use for high-level
            'high_level_basis': '6-311++G(d,p)',
            # Integral grid for Gaussian, only for the high-level calculations
            'integral': '',
            # for Gaussian irc: IRC(MaxPoints=n)
            'irc_maxpoints': 30,
            # for Gaussian irc, IRC(StepSize=n)
            'irc_stepsize': 20,
            # for Gaussian, request CalcAll for TS optimization
            'calcall_ts': 0,
            # for Gaussian, allow Guess=(Mix,Always) 
            'guessmix' : 0,
            # name of the single point code's name
            'single_point_qc': 'molpro',
            # Name of the template for the single-point calculation (L3)
            # If not specified, then the tpl/[single_point_qc].inp is used
            'single_point_template': '',
            # if there is a key (e.g., Molpro), what it is to read L3
            "single_point_key": "MYENERGY",
            # whether to use sella
            "sella": 0,

            # COMPUTATIONAL ENVIRONEMNT
            # Which queuing system to use
            'queuing': 'pbs',  # or slurm
            # Template for queue:
            'queue_template': '',
            # Name of the queue
            'queue_name': 'medium',
            # E.g. the type of node or anything that comes with -C in SLURM
            'slurm_feature': '',
            # Number of cores to run the L0-L2 qc jobs on
            'ppn': 1,
            # Number of cores to run the L3 qc jobs on
            'single_point_ppn': 1,
            # This many spaces can be used for numbering files, e.g., in ga
            'zf': 4,
            # Scratch directory
            'scratch': '/scratch/jzador',
            # User name
            'username': 'jzador',


            # MASTER EQUATION
            # Which ME code to use:
            'me_code': 'mess',  # or mesmer
            # MESS specific keywords
            'mess_command': 'mess',
            'TemperatureList': [500+100*i for i in range(16)],
            'PressureList': [760],
            'EnergyStepOverTemperature': .2,
            'ExcessEnergyOverTemperature': 30,
            'ModelEnergyLimit': 400,
            'CalculationMethod': 'direct',  # or low-eigenvalue
            'ChemicalEigenvalueMax': 0.2,
            'EnergyRelaxationFactor': 200,
            'EnergyRelaxationPower': .85,
            'EnergyRelaxationExponentCutoff': 15,
            'Epsilons': [7.08, 310.387],
            'Sigmas': [2.576, 6.000],
            'Masses': [4.0, 87.0],
            # MESMER specific keywords
            'mesmer_command': 'mesmer',
        }

        if self.input_file is not None:
            # Read the user input and overwrite the user-defined parameters
            self.read_user_input()

    def read_user_input(self):
        """
        Read the user input file and overwrite the default values
        """
        try:
            with open(self.input_file) as json_file:
                try:
                    user_data = json.load(json_file)
                except ValueError:
                    msg = 'There is an error in the json input file'
                    raise ValueError(msg)
        except IOError:
            msg = 'Input file {} does not exist'.format(self.input_file)
            raise IOError(msg)
        for key in user_data:
            if key in self.par:
                self.par[key] = user_data[key]
            else:
                err = 'KinBot does not recognize option {} with value {}'
                logging.error(err.format(key, user_data[key]))

    def print_parameters(self):
        """
        Make a string out of the parameters
        """
        s = ''
        for key in self.par:
            s += '{}\t{}\n'.format(key, self.par[key])
        return s

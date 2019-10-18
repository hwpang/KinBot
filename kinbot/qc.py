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
import os, sys
import subprocess
import logging
import numpy as np
import re
import time
import copy
import pkg_resources

from ase.db import connect

from kinbot import constants
from kinbot import geometry

class QuantumChemistry:
    """
    This class provides the link between KinBot and the qc code
    It choses the write options, submits the jobs and checks
    the jobs for success or failure.
    A lot of the functionalities are delegated to the templates.
    """
    
    def __init__(self, par, mult, charge):
        self.par = par
        self.qc = par.par['qc']
        self.method = par.par['method']
        self.basis = par.par['basis']
        self.high_level = par.par['high_level']
        self.high_level_method = par.par['high_level_method']
        self.high_level_basis = par.par['high_level_basis']
        self.integral = par.par['integral']
        self.ppn = par.par['ppn']
        self.queuing = par.par['queuing']
        self.queue_name = par.par['queue_name']
        self.slurm_feature = par.par['slurm_feature']
        self.zf = par.par['zf']
        self.db = connect('kinbot.db')
        self.job_ids = {}
        self.irc_maxpoints = par.par['irc_maxpoints']
        self.irc_stepsize = par.par['irc_stepsize']
        self.qc_command = par.par['qc_command']
        self.sella = par.par['sella']
        self.mult = mult
        self.charge = charge
        self.calcall_ts = par.par['calcall_ts']
        self.guessmix = par.par['guessmix']

        self.mem0 = par.par['mem0'].split()[0]
        self.mem0u = par.par['mem0'].split()[1]
        self.mem = par.par['mem'].split()[0]
        self.memu = par.par['mem'].split()[1]
        self.memmp2 = par.par['memmp2'].split()[0]
        self.memmp2u = par.par['memmp2'].split()[1]
        self.memhl = par.par['memhl'].split()[0]
        self.memhlu = par.par['memhl'].split()[1]

        # sometimes there is no slurm feature at all
        if par.par['slurm_feature'] == '':
            self.slurm_feature = ''
        else:
            self.slurm_feature = '#SBATCH -C ' + par.par['slurm_feature']


    def qc_opt(self, species, geom, high_level=0, mp2=0):
        """ 
        Creates a geometry optimization input and runs it. 
        """
       
        if species.wellorts == 0:
            job = str(species.chemid) + '_well'
            if mp2:
                job += '_mp2'
                self.assemble_ase_template(job, 'optmp2', species, geom, species.wellorts, self.sella, fix=[], change=[])
            elif high_level:
                job += '_high'
                self.assemble_ase_template(job, 'opthl', species, geom, species.wellorts, self.sella, fix=[], change=[])
            else:
                self.assemble_ase_template(job, 'opt', species, geom, species.wellorts, self.sella, fix=[], change=[])
        else:
            job = str(species.name)
            if mp2:
                job += '_mp2'
                self.assemble_ase_template(job, 'optmp2', species, geom, species.wellorts, self.sella, fix=[], change=[])
            elif high_level:
                job += '_high'
                self.assemble_ase_template(job, 'opthl', species, geom, species.wellorts, self.sella, fix=[], change=[])
            else:
                self.assemble_ase_template(job, 'opt', species, geom, species.wellorts, self.sella, fix=[], change=[])
        
        return 0


    def qc_freq(self, species, name, geom, wellorts, high_level=0, mp2=0):
        """ 
        Creates a frequency calculation and runs it. Always done with internal calculation of the qc code.
        """

        if wellorts == 0:
            job = str(species.chemid) + '_well'
        else:
            job = str(name)

        if mp2:
            job += '_mp2'
            self.assemble_ase_template(job, 'freqmp2', species, geom, wellorts, 0, fix=[], change=[])
        elif high_level:
            job += '_high'
            self.assemble_ase_template(job, 'freqhl', species, geom, wellorts, 0, fix=[], change=[])
        else:
            self.assemble_ase_template(job, 'freq', species, geom, wellorts, 0, fix=[], change=[])

        return 0


    def qc_hir(self, species, geom, rot_index, ang_index, fix):
        """ 
        Creates a constrained geometry optimization input and runs it. 
        wellorts: 0 for wells and 1 for saddle points
        rot_index: index of the rotor in the molecule
        ang_index: index for the current size of the angle
        fix: four atoms of the dihedral that is currently fixed
        """
        if species.wellorts:
            job = 'hir/' + species.name + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)
        else:
            job = 'hir/' + str(species.chemid) + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)
       
        self.assemble_ase_template(job, 'hir', species, geom, species.wellorts, self.sella, fix, change=[])

        return 0


    def qc_ring_conf(self, species, geom, fix, change, conf_nr, scan_nr):
        """ 
        Creates a constrained geometry optimization input for the
        conformational search of cyclic structures and runs it.
        Make use of the ASE optimizer PCOBFGS

        qc: 'gauss' or 'nwchem'
        scan: list of dihedrals to be scanned and their values
        wellorts: 0 for wells and 1 for saddle points
        conf_nr: number of the conformer in the conformer search
        scan_nr: number of the scan for this conformer
        """
        if species.wellorts:
            job = 'conf/' + species.name + '_r' + str(conf_nr).zfill(self.zf) + '_' + str(scan_nr).zfill(self.zf)
        else:
            job = 'conf/' + str(species.chemid) + '_r' + str(conf_nr).zfill(self.zf) + '_' + str(scan_nr).zfill(self.zf)
        
        self.assemble_ase_template(job, 'ringconf', species, geom, species.wellorts, self.sella, fix, change)

        return 0


    def qc_conf(self, species, geom, index=-1, freq=0, mp2=0):
        """ 
        Creates a geometry optimization input for the conformational search and runs it.
        wellorts: 0 for wells and 1 for saddle points
        index: >=0 for sampling, each job will get numbered with index
        freq: request a single freq calculation
        mp2: calculate at the mp2 level
        """

        if index == -1:
            job = 'conf/' + str(species.chemid) + '_well'
        else:
            if species.wellorts:
                job = 'conf/' + species.name + '_' + str(index).zfill(self.zf)
            else:
                job = 'conf/' + str(species.chemid) + '_' + str(index).zfill(self.zf)
       
        if freq:
            if mp2:
                self.assemble_ase_template(job, 'freqmp2', species, geom, species.wellorts, self.sella, fix=[], change=[])
            else:
                self.assemble_ase_template(job, 'freq', species, geom, species.wellorts, self.sella, fix=[], change=[])
        else:
            if mp2:
                self.assemble_ase_template(job, 'confmp2', species, geom, species.wellorts, self.sella, fix=[], change=[])
            else:
                self.assemble_ase_template(job, 'conf', species, geom, species.wellorts, self.sella, fix=[], change=[])

        return 0


    def submit_qc(self, job, mem, memu, freq, singlejob):
        """
        Submit a job to the queue, unless the job:
            * has finished with Normal termination
            * has finished with Error termination
            * is currently running
        It means that the job is only submitted if it's
            * not running
            * never started
            * was killed half way

        However, if the optional parameter singlejob is set to False, then 
        the job is run only if it has finished, which means that Normal or Error 
        termination are both acceptable finish states in this case.
        This is for continuations, when the continuing jobs overwrite each other.
        """

        check = self.check_qc(job) # returns 0, 'running', 'normal', 'normal freq', or 'error'
        if freq:
            if check == 'running' or check == 'normal freq' or check == 'error':
                return 0
        else:
            if singlejob == True:
                if check != 0: # running or finished, so no need to submit again 
                    return 0 
            else:
                if check == 'running': 
                    return 0  # running

        try: 
            if self.par.par['queue_template'] == '':
                template_head_file = pkg_resources.resource_filename('tpl', self.queuing + '.tpl')
            else:
                template_head_file = self.par.par['queue_template']
        except OSError:
            logging.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logging.error('Or no file is found at {}.'.format(self.par.par['queue_template']))
            logging.error('Exiting')
            sys.exit()

        template_file = pkg_resources.resource_filename('tpl', self.queuing + '_python.tpl')
        python_file = '{}.py'.format(job)
        
        python_template = open(template_head_file, 'r').read() + open(template_file, 'r').read()

        if self.queuing == 'pbs':
            python_template = python_template.format(name=job, 
                                                     ppn=self.ppn, 
                                                     queue_name=self.queue_name, 
                                                     dir='perm', 
                                                     python_file=python_file, 
                                                     arguments='',
                                                     mem='{}{}'.format(mem, memu))
        elif self.queuing == 'slurm':
            python_template = python_template.format(name=job, 
                                                     ppn=self.ppn, 
                                                     queue_name=self.queue_name, 
                                                     dir='perm', 
                                                     slurm_feature=self.slurm_feature, 
                                                     python_file=python_file, 
                                                     arguments='',
                                                     mem='{}{}'.format(mem, memu))
        else:
            logging.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logging.error('Exiting')
            sys.exit()

        qu_file = '{}{}'.format(job, constants.qext[self.queuing])
        with open(qu_file, 'w') as f_out_qu:
            f_out_qu.write(python_template)

        command = [constants.qsubmit[self.queuing], job + constants.qext[self.queuing]]
        process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = process.communicate()
        out = out.decode()
        if self.queuing == 'pbs':
            pid = out.split('\n')[0].split('.')[0]
        elif self.queuing == 'slurm':
            pid = out.split('\n')[0].split()[3]
        self.job_ids[job] = pid
        
        return 1  # important to keep it 1, this is the natural counter of jobs submitted


    def get_qc_geom(self, job, natom, wait=0, allow_error=0):
        """
        Get the geometry from the ase database file.
        Returns it, with the following conditions about the status of the job.
        If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
        If wait = 1, wait for the job to finish.
        If wait = 2, return the last geometry while the job is still running.
            This option is to monitor the formation of bimolecular products.
        If allow_error = 0, do not read the final geometry if the status is not "normal"
        if allow_error = 1, read the geometry even though there is an error in the output file
            This option is to read the final IRC geometry when it did not converge
        """ 
        geom = np.zeros((natom,3))    
        
        check = self.check_qc(job)
        if check == 'error' and not allow_error: return -1, geom
        status = 0
        while 1:
            check = self.check_qc(job)
            if check == 'running': 
                if wait == 1:
                    time.sleep(1)
                elif wait == 2:
                    status = 2
                    break
                else: 
                    return 1, geom 
            else:
                break
        if check != 'normal' and check != 'normal freq':
            if not allow_error:
                if wait != 2: return -1, geom

        #open the database
        rows = self.db.select(name = job)
        
        found_entry = 0
        #take the last entry
        for row in rows:
            mol = row.toatoms()
            geom = mol.positions
            found_entry = 1
        
        if found_entry:
            return status, geom
        else:
            return -1, geom

    def get_second_to_last_geom(self, job, natom, wait=0, allow_error = 0):
        """
        Get the geometry from the ase database file.
        Returns it, with the following conditions about the status of the job.
        If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
        If wait = 1, wait for the job to finish.
        If wait = 2, return the last geometry while the job is still running.
            This option is to monitor the formation of bimolecular products.
        If allow_error = 0, do not read the final geometry if the status is not "normal"
        if allow_error = 1, read the geometry even though there is an error in the output file
            This option is to read the final IRC geometry when it did not converge
        """ 
        
       
        geom = np.zeros((natom,3))    
        
        check = self.check_qc(job)
        if check == 'error' and not allow_error: return -1, geom
        status = 0
        while 1:
            check = self.check_qc(job)
            if check == 'running': 
                if wait == 1:
                    time.sleep(1)
                elif wait == 2:
                    status = 2
                    break
                else: 
                    return 1, geom 
            else:
                break
        if check != 'normal' and check != 'normal freq':
            if not allow_error:
                if wait != 2: return -1, geom

        #open the database
        rows = self.db.select(name = job)
        
        found_entry = 0
        geoms = []
        #take the last entry
        for row in rows:
            mol = row.toatoms()
            geoms.append(mol.positions)
            found_entry = 1
        
        if found_entry:
            return status, geoms[-2]
        else:
            return -1, geom

    def get_qc_freq(self, job, natom, wait=0, allow_error = 0):
        """
        Get the frequencies from the ase database file
        If wait is set to 1, it will wait for the job to finish.
        """ 
        
        check = self.check_qc(job)
        if check == 'error': return -1, [0]
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 1, []
            else:
                break
        
        if check != 'normal' and check != 'normal freq':
            return -1, [0]

        freq = []
        
        #open the database
        rows = self.db.select(name = job)
        
        #take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                if not row.data.get('frequencies') is None:
                    freq = list(row.data.get('frequencies'))
        
        if len(freq) == 0 and natom > 1:
            return -1,freq

        return 0, freq

    def get_qc_energy(self, job, wait=0):
        """
        Read the last energy from a job. 
        For Gaussian currently works for DFT and HF only.
        For NWChem it works for optimization jobs, using the @ handle.
        If wait is set to 1, it will wait for the job to finish, otherwise
        just reads the last energy in the file.
        Returns the error code and the energy
        Error codes:
        -1: error
         0: normal and done
         1: running
        """

        check = self.check_qc(job)
        if check == 'error': return -1, 0.
        while 1:
            check = self.check_qc(job)
            if check == 'running': 
                if wait == 1:
                    time.sleep(1)
                else:
                    return 1, 0.
            else:
                break

        #open the database
        rows = self.db.select(name = job)
        
        #take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                if not row.data.get('energy') is None:
                    energy = row.data.get('energy')
        
        #ase energies are always in ev, convert to hartree
        energy *= constants.EVtoHARTREE
        
        return 0, energy



    def get_qc_zpe(self, job, wait=1):
        """
        Read the zero point energy. 
        If wait is set to 1 (default), it will wait for the job to finish.
        """
        
        check = self.check_qc(job)
        if check == 'error': return -1, 0.
        while 1:
            check = self.check_qc(job)
            if check == 'running': 
                if wait == 1:
                    time.sleep(1)
                else:
                    return 0, 0.  # TODO is it not a problem for atoms??
            else:
                break

        #open the database
        rows = self.db.select(name = job)

        #take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                zpe = row.data.get('zpe')
        
        return 0, zpe

    def is_in_database(self, job):
        """
        Checks if the current job is in the database:
        """
        #open the database
        rows = self.db.select(name = job)
        
        mol = None
        
        #take the last entry
        for row in rows:
            mol = row.toatoms()
        
        if mol is None:
            return 0
        
        return 1

    def check_qc(self, job):
        """
        Checks the status of the qc job.
        0: not in database (yet)
        'running'
        data['status'] can be 'normal', 'normal freq', or 'error'
        """
        if self.qc == 'gauss':
            log_file = job + '.log'
        elif self.qc == 'nwchem':
            log_file = job + '.out'
        log_file_exists = os.path.exists(log_file)
        
        devnull = open(os.devnull, 'w')
        if self.queuing == 'pbs':
            command = 'qstat -f | grep ' + '"Job Id: ' + self.job_ids.get(job,'-1') + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell = True, stdout=devnull, stderr=devnull)) == 0: 
                return 'running'
        elif self.queuing == 'slurm':
            #command = 'scontrol show job ' + self.job_ids.get(job,'-1') + ' | grep "JobId=' + self.job_ids.get(job,'-1') + '"' + ' > /dev/null'
            command = 'squeue'
            process = subprocess.Popen(command,
                                       shell=True,
                                       stdout=subprocess.PIPE,
                                       stdin=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            out,err = process.communicate()
            out = out.decode()
            for line in out.split('\n'):
                if len(line) > 0:
                    while line.startswith(' '):
                        line = line[1:]
                    pid = line.split()[0]
                    if pid == self.job_ids.get(job,'-1'):
                        return 'running'
        else:
            logging.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logging.error('Exiting')
            sys.exit()

        if self.is_in_database(job) and log_file_exists: #by deleting a log file, you allow restarting a job
            #open the database
            rows = self.db.select(name = job)
            data = None
            #take the last entry
            for row in rows:
                if hasattr(row,'data'):
                    data = row.data
            if data is None:
                return 0
            else:
                return data['status']
        else: 
            return 0


    def add_dummy(self, atom, geom, bond):
            """
            Add a dummy atom if needed to linear substructures.
            """

            dummy = geometry.is_linear(geom, bond)

            if len(dummy) > 0: # add a dummy atom for each close to linear angle
                for d in dummy:
                    atom = np.append(atom,['X'])
                    geom = np.concatenate((geom, [d]), axis=0)
            dummy = [d.tolist() for d in dummy]

            return atom, geom, dummy
       

    def assemble_ase_template(self, job, task, species, geom, wellorts, sella, fix=[], change=[], release=[], 
            app_traj=None, tight=True, singlejob=True):
        """
        Assemble the template for an ASE.
        """

        atom = copy.deepcopy(species.atom)
        if not sella:
            atom, geom, dummy = self.add_dummy(atom, geom, species.bond) 
        else:
            dummy = []

        # frequency calculations are not done by sella
        nosella = [ 'freq', 'freqmp2', 'freqhigh']
        if task in nosella:
            sella = 0

        chk = True

        # TASKS AND OPTIONS
        # OPTIMIZATIONS
        if task == 'opt':
            method = self.method
            basis = self.basis 
            integral = ''
            order = wellorts
            freq = False
            guess = False
            maxattempt = 2
            mem = self.mem
            memu = self.memu

        elif task == 'optmp2':
            method = 'mp2'
            basis = self.basis 
            integral = ''
            order = wellorts
            freq = False
            guess = False
            maxattempt = 2
            mem = self.memmp2
            memu = self.memmp2u

        elif task == 'opthl':
            method = self.high_level_method
            basis = self.high_level_basis
            integral = self.integral
            order = wellorts
            freq = False
            guess = False
            maxattempt = 2
            mem = self.memhl
            memu = self.memhlu

        elif task == 'preopt0':
            method = 'am1'
            basis = ''
            integral = ''
            order = 0
            freq = False
            guess = False
            maxattempt = 2
            mem = self.mem0
            memu = self.mem0u

        elif task == 'preopt':
            method = 'am1'
            basis = ''
            integral = ''
            order = 0
            freq = False
            guess = True
            maxattempt = 2
            mem = self.mem0
            memu = self.mem0u

        # FREQUENCY

        if task == 'freq':
            method = self.method
            basis = self.basis 
            integral = ''
            order = -1
            freq = True
            guess = True
            maxattempt = 1
            mem = self.mem
            memu = self.memu

        elif task == 'freqmp2':
            method = 'mp2'
            basis = self.basis 
            integral = ''
            order = -1
            freq = True
            guess =  True
            maxattempt = 1
            mem = self.memmp2
            memu = self.memmp2u

        elif task == 'freqhl':
            method = self.high_level_method
            basis = self.high_level_basis
            integral = self.integral
            order = -1
            freq = True
            guess = True
            maxattempt = 1
            mem = self.memhl
            memu = self.memhlu

        # CONFORMERS

        elif task == 'conf':
            method = self.method
            basis = self.basis 
            integral = ''
            order = wellorts
            freq = False
            guess = False
            maxattempt = 2
            mem = self.mem
            memu = self.memu

        elif task == 'confmp2':
            method = 'mp2'
            basis = self.basis 
            integral = ''
            order = wellorts
            freq = False
            guess = False
            maxattempt = 2
            mem = self.mem
            memu = self.memu

        elif task == 'ringconf':
            method = 'am1'
            basis = ''
            integral = ''
            order = 0
            freq = False
            guess = False
            maxattempt = 1
            mem = self.mem0
            memu = self.mem0u

        # HINDERED ROTORS
        elif task == 'hir':
            method = self.high_level_method
            basis = self.high_level_basis 
            integral = self.integral
            order = wellorts
            freq = False
            guess = False
            maxattempt = 2
            mem = self.mem
            memu = self.memu
 
        # IRC
        elif task == 'ircf' or task == 'ircr':
            method = self.method
            basis = self.basis 
            integral = ''
            order = -1  # do not optimize
            freq = False
            guess = True
            maxattempt = 1
            mem = self.mem
            memu = self.memu

        elif task == 'ircfmp2' or task == 'ircrmp2':
            method = 'mp2'
            basis = self.basis 
            integral = ''
            order = -1  # do not optimize
            freq = False
            guess = True
            maxattempt = 1
            mem = self.memmp2
            memu = self.memmp2u
 
        # IRC PRODUCT
        elif task == 'prodirc':
            method = self.method
            basis = self.basis 
            integral = ''
            order = 0  
            freq = False
            guess = True
            maxattempt = 1
            mem = self.mem
            memu = self.memu

        elif task == 'prodircmp2':
            method = 'mp2'
            basis = self.basis 
            integral = ''
            order = 0  
            freq = False
            guess = True
            maxattempt = 1
            mem = self.memmp2
            memu = self.memmp2u

        # TEMPLATES 
        
        header = pkg_resources.resource_filename('tpl', 'ase_header.tpl.py')
        with open(header) as f:
            tpl_header = f.read()

        qc_translate = pkg_resources.resource_filename('tpl', 'ase_{qc}_translate.tpl.py'.format(qc = self.qc)) 
        with open(qc_translate) as f:
            tpl_translate = f.read()

        if not sella and self.qc == 'gauss':
            qc_translate_qc = pkg_resources.resource_filename('tpl', 'ase_{qc}_translate_{qc}.tpl.py'.format(qc = self.qc))
            with open(qc_translate_qc) as f:
                tpl_qc = f.read()
        else:
            qc_translate_qc = pkg_resources.resource_filename('tpl', 'ase_none.tpl.py'.format(qc = self.qc))
            with open(qc_translate_qc) as f:
                tpl_qc = f.read()

        qc_calc = pkg_resources.resource_filename('tpl', 'ase_calc.tpl.py')
        with open(qc_calc) as f:
            tpl_calc = f.read()

        if sella:
            task_sella = pkg_resources.resource_filename('tpl', 'ase_task_sella.tpl.py')
            with open(task_sella) as f:
                tpl_task = f.read()
        elif self.qc == 'gauss':
            task_qc = pkg_resources.resource_filename('tpl', 'ase_task_{qc}.tpl.py'.format(qc = self.qc))
            with open(task_qc) as f:
                tpl_task = f.read()

        write_db = pkg_resources.resource_filename('tpl', 'ase_write_db.tpl.py') 
        with open(write_db) as f:
            tpl_write_db = f.read()

        #ASSEMBLE TEMPLATES

        template = '{}{}{}{}{}{}'.format(tpl_header, tpl_translate, tpl_qc, tpl_calc, tpl_task, tpl_write_db)

        #SUBSTITUTE TEMPLATES 
        
        template = template.format(label=job, 
                                   atom=list(atom), 
                                   geom=[list(gi) for gi in geom], 
                                   ppn=self.ppn,
                                   method=method,
                                   basis=basis,
                                   mult=self.mult,
                                   charge=self.charge,
                                   chk=chk,
                                   guess=guess,
                                   integral=integral,
                                   dummy=dummy,
                                   sella=sella,
                                   order=order,
                                   freq=freq,
                                   task=task,
                                   irc_maxpoints=self.irc_maxpoints,
                                   irc_stepsize=self.irc_stepsize,
                                   qc=self.qc,
                                   fix=fix,
                                   change=change,
                                   release=release,
                                   maxattempt=maxattempt,
                                   qc_command=self.qc_command,
                                   guessmix=self.guessmix,
                                   calcall_ts=self.calcall_ts,
                                   mem=mem,
                                   memu=memu,
                                   app_traj=app_traj,
                                   tight=tight)


        f_out = open('{}.py'.format(job),'w')
        f_out.write(template)
        f_out.close()

        # this will return the same number as submit_qc
        return self.submit_qc(job, mem, memu, freq, singlejob)




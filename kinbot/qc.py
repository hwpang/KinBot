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
    
    def __init__(self,par):
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
        self.mult = par.par['mult']
        self.charge = par.par['charge']
        # sometimes there is no slurm feature at all
        if par.par['slurm_feature'] == '':
            self.slurm_feature = ''
        else:
            self.slurm_feature = '#SBATCH -C ' + par.par['slurm_feature']
        

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


    def qc_conf(self, species, geom, index=-1, ring=0):
        """ 
        Creates a geometry optimization input for the conformational search and runs it.
        wellorts: 0 for wells and 1 for saddle points
        index: >=0 for sampling, each job will get numbered with index
        """
        if index == -1:
            job = 'conf/' + str(species.chemid) + '_well'
        else:
            r = ''
            if ring: r = 'r'
            if species.wellorts:
                job = 'conf/' + species.name + '_' + r + str(index).zfill(self.zf)
            else:
                job = 'conf/' + str(species.chemid) + '_' + r + str(index).zfill(self.zf)
        
        self.assemble_ase_template(job, 'conf', species, geom, species.wellorts, self.sella, fix=[], change=[])

        return 0

    def qc_opt(self, species, geom, high_level=0, mp2=0):
        """ 
        Creates a geometry optimization input and runs it. 
        """
        
        job = str(species.chemid) + '_well'
        if high_level:
            job = str(species.chemid) + '_well_high'
        if mp2:
            job = str(species.chemid) + '_well_mp2'
        
        if mp2:
            self.assemble_ase_template(job, 'optmp2', species, geom, species.wellorts, self.sella, fix=[], change=[])
        elif high_level:
            self.assemble_ase_template(job, 'opthl', species, geom, species.wellorts, self.sella, fix=[], change=[])
        else:
            self.assemble_ase_template(job, 'opt', species, geom, species.wellorts, self.sella, fix=[], change=[])

        return 0


    def qc_opt_ts(self, species, geom, high_level=0):
        """
        Creates a ts optimization input and runs it
        """
        
        job = str(species.name)
        if high_level:
            job += '_high'

        if high_level:
            self.assemble_ase_template(job, 'opthl', species, geom, species.wellorts, self.sella, ts=0, fix=[], change=[])
        elif step == 0:
            self.assemble_ase_template(job, 'preopt0', species, geom, 0, self.sella, ts=0, fix=[], change=[])
        elif step < max_step:
            self.assemble_ase_template(job, 'preopt', species, geom, 0, self.sella, ts=0, fix=[], change=[])
        else:
            self.assemble_ase_template(job, 'opt', species, geom, species.wellorts, self.sella, ts=0, fix=[])

        return 0

    def submit_qc(self, job, singlejob=1):
        """
        Submit a job to the queue, unless the job:
            * has finished with Normal termination
            * has finished with Error termination
            * is currently running
        However, if the optional parameter singlejob is set to zero, then 
        the job is run only if it has finished earlier with normal termination.
        This is for continuations, when the continuing jobs overwrite each other.
        """

        check = self.check_qc(job)
        if singlejob == 1:
            if check != 0: return 0  # either normal or error termination, but is in database and done
        else:
            if check == 'running': return 0  # still running

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
        
        python_template = open(template_head_file, 'r').read() 
        python_template = open(template_head_file, 'r').read() + open(template_file, 'r').read()

        if self.queuing == 'pbs':
            python_template = python_template.format(   name=job, ppn=self.ppn, queue_name=self.queue_name, 
                                                        dir='perm', python_file=python_file, arguments='' )
        elif self.queuing == 'slurm':
            python_template = python_template.format(   name=job, ppn=self.ppn, queue_name=self.queue_name, dir='perm', 
                                                        slurm_feature=self.slurm_feature, python_file=python_file, arguments='' )
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


    def get_qc_geom(self, job, natom, wait=0, allow_error = 0):
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
        if check != 'normal':
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
        if check != 'normal':
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
        
        if check != 'normal': return -1, [0]

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
                    return 0, 0.
            else:
                break

        #open the database
        rows = self.db.select(name = job)

        #take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                zpe = row.data.get('zpe')
        
        return 0, zpe

    def read_qc_hess(self, job, natom):
        """
        Read the hessian of a gaussian chk file
        """
        
        check = self.check_qc(job)
        if check != 'normal': 
            return []
        
        #initialize hessian matrix
        hess = np.zeros((3*natom,3*natom))
        
        if self.qc == 'gauss':
            
            fchk = str(job) + '.fchk'
            chk = str(job) + '.chk'
            if os.path.exists(chk):
            #create the fchk file using formchk
                os.system('formchk ' + job + '.chk > /dev/null')
            
            with open(fchk) as f:
                lines = f.read().split('\n')
            
            nvals = 3 * natom * (3 * natom + 1) / 2

            for index, line in enumerate(reversed(lines)):
                if re.search('Cartesian Force Constants', line) != None:
                    hess_flat = []
                    n = 0
                    while len(hess_flat) < nvals:
                        hess_flat.extend([float(val) for val in lines[-index + n].split()])
                        n += 1
                    n = 0
                    for i in range(3*natom):
                        for j in range(i+1):
                            hess[i][j] = hess_flat[n]
                            hess[j][i] = hess_flat[n]
                            n += 1
                    break
        return hess

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
        data['status'] can be 'normal' or 'error'
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
        #if int(subprocess.call(command, shell = True, stdout=devnull, stderr=devnull)) == 0: 
        #    return 'running' 
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
       

    def assemble_ase_template(self, job, task, species, geom, wellorts, sella, fix, change):
        """
        Assemble the template for an ASE.
        """

        atom = copy.deepcopy(species.atom)
        atom, geom, dummy = self.add_dummy(atom, geom, species.bond) 
        wellorts = bool(wellorts)

        # TASKS AND OPTIONS
        if task == 'opt':
            method = self.method
            basis = self.basis 
            integral = ''
            opt = True
            order = wellorts
            freq = True
            guess = False
            chk = True
            maxattempt = 2
            singlejob = True

        elif task == 'optmp2':
            method = 'mp2'
            basis = self.basis 
            integral = ''
            opt = True
            order = wellorts
            freq = True
            guess = False
            chk = True
            maxattempt = 2
            singlejob = True

        elif task == 'opthl':
            method = self.high_level_method
            basis = self.high_level_basis
            integral = self.integral
            opt = True
            order = wellorts
            freq = True
            guess = False
            chk = True
            maxattempt = 2
            singlejob = True

        elif task == 'preopt0':
            method = 'am1'
            basis = ''
            integral = ''
            opt = True
            order = 0
            freq = False
            guess = False
            chk = True
            maxattempt = 2
            singlejob = True

        elif task == 'preopt':
            method = 'am1'
            basis = ''
            integral = ''
            opt = True
            order = 0
            freq = False
            guess = True
            chk = True
            maxattempt = 2
            singlejob = False

        #CONFORMERS
        elif task == 'conf':
            method = self.method
            basis = self.basis 
            integral = ''
            opt = True
            order = wellorts
            freq = True
            guess = False
            chk = False
            maxattempt = 2
            singlejob = True

        elif task == 'ringconf':
            method = 'am1'
            basis = ''
            integral = ''
            opt = False
            order = 0
            freq = False
            guess = False
            chk = False
            maxattempt = 1
            singlejob = True

        # HINDERED ROTORS
        elif task == 'hir':
            method = self.high_level_method
            basis = self.high_level_basis 
            integral = self.integral
            opt = True
            order = wellorts
            freq = False
            guess = False
            chk = False
            maxattempt = 2
            singlejob = True
 
        # IRC
        elif task == 'irc':
            method = self.method
            basis = self.basis 
            integral = ''
            opt = False
            order = 0
            freq = False
            guess = True
            chk = True
            maxattempt = 1
            singlejob = False


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

        if sella:
            task_sella = pkg_resources.resource_filename('tpl', 'ase_task_sella.tpl.py')
            with open(task_sella) as f:
                tpl_task = f.read()
        elif self.qc == 'gauss':
            task_qc = pkg_resources.resource_filename('tpl', 'ase_task_{qc}.tpl.py'.format(qc = self.qc))
            with open(task_qc) as f:
                tpl_task = f.read()

        #ASSEMBLE TEMPLATES

        template = tpl_header + tpl_translate + tpl_qc + tpl_task 

        #SUBSTITUTE TEMPLATES 
        #CalcAll TODO

        if 0:
            template = template.format(label=job, 
                                       atom=list(atom), 
                                       geom=list([list(gi) for gi in geom]),
                                       ppn=self.ppn,
                                       method=method,
                                       basis=basis,
                                       mult=self.mult,
                                       charge=self.charge,
                                       chk=chk,
                                       guess=guess,
                                       integral=integral,
                                       dummy=dummy,
                                       sella=self.sella,
                                       order=order,
                                       freq=freq,
                                       task=task,
                                       irc_maxpoints=self.irc_maxpoints,
                                       irc_stepsize=self.irc_stepsize,
                                       qc=self.qc,
                                       fix=fix,
                                       change=change,
                                       maxattempt=maxattempt,
                                       qc_command=self.qc_command)


        f_out = open('{}.py'.format(job),'w')
        f_out.write(template)
        f_out.close()

        # this will return the same number as submit_qc
        return self.submit_qc(job, singlejob)




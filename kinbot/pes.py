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
This is the main class to run KinBot to explore
a full PES instead of only the reactions of one well
"""
from __future__ import print_function
from __future__ import absolute_import
import sys
import os
import stat
import shutil
import logging
import datetime
import time
import subprocess
import json
from distutils.dir_util import copy_tree
import pkg_resources
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from ase.db import connect

from kinbot import constants
from kinbot import license_message
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.mess import MESS


def main():
    try:
        input_file = sys.argv[1]
    except IndexError:
        print('To use the pes script, supply one argument being the input file!')
        sys.exit(-1)

    # TODO: write information about the arguments
    # change this to nice argument parsers with
    # dashes etc.
    no_kinbot = 0
    task = 'all'
    names = []
    if len(sys.argv) > 2:
        if sys.argv[2] == 'no-kinbot':
            no_kinbot = 1
    if len(sys.argv) > 3:
        # possible tasks are:
        # 1. all: This is the default showing all pathways
        # 2. lowestpath: show the lowest path between the species
        # corresponding to the names
        # 3. allpaths: show all paths between the species
        # corresponding to the names
        # 4. wells: show all reactions of one wells
        # corresponding to the names
        task = sys.argv[3]
        names = sys.argv[4:]

    # print the license message to the console
    print(license_message.message)

    # initialize the parameters
    par = Parameters(input_file)

    # set up the logging environment
    logging.basicConfig(filename='pes.log', level=logging.INFO)

    logging.info(license_message.message)
    msg = 'Starting the PES search at {}'.format(datetime.datetime.now())
    logging.info(msg)

    if par.par['pes'] and par.par['specific_reaction']:
        logging.error('Specific reaction cannot be searched in PES mode.')
        return

    well0 = StationaryPoint('well0',
                            par.par['charge'],
                            par.par['mult'],
                            smiles=par.par['smiles'],
                            structure=par.par['structure'])
    well0.characterize(par.par['dimer'])
    write_input(par, well0, par.par['barrier_threshold'], os.getcwd())

    # add the initial well to the chemids
    with open('chemids', 'w') as f:
        f.write(str(well0.chemid) + '\n')

    # create a directory for the L3 single point calculations 
    # direcotry has the name of the code, e.g., molpro
    try:
        os.mkdir(par.par['single_point_qc'])
    except OSError:
        pass

    # maximum number of kinbot jobs that run simultaneously
    max_running = par.par['simultaneous_kinbot']
    # jobs that are running
    running = []
    # jobs that are finished
    finished = []
    # list of all jobs
    jobs = []
    # dict of the pid's for all jobs
    pids = {}
    while 1:
        j = len(jobs)
        with open('chemids', 'r') as f:
            jobs = f.read().split('\n')
            jobs = [ji for ji in jobs if ji != '']

        if len(jobs) > j:
            logging.info('\tPicked up new jobs: ' + ' '.join(jobs[j:]))

        if len(finished) == len(jobs):
            break

        while (len(running) < max_running and
               len(running) + len(finished) < len(jobs)):
            # start a new job
            job = jobs[len(running) + len(finished)]
            pid = 0
            if not no_kinbot:
                pid = submit_job(job, par)  # kinbot is submitted here
            else:
                get_wells(job)
            pids[job] = pid
            t = datetime.datetime.now()
            logging.info('\tStarted job {} at {}'.format(job, t))
            running.append(job)
        # check if a thread is done
        for job in running:
            if not check_status(job, pids[job]):
                t = datetime.datetime.now()
                logging.info('\tFinished job {} at {}'.format(job, t))
                finished.append(job)
                if not no_kinbot:
                    # write a temporary pes file
                    # remove old xval and im_extent files
                    try:
                        os.remove('{}_xval.txt'.format(par.par['title']))
                    except OSError:
                        pass
                    try:
                        os.remove('{}_im_extent.txt'.format(par.par['title']))
                    except OSError:
                        pass
                    postprocess(par, jobs, task, names)
        # remove the finished threads
        for job in finished:
            if job in running:
                running.remove(job)
        if not no_kinbot:
            # write a summary of what is running and finished
            summary_lines = []
            summary_lines.append('Total\t\t{}'.format(len(jobs)))
            summary_lines.append('Running\t\t{}'.format(len(running)))
            summary_lines.append('Finished\t{}'.format(len(finished)))
            summary_lines.append('')
            summary_lines.append('Running:')
            for job in running:
                summary_lines.append('\t{}'.format(job))
            summary_lines.append('')
            summary_lines.append('Finished:')
            for job in finished:
                summary_lines.append('\t{}'.format(job))
            with open('pes_summary.txt', 'w') as f:
                f.write('\n'.join(summary_lines))
            time.sleep(1)

    postprocess(par, jobs, task, names)
    # make molpro inputs for all keys above
    # place submission script in the directory for offline submission
    # read in the molpro energies for the keys in the above three dicts
    # for key in newdict.keys():
    #      print(key)
    # if all energies are there
    # do something like postprocess, but with new energies
    # postprocess_L3(saddle_zpe, well_zpe, prod_zpe, saddle_energy, well_energy, prod_energyi, conn)

    # Notify user the search is done
    logging.info('PES search done!')
    print('PES search done!')


def get_wells(job):
    """
    Read the summary file and add the wells to the chemid list
    """
    try:
        summary = open(job + '/summary_' + job + '.out', 'r').readlines()
    except:
        return 0
    with open('chemids', 'r') as f:
        jobs = f.read().split('\n')
    jobs = [ji for ji in jobs if ji != '']

    new_wells = []
    for line in summary:
        if line.startswith('SUCCESS'):
            pieces = line.split()
            prod = pieces[3:]
            if (len(prod) == 1 and
                    prod[0] not in jobs and
                    prod[0] not in new_wells):
                new_wells.append(prod[0])
    if len(new_wells) > 0:
        with open('chemids', 'a') as f:
            f.write('\n'.join(new_wells) + '\n')


def postprocess(par, jobs, task, names):
    """
    postprocess a pes search
    par: parameters of the search
    jobs: all of the jobs that were run
    temp: this is a temporary output file writing
    """

    l3done = 1  # flag for L3 calculations to be complete

    # base of the energy is the first well, these are L2 energies
    base_energy = get_energy(jobs[0], jobs[0], 0, par.par['high_level'])
    # L3 energies
    status, base_l3energy = get_l3energy(jobs[0], par)
    if not status:
        l3done = 0
    # L2 ZPE
    base_zpe = get_zpe(jobs[0], jobs[0], 0, par.par['high_level'])

    # list of lists with four elements
    # 1. reactant chemid
    # 2. reaction name
    # 3. products chemid list
    # 4. reaction barrier height
    reactions = []

    # list of the parents for each calculation
    # the key is the name of the calculation
    # the value is the parent directory,
    # i.e. the well kinbot started from to find
    # this calculation
    parent = {}

    wells = []
    failedwells = []
    products = []

    # read all the jobs
    for ji in jobs:
        try:
            summary = open(ji + '/summary_' + ji + '.out', 'r').readlines()
        except:
            failedwells.append(ji)
            continue
        # read the summary file
        for line in summary:
            if line.startswith('SUCCESS'):
                pieces = line.split()
                reactant = ji
                ts = pieces[2]  # this is the long specific name of the reaction
                prod = pieces[3:]  # this is the chemid of the product

                # calculate the barrier based on the new energy base
                barrier = 0. - base_energy - base_zpe
                # overwrite energies with mp2 energy if needed
                if ('R_Addition_MultipleBond' in ts and
                        not par.par['high_level']):
                    base_energy_mp2 = get_energy(jobs[0],
                                                 jobs[0],
                                                 0,
                                                 par.par['high_level'],
                                                 mp2=1)
                    base_zpe_mp2 = get_zpe(jobs[0],
                                           jobs[0],
                                           0,
                                           par.par['high_level'],
                                           mp2=1)
                    barrier = 0. - base_energy_mp2 - base_zpe_mp2
                ts_energy = get_energy(reactant, ts, 1, par.par['high_level'])
                ts_zpe = get_zpe(reactant, ts, 1, par.par['high_level'])
                barrier += ts_energy + ts_zpe
                barrier *= constants.AUtoKCAL

                if reactant not in wells:
                    wells.append(reactant)
                    parent[reactant] = reactant
                if len(prod) == 1:
                    if prod[0] not in wells:
                        if prod[0] not in parent:
                            parent[prod[0]] = reactant
                        wells.append(prod[0])
                else:
                    prod_name = '_'.join(sorted(prod))
                    if prod_name not in products:
                        if prod_name not in parent:
                            parent[prod_name] = reactant
                        products.append('_'.join(sorted(prod)))
                new = 1
                temp = None
                for i, rxn in enumerate(reactions):
                    rxn_prod_name = '_'.join(sorted(rxn[2]))
                    if (reactant == rxn[0] and
                            '_'.join(sorted(prod)) == rxn_prod_name):
                        new = 0
                        temp = i
                    if reactant == ''.join(rxn[2]) and ''.join(prod) == rxn[0]:
                        new = 0
                        temp = i
                if new:
                    reactions.append([reactant, ts, prod, barrier])
                else:
                    # check if the previous reaction has a lower energy or not
                    if reactions[temp][3] > barrier:
                        reactions.pop(temp)
                        reactions.append([reactant, ts, prod, barrier])
        # copy the xyz files
        copy_from_kinbot(ji, 'xyz')
        # copy the L3 calculations here, whatever was in those directories, inp, out, pbs, etc.
        copy_from_kinbot(ji, par.par['single_point_qc'])
    # create a connectivity matrix for all wells and products
    conn, bars = get_connectivity(wells, products, reactions)
    # create a batch submission for all L3 jobs
    # TODO slurm
    if par.par['queuing'] == 'pbs':
        batch = 'batch_L3_pbs.sub'
        with open(batch, 'w') as f:
            for well in wells:
                f.write('qsub molpro/' + well + '.pbs' + '\n')
            for prod in products:
                for frag in prod.split('_'):
                    f.write('qsub molpro/' + frag + '.pbs' + '\n')
            for reac in reactions:
                f.write('qsub molpro/' + reac[1] + '.pbs' + '\n')
        os.chmod(batch, stat.S_IRWXU)  # read, write, execute by owner


    well_energies = {}
    well_l3energies = {}
    for well in wells:
        energy = get_energy(parent[well], well, 0, par.par['high_level'])  # from the db
        zpe = get_zpe(parent[well], well, 0, par.par['high_level'])
        well_energies[well] = ((energy + zpe) - (base_energy + base_zpe)) * constants.AUtoKCAL
        status, l3energy = get_l3energy(well, par)
        if not status:
            l3done = 0  # not all L3 calculations are done
        else:
            well_l3energies[well] = ((l3energy + zpe) - (base_l3energy + base_zpe)) * constants.AUtoKCAL

    prod_energies = {}
    prod_l3energies = {}
    for prods in products:
        energy = 0. - base_energy - base_zpe
        l3energy = 0. - base_l3energy - base_zpe
        for pr in prods.split('_'):
            energy += get_energy(parent[prods], pr, 0, par.par['high_level'])
            zpe = get_zpe(parent[prods], pr, 0, par.par['high_level'])
            energy += zpe
            status, l3energy = get_l3energy(pr, par)
            if not status:
                l3done = 0  # not all L3 calculations are done
            else:
                l3energy = ((l3energy + zpe) - (base_l3energy + base_zpe)) 
        prod_energies[prods] = energy * constants.AUtoKCAL
        prod_l3energies[prods] = l3energy * constants.AUtoKCAL

    ts_l3energies = {}
    for reac in reactions:
        zpe = get_zpe(reac[0], reac[1], 1, par.par['high_level'])
        status, l3energy = get_l3energy(reac[1], par)
        if not status:
            l3done = 0
        else:
            ts_l3energies[reac[1]] =  ((l3energy + zpe) - (base_l3energy + base_zpe)) * constants.AUtoKCAL


    if l3done == 1 and len(reactions) > 1:
        well_energies = well_l3energies
        prod_energies = prod_l3energies
        for reac in reactions:  # swap out the barrier
            reac[3] = ts_l3energies[reac[1]]

    # if L3 was done, everything below is done with that
    # filter according to tasks
    wells, products, reactions, highlight = filter(wells,
                                                   products,
                                                   reactions,
                                                   conn,
                                                   bars,
                                                   well_energies,
                                                   task,
                                                   names)
    # write full pesviewer input
    create_pesviewer_input(par,
                           wells,
                           products,
                           reactions,
                           well_energies,
                           prod_energies,
                           highlight)
    # draw a graph of the network
#    create_graph(wells,
#                 products,
#                 reactions,
#                 well_energies,
#                 prod_energies,
 #                highlight)
    # write_mess
    create_mess_input(par,
                      wells,
                      products,
                      reactions,
                      well_energies,
                      prod_energies,
                      parent)


def filter(wells, products, reactions, conn, bars, well_energies, task, names):
    """
    Filter the wells, products and reactions according to the task
    and the names
    """
    # list of reactions to highlight
    highlight = []

    # 1. all: This is the default showing all pathways
    # 2. lowestpath: show the lowest path between the species
    # corresponding to the names
    # 3. allpaths: show all paths between the species
    # corresponding to the names
    # 4. wells: show all reactions of one wells
    # corresponding to the names
    # 5. temperature
    # 6. threshold_reapply: apply the barrier threshold 
    # cutoff at the highest level that was done

    # filter the reactions according to the task
    if task == 'all':
        filtered_reactions = reactions
        pass
    elif task == 'lowestpath':
        all_rxns = get_all_pathways(wells, products, reactions, names, conn)
        # this is the maximum energy along the minimun energy pathway
        min_energy = None
        min_rxn = None
        for rxn_list in all_rxns:
            barriers = [ri[3] for ri in rxn_list]
            if min_energy is None:
                min_energy = max(barriers)
                min_rxn = rxn_list
            else:
                if max(barriers) < min_energy:
                    min_energy = max(barriers)
                    min_rxn = rxn_list
        filtered_reactions = min_rxn
    elif task == 'allpaths':
        all_rxns = get_all_pathways(wells, products, reactions, names, conn)
        filtered_reactions = []
        for list in all_rxns:
            for rxn in list:
                new = 1
                for r in filtered_reactions:
                    if r[1] == rxn[1]:
                        new = 0
                if new:
                    filtered_reactions.append(rxn)
        # this is the maximum energy along the minimun energy pathway
        min_energy = None
        min_rxn = None
        for rxn_list in all_rxns:
            barriers = [ri[3] for ri in rxn_list]
            if min_energy is None:
                min_energy = max(barriers)
                min_rxn = rxn_list
            else:
                if max(barriers) < min_energy:
                    min_energy = max(barriers)
                    min_rxn = rxn_list
        for rxn in min_rxn:
            highlight.append(rxn[1])
    elif task == 'well':
        if len(names) == 1:
            filtered_reactions = []
            for rxn in reactions:
                prod_name = '_'.join(sorted(rxn[2]))
                if names[0] == rxn[0] or names[0] == prod_name:
                    filtered_reactions.append(rxn)
        else:
            logging.error('Only one name should be given for a well filter')
            logging.error('Received: ' + ' '.join(names))
            sys.exit(-1)
    elif task == 'temperature':
        if len(names) == 1:
            try:
                # read the temperature
                temperature = float(names[0])
            except ValueError:
                logging.error('A float is needed for a temperature filter')
                logging.error('Received: ' + ' '.join(names))
                sys.exit(-1)
            filtered_reactions = []
            # iterate the wells
            wells, filtered_reactions = filter_boltzmann(wells[0],
                                                         [wells[0]],
                                                         reactions,
                                                         filtered_reactions,
                                                         well_energies,
                                                         temperature)
        else:
            logging.error('Only one argument should be given for a temperature filter')
            logging.error('Received: ' + ' '.join(names))
            sys.exit(-1)
    elif task == 'l2threshold':
        filtered_reactions = []
        for rxn in reactions:
            if rxn[3] < par.par['barrier_threshold']: 
                filtered_reactions.append(rxn)
    else:
        logging.error('Could not recognize task ' + task)
        sys.exit(-1)

    # filter the wells
    filtered_wells = []
    for well in wells:
        for rxn in filtered_reactions:
            prod_name = '_'.join(sorted(rxn[2]))
            if well == rxn[0] or well == prod_name:
                if well not in filtered_wells:
                    filtered_wells.append(well)

    # filter the products
    filtered_products = []
    for prod in products:
        for rxn in filtered_reactions:
            prod_name = '_'.join(sorted(rxn[2]))
            if prod == prod_name:
                if prod not in filtered_products:
                    filtered_products.append(prod)

    return filtered_wells, filtered_products, filtered_reactions, highlight


def filter_boltzmann(well, wells, reactions, filtered_reactions,
                     well_energies, temperature):
    """
    Filter the reactions based on branching fractions at a given temperature
    """
    well_energy = well_energies[well]
    # all the reactions that include this well
    one_well_rxns = []
    # iterate all reactions
    for rxn in reactions:
        # check if this reactions belongs to the current well
        if rxn[0] == well or rxn[2][0] == well:
            one_well_rxns.append(rxn)
    # list containing the branching fractions for this well
    branching = []
    for rxn in one_well_rxns:
        # calculate the boltzmann factor
        value = np.exp(-(rxn[3] - well_energy) * 1000 / 1.9872036 / temperature)
        branching.append(value)
    # calculate the actual branching fractions
    br_sum = sum(branching)
    branching = np.array(branching) / br_sum
    for i, rxn in enumerate(one_well_rxns):
        # only add a reaction if the branching fraction is below 1%
        if branching[i] > 0.01:
            new = 1
            for r in filtered_reactions:
                if r[1] == rxn[1]:
                    new = 0
            if new:
                filtered_reactions.append(rxn)
                if not rxn[0] in wells:
                    wells.append(rxn[0])
                    wells, filtered_reactions = filter_boltzmann(rxn[0],
                                                                 wells,
                                                                 reactions,
                                                                 filtered_reactions,
                                                                 well_energies,
                                                                 temperature)
                if len(rxn[2]) == 1:
                    if not rxn[2][0] in wells:
                        wells.append(rxn[2][0])
                        wells, filtered_reactions = filter_boltzmann(rxn[2][0],
                                                                     wells,
                                                                     reactions,
                                                                     filtered_reactions,
                                                                     well_energies,
                                                                     temperature)
    return wells, filtered_reactions


def get_connectivity(wells, products, reactions):
    """
    Create two matrices:
    conn: connectivity (1 or 0) between each pair of stationary points
    bars: barrier height between each pair of stationary points
          0 if not connected
    """
    conn = np.zeros((len(wells) + len(products), len(wells) + len(products)), dtype=int)
    bars = np.zeros((len(wells) + len(products), len(wells) + len(products)))
    for rxn in reactions:
        reac_name = rxn[0]
        prod_name = '_'.join(sorted(rxn[2]))
        i = get_index(wells, products, reac_name)
        j = get_index(wells, products, prod_name)
        conn[i][j] = 1
        conn[j][i] = 1
        barrier = rxn[3]
        bars[i][j] = barrier
        bars[j][i] = barrier
    return conn, bars


def get_all_pathways(wells, products, reactions, names, conn):
    """
    Get all the pathways in which all intermediate species
    are wells and not bimolecular products
    """
    if len(names) == 2:
        # the maximum length between two stationary points
        # is the number of wells+2
        max_length = 5
        n_mol = len(wells) + len(products)
        start = get_index(wells, products, names[0])
        end = get_index(wells, products, names[1])
        # make a graph out of the connectivity
        # nodes of the graph
        nodes = [i for i in range(n_mol)]
        G = nx.Graph()
        G.add_nodes_from(nodes)
        # add the edges of the graph
        for i, ci in enumerate(conn):
            for j, cij in enumerate(ci):
                if cij > 0:
                    G.add_edge(i, j)
        # list of reaction lists for each pathway
        paths = nx.all_simple_paths(G, start, end, cutoff=max_length)
        rxns = []
        for path in paths:
            if is_pathway(wells, products, path, names):
                rxns.append(get_pathway(wells, products, reactions, path, names))
        return rxns
    else:
        logging.error('Cannot find a lowest path if the number of species is not 2')
        logging.error('Found species: ' + ' '.join(names))


def get_index(wells, products, name):
    try:
        i = wells.index(name)
    except ValueError:
        try:
            i = products.index(name) + len(wells)
        except ValueError:
            logging.error('Could not find reactant ' + name)
            sys.exit(-1)
    return i


def get_name(wells, products, i):
    if i < len(wells):
        name = wells[i]
    else:
        name = products[i - len(wells)]
    return name


def get_pathway(wells, products, reactions, ins, names):
    """
    Return the list of reactions between the species in
    the names, according to the instance ins
    """
    # list of reactions
    rxns = []
    for index, i in enumerate(ins[:-1]):
        j = ins[index + 1]
        rxns.append(get_reaction(wells, products, reactions, i, j))
    return rxns


def get_reaction(wells, products, reactions, i, j):
    """
    method to get a reaction in the reactions list
    according to the indices i and j which correspond
    to the index in wells or products
    """
    name_1 = get_name(wells, products, i)
    name_2 = get_name(wells, products, j)
    for rxn in reactions:
        reac_name = rxn[0]
        prod_name = '_'.join(sorted(rxn[2]))
        if ((name_1 == reac_name and name_2 == prod_name) or
                (name_2 == reac_name and name_1 == prod_name)):
            return rxn
    return None


def is_pathway(wells, products, ins, names):
    """
    Method to check if the instance ins
    corresponds to a pathway between the species
    in the names list
    """
    # check of all intermediate species are wells
    if all([insi < len(wells) for insi in ins[1:-1]]):
        name_1 = get_name(wells, products, ins[0])
        name_2 = get_name(wells, products, ins[-1])
        # check if the names correspond
        if ((name_1 == names[0] and name_2 == names[1]) or
                (name_2 == names[0] and name_1 == names[1])):
            return 1
    return 0


def copy_from_kinbot(well, dirname):
    dirname = dirname + '/'
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    copy_tree(well + '/' + dirname, dirname)


def get_rxn(prods, rxns):
    for rxn in rxns:
        if prods == '_'.join(sorted(rxn[2])):
            return rxn


def create_short_names(wells, products, reactions):
    """
    Create a short name for all the wells, all the
    bimolecular products and all the transition states
    """
    # short names of the wells
    # key: chemid
    # value: short name
    well_short = {}
    # short names of the products
    # key: sorted chemids with underscore
    # value: short name
    pr_short = {}
    # short names of the fragments
    # key: chemid
    # value: short name
    fr_short = {}
    # short name of the ts's
    # key: reaction name (chemid, type and instance)
    # value: short name
    ts_short = {}

    for well in wells:
        if well not in well_short:
            short_name = 'w_' + str(len(well_short) + 1)
            well_short[well] = short_name
    for prod in products:
        if prod not in pr_short:
            short_name = 'pr_' + str(len(pr_short) + 1)
            pr_short[prod] = short_name
        for frag in prod.split('_'):
            if frag not in fr_short:
                short_name = 'fr_' + str(len(fr_short) + 1)
                fr_short[frag] = short_name
    for rxn in reactions:
        if rxn[1] not in ts_short:
            short_name = 'rxn_' + str(len(ts_short) + 1)
            ts_short[rxn[1]] = short_name
    return well_short, pr_short, fr_short, ts_short


def write_header(par, well0):
    """
    Create the header block for MESS
    """
    # Read the header template
    header_file = pkg_resources.resource_filename('tpl', 'mess_header.tpl')
    with open(header_file) as f:
        tpl = f.read()
    header = tpl.format(TemperatureList=' '.join([str(ti) for ti in par.par['TemperatureList']]),
                        PressureList=' '.join([str(pi) for pi in par.par['PressureList']]),
                        EnergyStepOverTemperature=par.par['EnergyStepOverTemperature'],
                        ExcessEnergyOverTemperature=par.par['ExcessEnergyOverTemperature'],
                        ModelEnergyLimit=par.par['ModelEnergyLimit'],
                        CalculationMethod=par.par['CalculationMethod'],
                        ChemicalEigenvalueMax=par.par['ChemicalEigenvalueMax'],
                        Reactant=well0,
                        EnergyRelaxationFactor=par.par['EnergyRelaxationFactor'],
                        EnergyRelaxationPower=par.par['EnergyRelaxationPower'],
                        EnergyRelaxationExponentCutoff=par.par['EnergyRelaxationExponentCutoff'],
                        Epsilons=' '.join([str(ei) for ei in par.par['Epsilons']]),
                        Sigmas=' '.join([str(si) for si in par.par['Sigmas']]),
                        Masses=' '.join([str(mi) for mi in par.par['Masses']]))
    return header


def create_mess_input(par, wells, products, reactions,
                      well_energies, prod_energies, parent):
    """
    When calculating a full pes, the files from the separate wells
    are read and concatenated into one file
    Two things per file need to be updated
    1. the names of all the wells, bimolecular products and ts's
    2. all the zpe corrected energies
    """
    # generate short names for all startionary points
    well_short, pr_short, fr_short, ts_short = create_short_names(wells,
                                                                  products,
                                                                  reactions)
    # list of the strings to write to mess input file
    s = []
    # write the header
    s.append(write_header(par, well_short[wells[0]]))

    # write the wells
    s.append('######################')
    s.append('# WELLS')
    s.append('######################')
    for well in wells:
        name = well_short[well] + ' ! ' + well
        energy = well_energies[well]
        with open(parent[well] + '/' + well + '.mess') as f:
            s.append(f.read().format(name=name, zeroenergy=energy))
        s.append('!****************************************')

    # write the products
    s.append('######################')
    s.append('# BIMOLECULAR PRODUCTS')
    s.append('######################')
    for prod in products:
        name = pr_short[prod] + ' ! ' + prod
        energy = prod_energies[prod]
        fr_names = {}
        for fr in prod.split('_'):
            key = 'fr_name_{}'.format(fr)
            value = fr_short[fr] + ' ! ' + fr
            fr_names[key] = value
        with open(parent[prod] + '/' + prod + '.mess') as f:
            s.append(f.read().format(name=name,
                                     ground_energy=energy,
                                     **fr_names))
        s.append('!****************************************')

    # write the barrier
    s.append('######################')
    s.append('# BARRIERS')
    s.append('######################')
    for rxn in reactions:
        name = [ts_short[rxn[1]]]
        name.append(well_short[rxn[0]])
        if len(rxn[2]) == 1:
            name.append(well_short[rxn[2][0]])
        else:
            name.append(pr_short['_'.join(sorted(rxn[2]))])
        name.append('!')
        name.append(rxn[1])
        energy = rxn[3]
        with open(rxn[0] + '/' + rxn[1] + '.mess') as f:
            s.append(f.read().format(name=' '.join(name), zeroenergy=energy))
        s.append('!****************************************')

    # add last end statement
    s.append('!****************************************')
    s.append('End ! end kinetics\n')

    if not os.path.exists('me'):
        os.mkdir('me')

    # write everything to a file
    with open('me/mess.inp', 'w') as f:
        f.write('\n'.join(s))

    dummy = StationaryPoint('dummy',
                            par.par['charge'],
                            par.par['mult'],
                            smiles=par.par['smiles'],
                            structure=par.par['structure'])

    mess = MESS(par, dummy)
    if par.par['me']:
        mess.run()


def create_pesviewer_input(par, wells, products, reactions,
                           well_energies, prod_energies, highlight):
    """
    highlight: list of reaction names that need a red highlight
    """
    # delete the im_extent and xval files
    try:
        os.remove('{}_xval.txt'.format(par.par['title']))
    except OSError:
        pass
    try:
        os.remove('{}_im_extent.txt'.format(par.par['title']))
    except OSError:
        pass

    if highlight is None:
        highlight = []

    well_lines = []
    for well in wells:
        energy = well_energies[well]
        well_lines.append('{} {:.2f}'.format(well, energy))

    bimol_lines = []
    for prods in products:
        energy = prod_energies[prods]
        bimol_lines.append('{} {:.2f}'.format(prods, energy))

    ts_lines = []
    for rxn in reactions:
        high = ''
        if rxn[1] in highlight:
            high = 'red'
        prod_name = '_'.join(sorted(rxn[2]))
        ts_lines.append('{} {:.2f} {} {} {}'.format(rxn[1],
                                                    rxn[3],
                                                    rxn[0],
                                                    prod_name,
                                                    high))

    well_lines = '\n'.join(well_lines)
    bimol_lines = '\n'.join(bimol_lines)
    ts_lines = '\n'.join(ts_lines)

    # write everything to a file
    fname = 'pesviewer.inp'
    template_file_path = pkg_resources.resource_filename('tpl', fname + '.tpl')
    with open(template_file_path) as template_file:
        template = template_file.read()
    template = template.format(id=par.par['title'],
                               wells=well_lines,
                               bimolecs=bimol_lines,
                               ts=ts_lines,
                               barrierless='')
    with open(fname, 'w') as f:
        f.write(template)


def create_graph(wells, products, reactions,
                 well_energies, prod_energies, highlight):
    """
    highlight: list of reaction names that need a red highlight
    """
    if highlight is None:
        highlight = []
    # update the connectivity with the filtered wells, products and reactions
    conn, bars = get_connectivity(wells, products, reactions)

    # get the minimum and maximum well and product energy
    try:
        minimum = min(min(well_energies.values()),
                      min(prod_energies.values()))
        maximum = max(max(well_energies.values()),
                      max(prod_energies.values()))
    except ValueError:
        # list of products can be empty, but list of wells not
        minimum = min(well_energies.values())
        maximum = max(well_energies.values())
    # define the inveresly proportial weights function
    max_size = 400
    min_size = 100
    slope = (min_size - max_size) / (maximum - minimum)
    offset = max_size - minimum * slope
    # define the graph nodes
    nodes = [i for i, wi in enumerate(wells)]
    nodes += [len(wells) + i for i, pi in enumerate(products)]
    # size of the nodes from the weights
    node_size = [slope * well_energies[wi] + offset for wi in wells]
    node_size += [slope * prod_energies[pi] + offset for pi in products]
    # color nodes and wells differently
    node_color = ['lightskyblue' for wi in wells]
    node_color += ['lightcoral' for pi in products]
    # labels of the wells and products
    labels = {}
    name_dict = {}
    for i, wi in enumerate(wells):
        labels[i] = 'w{}'.format(i+1)
        name_dict[labels[i]] = wi
    for i, pi in enumerate(products):
        labels[i + len(wells)] = 'b{}'.format(i+1)
        name_dict[labels[i + len(wells)]] = pi
    # write the labels to a file
    with open('species_dict.txt', 'w') as f:
        lines = []
        for name in sorted(name_dict.keys()):
            lines.append('{}  {}'.format(name, name_dict[name]))
        f.write('\n'.join(lines))
    # make a graph object
    G = nx.Graph()
    # add the nodes
    for i, node in enumerate(nodes):
        G.add_node(node, weight=node_size[i])

    # define the inversely proportional weights for the lines
    minimum = min(rxn[3] for rxn in reactions)
    maximum = max(rxn[3] for rxn in reactions)
    max_size = 5
    min_size = 0.5
    try:
        slope = (min_size - max_size) / (maximum - minimum)
    except:
        slope = 1.
    offset = max_size - minimum * slope

    # add the edges
    for i, ci in enumerate(conn):
        for j, cij in enumerate(ci):
            if cij > 0:
                weight = slope * bars[i][j] + offset
                G.add_edge(i, j, weight=weight)
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]

    # position the nodes
    pos = nx.spring_layout(G, scale=1)

    # make the matplotlib figure
    plt.figure(figsize=(8, 8))
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights)
    nx.draw_networkx_nodes(G,
                           pos,
                           nodelist=G.nodes(),
                           node_size=node_size,
                           node_color=node_color)
    nx.draw_networkx_labels(G, pos, labels, font_size=8)
    plt.axis('off')
    plt.savefig('graph.png')


def get_energy(dir, job, ts, high_level, mp2=0):
    db = connect(dir + '/kinbot.db')
    if ts:
        j = job
    else:
        j = job + '_well'
    if mp2:
        j += '_mp2'
    if high_level:
        j += '_high'
    rows = db.select(name=j)
    for row in rows:
        if hasattr(row, 'data'):
            energy = row.data.get('energy')
    try:
        # ase energies are always in ev, convert to hartree
        energy *= constants.EVtoHARTREE
    except UnboundLocalError:
        # this happens when the job is not found in the database
        logging.error('Could not find {} in directory {}'.format(job, dir))
        logging.error('Exiting...')
        sys.exit(-1)
    return energy


def get_l3energy(job, par):
    """ 
    Get the L3, single-point energies. 
    This is not object oriented.
    """

    if par.par['single_point_qc'] == 'molpro':
        if os.path.exists('molpro/' + job + '.out'):
            with open('molpro/' + job + '.out') as f:
                lines = f.readlines()
                for index, line in enumerate(reversed(lines)):
                    if ('SETTING ' + par.par['single_point_key']) in line:
                        return 1, float(line.split()[2])  # energy was found
                    else:
                        return 0, -1  # the job not yet done
        else:
            return 0, -1  # job not yet started to run

    
def get_zpe(dir, job, ts, high_level, mp2=0):
    db = connect(dir + '/kinbot.db')
    if ts:
        j = job
    else:
        j = job + '_well'
    if mp2:
        j += '_mp2'
    if high_level:
        j += '_high'
    rows = db.select(name=j)
    for row in rows:
        if hasattr(row, 'data'):
            zpe = row.data.get('zpe')
    return zpe


def check_status(job, pid):
    command = ['ps', '-u', 'root', '-N', '-o', 'pid,s,user,%cpu,%mem,etime,args']
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = out.decode()
    lines = out.split('\n')
    for line in lines:
        if len(line) > 0:
            if str(pid) == line.split()[0]:
                return 1
    return 0


def submit_job(chemid, par):
    """
    Submit a kinbot run using subprocess and return the pid
    """
    command = ["kinbot", chemid + ".json", "&"]
    # purge previous summary and monitor files, so that pes doesn't think
    # everything is done
    # relevant if jobs are killed
    try:
        os.system('rm -f {dir}/summary_*.out'.format(dir=chemid))
    except OSError:
        pass
    try:
        os.system('rm -f {dir}/kinbot_monitor.out'.format(dir=chemid))
    except OSError:
        pass
 
    if par.par['queue_template'] != '':
        shutil.copyfile('{}'.format(par.par['queue_template']), '{}/{}'.format(chemid, par.par['queue_template']))
    outfile = open('{dir}/kinbot.out'.format(dir=chemid), 'w')
    errfile = open('{dir}/kinbot.err'.format(dir=chemid), 'w')
    process = subprocess.Popen(command,
                               cwd=chemid,
                               stdout=outfile,
                               stdin=subprocess.PIPE,
                               stderr=errfile)
    time.sleep(1)
    pid = process.pid
    return pid


def write_input(par, species, threshold, root):
    # directory for this particular species
    dir = root + '/' + str(species.chemid) + '/'
    if not os.path.exists(dir):
        os.makedirs(dir)

    # make a new parameters instance and overwrite some keys
    par2 = Parameters(par.input_file)
    # overwrite the title
    par2.par['title'] = str(species.chemid)
    # make a structure vector and overwrite the par structure
    structure = []
    for at in range(species.natom):
        pos = species.geom[at]
        sym = species.atom[at]
        structure += [sym, pos[0], pos[1], pos[2]]
    par2.par['structure'] = structure
    # delete the par smiles
    par2.par['smiles'] = ''
    # overwrite the barrier treshold
    par2.par['barrier_threshold'] = threshold
    # set the pes option to 1
    par2.par['pes'] = 1

    file_name = dir + str(species.chemid) + '.json'
    with open(file_name, 'w') as outfile:
        json.dump(par2.par, outfile, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()

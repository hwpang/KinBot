import numpy as np
import os
import sys
import logging
import subprocess

def write_cnn_input(reaction):
        """
        Write the input for the neural network which predicts the TS bond lengths.
        The input is a reaction instance from one of the reaction families.
        """
        cnn = 'cnn/{}/'.format(reaction.instance_name)
        if not os.path.exists(cnn):
            os.makedirs(cnn)
        with open(cnn + 'atom.inp', 'w') as f:
            f.write(str(reaction.species.natom) + '\n')
            f.write(' '.join(reaction.species.atom))
        with open(cnn + 'reactant_bond.inp', 'w') as f:
            f.write('\n'.join([' '.join(str(bij) for bij in bi) for bi in reaction.species.bond]))
        with open(cnn + 'reactant_dist.inp', 'w') as f:
            f.write('\n'.join([' '.join(str(bij) for bij in bi) for bi in reaction.species.dist]))
        with open(cnn + 'product_bond.inp', 'w') as f:
            f.write('\n'.join([' '.join(str(bij) for bij in bi) for bi in reaction.product_bonds]))

def predict():
    """
    This method write files the the necessary information to predict the 
    TS bond lenghts using a CNN
    It runs a separate python script with the input files that does the 
    actual prediction
    And it reads the output matrix. 
    
    """
    # do the calculation
    command = 'python cnn.py'
    process = subprocess.Popen(command,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    out,err = process.communicate()
    out = out.decode()
    err = err.decode()
    logging.debug('Standard output of command:' + out)
    logging.debug('Standard error of command:' + err)

def get_cnn_results(reaction):
    """
    Read the output of the neural network for a reaction
    """
    # read the output file
    cnn = 'cnn/{}/'.format(reaction.instance_name)
    with open(cnn + 'ts_bond.out') as f:
        lines = f.read().split('\n')
    ts_bonds = np.zeros((reaction.species.natom, reaction.species.natom))
    for i, line in enumerate(lines):
        for j, value in enumerate(line.split()):
            ts_bonds[i][j] = float(value)
    
    return ts_bonds
    



import numpy as np
import os
import sys
import logging
import subprocess

def predict(name, reactant, product_bond_mx):
    """
    This method write files the the necessary information to predict the 
    TS bond lenghts using a CNN
    It runs a separate python script with the input files that does the 
    actual prediction
    And it reads the output matrix. 
    
    """
    # write the input files
    cnn = 'cnn/{}/'.format(name)
    if not os.path.exists(cnn):
        os.makedirs(cnn)
    with open(cnn + 'atom.inp', 'w') as f:
        f.write(str(reactant.natom) + '\n')
        f.write(' '.join(reactant.atom))
    with open(cnn + 'reactant_bond.inp', 'w') as f:
        f.write('\n'.join([' '.join(str(bij) for bij in bi) for bi in reactant.bond]))
    with open(cnn + 'reactant_dist.inp', 'w') as f:
        f.write('\n'.join([' '.join(str(bij) for bij in bi) for bi in reactant.dist]))
    with open(cnn + 'product_bond.inp', 'w') as f:
        f.write('\n'.join([' '.join(str(bij) for bij in bi) for bi in product_bond_mx]))

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
    
    # read the output file
    with open(cnn + 'ts_bond.out') as f:
        lines = f.read().split('\n')
    ts_bonds = np.zeros((reactant.natom, reactant.natom))
    for i, line in enumerate(lines):
        for j, value in enumerate(line.split()):
            ts_bonds[i][j] = float(value)
    
    return ts_bonds
    



import numpy as np
import importlib
import os
import logging
import subprocess

try:
    from keras.models import load_model
except:
    logging.warn('Keras was not imported, neural network will not work.')

from kinbot import stationary_pt

class CNN():
    def __init__(self):
        # Positional (one-hot) encoding of bond types
        self.bond_types = {'CC':0, 'CO':1, 'CH':2, 'OO':3, 'HO':4, 'HH':5}

        # Positional (one-hot) encoding of bond changes:
        # Example: 1-0 = breaking of single bond; 1-2 = increase in bond order
        # doesn't exist in dataset, '3-2':7}
        self.bond_changes = {'1-0':0, '0-1':1, '2-1':2, '1-2':3, '2-0':4, '0-2':5, '2-3':6} 

        # Positional (one-hot) encoding of atom neighbors
        self.neighbors = {'H':0, 'C':1, 'O':2}  #no nitrogen atm, 'N':3}

        # Input normalization
        self.max_mol_size = 9.0
        self.max_num_changes = 9  #9 changes, 0 changes does not occur.
        self.max_size = 25
        self.number_of_permutations = 10
        
        # Read the model
        self.model_path = '/kyukon/data/gent/408/vsc40810/KinBot/kinbot/cnn/CNN_v6.model'
        self.model = load_model(self.model_path)

    def predict(self, reactant, product_bond_mx):
        """
        This method predicts the transition state bond distances of the bonds
        that change in connectivity throughout the reaction. 
        The stationary point object of the reactant and the bond matrix of the product
        are needed as input. 
        It is important to have the product bond matrix with the same atom order as the reactant.
        
        Q:
        1. One prediction or different predictions with different permutations to average?
        """

        natom = reactant.natom
        atom = np.array(reactant.atom)
        reactant_bond_mx = reactant.bond
        ts_bonds = product_bond_mx - reactant_bond_mx
        ts_bonds[ts_bonds!=0] = 1
        
        # Make a matrix with the bond types of the atoms
        reactant_bond_types = []
        for i, ati in enumerate(atom):
            bond_type_line = []
            for j, atj in enumerate(atom):
                bond_type_line.append(''.join(sorted([ati,atj])))
            bond_type_line = np.array(bond_type_line)
            reactant_bond_types.append(bond_type_line)
        
        # Make the distance matrix with the correct dimenstions
        reactant_dist_mx = np.zeros((self.max_size, self.max_size))
        for i in range(self.max_size):
            for j in range(self.max_size):
                if i < natom and j < natom:
                    reactant_dist_mx[i][j] = reactant.dist[i][j]
        reactant_dist_mx = reactant_dist_mx.reshape((self.max_size, self.max_size, 1))
        
        inputA = []
        inputB = []
        
        ret = np.zeros((natom, natom))    
        
        # Split transition state into different datapoints for each bond 
        # that has been identified as ts bond
        for i in range(natom - 1):
            for j in range(i+1, natom):
                if ts_bonds[i][j] != 0:
                    # Initialize all feature vectors
                    # Transition state bond types
                    tsbt = [0 for i in range(len(self.bond_types))]
                    # Transition state bond change
                    tsbc = [0 for i in range(len(self.bond_changes))]
                    # Total number of changes in transition state
                    tsbs = [0 for i in range(self.max_num_changes)]
                    # Valence of bond atom1 [1,2,3,4]
                    vals1= [0 for i in range(4)]
                    # Valence of bond atom2 [1,2,3,4]
                    vals2= [0 for i in range(4)]
                    # Neighbors bond atom1 [H, C, O]
                    nct1 = [0 for i in range(len(self.neighbors))]
                    # Neighbors bond atom2 [H, C, O]
                    nct2 = [0 for i in range(len(self.neighbors))]
                    
                    # Original (reactant) bond length of bond that became TS bond
                    ts_length = reactant_dist_mx[i][j]
                    # Bond type of ts bond
                    ts_bond_type = self.bond_types[reactant_bond_types[i][j]]
                    # Set one-hot encoded type to 1
                    tsbt[ts_bond_type] = 1.0
                    # Get valences of atom 1&2 (via number of connected bonds in adjacency matrix)
                    val1 = sum(reactant_bond_mx[i])
                    val2 = sum(reactant_bond_mx[j])
                    # One-hot valence position set to 1
                    vals1[val1-1]=1.0
                    vals2[val2-1]=1.0
                    # Keep track of the neighbors
                    for jj in range(natom):
                        # Atoms[jj] gives 'C', 'H' or 'O', which is a key in the neighbors dictionary,
                        # linked to the correct position in the nct1 one-hot feature vector for the
                        # number of neighbors
                        if reactant_bond_mx[i][jj]!=0:
                            nct1[self.neighbors[atom[jj]]]+=1.0
                        if reactant_bond_mx[j][jj]!=0:
                            nct2[self.neighbors[atom[jj]]]+=1.0

                    # Make bond change string
                    ts_bond_change = self.bond_changes["{}-{}".format(reactant_bond_mx[i][j],product_bond_mx[i][j])]
                    # One-hot bond change position set to 1
                    tsbc[ts_bond_change] = 1.0

                    ts_bond_sum = int(np.sum(ts_bonds)/2)
                    # One-hot encoding for the TOTAL number of changes in the transition state
                    # so input is not fully independent of other inputs for the same 
                    # reaction/transition state
                    tsbs[ts_bond_sum-1] = 1.0
                    
                    #Combine everything into a single feature vector
                    iA = [ts_length]
                    iA.extend(tsbt)
                    iA.extend(tsbc)
                    iA.extend(tsbs)
                    iA.extend(vals1)
                    iA.extend(vals2)
                    iA.extend(nct1)
                    iA.extend(nct2)
                    iA = np.array(iA)
                    iB = reactant_dist_mx.copy()
                    preds = self.model.predict([[iA], [iB]])[0][0]
                    #preds = 0
                    print(ts_length + preds)
                    
                    ret[i][j] = ts_length + preds
                    ret[j][i] = ts_length + preds
                    
        print(ret)
        return ret

if __name__ == "__main__":
    structure = ['C', 1.149496,   0.171148,   0.154969,
                 'H', 0.582461,  -0.501854,   0.792486,
                 'H', 0.623938,   1.012385,  -0.286640,
                 'C', 2.633431,   0.032105,   0.027857,
                 'C', 3.136464,  -1.409330,   0.237072,
                 'H', 2.960643,   0.399948,  -0.956684,
                 'H', 3.141534,   0.690537,   0.759810,
                 'H', 2.710277,  -2.084498,  -0.513710,
                 'H', 4.228576,  -1.463407,   0.163698,
                 'H', 2.846698,  -1.783635,   1.226773]

    reactant = stationary_pt.StationaryPoint('propy', 0, 2, structure=structure)
    reactant.characterize()

    product_bond_mx = [[0, 1, 1, 1, 0, 0, 1, 0, 0, 0],
                       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [1, 0, 0, 0, 1, 1, 0, 0, 0, 0],
                       [0, 0, 0, 1, 0, 0, 0, 1, 1, 1],
                       [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                       [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]]
    product_bond_mx = np.array(product_bond_mx)
    cnn = CNN()
    cnn.predict(reactant, product_bond_mx)
    print(cnn.model)
    print('TEST 4')
    arguments = ['lsc']
    print('TEST 5')
    print(arguments)
    print(subprocess.PIPE)
    print('aaa')
    importlib.reload(subprocess)
    process = subprocess.Popen(args = arguments, close_fds=True)
                               #shell=True,
                               #stdout=subprocess.PIPE,
                               #stdin=subprocess.PIPE,
                               #stderr=subprocess.PIPE)
    print('TEST 6')
    #out,err = process.communicate()
    #out = out.decode()
    #err = err.decode()
    print('TEST 7')
    #print(out)
    #print(err)
    print('TEST 8')
    
    
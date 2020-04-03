import numpy as np
import copy
import time

from kinbot import reac_family
from kinbot import geometry

class IntraRAddExocyclicF:
    
    def __init__(self,species,qc,par,instance,instance_name):
        #st_pt of the reactant
        self.species = species
        #st_pt of the ts
        self.ts = None
        #st_pt of the product(s)
        self.products = []
        #bond matrix of the products
        self.product_bonds = [] 
        
        #optimization objects
        self.ts_opt = None
        self.prod_opt = []
        
        self.qc = qc
        self.par = par
        
        #indices of the reactive atoms
        self.instance = instance
        #name of the reaction
        self.instance_name = instance_name
        
        #maximum number of steps for this reaction family
        self.max_step = 14
        #do a scan?
        self.scan = 0
        #skip the first 12 steps in case the instance has a length of 3?
        self.skip = 1

    def get_constraints(self,step, geom):
        """
        There are three types of constraints:
        1. fix a coordinate to the current value
        2. change a coordinate and fix is to the new value
        3. release a coordinate (only for gaussian)
        """
        fix = []
        change = []
        release = []
        if step < self.max_step:
            #fix all the bond lengths
            for i in range(self.species.natom - 1):
                for j in range(i+1, self.species.natom):
                    if self.species.bond[i][j] > 0:
                        fix.append([i+1,j+1])
        if step < 12:
            new_dihs = geometry.new_ring_dihedrals(self.species, self.instance, step, 12)
            for dih in range(len(self.instance)-4): #do not use the last atom
                constraint = []
                for i in range(4):
                    constraint.append(self.instance[dih+i] + 1)
                constraint.append(new_dihs[dih])
                change.append(constraint)
        elif step == 12:
            for dih in range(len(self.instance)-3):
                f = []
                for i in range(4):
                    f.append(self.instance[dih+i] + 1)
                fix.append(f)
            for angle in range(len(self.instance)-2):
                constraint = []
                for i in range(3):
                    constraint.append(self.instance[angle+i] + 1)
                constraint.append(180. * (len(self.instance)-2.) / len(self.instance))
                change.append(constraint)
        elif step == 13:
            for angle in range(len(self.instance)-2):
                constraint = []
                for i in range(3):
                    constraint.append(self.instance[angle+i] + 1)
                release.append(constraint)
            for dih in range(len(self.instance)-3):  
                constraint = []
                for i in range(4):
                    constraint.append(self.instance[dih+i] + 1)
                release.append(constraint)
                
            fval = 2.2
            constraint = [self.instance[0] + 1,self.instance[-2] + 1,fval]
            change.append(constraint)

        #remove the bonds from the fix if they are in another constaint
        for c in change:
            if len(c) == 3:
                index = -1
                for i,fi in enumerate(fix):
                    if len(fi) == 2:
                        if sorted(fi) == sorted(c[:2]):
                            index = i
                if index > -1:
                    del fix[index]
        
        return step, fix, change, release


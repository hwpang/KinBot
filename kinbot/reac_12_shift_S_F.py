from kinbot import constants
from kinbot import geometry
from kinbot.reac_General import GeneralReac

class S12ShiftF(GeneralReac):
    max_step = 12
    scan = 0
    skip = 0
        

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)
           
        if step < self.max_step:
            final_dist1 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[1]] + self.species.atom[self.instance[2]]))]
            val = geometry.new_bond_length(self.species, self.instance[1], self.instance[2], step, self.max_step, final_dist1, geom)
            self.set_bond(1, 2, val, change)
            
            final_dist2 = constants.st_bond[''.join(sorted(self.species.atom[self.instance[0]] + self.species.atom[self.instance[2]]))]
            val = geometry.new_bond_length(self.species,self.instance[0], self.instance[2], step, self.max_step, final_dist2, geom)
            self.set_bond(0, 2, val, change)

        self.clean_constraints(change, fix)        
        
        return step, fix, change, release

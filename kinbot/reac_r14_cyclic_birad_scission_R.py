from kinbot.reac_General import GeneralReac
from kinbot import geometry

class R14CyclicBiradScission(GeneralReac):
    max_step = 22
    scan = 0
    skip = 0
    dihstep = 12
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []

        self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change)

        elif step == self.dihstep:  # originally was < 12, must be a typo
            self.release_dihedrals(release)
                
            fval = 1.8
            val = geometry.new_bond_length(self.species,self.instance[0],self.instance[-1],step - 11,10,fval,geom)
            self.set_bond(0, -1, val, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release

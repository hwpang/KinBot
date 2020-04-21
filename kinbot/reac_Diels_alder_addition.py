class DielsAlder(GeneralReac):
    max_step = 1
    scan = 0
    skip = 0


    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step ==0 :
            fval = 2.2
            self.set_bond(2, 3, fval)
            self.set_bond(4, 5, fval)
        
        self.clean_constraints(change, fix)
        
        return step, fix, change, release

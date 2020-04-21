class IntraOHMigration(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    dihstep = max_step - 2
    

    def get_constraints(self,step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            self.fix_bonds(fix)

        if step < self.dihstep:
            self.set_dihedrals(change, step)

        elif step == self.dihstep:
            self.fix_dihedrals(fix)
            self.set_angles(change) 

        elif step == self.dihstep + 1:
            self.release_angles(release)
            self.release_dihedrals(release)
            
            if self.species.atom[self.instance[0]] == 'C':
                fval1 = 2.0
                fval2 = 1.7
            else:
                fval1 = 1.7 
                fval2 = 2.0

            self.set_bond(0, -1, fval1, change)
            self.set_bond(-2, -1, fval2, change)

        self.clean_constraints(change, fix)
        
        return step, fix, change, release

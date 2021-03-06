from kinbot.reac_General import GeneralReac


class CpdHMigration(GeneralReac):
    max_step = 14
    scan = 0
    skip = 1
    

    def get_constraints(self, step, geom):
        fix = []
        change = []
        release = []
        if step < self.max_step:
            fix_bonds(fix)

        if step == 0:
            self.set_angle_single(-2, -1, 0, 70., change)

        if step == 1:
            fval = 1.35
            self.set_bond(-2, -1, fval, change)
            self.set_bond(0, -1, fval, change)
            
        self.clean_constraints(change, fix)

        return step, fix, change, release

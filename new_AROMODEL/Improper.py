# Improper Class

class Improper(object):
    """
    Class defining an improper interactions between 4 atoms
    This is mainly just for aromatic carbons to enforce planarity
    instance variables
        Type = int
        atom_1 = Atom object  # was master
        atom_2 = Atom object  # was slave1
        atom_3 = Atom object  # was slave2
        atom_4 = Atom object  # was slave3 ( This one is the "odd one out")
        Ki = float (Force/degree)
        improper_Eq = float (degrees)  
    """

    def __init__(self, a1, a2, a3, a4, Ki, improper_eq, improper_id, d=-1, n=2):
        self.atom_1 = a1
        self.atom_2 = a2
        self.atom_3 = a3
        self.atom_4 = a4 # outside the ring
        self.Ki = Ki
        self.improper_eq = improper_eq
        self.improper_id = improper_id
        self.lammps_type = 1
        self.system_id = 0
        self.d = d
        self.n = n
        return

    def Compare_Impropers(self,atom_2,atom_1,atom_3,atom_4):
        if (self.atom_1 == atom_1 and self.atom_2 == atom_2 and self.atom_3 == atom_3 and self.atom_4 == atom_4) or (self.atom_1 == atom_1 and self.atom_2 == atom_3 and self.atom_3 == atom_2 and self.atom_4 == atom_4):
            return True
        else:
            return False
# Improper Class

class Improper(object):
    """
    Class defining an improper interactions between 4 atoms
    This is mainly just for aromatic carbons to enforce planarity
    instance variables
        Type = int
        Improper_Main = Atom object  # was master
        Improper_Node1 = Atom object  # was slave1
        Improper_Node2 = Atom object  # was slave2
        Improper_Node3 = Atom object  # was slave3 ( This one is the "odd one out")
        Ki = float (Force/degree)
        improper_Eq = float (degrees)  
    """

    def __init__(self, a1, a2, a3, a4, Ki, improper_eq, improper_id, d=-1, n=2):
        self.Improper_Main = a1
        self.Improper_Node1 = a2
        self.Improper_Node2 = a3
        self.Improper_Node3 = a4 # outside the ring
        self.Ki = Ki
        self.improper_eq = improper_eq
        self.improper_id = improper_id
        self.lammps_type = 1
        self.system_id = 0
        self.d = d
        self.n = n
        return

    def Compare_Impropers(self,atom_2,atom_1,atom_3,atom_4):
        if (self.Improper_Main == atom_1 and self.Improper_Node1 == atom_2 and self.Improper_Node2 == atom_3 and self.Improper_Node3 == atom_4) or (self.Improper_Main == atom_1 and self.Improper_Node1 == atom_3 and self.Improper_Node2 == atom_2 and self.Improper_Node3 == atom_4):
            return True
        else:
            return False
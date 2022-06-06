# Bond Class
import numpy as np

class Bond(object):
    """
    Class defining a bond between atom objects
    instance variables:
        Type = int
        atom_1 = Atom object
        atom_2 = Atom object
        Bond_Length = float
        Bond_Vector = np[3, float]
        req = float (equilibrium radius)
        kb = float
    """
    def __init__(self, atom_1, atom_2, req):
        self.atom_1 = atom_1 # Atom object
        self.atom_2 = atom_2 # Atom object
        self.req = req # distance
        self.kb = 0.0
        self.bond_id = 0
        self.lammps_type = 0
        self.system_id = 0
        self.kb = 0
        self.bond_length = 0
        self.bond_vector = np.zeros(3)
        self.type = 1
        return





# Bond Class
import numpy as np

class Bond(object):
    """
    Class defining a bond between atom objects
    instance variables:
        Type = int
        Bond_Main = Atom object
        Bond_Node = Atom object
        Bond_Length = float
        Bond_Vector = np[3, float]
        req = float (equilibrium radius)
        kb = float
    """
    def __init__(self, Bond_Main, Bond_Node, req):
        self.Bond_Main = Bond_Main # Atom object
        self.Bond_Node = Bond_Node # Atom object
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





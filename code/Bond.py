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
        self.Bond_id = 0
        self.Lammps_Type = 0
        self.System_ID = 0
        self.Bond_Length = 0
        self.Bond_Vector = np.zeros(3)
        self.Type = 1
        return





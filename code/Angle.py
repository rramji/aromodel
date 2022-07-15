class Angle(object):
    """ 
    Class defining an angle between 3-atom objects
    Instance Variables:
        Type = int
        Angle_Main = Atom object
        Angle_Node1 = Atom object 
        Angle_Node2 = Atom object
        angle_Normal = np[3,float]
        angle_Th = float (radians0
        angle_Eq = float (degrees)
        ka = float
    """
    def __init__(self, Angle_Main, Angle_Node1, Angle_Node2, angle_eq):
        self.Angle_Main = Angle_Main #Atom Object (Atom in the middle)
        self.Angle_Node1 = Angle_Node1 #Atom Object
        self.Angle_Node2 = Angle_Node2 #Atom Object
        self.angle_eq = angle_eq
        self.ka = 0.0 # Defined in OPLS.Assign_OPLS
        self.angle_id = 0 # Defined in OPLS.Assign_OPLS
        self.lammps_type = 0 # Defined in Molecule.Assign_Lammps
        self.system_id = 0 # Defined System.Write_LAMMPS_Data()
        return

    def compare_angles(self, node1, main, node2):
    #returns true if the given objects form an identical angle object to the existing object, False if not
        if (self.Angle_Node1 == node1 and self.Angle_Main == main and self.Angle_Node2 == node2) or (self.Angle_Main == main and self.Angle_Node1 == node2 and self.Angle_Node2 == node1):
            return True
        else:
            return False

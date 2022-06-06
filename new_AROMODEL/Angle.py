class Angle(object):
    """ 
    Class defining an angle between 3-atom objects
    Instance Variables:
        Type = int
        atom_1 = Atom object  # was master
        atom_2 = Atom object  # was slave1
        atom_3 = Atom object  # was slave2
        angle_Normal = np[3,float]
        angle_Th = float (radians0
        angle_Eq = float (degrees)
        ka = float
    """
    def __init__(self, atom_2, atom_1, atom_3, angle_eq):
        self.atom_2 = atom_2 #Atom Object (Atom in the middle)
        self.atom_1 = atom_1 #Atom Object
        self.atom_3 = atom_3 #Atom Object
        self.angle_eq = angle_eq
        self.ka = 0.0 # Defined in OPLS.Assign_OPLS
        self.angle_id = 0 # Defined in OPLS.Assign_OPLS
        self.lammps_type = 0 # Defined in Molecule.Assign_Lammps
        self.system_id = 0 # Defined System.Write_LAMMPS_Data()
        return

    def compare_angles(self,atom_1,atom_2,atom_3):
    #returns true if the given objects form an identical angle object to the existing object, False if not
        if (self.atom_1 == atom_1 and self.atom_2 == atom_2 and self.atom_3 == atom_3) or (self.atom_2 == atom_2 and self.atom_1 == atom_3 and self.atom_3 == atom_1):
            return True
        else:
            return False

import numpy as np
from numpy import cos
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

def Fourier_Series(Phi, A1, A2, A3, A4):
    return 0.5*(A1*(1+cos(Phi)) + A2*(1-cos(2.0*Phi)) + A3*(1+cos(3.0*Phi)) + A4*(1-cos(4.0*Phi)))

def Fourier_Series_Val(Phi, Pop):
    return 0.5*(Pop[0]*(1+cos(Phi)) + Pop[1]*(1-cos(2.0*Phi)) + Pop[2]*(1+cos(3.0*Phi)) + Pop[3]*(1-cos(4.0*Phi)))

def Multi_Harmonic(Phi, A1, A2, A3, A4, A5):
    return A1 + A2*cos(Phi) + A3*cos(Phi)**2 + A4*cos(Phi)**3 + A5*cos(Phi)**4

def Multi_Harmonic_Val(Phi, Pop):
    return Pop[0] + Pop[1]*cos(Phi) + Pop[2]*cos(Phi)**2 + Pop[3]*cos(Phi)**3 + Pop[4]*cos(Phi)**4

class Dihedral(object):
    """
    Class defining a dihedral angle between 4 atom objects
    Instance Variables
        Type = int
        Dihedral_Main1 = Atom object  # was Main1
        Dihedral_Main2 = Atom object  # was Main2
        Dihedral_Node1 = Atom object  # was Node1
        Dihedral_Node2= Atom object   # was Node2
        dihedral_eq = float
        Coeffs = [list]

        Dihedral 1 ---- Dihedral 2 ---- Dihedral 3 ---- Dihedral 4
                    (primary atom #1) (primary atom #2)
        

    """

    def __init__(self, Dihedral_Main1, Dihedral_Main2, Dihedral_Node1, Dihedral_Node2, dihedral_eq):
        self.Dihedral_Main1 = Dihedral_Main1 # Atom object
        self.Dihedral_Main2 = Dihedral_Main2 # Atom object
        self.Dihedral_Node1 = Dihedral_Node1
        self.Dihedral_Node2 = Dihedral_Node2
        self.Dihedral_Eq = dihedral_eq # Assigned in OPLS.Assign_OPLS
        self.Coeffs = [] # Assigned in OPLS.Assign_OPLS
        self.dihedral_ID = 0 # Assigned in OPLS.Assign_OPLS
        self.Lammps_Type = 0 # Assigned in Molecule.Assign_Lammps
        self.System_ID = 0 # Assigned in System.Write_LAMMPS_Data()
        self.Style = "" #Assigned in OPLS.Assign_OPLS
        return
        

    def Fit_Parameters(self, dihedral_energy, dihedral_angles):
        Pop, Pco = curve_fit(Fourier_Series, dihedral_angles, dihedral_energy)
        plt.plot(dihedral_angles, dihedral_energy, 'o', label= "MP2 Derived")
        plt.plot(dihedral_angles, Fourier_Series_Val(dihedral_angles, Pop), '-', label="Fourier Fit")
        plt.ylabel('Relative Energy (kcal/mol)', fontsize=25)
        plt.xlabel('dihedral angle (Radians)', fontsize=25)
        plt.xlim((dihedral_angles[0],dihedral_angles[-1]))
        plt.axhline(y=0.0, linestyle='--')
        plt.legend(frameon=False)
        plt.savefig('MP2_dihedral_%d.png' % self.dihedral_id)
        plt.close()
        self.coeffs = Pop
        print(self.coeffs)
        return

    def Compare_Dihedrals(self,atom_1,atom_2,atom_3,atom_4):
        if (self.Dihedral_Node1 == atom_1 and self.Dihedral_Main1 == atom_2 and self.Dihedral_Main2 == atom_3 and self.Dihedral_Node2 == atom_4) or (self.Dihedral_Node1 == atom_4 and self.Dihedral_Main1 == atom_3 and self.Dihedral_Main2 == atom_2 and self.Dihedral_Node2 == atom_1):
            return True
        else:
            return False






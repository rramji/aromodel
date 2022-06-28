# Import relevant modules
import numpy as np
from matplotlib import pyplot as plt
import math
import pickle
import Polymer
import matplotlib as mpl
from pylab import rcParams
import time
import copy


rcParams['figure.figsize'] = 14, 12

mpl.rcParams['axes.color_cycle'] = ['b', 'k','r',  'c', 'y', 'm']

# Define Classes
class Atom(object):
    """
    class defining an atom
    instance variables: Type, id,  Mol_ID, position[3], image_flags[3]
    """
    
    def __init__(self, id, Type, Mol_ID, position, image_flags):
        self.Type = Type
        self.id = id
        self.position = np.asarray(position, dtype= float)
        self.image_flags = np.asarray(image_flags, dtype= int)
        self.Mol_ID = Mol_ID
        self.unwrapped_position = np.zeros(3, dtype =float)
        self.Mass = Polymer.Mass[Type]
        return

    def Print_Info(self):
        print("Atom ID = %d, Type = %d\n, Mol_ID = %d" % (self.id, self.Type, self.Mol_ID))
        print("Position")
        print(self.position)
        print(self.image_flags)
        print("-----------------------------------------")
        return

class CG_Bead(object):
    """
        Class defining a coarse-grained grouping of atoms
        instance variables: CG_ID, CG_Type, Atom_list,  CG_Mass, COM
    """
    def __init__(self, CG_ID, CG_Type, CG_Num, Atom_List, Mol_Num):
        self.CG_ID = CG_ID
        self.CG_Type = CG_Type
        self.CG_Num = CG_Num
        self.Atom_List = Atom_List
        self.COM = np.zeros(3, dtype=float)
        self.COM_PBC = np.zeros(3,dtype=float)
        self.IF = np.zeros(3,dtype=int)
        self.CG_Mass = 0.0
        Mass_Weighted_Sum = np.zeros(3, dtype=float)
        self.Orientation_Vector = np.zeros(3,dtype=float)
        self.basis = np.zeros((3,3),dtype=float)
        self.End_Group = False
        self.Bond_List = []
        self.Mol_Num = Mol_Num
        try:
            for Atom_Obj in self.Atom_List:
                Mass_Weighted_Sum  += Polymer.Mass[Atom_Obj.Type]*Atom_Obj.unwrapped_position
                self.CG_Mass += Polymer.Mass[Atom_Obj.Type]
            self.COM = Mass_Weighted_Sum/self.CG_Mass
            #print self.CG_ID, self.CG_Type, self.COM, self.CG_Mass
        except:
            self.CG_Mass = Polymer.Mass[Atom_List.Type]
            self.COM = Atom_List.unwrapped_position
            #print "Exception Handling Fuck Yah"
            #print self.CG_ID, self.CG_Type, self.COM, self.CG_Mass
    
        return


    def Assign_Basis( self, A, B):
        Vec1 = self.Atom_List[A[0]].unwrapped_position - self.Atom_List[A[1]].unwrapped_position
        Vec2 = self.Atom_List[B[0]].unwrapped_position - self.Atom_List[B[1]].unwrapped_position
        Vec1 = Vec1/ np.linalg.norm(Vec1)
        Vec2 = Vec2/ np.linalg.norm(Vec2)
        Vec3 = np.zeros(3,dtype=float)
        Vec3[0] = Vec1[1]*Vec2[2] - Vec1[2]*Vec2[1]
        Vec3[1] = Vec1[2]*Vec2[0] - Vec1[0]*Vec2[2]
        Vec3[2] = Vec1[0]*Vec2[1] - Vec1[1]*Vec2[0]
        print(np.dot(Vec1,Vec1), np.dot(Vec1,Vec2), np.dot(Vec1,Vec3))
        self.basis[0] = Vec1
        self.basis[1] = Vec2
        self.basis[2] = Vec3
        return
    
    def Check_Elements(self,list1,list2):
        count = 0
        for element1 in list1:
            for element2 in list2:
                if element1 == element2:
                    count += 1
        if count == len(list1):
            return False
        else:
            return True
    
    def Find_CG_Angles(self,SS, Max_Num_Bonds = 3,Bead_List = []):
        Temp_List = copy.deepcopy(Bead_List)
        Temp_List.append(self.CG_ID)
        for Bead_ID_2 in self.Bond_List:
            if Bead_ID_2 not in Temp_List and len(Temp_List) < Max_Num_Bonds:
                SS.Find_Bead(Bead_ID_2).Find_CG_Angles(SS,Max_Num_Bonds,Temp_List)
            elif len(Temp_List) == Max_Num_Bonds:
                if Max_Num_Bonds == 3:
                    Check_List = np.zeros(len(SS.Angle_List))
                    for i in range(len(SS.Angle_List)):
                        Check_List[i] = self.Check_Elements(Temp_List,SS.Angle_List[i])
                    if all(Check_List):
                        SS.Angle_List.append(Temp_List)
                if Max_Num_Bonds == 4:
                    Check_List = np.zeros(len(SS.Dihedral_List))
                    for i in range(len(SS.Dihedral_List)):
                        Check_List[i] = self.Check_Elements(Temp_List,SS.Dihedral_List[i])
                    if all(Check_List):
                        SS.Dihedral_List.append(Temp_List)

class CG_Bond(object):
    def __init__(self, Bond_Type, Bond_Main, Bond_Node):
        self.Bond_Type = Bond_Type # integer
        self.Bond_Main = Bond_Main # CG_Bead object
        self.Bond_Node = Bond_Node
        self.Bond_Vector = (Bond_Main.COM - Bond_Node.COM)
        self.Bond_Length = math.sqrt(self.Bond_Vector[0]**2 + self.Bond_Vector[1]**2 + self.Bond_Vector[2]**2)
        self.Bond_Vector /= self.Bond_Length
        return
    def __repr__(self):
        return "Bond Type %d with length %.2f angstroms\n" % (self.Bond_Type,self.Bond_Length)
    
    def __str__(self):
        return "Bond Type %d with length %.2f angstroms\n" % (self.Bond_Type,self.Bond_Length)

class CG_Angle(object):
    def __init__(self, Angle_Type, Angle):
        self.Angle_Type = Angle_Type # integer
        self.Angle = Angle
        return

    def __repr__(self):
        return "Angle Type %d with angle %.2f degrees\n" % (self.Angle_Type,self.Angle)

    def __str__(self):
        return "Angle Type %d with angle %.2f degrees\n" % (self.Angle_Type,self.Angle)


class CG_Mol(object):
    """
        Class defining a coarse-grained molecule
    """
    def __init__(self, Mol_ID, CG_Bead_List):
        self.Mol_ID = Mol_ID
        self.CG_Bead_List = CG_Bead_List
        self.RG = 0.0
        self.COM = np.zeros(3, dtype= float)
        self.Orientation_Vector_List = []
        self.Q = np.zeros([3,3], dtype=float)
        self.Eig = np.zeros(1,dtype=float)
        return

class CG_Monomer(object):
    """
        Class defining a coarse-grained monomer
    """
    def __init__(self, Mono_ID, CG_Bead_List):
        self.Mono_ID = Mono_ID
        self.CG_Bead_List = CG_Bead_List
        self.COM = np.zeros(3, dtype = float)
        return



class Molecule(object):
    """
    class defining a molecule
    instance variable: Mol_ID, N, MW, COM, RG
    Note: COM and RG are computed as mass averages
    """
        
    def __init__(self, Mol_ID):
        self.Mol_ID = Mol_ID
        self.N = 0
        self.MW = 0.0
        self.COM = np.empty(3, dtype = float)
        self.RG = 0.0
        self.RG_x = 0.0
        self.RG_y = 0.0
        self.E2E = 0.0
        self.Atom_List = []
        self.Orientation_Vector = np.zeros([148,3], dtype =float)
        return

    def Add_Atom(self, Atom):
        self.N += 1
        self.Atom_List.append(Atom)
        self.MW += Polymer.Mass[Atom.Type]
        return

    def NP_Convert(self):
        """
            This method converts the list of atoms attribute into a numpy array to improve efficiency
            
            It is originally a python list because the length of the array isn't known a priori
            
            Note: I think the numpy array stores the elements in a closer vicinity for more convenient access.
            Not really sure about this though.
        """
        self.Atom_List = np.asarray(self.Atom_List)
        return





class Snap_Shot(object):
    """
    Class defining a snap shot of the trajectory 
    instance variables: Time_step, Mol_List, Box_Dim, Rg_Dist, E2E_Dist, RDF_Dist
    """
    def __init__(self, Time_Step, Mol_List, Box_Dim):
        self.Time_Step = Time_Step
        self.Mol_List = Mol_List
        self.Box_Dim = Box_Dim
        self.RG_Dist = np.zeros(len(Mol_List), dtype = float)
        self.RG_x_Dist = np.zeros(len(Mol_List), dtype = float)
        self.RG_y_Dist = np.zeros(len(Mol_List), dtype = float)
        self.E2E_Dist = np.zeros(len(Mol_List), dtype = float)
        self.C_Ratio = np.zeros(len(Mol_List), dtype = float)
        self.Order_Hist = []
        self.Strain = 0.0
        self.Alignment_x = 0.0
        self.Alignment_y = 0.0
        self.RDF_Dist = []
        self.Tangent_Correlation = np.zeros(24, dtype=float)
        self.Binormal_Correlation = np.zeros(24, dtype=float)
        self.Contact_Map = np.zeros([24,24], dtype=float)
        self.CG_List = []
        self.CG_Mol_List = []
        self.CG_Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Dih_Values = []
        self.Angle_Values = []
        self.Angle_Map = {}
        self.Bond_Groups = {}
        self.Dih_Map = {}
        return

    def Compute_COM(self):
        """
            Compute the centers of mass of all the polymer chains
        """
        for Mol in self.Mol_List:
            Mass_Weighted_Sum = np.zeros(3,dtype=float)
            for Atom in Mol.Atom_List:
                Atom.unwrapped_position = Atom.position + np.multiply(self.Box_Dim, Atom.image_flags)
                Mass_Weighted_Sum += Polymer.Mass[Atom.Type]*Atom.unwrapped_position
            Mol.COM = Mass_Weighted_Sum/Mol.MW
        return

    def Compute_RG(self):
        """ 
            Compute the Radius of gyration for all the polymer chains
        """
        
        for Mol in self.Mol_List:
            MSD = 0.0
            MSD_x = 0.0
            MSD_y = 0.0
            for atom in Mol.Atom_List:
                MSD_x += Polymer.Mass[atom.Type]*((atom.unwrapped_position[0] - Mol.COM[0])**2)
                MSD_y += Polymer.Mass[atom.Type]*((atom.unwrapped_position[1] - Mol.COM[1])**2)
                MSD += Polymer.Mass[atom.Type]*((atom.unwrapped_position[2] - Mol.COM[2])**2 + (atom.unwrapped_position[0] - Mol.COM[0])**2 + (atom.unwrapped_position[1] - Mol.COM[1])**2)
            Mol.RG = math.sqrt(MSD/Mol.MW)
            Mol.RG_x = math.sqrt(MSD_x/Mol.MW)
            Mol.RG_y = math.sqrt(MSD_y/Mol.MW)
            self.RG_Dist[Mol.Mol_ID-1] = Mol.RG
            self.RG_x_Dist[Mol.Mol_ID-1] = Mol.RG_x
            self.RG_y_Dist[Mol.Mol_ID-1] = Mol.RG_y
        return

    def Compute_E2E(self):
        """
            Compute the End to end vector for all the polymer chains
        """
        for Mol in self.Mol_List:
            Mol.E2E = np.linalg.norm(Mol.Atom_List[-1].unwrapped_position - Mol.Atom_List[-2].unwrapped_position)
            self.E2E_Dist[Mol.Mol_ID-1] = Mol.E2E
            self.C_Ratio[Mol.Mol_ID-1] = (self.E2E_Dist[Mol.Mol_ID-1]/self.RG_Dist[Mol.Mol_ID - 1])**2
        #print Mol.Mol_ID, self.RG_Dist[Mol.Mol_ID - 1], self.E2E_Dist[Mol.Mol_ID-1], (self.E2E_Dist[Mol.Mol_ID-1]/self.RG_Dist[Mol.Mol_ID - 1])**2
        return


    
    def Compute_Tangent_Correlation(self):
        """
            Compute the Tangent-Tangent correlation function decay to estimate the persistence length
        """
        # Sum over all chains
        for Mol in self.CG_Mol_List:
            for i in range(len(Mol.Orientation_Vector_List)):
                self.Tangent_Correlation[i] += np.dot(Mol.Orientation_Vector_List[0], Mol.Orientation_Vector_List[i])
        # Normalize
        self.Tangent_Correlation /= float(len(self.Mol_List))
       
        """
        plt.plot(self.Tangent_Correlation, linewidth=3)
        plt.ylim((-.2, 1.0))
        plt.xlim ((0, 20))
        plt.axhline(y=-.2,linewidth=4, color='k');
        plt.axhline(y=1.0,linewidth=4, color='k');
        plt.axvline(linewidth=4, color='k');
        plt.axvline( x=20 ,linewidth=4, color='k');
        plt.tick_params( labelsize = 20, width=2, length=7)
    
        plt.ylabel("Tangent Correlation Function", fontsize=25)
        plt.xlabel("Distance along chain (N)", fontsize=25)
        plt.show()
        """
        return
    
                
        # Re-wrap coordinates corresponding to periodic boundary conditions
        for CG_Mol_Obj in self.CG_Mol_List:
            for CG_Bead_Obj in CG_Mol_Obj.CG_Bead_List:
                for i in range(3):
                    if CG_Bead_Obj.COM[i] > self.Box_Dim[i]:
                        CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i] - self.Box_Dim[i]
                        CG_Bead_Obj.IF[i] += 1
                    elif CG_Bead_Obj.COM[i] < 0.0:
                        CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i] + self.Box_Dim[i]
                        CG_Bead_Obj.IF[i] -= 1
                    else:
                        CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i]
                        CG_Bead_Obj.IF[i] = 0

            
        return

    def Print_CG(self, Filename):
        

        
        File_Obj = open(Filename, 'a')
        File_Obj.write('ITEM: TIMESTEP\n')
        File_Obj.write('%d\n' % self.Time_Step)
        File_Obj.write('ITEM: NUMBER OF ATOMS\n')
        File_Obj.write('%d\n' % Polymer.Num_CG)
        File_Obj.write('ITEM: BOX BOUNDS pp pp pp\n')
        File_Obj.write('%.3f %.3f\n' % (0.000, self.Box_Dim[0]))
        File_Obj.write('%.3f %.3f\n' % (0.000, self.Box_Dim[1]))
        File_Obj.write('%.3f %.3f\n' % (0.000, self.Box_Dim[2]))
        File_Obj.write( 'ITEM: ATOMS id type mol x y z\n')
        i=1
        j=1
        for Mol in self.CG_Mol_List:
            for CG_Bead in Mol.CG_Bead_List:
                File_Obj.write('%d %d %d %.3f %.3f %.3f\n' % (i, CG_Bead.CG_Num, j, CG_Bead.COM_PBC[0], CG_Bead.COM_PBC[1],CG_Bead.COM_PBC[2]))
                i += 1
            j += 1

        return
        
    def Print_Unwrapped(self, Filename):
        File_Obj = open(Filename, 'a')
        File_Obj.write('ITEM: TIMESTEP\n')
        File_Obj.write('%d\n' % self.Time_Step)
        File_Obj.write('ITEM: NUMBER OF ATOMS\n')
        Num_Atoms = len(self.Mol_List)*len(self.Mol_List[0].Atom_List)
        File_Obj.write('%d\n' % Num_Atoms)
        File_Obj.write('ITEM: BOX BOUNDS pp pp pp\n')
        File_Obj.write('%.3f %.3f\n' % (0.000, self.Box_Dim[0]))
        File_Obj.write('%.3f %.3f\n' % (0.000, self.Box_Dim[1]))
        File_Obj.write('%.3f %.3f\n' % (0.000, self.Box_Dim[2]))
        File_Obj.write( 'ITEM: ATOMS id type mol x y z\n')
        i=1
        j=1
        for Mol in self.Mol_List:
            for Atom in Mol.Atom_List:
                File_Obj.write('%d %d %d %.3f %.3f %.3f\n' % (Atom.id , Atom.Type , j, Atom.unwrapped_position[0], Atom.unwrapped_position[1], Atom.unwrapped_position[2]))
                i += 1
            j += 1
        
        return
    
    def Find_Atom_Types(self):
        Atom_Types = []
        for Bead in self.CG_List:
            if Bead.CG_Type not in Atom_Types:
                Atom_Types.append(Bead.CG_Type)
        return Atom_Types

## Appears to be unused
    # def Define_CG_Bonds(self, Bond_Type, Type1, Type2):
    #     for Mol in self.CG_Mol_List:
    #         i= 0
    #         Temp_Bond_List = []
    #         for CG_Bead in Mol.CG_Bead_List:
    #             if CG_Bead.CG_Type == "D" and i == 0:
    #                 Bond_Main = CG_Bead
    #             if CG_Bead.CG_Type == "D" and i != 0:
    #                 Bond_Main = CG_Bead
    #                 Temp_Bond = CG_Bond("AD", Bond_Node, Bond_Main)
    #                 Temp_Bond_List.append(Temp_Bond)
    #                 Bond_List_AD.append(Temp_Bond.Bond_Length)
    #             if CG_Bead.CG_Type == "A":
    #                 Bond_Node = CG_Bead
    #                 Temp_Bond = CG_Bond("DA", Bond_Main, Bond_Node)
    #                 Temp_Bond_List.append(Temp_Bond)
    #                 Bond_List_DA.append(Temp_Bond.Bond_Length)
    #                 i += 1
    #         self.CG_Bond_List.append(Temp_Bond_List)
        
    #     """
    #     plt.hist([Bond_List_DA, Bond_List_AD], bins = 40, normed=True, histtype = 'bar', label = ["DA Bond", "AD Bond"])
    #     plt.legend( loc = 'upper right', frameon = False, fontsize= 25)
    #     plt.tick_params( labelsize = 20, width=2, length=7)
    #     plt.xlabel("Distance ($\AA$)", fontsize=25)
    #     plt.ylabel("Probability", fontsize=25)
    #     plt.xlim((6.5,8))
    #     plt.show()
    #     """


    #     """
    #     plt.hist([Angle_List_DAD, Angle_List_ADA], bins = 40, normed=True, histtype = 'bar', label = ["DAD Angle", "ADA Angle"])
    #     plt.legend( loc = 'upper right', frameon = False, fontsize= 25)
    #     plt.tick_params( labelsize = 20, width=2, length=7)
    #     plt.xlabel("Distance ($\AA$)", fontsize=25)
    #     plt.ylabel("Probability", fontsize=25)
    #     plt.xlim((0,120))
    #     plt.show()
    #     """
            
    def Group_CG_Bonds(self):
        Grouped_Bonds = []
        for i in range(0,max(self.Bond_Groups.values())):
            temp_list = np.zeros(1000)
            for bond in self.CG_Bond_List:
                if bond.Bond_Type == (i+1):
                    temp_list[int(bond.Bond_Length*100)]+=1
            temp_list = temp_list/sum(temp_list)
            Grouped_Bonds.append(temp_list)
        return Grouped_Bonds

    def Compute_CG_Angles(self):
        angle_plots = []
        for Angle in self.Angle_List:
            Bead1 = self.Find_Bead(Angle[0])
            Bead2 = self.Find_Bead(Angle[1])
            Bead3 = self.Find_Bead(Angle[2])
            Bond1 = (Bead1.COM_PBC - Bead2.COM_PBC) / math.sqrt((Bead1.COM_PBC[0]-Bead2.COM_PBC[0])**2+(Bead1.COM_PBC[1]-Bead2.COM_PBC[1])**2+(Bead1.COM_PBC[2]-Bead2.COM_PBC[2])**2)
            Bond2 = (Bead2.COM_PBC - Bead3.COM_PBC) / math.sqrt((Bead2.COM_PBC[0]-Bead3.COM_PBC[0])**2+(Bead2.COM_PBC[1]-Bead3.COM_PBC[1])**2+(Bead2.COM_PBC[2]-Bead3.COM_PBC[2])**2)
            Cos = np.dot(Bond1, Bond2)
            angle = np.arccos(Cos)
            angle *= (180./3.14)
            try:
                self.Angle_Values.append(CG_Angle((self.Angle_Map[Bead1.CG_Type + Bead2.CG_Type + Bead3.CG_Type]),angle))
            except:
                if len(self.Angle_Map) == 0:
                    self.Angle_Map[Bead1.CG_Type + Bead2.CG_Type + Bead3.CG_Type] = 1
                    self.Angle_Map[Bead3.CG_Type + Bead2.CG_Type + Bead1.CG_Type] = 1
                else:
                    self.Angle_Map[Bead1.CG_Type + Bead2.CG_Type + Bead3.CG_Type] = max(self.Angle_Map.values())+1
                    self.Angle_Map[Bead3.CG_Type + Bead2.CG_Type + Bead1.CG_Type] = max(self.Angle_Map.values())
                self.Angle_Values.append(CG_Angle((self.Angle_Map[Bead1.CG_Type + Bead2.CG_Type + Bead3.CG_Type]),angle))
        print(self.Angle_Map)
    #print self.Angle_Values
        for i in range(max(self.Angle_Map.values())):
            print(i)
            bin_list = list(range(60))
            angle_bins = np.zeros(len(bin_list))
            for Angle in self.Angle_Values:
                if Angle.Angle_Type == i+1:
                    angle_bins[int(Angle.Angle/3)] += 1
            angle_plots.append((angle_bins))
        return angle_plots

    # def Harmonic(K,theta0):
    #     P = K*(theta-theta0)**2

    def Compute_CG_Dihedrals(self):
        angle_plots = []
        for Dih in self.Dihedral_List:
            vector1 = self.Find_Bead(Dih[0]).COM_PBC - self.Find_Bead(Dih[1]).COM_PBC
            vector2 = self.Find_Bead(Dih[1]).COM_PBC - self.Find_Bead(Dih[2]).COM_PBC
            plane1 = np.cross(vector1,vector2)
            print(vector1)
            print(vector2)
            plane1 = plane1/math.sqrt(plane1[0]**2+plane1[1]**2+plane1[2]**2)
            print("Plane1")
            print(plane1)
            vector1 = self.Find_Bead(Dih[1]).COM_PBC - self.Find_Bead(Dih[2]).COM_PBC
            vector2 = self.Find_Bead(Dih[2]).COM_PBC - self.Find_Bead(Dih[3]).COM_PBC
            print(vector1)
            print(vector2)
            plane2 = np.cross(vector1,vector2)
            plane2 = plane2/math.sqrt(plane2[0]**2+plane2[1]**2+plane2[2]**2)
            Cos = np.dot(plane1, plane2)
            angle = np.arccos(Cos)
            angle *= (180./3.14)
            try:
                self.Dih_Values.append(CG_Angle((self.Dih_Map[self.Find_Bead(Dih[0]).CG_Type + self.Find_Bead(Dih[1]).CG_Type + self.Find_Bead(Dih[2]).CG_Type + self.Find_Bead(Dih[3]).CG_Type]),angle))
            except:
                if len(self.Dih_Map) == 0:
                    self.Dih_Map[self.Find_Bead(Dih[0]).CG_Type + self.Find_Bead(Dih[1]).CG_Type + self.Find_Bead(Dih[2]).CG_Type + self.Find_Bead(Dih[3]).CG_Type] = 1
                    self.Dih_Map[self.Find_Bead(Dih[3]).CG_Type + self.Find_Bead(Dih[2]).CG_Type + self.Find_Bead(Dih[1]).CG_Type + self.Find_Bead(Dih[0]).CG_Type] = 1
                else:
                    self.Dih_Map[self.Find_Bead(Dih[0]).CG_Type + self.Find_Bead(Dih[1]).CG_Type + self.Find_Bead(Dih[2]).CG_Type + self.Find_Bead(Dih[3]).CG_Type] = max(self.Dih_Map.values())+1
                    self.Dih_Map[self.Find_Bead(Dih[3]).CG_Type + self.Find_Bead(Dih[2]).CG_Type + self.Find_Bead(Dih[1]).CG_Type + self.Find_Bead(Dih[0]).CG_Type] = max(self.Dih_Map.values())
                self.Dih_Values.append(CG_Angle((self.Dih_Map[self.Find_Bead(Dih[0]).CG_Type + self.Find_Bead(Dih[1]).CG_Type + self.Find_Bead(Dih[2]).CG_Type + self.Find_Bead(Dih[3]).CG_Type]),angle))
        print(self.Dih_Map)
        print(self.Dih_Values)
        print(plane1)
        print(plane2)
        for i in range(max(self.Angle_Map.values())):
           print(i)
           bin_list = list(range(60))
           angle_bins = np.zeros(len(bin_list))
           for Angle in self.Dih_Values:
               if Angle.Angle_Type == i+1:
                   print(int(Angle.Angle/3))
                   angle_bins[int(Angle.Angle/3)-1] += 1
           angle_plots.append((angle_bins)/sum(angle_bins))
        return angle_plots



    def Compute_CG_Dihedrals_2(self):
        Dihedral_List_AD = []
        Dihedral_List_DA = []
        Conjugation_Dist = [1]
        i = 0
        for CG_Mol in self.CG_Mol_List:
            for CG_Bead in CG_Mol.CG_Bead_List:
                if CG_Bead.CG_Type == "D":
                    Basis_D = CG_Bead.basis
                    try:
                        N1 = Basis_A[1]
                        N2 = Basis_D[1]
                        M1 = Basis_A[2]
                        x = np.dot(N1,N2)
                        y = np.dot(M1,N2)
                        angle = np.arctan2(y,x)*(180./3.14)
                        Dihedral_List_AD.append(angle)
                    except:
                        i+=1
                        print("exception %d" % i)
                        continue
                                
                            
                elif CG_Bead.CG_Type == "A":
                    Basis_A = CG_Bead.basis
                    try:
                        N1 = Basis_D[1]
                        N2 = Basis_A[1]
                        M1 = Basis_D[2]
                        x = np.dot(N1,N2)
                        y = np.dot(M1,N2)
                        angle = np.arctan2(y,x)*(180./3.14)
                        Dihedral_List_DA.append(angle)
                    except:
                        i+=1
                        print("exception %d" % i)
                        continue
                
                if CG_Bead.CG_Type == "A" or CG_Bead.CG_Type == "D":
                    try:
                        if abs(angle) < 40.0 or abs(angle) > 140.0:
                            print("Conjugated %.2f" % abs(angle))
                            Conjugation_Dist[-1] += 1
                        else:
                            print("Not Conjugated %.2f" % abs(angle))
                            Conjugation_Dist.append(1)
                    except:
                            continue

            Conjugation_Dist.append(1)
                

    
        return Dihedral_List_DA, Dihedral_List_AD, Conjugation_Dist
    
    def Export2Z(self, Filename):
        File_Obj = open(Filename, 'w')
        Num_Chains = 60
        Chain_Length = 24
        File_Obj.write('%d\n' % Num_Chains)
        File_Obj.write('%.2f %.2f %.2f\n' % (self.Box_Dim[0], self.Box_Dim[1], self.Box_Dim[2]))
        for i in range(Num_Chains):
            File_Obj.write('%d ' % Chain_Length)
        File_Obj.write('\n')
        for Mol in self.CG_Mol_List:
            j = 0
            for bead in Mol.CG_Bead_List:
                if bead.CG_Type == "D" or bead.CG_Type == "A":
                    j += 1
                    File_Obj.write('%.2f %.2f %.2f\n' % (bead.COM[0], bead.COM[1], bead.COM[2]))
        File_Obj.close()
        return


    def Export2Z_SC(self, Filename):
        File_Obj = open(Filename, 'w')
        Num_Chains = 60
        Chain_Length = 24
        File_Obj.write('%d\n' % Num_Chains)
        File_Obj.write('%.2f %.2f %.2f\n' % (self.Box_Dim[0], self.Box_Dim[1], self.Box_Dim[2]))
        for i in range(Num_Chains):
            File_Obj.write('%d ' % Chain_Length)
            File_Obj.write('\n')
        for Mol in self.CG_Mol_List:
            j = 0
            for bead in Mol.CG_Bead_List:
                if bead.CG_Type == "D" or bead.CG_Type == "A":
                    j += 1
                    File_Obj.write('%.2f %.2f %.2f\n' % (bead.COM[0], bead.COM[1], bead.COM[2]))
        File_Obj.close()
        return


    
    def Compute_Tangent_Correlation2(self):
        """
            Compute the Tangent-Tangent correlation function decay to estimate the persistence length
            """
        # Sum over all chains
        i = 0
        j = 0
        Vector_List = np.zeros([3, 24, 60], dtype=float)
        Vector_List2 = np.zeros([3, 24, 60], dtype=float)
        for Mol in self.CG_Mol_List:
            i = 0
            for CG_Bead in Mol.CG_Bead_List:
                if CG_Bead.CG_Type == "D" or CG_Bead.CG_Type == "A":
                    Vector_List[:,i,j] = CG_Bead.basis[0]
                    Vector_List2[:,i,j] = CG_Bead.basis[1]
                    i += 1
            j += 1
        
        Tangent_Correlation = np.zeros(24)
        Binormal_Correlation = np.zeros(24)
        for i in range(60):
            for j in range(24):
                for k in range(24-j):
                    Tangent_Correlation[k] += np.dot(Vector_List[:,j,i], Vector_List[:,k+j,i])
                    Binormal_Correlation[k] += np.dot(Vector_List2[:,j,i], Vector_List2[:,k+j,i])
        
        for i in range(24):
            Tangent_Correlation[i] /= (60.0*24.-float(i))
            Binormal_Correlation[i] /= (60.0*24.-float(i))
            
        self.Tangent_Correlation = Tangent_Correlation
        self.Binormal_Correlation = Binormal_Correlation
        print(Binormal_Correlation)
        
        """
        plt.plot(Tangent_Correlation, linewidth=3)
        plt.ylim((-.2, 1.0))
        plt.xlim ((0, 23))
        plt.axhline(y=-.2,linewidth=4, color='k');
        plt.axhline(y=1.0,linewidth=4, color='k');
        plt.axvline(linewidth=4, color='k');
        plt.axvline( x=23 ,linewidth=4, color='k');
        plt.tick_params( labelsize = 20, width=2, length=7)
            
        plt.ylabel("Tangent Correlation Function", fontsize=25)
        plt.xlabel("Distance along chain (N)", fontsize=25)
        plt.show()
        """
        
        return
    
    def Compute_RDF(self, type1, type2):
        RDF = []
        # Loop over all molecules
        dr = .3
        Rmax = 94.0
        MaxBin = int(Rmax/dr)
        Hist = np.zeros(MaxBin, dtype= int)
        R = np.arange(dr, Rmax, dr)
        NA = 0
        NB = 0
        for CG_Bead_Obj in self.CG_List:
            for i in range(3):
                if CG_Bead_Obj.COM[i] > self.Box_Dim[i]:
                    CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i] - self.Box_Dim[i]
                    CG_Bead_Obj.IF[i] += 1
                elif CG_Bead_Obj.COM[i] < 0.0:
                    CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i] + self.Box_Dim[i]
                    CG_Bead_Obj.IF[i] -= 1
                else:
                    CG_Bead_Obj.COM_PBC[i] = CG_Bead_Obj.COM[i]
                    CG_Bead_Obj.IF[i] = 0

        for CG_Bead1 in self.CG_List:
            if CG_Bead1.CG_Type == type1 and CG_Bead1.End_Group == False:
                NA += 1
                NB = 0
                for CG_Bead2 in self.CG_List:
                    if CG_Bead2.CG_Type == type2 and CG_Bead1.CG_ID != CG_Bead2.CG_ID:
                        NB += 1
                        if CG_Bead2.End_Group == False and CG_Bead2.CG_ID not in CG_Bead1.Bond_List:
                            R12 = CG_Bead1.COM_PBC - CG_Bead2.COM_PBC
                            # Apply minimum image convention
                            for i in range(3):
                                if R12[i] > (0.5*self.Box_Dim[i]):
                                    R12[i] -= self.Box_Dim[i]
                                if R12[i] <= (-0.5*self.Box_Dim[i]):
                                    R12[i] += self.Box_Dim[i]
                            R12SQ = R12[0]**2 + R12[1]**2 + R12[2]**2
                            Distance = math.sqrt(R12SQ)
                            if Distance < 3.5:
                                print(CG_Bead1.COM_PBC)
                                print(CG_Bead2.COM_PBC)
                                print(CG_Bead1.COM)
                                print(CG_Bead2.COM)
                                print("Distance smaller than 3.5", Distance, CG_Bead1.CG_ID, CG_Bead2.CG_ID)
                            RDF.append(Distance)
                            Bin = int(Distance/dr)
                            # Compute a histogram of pair seperations
                            if Bin <  MaxBin:
                                Hist[Bin] += 1
        print("Max Bin = ", MaxBin)

        rhoB = NB/(self.Box_Dim[0]**3)
        rdf = np.zeros(MaxBin, dtype=float)
        for i in range(len(Hist)):
            if float((NA*rhoB*4*3.14*((i+1)**3 - i**3)*dr**3)) != 0:
                rdf[i] = 3*float(Hist[i])/float((NA*rhoB*4*3.14*((i+1)**3 - i**3)*(dr**3)))
        
        return rdf


    def Compute_Order(self):
        """
        Compute the orientation vector for each monomer, as well as chain alignment order parameter for the system
        """
        Num_Beads = Polymer.Num_Beads
        
        print("Computing order parameter")
        for Mol in self.CG_Mol_List:
            j=0
            Donor_Index = Polymer.Donor_Index
            Acceptor_Index = Polymer.Acceptor_Index
            
            for CG_Bead in Mol.CG_Bead_List:
                j += 1
                if j > Num_Beads and j < Num_Beads*11 and CG_Bead.CG_Type == "D":
                    #print j, Mol.CG_Bead_List[ j - Num_Beads + Acceptor_Index -1 ].CG_Type, Mol.CG_Bead_List[j + Acceptor_Index - 1].CG_Type
                    Bond = Mol.CG_Bead_List[ j - Num_Beads + Acceptor_Index - 1].COM - Mol.CG_Bead_List[j + Acceptor_Index -1].COM
                    Bond = Bond/np.linalg.norm(Bond)
                    CG_Bead.Orientation_Vector = Bond
                    Mol.Orientation_Vector_List.append(Bond)
                if j > Num_Beads and j < Num_Beads*11 and CG_Bead.CG_Type == "A":
                    #print j, Mol.CG_Bead_List[ j - Acceptor_Index - 1 ].CG_Type, Mol.CG_Bead_List[j + Num_Beads - Acceptor_Index  -1].CG_Type
                    Bond = Mol.CG_Bead_List[ j - Acceptor_Index - 1].COM - Mol.CG_Bead_List[j + Num_Beads - Acceptor_Index -1].COM
                    Bond = Bond/np.linalg.norm(Bond)
                    CG_Bead.Orientation_Vector = Bond
                    Mol.Orientation_Vector_List.append(Bond)
        return

    def Compute_Order2(self):
        for Mol in self.CG_Mol_List:
            for CG_Bead in Mol.CG_Bead_List:
                if CG_Bead.CG_Type == "D" or CG_Bead.CG_Type == "A":
                    Mol.Orientation_Vector_List.append(CG_Bead.basis[0])
        return
    
    


    def Compute_Contact_Map(self):
        
        # Sum over all chains
        f = 0
        Eig_List = []
        print("Computing Order")
        for Mol in self.CG_Mol_List:
            f += 1
            for i in range(len(Mol.Orientation_Vector_List)):
                for j in range(len(Mol.Orientation_Vector_List)):
                    self.Contact_Map[i,j] = np.dot(Mol.Orientation_Vector_List[i], Mol.Orientation_Vector_List[j])
            
            for i in range(len(Mol.Orientation_Vector_List)):
                for j in range(3):
                    for k in range(3):
                        if j == k:
                            Mol.Q[(j,k)] += (3.0/2.0)*(Mol.Orientation_Vector_List[i][j]*Mol.Orientation_Vector_List[i][k] - 1.0/2.0)
                        if j != k:
                            Mol.Q[(j,k)] += (3.0/2.0)*Mol.Orientation_Vector_List[i][j]*Mol.Orientation_Vector_List[i][k]
            Mol.Q = Mol.Q / float(len(Mol.Orientation_Vector_List))
            print(np.sort(np.linalg.eig(Mol.Q)[0]), np.max(np.linalg.eig(Mol.Q)[0]))
            fig, ax = plt.subplots()
            heatmap = ax.pcolor(self.Contact_Map, cmap=plt.cm.hot)
            ax.set_xlim((0,24))
            ax.set_ylim((0,24))
            #fig.colorbar()
            plt.pcolormesh(self.Contact_Map, cmap=plt.cm.hot)
            plt.colorbar()
            #plt.show()
            
            Mol.Eig = np.max(np.linalg.eig(Mol.Q)[0])
            Eig_List.append(Mol.Eig)
            
        Eig_List = np.asarray(Eig_List)
        Average_Order = Eig_List.mean()
        Stdv = Eig_List.std()
        print(" Average intrachain order is ", Average_Order, "+-", Stdv)

        """
        fig, ax = plt.subplots()
        heatmap = ax.pcolor(self.Contact_Map, cmap=plt.cm.hot)
        ax.set_xlim((0,24))
        ax.set_ylim((0,24))
        fig.colorbar()
        plt.pcolormesh(self.Contact_Map, cmap=plt.cm.hot)
        plt.colorbar()
        plt.show()
        """
        # Normalize
        self.Contact_Map /= float(len(self.Mol_List))

    def Find_Bead(self,CG_ID):
        for Bead in self.CG_List:
            if Bead.CG_ID == CG_ID:
                return Bead

    def Find_Angles(self):
        for Bead in self.CG_List:
            Bead.Find_CG_Angles(self)
        return self.Angle_List

    def Find_Dihedrals(self):
        for Bead in self.CG_List:
            Bead.Find_CG_Angles(self, Max_Num_Bonds = 4)
        print(self.Dihedral_List)
        return self.Dihedral_List



# Classes for storing output data



class RDF_OBJ:
    def __init__(self, Type1, Type2):
        self.Radius = np.zeros(312,dtype=float)
        self.RDF = np.zeros(312, dtype=float)
        self.Num_Snaps = 0.0
        self.Type1 = Type1
        self.Type2 = Type2
    def Add_Snap(self, RDF, radius):
        self.RDF += RDF
        self.Radius = radius
        self.Num_Snaps += 1.0
    def Normalize(self):
        self.RDF = self.RDF/self.Num_Snaps


class Output(object):
    def __init__(self, Rg, Tangent_Correlation, Binormal_Correlation, RDF, Conjugation_Dist, Dihedral):
        self.Rg = Rg
        self.Tangent_Correlation = Tangent_Correlation
        self.Binormal_Correlation = Binormal_Correlation
        self.RDF = RDF
        self.Conjugation_Dist = Conjugation_Dist
        self.Dihedral = Dihedral
        return


class Trajectory(object):
    """
    Class defining a complete trajectory outputted from an MD simulation
    instance variables: Num_Snaps, Snap_List
    """

    def __init__(self, File_Name,Num_Poly,Bond_List = [],Angle_List = [],Dihedral_List = [],Atom_Types = {}):
            File = open(File_Name,'r')
            File_Lines = File.readlines()
            self.Temp_Snap_List = []
            self.Order_Evolution_x = []
            self.Order_Evolution_y = []
            self.RG_x_Ave = []
            self.RG_y_Ave = []
            self.Strain= []
            self.RG_Dist = []
            self.E2E_Dist = []
            self.C_Ratio_Dist = []
            self.Bond_List = Bond_List
            self.Bond_Values = []
            self.Angle_List = Angle_List
            self.Angle_Values = []
            self.Dihedral_List = Dihedral_List
            self.Dih_Values = []
            self.Tangent_Correlation = np.zeros(24,dtype=float)
            self.Binormal_Correlation = np.zeros(24,dtype=float)
            self.Conjugation_Dist = []
            self.Atom_Types = Atom_Types
            iter = 0
            # Initialize RDF_OBJ's
            self.rdf = {}
            N_MOL = Num_Poly
            # Parse through file lines and extract the trajectory data
            for i in range(len(File_Lines)):
                line = File_Lines[i]
                if line == "ITEM: TIMESTEP\n":
                    Time_Step = int(File_Lines[i+1])
                    N = int(File_Lines[i+3])
                    Atom_List = []
                    XBounds = [ float(File_Lines[i+5].split()[0]), float(File_Lines[i+5].split()[1])]
                    YBounds = [ float(File_Lines[i+6].split()[0]), float(File_Lines[i+6].split()[1])]
                    ZBounds = [ float(File_Lines[i+7].split()[0]), float(File_Lines[i+7].split()[1])]
                    Box_Dim = [ abs(XBounds[1]- XBounds[0]), abs(YBounds[1] - YBounds[0]), abs(ZBounds[1]- ZBounds[0])]
                    
                    # Define the unstrained Box Length
                    if (Box_Dim[0] == Box_Dim[1]) and (Box_Dim[1] == Box_Dim[2]):
                        self.L0 = Box_Dim[1]
                        print("Setting L0 to ", self.L0)
                    
                    
                    print("Timestep = %d, N = %d\n" % (Time_Step, N))
                    print(Box_Dim)
                    # Extract Snapshot as list of atoms
                    print("Extracting Atoms...\n")
                    for j in range(N):
                        Atom_Line = File_Lines[i + 9 + j].split(' ')
                        ID = int(Atom_Line[0])
                        TYPE = int(Atom_Line[1])
                        MOL = int(Atom_Line[2])
                        POS= [ float(Atom_Line[3]), float(Atom_Line[4]), float(Atom_Line[5])]
                        IMAGE = [ int(Atom_Line[6]), int(Atom_Line[7]), int(Atom_Line[8])]
                        if MOL < (N_MOL + 1):
                            Atom_List.append(Atom(ID, TYPE, MOL, POS, IMAGE))
                    print("Finished extracting atoms, now sorting them")
                    Sorted_Atom_List = sorted(Atom_List, key=lambda Atom: Atom.id)
                    # Extract Molecule Objects from Atom List
                    print("Instantiating Molecule Objects")
                    print("N_MOL = %d " % N_MOL)
                    Mol_List = np.empty(N_MOL, dtype=object)
                    for i in range(N_MOL):
                        Mol_List[i] = Molecule(i+1)

                    for Atom_Obj in Sorted_Atom_List:
                        #Atom_Obj.Print_Info()
                        MOLID = Atom_Obj.Mol_ID
                        Mol_List[MOLID-1].Add_Atom(Atom_Obj)
                        Atom_Obj.position -= [ XBounds[0], YBounds[0], ZBounds[0]]
                    
                    print("Converting to Numpy array")
                    for Mol_Obj in Mol_List:
                        Mol_Obj.NP_Convert()
                    print("Instantiating Snap Shot Object, Computing properties")
                    Snap_Shot_Obj = Snap_Shot(Time_Step, Mol_List, Box_Dim)
                    
                    # Compute center of mass, Radius of Gyration, E2E, order parameter, Backbone Tangent Correlation, Contact Map, etc.
                    Snap_Shot_Obj.Compute_COM()
                    Snap_Shot_Obj.Compute_RG()
                    """Snap_Shot_Obj.Compute_E2E()"""
                    #Snap_Shot_Obj.Map2CG_PTB7()
                    tic = time.clock()
                    Polymer.Map2CG(Snap_Shot_Obj)
                    toc = time.clock()
                    print("Mapping Time:")
                    print((toc-tic))
                    tic = time.clock()
                    if iter == 0:
                        self.Atom_Types = Snap_Shot_Obj.Find_Atom_Types()
                    for i in range(len(self.Atom_Types)):
                        for j in range(i,len(self.Atom_Types)):
                            if iter == 0:
                                self.rdf[self.Atom_Types[i] + self.Atom_Types[j]]=(Snap_Shot_Obj.Compute_RDF(self.Atom_Types[i],self.Atom_Types[j]))
                            else:
                                self.rdf[self.Atom_Types[i] + self.Atom_Types[j]]=(self.rdf[self.Atom_Types[i] + self.Atom_Types[j]]*iter + Snap_Shot_Obj.Compute_RDF(self.Atom_Types[i],self.Atom_Types[j]))
                                self.rdf[self.Atom_Types[i] + self.Atom_Types[j]]/=(iter+1)
                    toc = time.clock()
                    print("RDF Time:")
                    print((toc-tic))
                    tic = toc
                    if iter == 0:
                        self.Bond_List = Snap_Shot_Obj.CG_Bond_List
                        self.Bond_Values = np.asarray(Snap_Shot_Obj.Group_CG_Bonds())
                    else:
                        self.Bond_Values += np.asarray(Snap_Shot_Obj.Group_CG_Bonds())
                    toc = time.clock()
                    print("Bond Time:")
                    print((toc-tic))
                    tic = toc
                    if self.Angle_List == []:
                        print("Empty Angles")
                        self.Angle_List = Snap_Shot_Obj.Find_Angles()
                    else:
                        Snap_Shot_Obj.Angle_List = self.Angle_List
                    if iter == 0:
                        self.Angle_Values = np.asarray(Snap_Shot_Obj.Compute_CG_Angles())
                    else:
                        print(sum(self.Angle_Values[0]))
                        self.Angle_Values += np.asarray(Snap_Shot_Obj.Compute_CG_Angles())
                        print(sum(self.Angle_Values[0]))
                    toc = time.clock()
                    print("Angle Time:")
                    print((toc-tic))
                    tic = toc
                    if self.Dihedral_List == []:
                        self.Dihedral_List = Snap_Shot_Obj.Find_Dihedrals()
                    else:
                        Snap_Shot_Obj.Dihedral_List = self.Dihedral_List
                    if iter == 0:
                        self.Dih_Values = np.asarray(Snap_Shot_Obj.Compute_CG_Dihedrals())
                    else:
                        self.Dih_Values += np.asarray(Snap_Shot_Obj.Compute_CG_Dihedrals())
                    toc = time.clock()
                    print("Dihedral Time:")
                    print((toc-tic))
                    """if self.Dihedral_List == []:
                        self.Dihedral_List = np.asarray(Snap_Shot_Obj.Find_Dihedrals())
                    else:
                        self.Dihedral_List += np.asarray(Snap_Shot_Obj.Find_Dihedrals())"""
                    """DADA, ADAD, Conjugation_Dist = Snap_Shot_Obj.Compute_CG_Dihedrals_2()
                    self.Dihedral_List_DADA.extend(DADA)
                    self.Dihedral_List_ADAD.extend(ADAD)
                    self.Conjugation_Dist.extend(Conjugation_Dist)
                    Snap_Shot_Obj.Print_CG("CG_PBC_420.lammpstrj")
                    Snap_Shot_Obj.Print_Unwrapped("Atomistic_Unwrapped_600.lammpstrj")"""
                    """RDF, RADIUS = Snap_Shot_Obj.Compute_RDF("D", "D")
                    self.rdf[0].Add_Snap(RDF, RADIUS)
                    RDF, RADIUS = Snap_Shot_Obj.Compute_RDF("D", "A")
                    self.rdf[1].Add_Snap(RDF, RADIUS)
                    #RDF, RADIUS = Snap_Shot_Obj.Compute_RDF("A", "D")
                    #self.rdf[2].Add_Snap(RDF, RADIUS)
                    RDF, RADIUS = Snap_Shot_Obj.Compute_RDF("A", "A")
                    self.rdf[2].Add_Snap(RDF, RADIUS)"""
                    """Snap_Shot_Obj.Compute_Order2()
                    Snap_Shot_Obj.Compute_Tangent_Correlation2()
                    Snap_Shot_Obj.Export2Z( "DA_Z1.txt")
                    Snap_Shot_Obj.Compute_Contact_Map()"""
                    #Snap_Shot_Obj.Plot_Rg_Dist()
                    #Snap_Shot_Obj.Strain = (Snap_Shot_Obj.Box_Dim[0] - self.L0)/self.L0
                    #self.Order_Evolution_x.append(Snap_Shot_Obj.Alignment_x)
                    #self.Order_Evolution_y.append(Snap_Shot_Obj.Alignment_y)
                    #self.RG_x_Ave.append(Snap_Shot_Obj.RG_x_Dist.mean())
                    #self.RG_y_Ave.append(Snap_Shot_Obj.RG_y_Dist.mean())
                    #self.Strain.append(Snap_Shot_Obj.Strain)
                    """print "Strain = ", Snap_Shot_Obj.Strain
                    print "(E2E/RG)$^2$ = ", (Snap_Shot_Obj.E2E_Dist.mean()/Snap_Shot_Obj.RG_Dist.mean())**2"""
                    self.Temp_Snap_List.append(Snap_Shot_Obj)
                    iter += 1
            """self.Plot_Angles()
            # Set instance variables
            self.Num_Snaps = len(Temp_Snap_List)
            self.Snap_Shot_List = Temp_Snap_List
            for RDF in self.rdf:
                RDF.Normalize()
            
            # Compile and print output object
            for Snap_Obj in self.Snap_Shot_List:
                self.C_Ratio_Dist.extend(Snap_Obj.C_Ratio)
                self.Tangent_Correlation += Snap_Obj.Tangent_Correlation
                self.Binormal_Correlation += Snap_Obj.Binormal_Correlation
                self.RG_Dist.extend(Snap_Obj.RG_Dist)
            
            self.Tangent_Correlation = self.Tangent_Correlation/ float(self.Num_Snaps)
            self.Binormal_Correlation = self.Binormal_Correlation/ float(self.Num_Snaps)
            
            output = Output(self.RG_Dist, self.Tangent_Correlation, self.Binormal_Correlation, self.rdf, self.Conjugation_Dist, [self.Dihedral_List_ADAD, self.Dihedral_List_DADA])
            Pickle_File = "Output_420.pickle"
            FileObj1 = open(Pickle_File, 'wb')
            pickle.dump(output, FileObj1)
            FileObj1.close()"""
            return

    def Plot_Angles(self):
        for plot in self.Angle_Values:
            plt.figure()
            plt.plot(list(range(len(plot))),plot)
        plt.show()


def main():
    Traj_Object1 = Trajectory('Traj_420.lammpstrj')
    #Traj_Object1.Plot_Bond()
    #Traj_Object1.Plot_Angle()
    #Traj_Object1.Plot_Dihedral()
    #Traj_Object1.Plot_RG_Dist()
    #Traj_Object1.Plot_RDF()
    #Traj_Object2 = Trajectory('NPT_600_2.lammpstrj')


if __name__=='__main__': main()












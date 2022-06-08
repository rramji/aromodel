#! usr/bin/python
import Molecule
import copy
import numpy as np
import scipy
import math
import OPLS
import Bond
import Angle
import Dihedral
import Improper
import Write_Inputs
import os
import Aux_Ring
import time
import Atom
import matplotlib.pyplot as plt
import random

class Conjugated_Polymer(Molecule.Molecule):
    def __init__(self,Ring_List,Scheduler = "Torque",Cluster_Location='/oasis/tscc/scratch/andrewk/Optimized_Monomers',Flip=10000000):
        self.Ring_List = copy.deepcopy(Ring_List)
        for ring in self.Ring_List:
            ring.Relink()
        self.Aux_Ring_List = {}
        self.Aux_Ring_OOP = {}
        self.Aux_Ring_Dih = {}
        Num_Aux_Rings = 0
        for ring in self.Ring_List:
            if ring.Aux_Ring_List != []:
                self.Aux_Ring_List[ring] = ring.Aux_Ring_List
                self.Aux_Ring_OOP[ring] = np.zeros(len(self.Aux_Ring_List[ring]))
                self.Aux_Ring_Dih[ring] = np.zeros(len(self.Aux_Ring_List[ring]))
                Num_Aux_Rings += len(ring.Aux_Ring_List)
        self.Num_Rings = len(self.Ring_List) + Num_Aux_Rings
        self.Name = ""
        self.OOP_Rotation = np.zeros(len(Ring_List)-1)
        self.Dih_Rotation = np.zeros(len(Ring_List)-1)

        self.Atom_List = []
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.Interring_Bond_List = []
        self.UnConverged = True
        self.Short_Name = ""
        for ring in self.Ring_List[:-1]:
            if ring.Aux_Ring_List == []:
                self.Name = self.Name + ring.Name + "_Phi_0_Theta_0_"
            else:
                self.Name = self.Name + ring.Name
                for i,aux_ring in enumerate(self.Aux_Ring_List[ring]):
                    self.Name = self.Name + "_Phi_0_Theta_0_Aux_%s_%d_" % (aux_ring.Name,i+1)
                self.Name = self.Name + "_Phi_0_Theta_0"
        self.Name = self.Name + self.Ring_List[-1].Name
        if self.Ring_List[-1].Aux_Ring_List != []:
            for i,aux_ring in enumerate(self.Aux_Ring_List[self.Ring_List[-1]]):
                self.Name = self.Name + "_Phi_0_Theta_0_Aux_%s_%d_" % (aux_ring.Name,i+1)

        for ring in self.Ring_List[:-1]:
            if ring.Aux_Ring_List == []:
                self.Short_Name = self.Short_Name + ring.Nickname + "_Phi_0_Theta_0_"
            else:
                self.Short_Name = self.Short_Name + ring.Nickname
                for i,aux_ring in enumerate(self.Aux_Ring_List[ring]):
                    self.Nickname = self.Nickname + "_Phi_0_Theta_0_Aux_%s_%d_" % (aux_ring.Nickname,i+1)
                self.Nickname = self.Nickname + "_Phi_0_Theta_0"
        self.Short_Name = self.Short_Name+ self.Ring_List[-1].Nickname
        if self.Ring_List[-1] in self.Aux_Ring_List.keys():
            for i,aux_ring in enumerate(self.Aux_Ring_List[self.Ring_List[-1]]):
                self.Short_Name = self.Short_Name + "_Phi_0_Theta_0_Aux_%s_%d_" % (aux_ring.Nickname,i+1)

        for i,ring in enumerate(self.Ring_List):
            ring.Ring_ID=i+1
        info = []
        for i in range(len(self.Ring_List)-1):
            if not self.Ring_List[i].Bonded_Atoms[-1].Is_Linked:
                center_position,bond_vector,new_ring_bonded_atom = self.Ring_List[i].Link_Rings(self.Ring_List[i+1])
                info.append((center_position,bond_vector,new_ring_bonded_atom))
                self.Ring_List[i].Link_Bonded_Atoms()
                self.Ring_List[i+1].Link_Bonded_Atoms()
        """for ring in self.Ring_List:
            center_position = info[i][0]
            bond_vector = info[i][1]
            new_ring_bonded_atom = info[i][2]
            ring.Update_Normal_Vector()
            New_X = self.Ring_List[i].Normal_Vector
            New_Y = bond_vector/np.linalg.norm(bond_vector)
            New_X = (New_X - np.dot(New_X,New_Y)*New_Y)/np.linalg.norm(New_X - np.dot(New_X,New_Y)*New_Y)
            New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
            Rotation_Matrix = [New_X,New_Y,New_Z]
            ring.Rotate_Ring(Rotation_Matrix)"""

        #Added 11/21
        """for i in range(len(self.Ring_List) - 1):
            if (i + 1) % Flip == 0:
                print("Flipped")
                self.Rotate_Ring("Dih", 180, self.Ring_List[i + 1], self.Ring_List[i])"""
        #End Add

        #for i in range(len(self.Ring_List)-1):
        #    for b_atom in self.Ring_List[i + 1].Bonded_Atoms:
        FlipIt = False
        for i in range(len(self.Ring_List)-1):
            """for j in range(i+1):
                    self.Ring_List[j].Translate_Ring(center_position)"""
            for b_atom in self.Ring_List[i+1].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Bonded_Ring == self.Ring_List[i] and b_atom.Interring_Bond_Atom.Bonded_Ring == self.Ring_List[i+1]:
                    #center = copy.deepcopy(b_atom.Central_Atom.Position)
                    #center_atom = b_atom
                    if (i+1) % Flip == 0:
                        if FlipIt:
                            FlipIt = False
                        else:
                            FlipIt = True
                    #self.Ring_List[i+1].Translate_Ring(center)
                    #self.Ring_List[i+1].Translate_Ring(-1*center_atom.Interring_Bond_Atom.Central_Atom.Position)
                    #self.Ring_List[i+1].Translate_Ring(-1*center_atom.Interring_Bond_Atom.Bonded_Vector)
                    self.Ring_List[i+1].Translate_Ring(b_atom.Central_Atom.Position)
                    self.Ring_List[i+1].Translate_Ring(-1*b_atom.Interring_Bond_Atom.Central_Atom.Position)
                    #self.Ring_List[i+1].Translate_Ring(-1*b_atom.Interring_Bond_Atom.Bonded_Vector)
                    self.Ring_List[i + 1].Translate_Ring(b_atom.Bonded_Vector)

                    New_X = self.Ring_List[i+1].Normal_Vector
                    New_Y = b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector)
                    New_X = (New_X - np.dot(New_X,New_Y)*New_Y)/np.linalg.norm(New_X - np.dot(New_X,New_Y)*New_Y)
                    New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
                    Rotation_Matrix = [New_X,New_Y,New_Z]
                    #self.Rotate_Polymer(Rotation_Matrix) #Commented 11/22
                    center = copy.deepcopy(b_atom.Central_Atom.Position)
                    self.Translate_Polymer(center)
                    #rotation_angle = -1*math.acos(np.dot(b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector),-1*b_atom.Interring_Bond_Atom.Bonded_Vector/np.linalg.norm(b_atom.Interring_Bond_Atom.Bonded_Vector)))
                    #possibly comment 11/20
                    #rotation_angle = math.acos(np.dot(b_atom.Original_Bonded_Vector / np.linalg.norm(b_atom.Original_Bonded_Vector),-1 * b_atom.Interring_Bond_Atom.Original_Bonded_Vector / np.linalg.norm(b_atom.Interring_Bond_Atom.Original_Bonded_Vector)))
                    #for j in range(i+1,len(self.Ring_List)):
                    #    self.Ring_List[j].Rotate_Ring([[1,0,0],[0,math.cos(rotation_angle),-1*math.sin(rotation_angle)],[0,math.sin(rotation_angle),math.cos(rotation_angle)]])
                    #end possible comment
            """deviation_angle = new_ring_bonded_atom.Check_Alignment()
            self.Ring_List[i+1].Rotate_Ring_Not_Interring([[1,0,0],[0,math.cos(deviation_angle),-math.sin(deviation_angle)],[0,math.sin(deviation_angle),math.cos(deviation_angle)]])"""
            #self.Ring_List[i+1].Translate_Ring(new_ring_bonded_atom.Bonded_Vector)
        #end delete
        self.N = 0
        for ring in self.Ring_List:
            self.N += len(ring.Atom_List)
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    self.N += 1
        atom_id = 1
        for ring in self.Ring_List:
            for atom in ring.Atom_List:
                atom.Atom_ID = atom_id
                atom_id += 1
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    b_atom.H_Atom.Atom_ID = atom_id
                    atom_id += 1
        #self.Write_XYZ()
        for ring in self.Ring_List:
            ring.Add_Bond_List()
        for ring in self.Ring_List:
            for atom in ring.Atom_List:
                if atom.Element == "C":
                    #atom.Assign_OPLS_ID()
                    self.Atom_List.append(atom)
            for atom in ring.Atom_List:
                if atom.Element != "C":
                    #atom.Assign_OPLS_ID()
                    self.Atom_List.append(atom)
            for b_atom in ring.Bonded_Atoms:
                if not b_atom.Is_Linked:
                    #b_atom.H_Atom.Assign_OPLS_ID()
                    self.Atom_List.append(b_atom.H_Atom)
        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
        for i in range(len(self.Ring_List) - 1):
            for b_atom in self.Ring_List[i].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Bonded_Ring == self.Ring_List[i+1]:
                    bond_master = b_atom.Central_Atom
            for b_atom in self.Ring_List[i+1].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Bonded_Ring == self.Ring_List[i]:
                    bond_slave = b_atom.Central_Atom
            self.Interring_Bond_List.append(Bond.Bond(bond_master,bond_slave,np.linalg.norm(bond_master.Position - bond_slave.Position)))
        #self.Map_From_Bonds()

        for ring in self.Ring_List:
            if ring.Bond_List != "":
                self.Bond_List = self.Bond_List + ring.Bond_List

        for ring in self.Ring_List:
            if ring.Angle_List != "":
                self.Angle_List = self.Angle_List + ring.Angle_List

        for ring in self.Ring_List:
            if ring.Dihedral_List != "":
                self.Dihedral_List = self.Dihedral_List + ring.Dihedral_List

        for ring in self.Ring_List:
            if ring.Improper_List != "":
                self.Improper_List = self.Improper_List + ring.Improper_List

        #OPLS.Assign_OPLS(self,ChelpG=False)

        for ring in self.Ring_List:
            for b_atom in ring.Bonded_Atoms:
                if b_atom.Is_Linked:
                    Existing_Improper = False
                    for imp in self.Improper_List:
                        if imp.Improper_Master == b_atom.Central_Atom and imp.Improper_Slave3 == b_atom.Interring_Bond_Atom:
                            Existing_Improper = True
                    if not Existing_Improper:
                        self.Improper_List.append(Improper.Improper(b_atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[0],b_atom.Same_Ring_Bonded_Atom_List[1],b_atom.Interring_Bond_Atom.Central_Atom,b_atom.K,0.0,len(self.Improper_List)+1,d=b_atom.d,n=b_atom.n))

        total_mass = 0.0
        self.COM = np.array([0.0,0.0,0.0])
        self.MW = 0.0
        for atom in self.Atom_List:
            #print("Atom ID: %d OPLS Type: %d OPLS Class: %d" % (atom.Atom_ID,atom.OPLS_Type,atom.OPLS_Class))
            self.COM = total_mass / (atom.Mass + total_mass) * self.COM + atom.Mass / (atom.Mass + total_mass) * atom.Position
            total_mass += atom.Mass
        self.MW = total_mass
        self.Update_All_Normal_Vectors()
        bond_ID = 1
        for bond in self.Bond_List:
            bond.Bond_ID = bond_ID
            bond_ID += 1
        angle_ID = 1
        for angle in self.Angle_List:
            angle.Angle_ID = angle_ID
            angle_ID += 1
        dih_ID = 1
        for dih in self.Dihedral_List:
            dih.Dihedral_ID = dih_ID
            dih_ID += 1
        imp_ID = 1
        for imp in self.Improper_List:
            imp.Improper_ID = imp_ID
            imp_ID += 1
        # Attempt 11/19

        for i in range(len(self.Ring_List)-1):
            if (i + 1) % Flip == 0:
                Original_Position = self.Ring_List[i+1].Bonded_Atoms[0].Central_Atom.Position
                self.Rotate_Ring("Dih", 180, self.Ring_List[i + 1], self.Ring_List[i],No_Bonded_Vector_Update=True)
                #Rotation_Angle = math.acos(np.dot(self.Ring_List[i].Bonded_Atoms[-1].Bonded_Vector/np.linalg.norm(self.Ring_List[i].Bonded_Atoms[-1].Bonded_Vector), (self.Ring_List[i].Bonded_Atoms[-1].Central_Atom.Position - self.Ring_List[i+1].Bonded_Atoms[0].Interring_Bond_Atom.Central_Atom.Position) / np.linalg.norm(b_atom.Central_Atom.Position - b_atom.Interring_Bond_Atom.Central_Atom.Position)))
                Rotation_Angle = math.acos(np.dot(-1*self.Ring_List[i].Bonded_Atoms[0].Bonded_Vector/np.linalg.norm(self.Ring_List[i].Bonded_Atoms[0].Bonded_Vector),self.Ring_List[i].Bonded_Atoms[-1].Bonded_Vector/np.linalg.norm(self.Ring_List[i].Bonded_Atoms[-1].Bonded_Vector)))*180/math.pi
                Rotation_Matrix = [[1, 0, 0], [0, math.cos(Rotation_Angle), math.sin(Rotation_Angle)],[0, -math.sin(Rotation_Angle), math.cos(Rotation_Angle)]]
                print("Rotate")
                print(Rotation_Angle)
                self.Rotate_Ring("OOP", Rotation_Angle, self.Ring_List[i + 1], self.Ring_List[i]) #This actually rotates in plane. I don't know why
                Rotation_Matrix = [[1, 0, 0], [0, math.cos(Rotation_Angle), math.sin(Rotation_Angle)],[0, -math.sin(Rotation_Angle), math.cos(Rotation_Angle)]]
                #self.Ring_List[i + 1].Rotate_Ring(Rotation_Matrix)
                self.Ring_List[i+1].Translate_Ring(self.Ring_List[i+1].Bonded_Atoms[0].Central_Atom.Position)
                self.Ring_List[i+1].Translate_Ring(-1*self.Ring_List[i].Bonded_Atoms[-1].Bonded_Vector)
            """New_Position = self.Ring_List[i + 1].Bonded_Atoms[0].Central_Atom.Position
            Move_Distance = New_Position - Original_Position"""
            #self.Ring_List[i + 1].Translate_Ring(Move_Distance)
            #self.Ring_List[i+1].Translate_Ring(self.Ring_List[i].Bonded_Atoms[-1].Bonded_Vector)
        """
            #Additional fix
            for b_atom in self.Ring_List[i].Bonded_Atoms:
                if b_atom.Is_Linked:
                    print(b_atom.Interring_Bond_Atom.Self_Ring.Ring_ID)
                    print(i+1)
                if b_atom.Is_Linked and b_atom.Interring_Bond_Atom.Bonded_Ring.Ring_ID == i + 1:
                    if np.dot(b_atom.Bonded_Vector / np.linalg.norm(b_atom.Bonded_Vector),-1 * b_atom.Interring_Bond_Atom.Bonded_Vector / np.linalg.norm(b_atom.Interring_Bond_Atom.Bonded_Vector)) > 1:
                        rotation_angle = 0
                    else:
                        rotation_angle = -1 * math.acos(np.dot(b_atom.Bonded_Vector / np.linalg.norm(b_atom.Bonded_Vector),-1 * b_atom.Interring_Bond_Atom.Bonded_Vector / np.linalg.norm(b_atom.Interring_Bond_Atom.Bonded_Vector)))
                    print("Check")
                    print(b_atom.Bonded_Vector)
                    print(b_atom.Interring_Bond_Atom.Bonded_Vector )
                    print(rotation_angle)
                    self.Rotate_Ring("In-Plane",rotation_angle,self.Ring_List[i + 1], self.Ring_List[i])"""
        #end attempt
        #print("END OF CP INITIALIZATION")

    def exponential(L,P):
        return np.exp(L/P)

    def Rotate_Polymer(self,Rotation_Matrix):
        #rotates every ring in the polymer about a specified rotation matrix
        for ring in self.Ring_List:
            ring.Rotate_Ring(Rotation_Matrix)
            ring.Update_Normal_Vector()
            for b_atom in ring.Bonded_Atoms:
                b_atom.Bonded_Vector = np.matmul(Rotation_Matrix,b_atom.Bonded_Vector)

    def Translate_Polymer(self,Center_Point):
        #translates every ring in the polymer along a specified vector
        for ring in self.Ring_List:
            ring.Translate_Ring(Center_Point)

    def Recursive_Rotate(self,Last_Ring,Current_Ring,Rotation_Matrix):
        for b_atom in Current_Ring.Bonded_Atoms:
            if b_atom.Is_Linked and b_atom.Bonded_Ring != Last_Ring:
                b_atom.Bonded_Ring.Rotate_Ring(Rotation_Matrix)
                self.Recursive_Rotate(Current_Ring,b_atom.Bonded_Ring,Rotation_Matrix)

    def Refresh_Name(self):
        self.Name = ""
        for ring,phi,theta in zip(self.Ring_List[:-1],self.OOP_Rotation,self.Dih_Rotation):
            if ring not in self.Aux_Ring_List.keys():
                self.Name = self.Name + ring.Name + "_Phi_%d_Theta_%d_" % (phi,theta)
            else:
                self.Name = self.Name + ring.Name
                for i,aux_ring,aux_oop,aux_dih in zip(range(len(self.Aux_Ring_List[ring])),self.Aux_Ring_List[ring],self.Aux_Ring_OOP[ring],self.Aux_Ring_Dih[ring]):
                    self.Name = self.Name + "_Phi_%d_Theta_%d_Aux_%s_%d_" % (aux_oop,aux_dih,aux_ring.Name,i+1)
                self.Name = self.Name + "_Phi_%d_Theta_%d" % (phi,theta)
        self.Name = self.Name + self.Ring_List[-1].Name
        if self.Ring_List[-1] in self.Aux_Ring_List.keys():
            for i,aux_ring,aux_oop,aux_dih in zip(range(len(self.Aux_Ring_List[self.Ring_List[-1]])),self.Aux_Ring_List[self.Ring_List[-1]],self.Aux_Ring_OOP[self.Ring_List[-1]],self.Aux_Ring_Dih[self.Ring_List[-1]]):
                    self.Name = self.Name + "_Phi_%d_Theta_%d_Aux_%s_%d_" % (aux_oop,aux_dih,aux_ring.Name,i+1)

        self.Short_Name = ""
        for ring,phi,theta in zip(self.Ring_List[:-1],self.OOP_Rotation,self.Dih_Rotation):
            if ring not in self.Aux_Ring_List.keys():
                self.Short_Name = self.Short_Name + ring.Nickname + "_Phi_%d_Theta_%d_" % (phi,theta)
            else:
                self.Short_Name = self.Short_Name + ring.Nickname
                for i,aux_ring,aux_oop,aux_dih in zip(range(len(self.Aux_Ring_List[ring])),self.Aux_Ring_List[ring],self.Aux_Ring_OOP[ring],self.Aux_Ring_Dih[ring]):
                    self.Short_Name = self.Short_Name + "_Phi_%d_Theta_%d_Aux_%s_%d_" % (aux_oop,aux_dih,aux_ring.Nickname,i+1)
                self.Short_Name = self.Short_Name + "_Phi_%d_Theta_%d" % (phi,theta)
        self.Short_Name = self.Short_Name + self.Ring_List[-1].Nickname
        if self.Ring_List[-1] in self.Aux_Ring_List.keys():
            for i,aux_ring,aux_oop,aux_dih in zip(range(len(self.Aux_Ring_List[self.Ring_List[-1]])),self.Aux_Ring_List[self.Ring_List[-1]],self.Aux_Ring_OOP[self.Ring_List[-1]],self.Aux_Ring_Dih[self.Ring_List[-1]]):
                self.Short_Name = self.Short_Name + "_Phi_%d_Theta_%d_Aux_%s_%d_" % (aux_oop,aux_dih,aux_ring.Nickname,i+1)


    def Rotate_Ring(self,Rotation_Type,Rotation_Angle,Rotation_Ring,Linked_Rotation_Ring,No_Bonded_Vector_Update=False):
        #Rotates the Rotation_Ring (and all rings along the chain in the direction away from the Linked_Rotation_Ring) about the Linked_Rotation_Ring. Rotation_Type can be Dih or OOP, for a dihedral rotation or improper rotation respectively. Rotation_Angle is given in degrees and converted to radians
        print(Rotation_Angle)
        for b_atom in Rotation_Ring.Bonded_Atoms:
            if b_atom.Is_Linked and b_atom.Bonded_Ring == Linked_Rotation_Ring:
                rotation_bonded_atom = b_atom
        self.Translate_Polymer(rotation_bonded_atom.Interring_Bond_Atom.Central_Atom.Position)
        Linked_Rotation_Ring.Update_Normal_Vector()
        New_X = Linked_Rotation_Ring.Normal_Vector
        New_Y = rotation_bonded_atom.Bonded_Vector/np.linalg.norm(rotation_bonded_atom.Bonded_Vector)
        New_X = (New_X - np.dot(New_X,New_Y)*New_Y)/np.linalg.norm(New_X - np.dot(New_X,New_Y)*New_Y)
        New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
        self.Rotate_Polymer([New_X,New_Y,New_Z])

        Rotation_Angle_Degrees = Rotation_Angle
        Rotation_Angle = Rotation_Angle * math.pi/180
        if Rotation_Type == "Dih":
            Rotation_Matrix = [[math.cos(Rotation_Angle),0,-math.sin(Rotation_Angle)],[0,1,0],[math.sin(Rotation_Angle),0,math.cos(Rotation_Angle)]]
        elif Rotation_Type == "OOP":
            Rotation_Matrix = [[math.cos(Rotation_Angle),-math.sin(Rotation_Angle),0],[math.sin(Rotation_Angle),math.cos(Rotation_Angle),0],[0,0,1]]
        elif Rotation_Type == "In-Plane":
            Rotation_Matrix = [[1,0,0],[0,math.cos(Rotation_Angle),math.sin(Rotation_Angle)],[0,-math.sin(Rotation_Angle),math.cos(Rotation_Angle)]]
        else:
            raise Exception("Unidentified Rotation_Type"())
        print(Rotation_Angle)
        print('hola')
        print(Rotation_Matrix)
        Rotation_Ring.Rotate_Ring(Rotation_Matrix)
        Current_Ring = Rotation_Ring
        Last_Ring = Linked_Rotation_Ring
        self.Recursive_Rotate(Last_Ring,Current_Ring,Rotation_Matrix)
        if not No_Bonded_Vector_Update:
            self.Update_All_Bonded_Vectors()
        self.Update_All_Normal_Vectors()
        for i,ring in enumerate(self.Ring_List):
            if ring == Rotation_Ring or ring == Linked_Rotation_Ring:
                Rotation_Index = i
                if Rotation_Type == "Dih":
                    self.Dih_Rotation[i] += Rotation_Angle_Degrees
                    if self.Dih_Rotation[i] >= 360:
                        self.Dih_Rotation[i] -= 360
                    elif self.Dih_Rotation[i] < 0:
                        self.Dih_Rotation[i] += 360
                if Rotation_Type == "OOP":
                    self.OOP_Rotation[i] += Rotation_Angle_Degrees
                break

        self.Refresh_Name()

    def Update_All_Normal_Vectors(self):
        for ring in self.Ring_List:
            ring.Update_Normal_Vector()

    def Rotate_Aux_Ring(self,Rotation_Type,Rotation_Angle,Rotation_Ring,Aux_Ring_Index):
        #Rotates the Rotation_Ring (and all rings along the chain in the direction away from the Linked_Rotation_Ring) about the Linked_Rotation_Ring. Rotation_Type can be Dih or OOP, for a dihedral rotation or improper rotation respectively. Rotation_Angle is given in degrees and converted to radians
        self.Translate_Polymer(Rotation_Ring.Aux_Ring_List[Aux_Ring_Index].Other_Bond_Atom.Position)
        otation_Ring.Update_Normal_Vector()
        New_X = Rotation_Ring.Normal_Vector
        New_Y = (Rotation_Ring.Aux_Ring_List[Aux_Ring_Index].Other_Bond_Atom.Position - Rotation_Ring.Aux_Ring_List[Aux_Ring_Index].Self_Bond_Atom.Position)/np.linalg.norm(Rotation_Ring.Aux_Ring_List[Aux_Ring_Index].Other_Bond_Atom.Position - Rotation_Ring.Aux_Ring_List[Aux_Ring_Index].Self_Bond_Atom.Position)
        New_X = (New_X - np.dot(New_X,New_Y)*New_Y)/np.linalg.norm(New_X - np.dot(New_X,New_Y)*New_Y)
        New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
        self.Rotate_Polymer([New_X,New_Y,New_Z])

        Rotation_Angle_Degrees = Rotation_Angle
        Rotation_Angle = Rotation_Angle * math.pi/180
        if Rotation_Type == "Dih":
            Rotation_Matrix = [[math.cos(Rotation_Angle),0,-math.sin(Rotation_Angle)],[0,1,0],[math.sin(Rotation_Angle),0,math.cos(Rotation_Angle)]]
        elif Rotation_Type == "OOP":
            Rotation_Matrix = [[math.cos(Rotation_Angle),-math.sin(Rotation_Angle),0],[math.sin(Rotation_Angle),math.cos(Rotation_Angle),0],[0,0,1]]
        else:
            raise Exception("Unidentified Rotation_Type")

        Rotation_Ring.Aux_Ring_List[Aux_Ring_Index].Rotate_Ring(Rotation_Matrix)
        self.Update_All_Bonded_Vectors()
        for ring in self.Ring_List:
            if ring == Rotation_Ring:
                if Rotation_Type == "Dih":
                    self.Aux_Ring_Dih[Rotation_Ring][Aux_Ring_Index] += Rotation_Angle_Degrees
                    if self.Aux_Ring_Dih[Rotation_Ring][Aux_Ring_Index] >= 360:
                        self.Aux_Ring_Dih[Rotation_Ring][Aux_Ring_Index] -= 360
                    elif self.Aux_Ring_Dih[Rotation_Ring][Aux_Ring_Index] < 0:
                        self.Aux_Ring_Dih[Rotation_Ring][Aux_Ring_Index] += 360
                if Rotation_Type == "OOP":
                    self.Aux_Ring_OOP[Rotation_Ring][Aux_Ring_Index] += Rotation_Angle_Degrees
                break

        self.Refresh_Name()

    def Rotate_Ring_Isolated(self,Rotation_Type,Rotation_Angle,Rotation_Ring,Linked_Rotation_Ring):
        #Rotates ONLY the Rotation_Ring  about the Linked_Rotation_Ring. Rotation_Type can be Dih or OOP, for a dihedral rotation or improper rotation respectively. Rotation_Angle is given in degrees and converted to radians
        for b_atom in Rotation_Ring.Bonded_Atoms:
            if b_atom.Is_Linked and b_atom.Bonded_Ring == Linked_Rotation_Ring:
                rotation_bonded_atom = b_atom
        self.Translate_Polymer(b_atom.Interring_Bond_Atom.Central_Atom.Position)
        Linked_Rotation_Ring.Update_Normal_Vector()
        New_X = Linked_Rotation_Ring.Normal_Vector
        New_Y = rotation_bonded_atom.Bonded_Vector/np.linalg.norm(rotation_bonded_atom.Bonded_Vector)
        New_X = (New_X - np.dot(New_X,New_Y)*New_Y)/np.linalg.norm(New_X - np.dot(New_X,New_Y)*New_Y)
        New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
        self.Rotate_Polymer([New_X,New_Y,New_Z])
        Rotation_Angle_Degrees = Rotation_Angle
        Rotation_Angle = Rotation_Angle * math.pi/180
        if Rotation_Type == "Dih":
            Rotation_Matrix = [[math.cos(Rotation_Angle),0,-math.sin(Rotation_Angle)],[0,1,0],[math.sin(Rotation_Angle),0,math.cos(Rotation_Angle)]]
        elif Rotation_Type == "OOP":
            Rotation_Matrix = [[math.cos(Rotation_Angle),-math.sin(Rotation_Angle),0],[math.sin(Rotation_Angle),math.cos(Rotation_Angle),0],[0,0,1]]
        else:
            raise Exception("Unidentified Rotation_Type")
        Rotation_Ring.Rotate_Ring(Rotation_Matrix)
        self.Update_All_Bonded_Vectors()
        for i,ring in enumerate(self.Ring_List):
            if ring == Rotation_Ring or ring == Linked_Rotation_Ring:
                Rotation_Index = i
                if Rotation_Type == "Dih":
                    self.Dih_Rotation[i] += Rotation_Angle_Degrees
                    if self.Dih_Rotation[i] >= 360:
                        self.Dih_Rotation -= 360
                    elif self.Dih_Rotation[i] < 0:
                        self.Dih_Rotation += 360
                    if self.Ring_List[i+1] == Linked_Rotation_Ring:
                        self.Dih_Rotation[i-1] -= Rotation_Angle_Degrees
                        if self.Dih_Rotation[i-1] >= 360:
                            self.Dih_Rotation[i-1] -= 360
                        elif self.Dih_Rotation[i-1] < 0:
                            self.Dih_Rotation[i-1] += 360
                    else:
                        self.Dih_Rotation[i+1] -= Rotation_Angle_Degrees
                        if self.Dih_Rotation[i+1] >= 360:
                            self.Dih_Rotation[i+1] -= 360
                        elif self.Dih_Rotation[i+1] < 0:
                            self.Dih_Rotation[i+1] += 360
                if Rotation_Type == "OOP":
                    self.OOP_Rotation[i] += Rotation_Angle_Degrees
                    if self.Ring_List[i+1] == Linked_Rotation_Ring:
                        self.OOP_Rotation[i-1] -= Rotation_Angle_Degrees
                    else:
                        self.OOP_Rotation[i+1] -= Rotation_Angle_Degrees
                break

        self.Refresh_Name()

    def Update_All_Bonded_Vectors(self):
        for ring in self.Ring_List:
            for b_atom in ring.Bonded_Atoms:
                b_atom.Update_Bonded_Vector()

    def Write_XYZ(self,Name = ""):
        if self.Num_Rings <= 4 and Name == "":
            XYZ_Filename = "%s.xyz" % (self.Name)
            f = open("%s.xyz" % (self.Name),'w')
        elif Name == "":
            XYZ_Filename = "%s.xyz" % (self.Short_Name)
            f = open("%s.xyz" % (self.Short_Name),'w')
        if Name != "":
            XYZ_Filename = "%s.xyz" % Name
            f = open("%s.xyz" % Name, 'w')
        f.write("%d\n\n" % (len(self.Atom_List)))
        i = 1
        for atom in self.Atom_List:
            f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        f.close()
        return XYZ_Filename

    def Write_Methylated_Dimer(self,Ring1_Num,Atoms_Per_Ring,Print_List,H_Atom_f,H_Atom_b,HB_Atom_f,HB_Atom_b,Methyl_Atoms_Ring1,Methyl_Atoms_Ring2,Name=""):
        Ring2_Num = Ring1_Num + 1
        if self.Num_Rings <= 4 and Name == "":
            XYZ_Filename = "%s.xyz" % (self.Name)
            f = open("%s.xyz" % (self.Name),'w')
        elif Name == "":
            XYZ_Filename = "%s.xyz" % (self.Short_Name)
            f = open("%s.xyz" % (self.Short_Name),'w')
        if Name != "":
            XYZ_Filename = "%s.xyz" % Name
            f = open("%s.xyz" % Name, 'w')
        f.write("%d\n\n" % (len(Print_List)+len(Print_List) + 8))
        #write only atoms in Print_List
        for pl in Print_List:
            for atom in self.Atom_List:
                if atom.Atom_ID == Ring1_Num*Atoms_Per_Ring+pl+1:
                    f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
                    break
        for pl in Print_List:
            for atom in self.Atom_List:
                if atom.Atom_ID == Ring2_Num*Atoms_Per_Ring+pl+1:
                    f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
                    break
        extra = 1
        for atom in self.Atom_List:
            if atom.Atom_ID == (Ring1_Num) * Atoms_Per_Ring + HB_Atom_b + 1 + extra:
                H_Bonded_Atom_b = atom.Position
            if atom.Atom_ID == (Ring2_Num) * Atoms_Per_Ring + HB_Atom_f + 1 + extra:
                H_Bonded_Atom_f = atom.Position
        #Change bonded points to outside bonds to Hs and scale
        for atom in self.Atom_List:
            if atom.Atom_ID == (Ring1_Num-1) * Atoms_Per_Ring + H_Atom_b + 1 + extra:
                back_atom_pos = (atom.Position - H_Bonded_Atom_b)/np.linalg.norm(atom.Position - H_Bonded_Atom_b)*1.08+H_Bonded_Atom_b
                f.write("%s\t%f\t%f\t%f\n" % ("H", back_atom_pos[0], back_atom_pos[1], back_atom_pos[2]))
            if atom.Atom_ID == (Ring2_Num+1)*Atoms_Per_Ring+H_Atom_f+1+extra:
                front_atom_pos = (atom.Position - H_Bonded_Atom_f) / np.linalg.norm(atom.Position - H_Bonded_Atom_f) * 1.08 + H_Bonded_Atom_f
                f.write("%s\t%f\t%f\t%f\n" % ("H", front_atom_pos[0], front_atom_pos[1], front_atom_pos[2]))
        #Add hydrogens to methyl group for first ring
        for atom in self.Atom_List:
            if Ring1_Num % 2 == 0:
                if atom.Atom_ID == (Ring1_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring1[0] + 1 + extra:
                    print(atom.Atom_ID)
                    Internal_Methyl = atom.Position
                if atom.Atom_ID == (Ring1_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring1[1] + 1 + extra:
                    External_Methyl = atom.Position
            else:
                if atom.Atom_ID == (Ring1_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring2[0] + 1 + extra:
                    print(atom.Atom_ID)
                    Internal_Methyl = atom.Position
                if atom.Atom_ID == (Ring1_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring2[1] + 1 + extra:
                    External_Methyl = atom.Position
        New_X = External_Methyl - Internal_Methyl
        New_X = New_X / np.linalg.norm(New_X)
        New_Y = np.cross(self.Atom_List[0].Position - self.Atom_List[1].Position,self.Atom_List[1].Position - self.Atom_List[2].Position)
        New_Y = New_Y/np.linalg.norm(New_Y)
        New_Y = (New_Y - np.dot(New_X,New_Y)*New_X)/np.linalg.norm(New_Y - np.dot(New_X,New_Y)*New_X)
        New_Z = np.cross(New_X,New_Y)
        #Set Position near origin, scale to 1.08 A
        New_H_Position_1 = New_X*1.08
        #Transform to new coordinate system
        New_H_Position_1 = np.matmul([New_X,New_Y,New_Z],New_H_Position_1)
        # Rotate to tetrahedral formation
        print(New_H_Position_1)
        print([[np.cos(72.0/180.0*np.pi),np.sin(72.0/180.0*np.pi),0], [np.sin(72.0/180.0*np.pi),-1*np.cos(72.0/180.0*np.pi),0], [0,0,1]])
        New_H_Position_1 = np.matmul([[np.cos(72.0/180.0*np.pi),np.sin(72.0/180.0*np.pi),0], [np.sin(72.0/180.0*np.pi),-1*np.cos(72.0/180.0*np.pi),0], [0,0,1]],New_H_Position_1)
        print(New_H_Position_1)
        #create other two H positions
        New_H_Position_2 = np.matmul([[1,0,0],[0,np.cos(-2*np.pi/3),np.sin(-2*np.pi/3)],[0,-np.sin(-2*np.pi/3),np.cos(-2*np.pi/3)]],New_H_Position_1)
        New_H_Position_3 = np.matmul([[1,0,0],[0,np.cos(2*np.pi/3),np.sin(2*np.pi/3)],[0,-np.sin(2*np.pi/3),np.cos(2*np.pi/3)]],New_H_Position_1)
        #conver back to original coordinate system
        New_H_Position_1 = np.matmul(np.linalg.inv([New_X,New_Y,New_Z]),New_H_Position_1)
        New_H_Position_2 = np.matmul(np.linalg.inv([New_X, New_Y, New_Z]), New_H_Position_2)
        New_H_Position_3 = np.matmul(np.linalg.inv([New_X, New_Y, New_Z]), New_H_Position_3)
        #Place by External_Methyl
        New_H_Position_1 += External_Methyl
        New_H_Position_2 += External_Methyl
        New_H_Position_3 += External_Methyl
        f.write("%s\t%f\t%f\t%f\n" % ("H", New_H_Position_1[0], New_H_Position_1[1], New_H_Position_1[2]))
        f.write("%s\t%f\t%f\t%f\n" % ("H", New_H_Position_2[0], New_H_Position_2[1], New_H_Position_2[2]))
        f.write("%s\t%f\t%f\t%f\n" % ("H", New_H_Position_3[0], New_H_Position_3[1], New_H_Position_3[2]))
        # Add hydrogens to methyl group for second ring
        for atom in self.Atom_List:
            if Ring2_Num % 2 == 0:
                if atom.Atom_ID == (Ring2_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring1[0] + 1 + extra:
                    Internal_Methyl = atom.Position
                if atom.Atom_ID == (Ring2_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring1[1] + 1 + extra:
                    External_Methyl = atom.Position
            else:
                if atom.Atom_ID == (Ring2_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring2[0] + 1 + extra:
                    Internal_Methyl = atom.Position
                if atom.Atom_ID == (Ring2_Num) * Atoms_Per_Ring + Methyl_Atoms_Ring2[1] + 1 + extra:
                    External_Methyl = atom.Position
        New_X = External_Methyl - Internal_Methyl
        New_X = New_X / np.linalg.norm(New_X)
        New_Y = np.cross(self.Atom_List[0].Position - self.Atom_List[1].Position,self.Atom_List[1].Position - self.Atom_List[2].Position)
        New_Y = New_Y/np.linalg.norm(New_Y)
        New_Y = (New_Y - np.dot(New_X,New_Y)*New_X)/np.linalg.norm(New_Y - np.dot(New_X,New_Y)*New_X)
        New_Z = np.cross(New_X,New_Y)
        #Set Position near origin, scale to 1.08 A
        New_H_Position_1 = New_X*1.08
        #Transform to new coordinate system
        New_H_Position_1 = np.matmul([New_X,New_Y,New_Z],New_H_Position_1)
        # Rotate to tetrahedral formation
        print(New_H_Position_1)
        print([[np.cos(72.0/180.0*np.pi),np.sin(72.0/180.0*np.pi),0], [np.sin(72.0/180.0*np.pi),-1*np.cos(72.0/180.0*np.pi),0], [0,0,1]])
        New_H_Position_1 = np.matmul([[np.cos(72.0/180.0*np.pi),np.sin(72.0/180.0*np.pi),0], [np.sin(72.0/180.0*np.pi),-1*np.cos(72.0/180.0*np.pi),0], [0,0,1]],New_H_Position_1)
        print(New_H_Position_1)
        #create other two H positions
        New_H_Position_2 = np.matmul([[1,0,0],[0,np.cos(-2*np.pi/3),np.sin(-2*np.pi/3)],[0,-np.sin(-2*np.pi/3),np.cos(-2*np.pi/3)]],New_H_Position_1)
        New_H_Position_3 = np.matmul([[1,0,0],[0,np.cos(2*np.pi/3),np.sin(2*np.pi/3)],[0,-np.sin(2*np.pi/3),np.cos(2*np.pi/3)]],New_H_Position_1)
        #conver back to original coordinate system
        New_H_Position_1 = np.matmul(np.linalg.inv([New_X,New_Y,New_Z]),New_H_Position_1)
        New_H_Position_2 = np.matmul(np.linalg.inv([New_X, New_Y, New_Z]), New_H_Position_2)
        New_H_Position_3 = np.matmul(np.linalg.inv([New_X, New_Y, New_Z]), New_H_Position_3)
        #Place by External_Methyl
        New_H_Position_1 += External_Methyl
        New_H_Position_2 += External_Methyl
        New_H_Position_3 += External_Methyl
        f.write("%s\t%f\t%f\t%f\n" % ("H", New_H_Position_1[0], New_H_Position_1[1], New_H_Position_1[2]))
        f.write("%s\t%f\t%f\t%f\n" % ("H", New_H_Position_2[0], New_H_Position_2[1], New_H_Position_2[2]))
        f.write("%s\t%f\t%f\t%f\n" % ("H", New_H_Position_3[0], New_H_Position_3[1], New_H_Position_3[2]))
        f.close()

    def Return_Symmetry_Name(self):
        Sym_Name = ""
        for ring,phi,theta in zip(self.Ring_List[len(self.Ring_List)-1::-1],self.OOP_Rotation[::-1],self.Dih_Rotation[::-1]):
            if ring.Aux_Ring_List == []:
                Sym_Name = Sym_Name+ ring.Name + "_Phi_%d_Theta_%d_" % (phi,theta)
            else:
                Sym_Name = Sym_Name + ring.Name
                for i,aux_ring,aux_phi,aux_theta in zip(range(len(ring.Aux_Ring_List)),ring.Aux_Ring_List,self.Aux_Ring_OOP[ring],self.Aux_Ring_Dih[ring]):
                    Sym_Name = Sym_Name + "_Phi_%d_Theta_%d_Aux_%s_%d" % (aux_phi,aux_theta,aux_ring.Name,i+1)
                Sym_Name = Sym_Name + "_Phi_%d_Theta_%d_" % (phi,theta)  
        Sym_Name = Sym_Name + self.Ring_List[0].Name
        if self.Ring_List[0].Aux_Ring_List != []:
            for i,aux_ring,aux_phi,aux_theta in zip(range(len(ring.Aux_Ring_List)),ring.Aux_Ring_List,self.Aux_Ring_OOP[ring],self.Aux_Ring_Dih[ring]):
                Sym_Name = Sym_Name + "_Phi_%d_Theta_%d_Aux_%s_%d" % (aux_phi,aux_theta,aux_ring.Name,i+1)
        return Sym_Name

    def Calculate_Total_Charge(self):
        Q = 0.0
        for atom in self.Atom_List:
            Q += atom.Charge
        return Q

    def Zero_Total_Charge(self):
        Q = 0.0
        for atom in self.Atom_List:
            Q += atom.Charge
        dQ = Q/len(self.Atom_List)
        Q = 0.0
        for atom in self.Atom_List:
            atom.Charge -= dQ
            Q += atom.Charge
        return Q

    def Calculate_Internal_Energy(self,Polymer_Name,lammps_bin,Exclude_All_Interring=False,Exclude_Interring_Bonds=False,Exclude_Interring_Angles=False,Exclude_Interring_Torsions=False,Interring_Angles_Only=False,dielectric=1.0,Symmetric = False,Adjusted_Angles=False,Plumed=""):
        LAMMPS_Filename = Polymer_Name + "_" + self.Name + "_Internal_Energy"
        if Adjusted_Angles:
            LAMMPS_Filename = LAMMPS_Filename + "_Adjusted_Angles"
        elif Exclude_All_Interring:
            LAMMPS_Filename = LAMMPS_Filename + "_No_Interring"
        elif Exclude_Interring_Bonds:
            LAMMPS_Filename = LAMMPS_Filename + "_No_Interring_Bonds"
        elif Exclude_Interring_Angles:
            LAMMPS_Filename = LAMMPS_Filename + "_No_Interring_Angles"
        elif Exclude_Interring_Torsions:
            LAMMPS_Filename = LAMMPS_Filename + "_No_Interring_Torsions"
        elif Interring_Angles_Only:
            LAMMPS_Filename = LAMMPS_Filename + "_Interring_Angles_Only"
        #end delete
        if not os.path.isfile("./Nontorsional_Outputs/log.%s" % LAMMPS_Filename):
            LAMMPS_Filename = LAMMPS_Filename + ".data"
            if Exclude_Interring_Bonds or Exclude_Interring_Torsions or Exclude_Interring_Angles or Exclude_All_Interring or Interring_Angles_Only:
                Bonded_Atom_List = []
                for ring in self.Ring_List:
                    for b_atom in ring.Bonded_Atoms:
                        if b_atom.Is_Linked:
                            Bonded_Atom_List.append(b_atom.Central_Atom)
            else:
                Bonded_Atom_List = []
            self.Write_Data_File(LAMMPS_Filename,Exclude_All_Interring=Exclude_All_Interring,Exclude_Interring_Bonds=Exclude_Interring_Bonds,Exclude_Interring_Angles=Exclude_Interring_Angles,Exclude_Interring_Torsions=Exclude_Interring_Torsions,Interring_Angles_Only=Interring_Angles_Only,Bond_Atoms=Bonded_Atom_List)
            Write_Inputs.Write_LAMMPS_Single_Point_Energy("in.%s" % LAMMPS_Filename.split('.data')[0],LAMMPS_Filename,LAMMPS_Filename.split('.data')[0],self.Find_Coul_Cutoff(),dielectric=dielectric,Plumed=Plumed)
            os.system("%s -in %s -log log.%s" % (lammps_bin,"in.%s" % LAMMPS_Filename.split('.data')[0],LAMMPS_Filename.split('.data')[0]))
            f = open("log.%s" % LAMMPS_Filename.split('.data')[0],'r')
            lines = f.readlines()
            Read_Energy = False
            Internal_Energy = 0.0
            for line in lines:
                if Read_Energy:
                    Internal_Energy = float(line.split()[-1].strip())
                    Read_Energy = False
                if len(line.split()) > 0 and line.split()[0].strip() == "Step":
                    Read_Energy = True
            if Internal_Energy == 0.0:
                print(LAMMPS_Filename)
                #raise Exception("Energy Not Found")
            os.system("scp %s ./Nontorsional_Inputs" % (LAMMPS_Filename))
            os.system("scp %s ./Nontorsional_Inputs" % ("in.%s" % LAMMPS_Filename.split('.data')[0]))
            os.system("scp %s ./Nontorsional_Outputs" % ("log.%s" % LAMMPS_Filename.split('.data')[0]))
            if Symmetric:
                Symmetry_LAMMPS_Name = Polymer_Name + "_" + self.Return_Symmetry_Name() + "_Internal_Energy"
                if Exclude_All_Interring:
                    Symmetry_LAMMPS_Name = Symmetry_LAMMPS_Name + "_No_Interring"
                elif Exclude_Interring_Bonds:
                    Symmetry_LAMMPS_Name = Symmetry_LAMMPS_Name + "_No_Interring_Bonds"
                elif Exclude_Interring_Angles:
                    Symmetry_LAMMPS_Name = Symmetry_LAMMPS_Name + "_No_Interring_Angles"
                elif Exclude_Interring_Torsions:
                    Symmetry_LAMMPS_Name = Symmetry_LAMMPS_Name + "_No_Interring_Torsions"
                elif Interring_Angles_Only:
                    LAMMPS_Filename = LAMMPS_Filename + "_Interring_Angles_Only"
                os.system("scp %s ./Nontorsional_Outputs/log.%s" % ("log.%s" % LAMMPS_Filename.split('.data')[0],Symmetry_LAMMPS_Name))
            os.system("rm -f %s" % (LAMMPS_Filename))
            os.system("rm -f %s" % ("in.%s" % LAMMPS_Filename.split('.data')[0]))
            os.system("rm -f %s" % ("log.%s" % LAMMPS_Filename.split('.data')[0]))
        else:
            f = open("./Nontorsional_Outputs/log.%s" % LAMMPS_Filename.split('.data')[0],'r')
            lines = f.readlines()
            Read_Energy = False
            Internal_Energy = 0.0
            for line in lines:
                if Read_Energy:
                    Internal_Energy = float(line.split()[-1].strip())
                    Read_Energy = False
                if len(line.split()) > 0 and line.split()[0].strip() == "Step":
                    Read_Energy = True
            if Internal_Energy == 0.0:
                print("./Nontorsional_Outputs/log.%s" % LAMMPS_Filename.split('.data')[0])
                #raise Exception("Energy Not Found")
        if Plumed != "":
            os.system("scp ./%s.txt ./Nontorsional_Outputs" % self.Name)
            f = open('./Nontorsional_Outputs/%s.txt' % self.Name)
            Lines = f.readlines()
            f.close()
            Internal_Energy += float(Lines[-1].split()[-1])/4.184

        return Internal_Energy

    def Write_Plumed_File(self,Base_File_Name,Fit_Energies,Force_y,Force_x,Delocalization_Modifier=1.0,Nonbonded_Modifier=1.0,Non_Interacting=False,Reverse=False,Folder_Name="./Full_Polymer_Plumed",Bias=True):
        #Write_Flags is a list of Booleans length of self.Molecule_List, representing whether to write a Plumed script for it
        #Molecules in self.Molecule_List must be Conjugated_Polymer objects
        Base_Calculate_Torsions_String = Base_File_Name + "_Calculate_Torsions"
        index = 1
        Num_Atoms = 1
        Num_Rings = 0
        if Delocalization_Modifier != 1:
            if Delocalization_Modifier == 0.0:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_No_Deloc"
            elif Delocalization_Modifier == 0.5:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_Half_Deloc"
            elif Delocalization_Modifier == 2.0:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_Double_Deloc"
            else:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_%.1f_Deloc" % Delocalization_Modifier
        if Nonbonded_Modifier != 1:
            if Nonbonded_Modifier == 0.0:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_No_Nonbonded"
            elif Nonbonded_Modifier == 0.5:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_Half_Nonbonded"
            elif Nonbonded_Modifier == 2.0:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_Double_Nonbonded"
            else:
                Base_Calculate_Torsions_String = Base_Calculate_Torsions_String + "_%.1f_Nonbonded" % Delocalization_Modifier
        if Bias:
            torsion_file = open(Base_Calculate_Torsions_String + ".dat", 'w')
        else:
            torsion_file = open(Base_Calculate_Torsions_String + "_Conventional.dat", 'w')

        plumed_index = 0
        oop_blocks = np.linspace(0,math.pi,200)
        dih_blocks = np.linspace(-math.pi,math.pi,200)

        if Non_Interacting == True:
            for atom in self.Atom_List:
                atom.Charge = 0.000000
                atom.Sigma = 0.0
                atom.Epsilon = 0.0
        torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (index*2-1,self.Ring_List[0].Plumed_Rings[plumed_index][0].Atom_ID))
        for atom in self.Ring_List[0].Core_Atom_List[1:]:
            torsion_file.write(",%d" % atom.Atom_ID)
        torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c1,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (index*2,self.Ring_List[0].Plumed_Rings[plumed_index][0].Atom_ID,self.Ring_List[0].Plumed_Rings[plumed_index][1].Atom_ID))

        Parameter_Index = 0
        for i in range(0,len(self.Ring_List)-1):
            if Bias:
                Bias_File_Name = "%s_Bias_File_%d" % (Base_File_Name,index)
                if Delocalization_Modifier != 1:
                    if Delocalization_Modifier == 0.0:
                        Bias_File_Name = Bias_File_Name + "_No_Deloc"
                    elif Delocalization_Modifier == 0.5:
                        Bias_File_Name = Bias_File_Name + "_Half_Deloc"
                    elif Delocalization_Modifier == 2.0:
                        Bias_File_Name = Bias_File_Name + "_Double_Deloc"
                    else:
                        Bias_File_Name = Bias_File_Name + "_%.1f_Deloc" % Delocalization_Modifier
                if Nonbonded_Modifier != 1:
                    if Nonbonded_Modifier == 0.0:
                        Bias_File_Name = Bias_File_Name + "_No_Nonbonded"
                    elif Nonbonded_Modifier == 0.5:
                        Bias_File_Name = Bias_File_Name + "_Half_Nonbonded"
                    elif Nonbonded_Modifier == 2.0:
                        Bias_File_Name = Bias_File_Name + "_Double_Nonbonded"
                    else:
                        Bias_File_Name = Bias_File_Name + "_%.1f_Nonbonded" % Delocalization_Modifier
            index += 1
            if Bias:
                f = open(Bias_File_Name + ".dat", 'w')
                f.write("#! FIELDS DIH_%d OOP_%d ext_%d.bias der_DIH_%d der_OOP_%d\n#! SET min_DIH_%d -pi\n#! SET max_DIH_%d pi\n#! SET nbins_DIH_%d 199\n#! SET periodic_DIH_%d true\n#! SET min_OOP_%d 0\n#! SET max_OOP_%d pi\n#! SET nbins_OOP_%d 199\n#! SET periodic_OOP_%d false\n" % (index,index,index,index,index,index,index,index,index,index,index,index,index))
                for q,energy_list in enumerate(Fit_Energies[Parameter_Index][:-1]):
                    for j,e in enumerate(energy_list[:-1]):
                        f.write("%.6f %.6f %.6f %.6f %.6f\n" % (dih_blocks[j],oop_blocks[q],e*Delocalization_Modifier,Force_x[Parameter_Index][q][j]*Delocalization_Modifier,Force_y[Parameter_Index][q][j]*Delocalization_Modifier))
                    f.write("\n")
                f.close()

                if Non_Interacting:
                    f = open(Bias_File_Name + "_Nonbonded.dat", 'w')
                    f.write("#! FIELDS DIH_%d OOP_%d ext_%d_nb.bias der_DIH_%d der_OOP_%d\n#! SET min_DIH_%d -pi\n#! SET max_DIH_%d pi\n#! SET nbins_DIH_%d 199\n#! SET periodic_DIH_%d true\n#! SET min_OOP_%d 0\n#! SET max_OOP_%d 2*pi\n#! SET nbins_OOP_%d 199\n#! SET periodic_OOP_%d false\n" % (index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1))
                    for q,energy_list in enumerate(Nonbonded_Fit_Energies[Parameter_Index][:-1]):
                        for j,e in enumerate(energy_list[:-1]):
                            f.write("%.6f %.6f %.6f %.6f %.6f\n" % (dih_blocks[j],oop_blocks[q],e*Nonbonded_Modifier,Nonbonded_Force_x[Parameter_Index][q][j]*Nonbonded_Modifier,Nonbonded_Force_y[Parameter_Index][q][j]*Nonbonded_Modifier))
                        f.write("\n")
                    f.close()
                os.system("scp %s.dat ./Nontorsional_Inputs" % (Bias_File_Name))
                #os.system("rm -f %s.dat" % (Bias_File_Name))
            index += 1
            torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (index*2-3, self.Ring_List[i].Plumed_Rings[0][0].Atom_ID)) #Writes the first atom in the first plumed ring
            for atom in self.Ring_List[i].Plumed_Rings[0][1:]:
                torsion_file.write(",%d" % atom.Atom_ID)
            torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (index*2-2, self.Ring_List[i+1].Plumed_Rings[1][0].Atom_ID)) #Writes the first atom in the second plumed ring
            for atom in self.Ring_List[i+1].Plumed_Rings[1][1:]:
                torsion_file.write(",%d" % atom.Atom_ID)
            for b_atom in self.Ring_List[i].Bonded_Atoms:
                if b_atom.Is_Linked and b_atom.Interring_Bond_Atom.Self_Ring.Ring_ID == i+2:
                    atom1 = b_atom.Interring_Bond_Atom.Central_Atom.Atom_ID
                    atom2 = b_atom.Central_Atom.Atom_ID
                    ring1_batom1 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Atom_ID
                    ring1_batom2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Atom_ID
                    ring2_batom1 = b_atom.Same_Ring_Bonded_Atom_List[0].Atom_ID
                    ring2_batom2 = b_atom.Same_Ring_Bonded_Atom_List[1].Atom_ID
            torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c%d,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (index*2-3,index*2-3,self.Ring_List[i].Plumed_Rings[0][0].Atom_ID,self.Ring_List[i].Plumed_Rings[0][1].Atom_ID))
            torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c%d,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (index*2-2,index*2-2,self.Ring_List[i+1].Plumed_Rings[1][0].Atom_ID,self.Ring_List[i+1].Plumed_Rings[1][1].Atom_ID))
            plumed_index += 1
            if plumed_index >= 2:
                plumed_index = 0
            if Reverse:
                try:
                    torsion_file.write("DIH_%d: TORSION VECTOR1=c%d_normal,c%d AXIS=%d,%d VECTOR2=c%d_normal,c%d\n\n" % (index-1,index*2-3,index*2-3,atom2,atom1,index*2-2,index*2-2))
                except:
                    torsion_file.write("DIH_%d: TORSION VECTOR1=c%d_normal,c%d AXIS=%d,%d VECTOR2=c%d_normal,c%d\n\n" % (index - 1, index*2-3, index*2-3, 0, 0, index*2-2, index*2-2))
            else:
                try:
                    torsion_file.write("DIH_%d: TORSION VECTOR1=c%d,c%d_normal AXIS=%d,%d VECTOR2=c%d,c%d_normal\n\n" % (index-1,index*2-3,index*2-3,atom1,atom2,index*2-2,index*2-2))
                except:
                    torsion_file.write("DIH_%d: TORSION VECTOR1=c%d,c%d_normal AXIS=%d,%d VECTOR2=c%d,c%d_normal\n\n" % (index - 1, index*2-3,index*2-3, 0, 1, index*2-2, index*2-2))
            torsion_file.write("OOP_%d_1: TORSION ATOMS=%d,%d,%d,%d\n\n" % (index-1,atom1,ring1_batom1,ring1_batom2,atom2))
            torsion_file.write("OOP_%d_2: TORSION ATOMS=%d,%d,%d,%d\n\n" % (index-1,atom2,ring2_batom1,ring2_batom2,atom1))
            torsion_file.write("OOP_%d: MATHEVAL ARG=OOP_%d_1,OOP_%d_2 FUNC=abs(x)+abs(y) PERIODIC=NO\n\n" % (index-1,index-1,index-1))
            if not Bias:
                Conventional_Add = "conventional_"
            else:
                Conventional_Add = ""
            if Bias:
                torsion_file.write("EXTERNAL ARG=DIH_%d,OOP_%d FILE=%s.dat LABEL=ext_%d\n\n" % (index-1,index-1,Bias_File_Name,index-1))
            if Non_Interacting:
                torsion_file.write("EXTERNAL ARG=DIH_%d,OOP_%d FILE=/scratch/andrewk/job_SLURM_JOBID/%s_Nonbonded.dat LABEL=ext_%d_nb\n\n" % (i,i,Bias_File_Name,i))
            Parameter_Index += 1
            if Parameter_Index == len(Fit_Energies)/2:
                Parameter_Index = 0

        os.system("mkdir ./%s" % (Folder_Name))

        index += 1


        torsion_file.write("WHOLEMOLECULES")
        torsion_file.write(" ENTITY0=1-%d,c3,c3_normal,c4,c4_normal" % len(self.Atom_List))
        torsion_file.write("\n\nPRINT ARG=ext_2.bias FILE=%s.txt STRIDE=1\n\n" % self.Name)
        torsion_file.close()
        os.system("scp %s ./Nontorsional_Inputs" % (Base_Calculate_Torsions_String + ".dat"))

    def Map_From_Bonds(self):
        for atom1 in self.Atom_List:
            for atom2 in atom1.Bond_List:
                for atom3 in atom2.Bond_List:
                    if atom1 != atom3:
                        Identical_Angle = False
                        for angle in self.Angle_List:
                            if angle.Compare_Angles(atom1,atom2,atom3):
                                Identical_Angle = True
                                break
                        if not Identical_Angle:
                            self.Angle_List.append(Angle.Angle(atom2,atom1,atom3,0.0))
                        for atom4 in atom3.Bond_List:
                            if atom1 != atom4 and atom2 != atom4:
                                Identical_Dihedral = False
                                for dih in self.Dihedral_List:
                                    if dih.Compare_Dihedrals(atom1,atom2,atom3,atom4):
                                        Identical_Dihedral = True
                                        break
                                if not Identical_Dihedral:
                                    self.Dihedral_List.append(Dihedral.Dihedral(atom2,atom3,atom1,atom4,0.0))
        for ring in self.Ring_List:
            for atom1 in ring.Core_Atom_List:
                for atom2 in atom1.Bond_List:
                    if atom2 in ring.Core_Atom_List:
                        for atom3 in atom2.Bond_List:
                            if atom1 != atom3 and atom3 in ring.Core_Atom_List:
                                for atom4 in atom2.Bond_List:
                                    if atom4 != atom1 and atom4 != atom3 and atom4 not in ring.Core_Atom_List and atom4 in ring.Atom_List:
                                        Identical_Improper = False
                                        for imp in self.Improper_List:
                                            if imp.Compare_Impropers(atom1,atom2,atom3,atom4):
                                                Identical_Improper = True
                                                break
                                        if not Identical_Improper:
                                            self.Improper_List.append(Improper.Improper(atom2,atom3,atom1,atom4,3.5,2,len(self.Improper_List)+1))

    def Find_Coul_Cutoff(self):
        coul_cutoff = 0
        for atom1 in self.Atom_List:
            for atom2 in self.Atom_List:
                if atom1 != atom2 and np.linalg.norm(atom1.Position - atom2.Position) > coul_cutoff:
                    coul_cutoff = np.linalg.norm(atom1.Position - atom2.Position)
        coul_cutoff += 2.5
        return coul_cutoff

    def Update_From_XYZ(self,XYZ_File):
        f = open(XYZ_File,'r')
        lines = f.readlines()
        if len(lines[len(lines)-len(self.Atom_List):]) != len(self.Atom_List):
            raise Exception("Mismatched Lengths")

        for line,atom in zip(lines[len(lines)-len(self.Atom_List):],self.Atom_List):
            if len(line.split()) == 4:
                atom.Position = np.array([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    # def LigParGen_Polymer(self):
    #     self.Parameterize_Interring_Bonds()
    #     self.Simplify_Interactions()

    def Create_Hydrogenated_Copy(self,Ring_Index,Ring_Normal):
        self.Update_All_Normal_Vectors()
        H_Atom_ID_List = []
        B_Atom_H_List = []
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            if not b_atom.Is_Linked:
                B_Atom_H_List.append(b_atom.H_Atom)
                self.Atom_List.remove(b_atom.H_Atom)
                H_Atom_ID_List.append(b_atom.H_Atom)
        i = 0
        for h_atom,c_atom in zip(self.Ring_List[Ring_Index].Hydrogenated_H_List,self.Ring_List[Ring_Index].Hydrogenated_C_List):
            if i < len(H_Atom_ID_List):
                h_atom.Atom_ID = H_Atom_ID_List[i]
                i+=1
            else:
                h_atom.Atom_ID = len(self.Atom_List)
            Linked_Atom_Flag = False
            for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
                if b_atom.Central_Atom == c_atom and b_atom.Is_Linked:
                    Linked_Atom_Flag = True
            if not Linked_Atom_Flag:
                self.Atom_List.append(h_atom)

        Directions = []
        b_atom_h_added = 0
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            if b_atom.Is_Linked:
                try1 = b_atom.Central_Atom.Position - 1.08*self.Ring_List[Ring_Index].Normal_Vector
                distance1 = (np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                try2 = b_atom.Central_Atom.Position + 1.08*self.Ring_List[Ring_Index].Normal_Vector
                distance2 = (np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                if distance1 > distance2:
                    New_H_Atom_Position = try1
                    Directions.append(0)
                else:
                    New_H_Atom_Position = try2
                    Directions.append(1)
                new_h_atom = Atom.Atom(New_H_Atom_Position,'H',len(self.Atom_List))
                self.Atom_List.append(new_h_atom)
                b_atom_h_added += 1
        Hydrogenated_Version = copy.deepcopy(self)

        if len(Directions) == 2:
            if Directions[0] == Directions[1]:
                Syn = 0
            elif self.Ring_List[Ring_Index].Symmetric:
                Syn = 1
            else:
                Syn = 2
        else:
            Syn = 0

        #remove all added hydrogens

        self.Atom_List = self.Atom_List[:-1*b_atom_h_added]
        for h_atom,c_atom in zip(self.Ring_List[Ring_Index].Hydrogenated_H_List,self.Ring_List[Ring_Index].Hydrogenated_C_List):
            Linked_Atom_Flag = False
            for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
                if b_atom.Central_Atom == c_atom and b_atom.Is_Linked:
                    Linked_Atom_Flag = True
            if not Linked_Atom_Flag:
                self.Atom_List.remove(h_atom)
        for h_atom in B_Atom_H_List:
            self.Atom_List.append(h_atom)
        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)

        return Hydrogenated_Version,Syn

    def Create_Hydrogenated_Copy_Alternate(self,Ring_Index,Ring_Normal):
        self.Update_All_Normal_Vectors()
        H_Atom_ID_List = []
        B_Atom_H_List = []
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            if not b_atom.Is_Linked:
                B_Atom_H_List.append(b_atom.H_Atom)
                self.Atom_List.remove(b_atom.H_Atom)
                H_Atom_ID_List.append(b_atom.H_Atom)
        i = 0
        for h_atom,c_atom in zip(self.Ring_List[Ring_Index].Hydrogenated_H_List,self.Ring_List[Ring_Index].Hydrogenated_C_List):
            if i < len(H_Atom_ID_List):
                h_atom.Atom_ID = H_Atom_ID_List[i]
                i+=1
            else:
                h_atom.Atom_ID = len(self.Atom_List)
            Linked_Atom_Flag = False
            for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
                if b_atom.Central_Atom == c_atom and b_atom.Is_Linked:
                    Linked_Atom_Flag = True
            if not Linked_Atom_Flag:
                self.Atom_List.append(h_atom)
        b_atom_h_added = 0
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            if b_atom.Is_Linked:
                try1 = b_atom.Central_Atom.Position - 1.08*self.Ring_List[Ring_Index].Normal_Vector
                distance1 = (np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                try2 = b_atom.Central_Atom.Position + 1.08*self.Ring_List[Ring_Index].Normal_Vector
                distance2 = (np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                if distance1 <= distance2:
                    New_H_Atom_Position = try1
                else:
                    New_H_Atom_Position = try2
                new_h_atom = Atom.Atom(New_H_Atom_Position,'H',len(self.Atom_List))
                self.Atom_List.append(new_h_atom)
                b_atom_h_added += 1
        Hydrogenated_Version = copy.deepcopy(self)

        #remove all added hydrogens

        self.Atom_List = self.Atom_List[:-1*b_atom_h_added]
        for h_atom,c_atom in zip(self.Ring_List[Ring_Index].Hydrogenated_H_List,self.Ring_List[Ring_Index].Hydrogenated_C_List):
            Linked_Atom_Flag = False
            for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
                if b_atom.Central_Atom == c_atom and b_atom.Is_Linked:
                    Linked_Atom_Flag = True
            if not Linked_Atom_Flag:
                self.Atom_List.remove(h_atom)
        for h_atom in B_Atom_H_List:
            self.Atom_List.append(h_atom)
        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)

        return Hydrogenated_Version

    def Create_Dihydrogenated_Monomer(self,Ring_Index):
        self.Update_All_Normal_Vectors()
        Name_List = []

        f = open("%s_Dihydrogenated_Syn.xyz" % self.Ring_List[Ring_Index].Name,'w')
        f.write("%d\n\n" % (len(self.Ring_List[Ring_Index].Atom_List)+10))
        for atom in self.Atom_List:
            if atom in self.Ring_List[Ring_Index].Atom_List:
                f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            H_Position = b_atom.Central_Atom.Position + 1.08*self.Ring_List[Ring_Index].Normal_Vector
            f.write("H\t%f\t%f\t%f\n" % (H_Position[0],H_Position[1],H_Position[2]))
            C_Position = b_atom.Interring_Bond_Atom.Central_Atom.Position
            f.write("C\t%f\t%f\t%f\n" % (C_Position[0],C_Position[1],C_Position[2]))

            Methyl_Bond = -1.08*self.Ring_List[Ring_Index].Normal_Vector
            rotations = [math.radians(0),math.radians(120),math.radians(240)]
            oop_rot = -1*math.radians(21)
            New_X = self.Ring_List[Ring_Index].Normal_Vector
            New_Y = b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector)
            New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
            Transform_Matrix = np.array([New_X,New_Y,New_Z])
            Untransform_Matrix = np.linalg.inv(Transform_Matrix)
            OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
            Methyl_Bond = np.matmul(Untransform_Matrix,np.matmul(OOP_Matrix,np.matmul(Transform_Matrix,Methyl_Bond)))
            for rot in rotations:
                Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                New_H_Atom_Position = C_Position + np.matmul(Untransform_Matrix,np.matmul(Dih_Matrix,np.matmul(Transform_Matrix,Methyl_Bond)))
                f.write("%s\t%f\t%f\t%f\n" % ("H",New_H_Atom_Position[0],New_H_Atom_Position[1],New_H_Atom_Position[2]))
        Name_List.append("%s_Dihydrogenated_Syn.xyz" % self.Ring_List[Ring_Index].Name)
        f.close()

        f = open("%s_Dihydrogenated_Anti.xyz" % self.Ring_List[Ring_Index].Name,'w')
        f.write("%d\n\n" % (len(self.Ring_List[Ring_Index].Atom_List)+10))
        for atom in self.Atom_List:
            if atom in self.Ring_List[Ring_Index].Atom_List:
                f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))

        Flip = 1
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            H_Position = b_atom.Central_Atom.Position + Flip*1.08*self.Ring_List[Ring_Index].Normal_Vector
            f.write("H\t%f\t%f\t%f\n" % (H_Position[0],H_Position[1],H_Position[2]))
            C_Position = b_atom.Interring_Bond_Atom.Central_Atom.Position
            f.write("C\t%f\t%f\t%f\n" % (C_Position[0],C_Position[1],C_Position[2]))

            Methyl_Bond = -1*Flip*1.08*self.Ring_List[Ring_Index].Normal_Vector
            rotations = [math.radians(0),math.radians(120),math.radians(240)]
            oop_rot = Flip*-1*math.radians(21)
            New_X = self.Ring_List[Ring_Index].Normal_Vector
            New_Y = b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector)
            New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
            Transform_Matrix = np.array([New_X,New_Y,New_Z])
            Untransform_Matrix = np.linalg.inv(Transform_Matrix)
            OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
            Methyl_Bond = np.matmul(Untransform_Matrix,np.matmul(OOP_Matrix,np.matmul(Transform_Matrix,Methyl_Bond)))
            for rot in rotations:
                Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                New_H_Atom_Position = C_Position + np.matmul(Untransform_Matrix,np.matmul(Dih_Matrix,np.matmul(Transform_Matrix,Methyl_Bond)))
                f.write("%s\t%f\t%f\t%f\n" % ("H",New_H_Atom_Position[0],New_H_Atom_Position[1],New_H_Atom_Position[2]))
            Flip = Flip*-1
        Name_List.append("%s_Dihydrogenated_Anti.xyz" % self.Ring_List[Ring_Index].Name)

        f.close()

        if self.Ring_List[Ring_Index].Symmetric == False:
            f = open("%s_Dihydrogenated_Anti_2.xyz" % self.Ring_List[Ring_Index].Name,'w')
            f.write("%d\n\n" % (len(self.Ring_List[Ring_Index].Atom_List)+10))
            for atom in self.Atom_List:
                if atom in self.Ring_List[Ring_Index].Atom_List:
                    f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
            Flip = -1
            for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
                H_Position = b_atom.Central_Atom.Position + Flip*1.08*self.Ring_List[Ring_Index].Normal_Vector
                f.write("H\t%f\t%f\t%f\n" % (H_Position[0],H_Position[1],H_Position[2]))
                C_Position = b_atom.Interring_Bond_Atom.Central_Atom.Position
                f.write("C\t%f\t%f\t%f\n" % (C_Position[0],C_Position[1],C_Position[2]))

                Methyl_Bond = -1*Flip*1.08*self.Ring_List[Ring_Index].Normal_Vector
                rotations = [math.radians(0),math.radians(120),math.radians(240)]
                oop_rot = Flip*-1*math.radians(21)
                New_X = self.Ring_List[Ring_Index].Normal_Vector
                New_Y = b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector)
                New_Z = np.cross(New_X,New_Y)/np.linalg.norm(np.cross(New_X,New_Y))
                Transform_Matrix = np.array([New_X,New_Y,New_Z])
                Untransform_Matrix = np.linalg.inv(Transform_Matrix)
                OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
                Methyl_Bond = np.matmul(Untransform_Matrix,np.matmul(OOP_Matrix,np.matmul(Transform_Matrix,Methyl_Bond)))
                for rot in rotations:
                    Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                    New_H_Atom_Position = C_Position + np.matmul(Untransform_Matrix,np.matmul(Dih_Matrix,np.matmul(Transform_Matrix,Methyl_Bond)))
                    f.write("%s\t%f\t%f\t%f\n" % ("H",New_H_Atom_Position[0],New_H_Atom_Position[1],New_H_Atom_Position[2]))
                Flip = Flip*-1
            Name_List.append("%s_Dihydrogenated_Anti_2.xyz" % self.Ring_List[Ring_Index].Name)
        f.close()

        return Name_List

    def Hydrogenated_Improper(self,Ring_Index,Ring_Normal,XYZ_Filename):
        self.Update_All_Normal_Vectors()
        f = open(XYZ_Filename,'w')
        f.write("%d\n\n" % (len(self.Ring_List[Ring_Index].Atom_List) + 5))
        for atom in self.Ring_List[Ring_Index].Atom_List:
            f.write("%s %.4f %.4f %.4f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        for b_atom in self.Ring_List[Ring_Normal].Bonded_Atoms:
            if b_atom.Is_Linked:
                f.write("%s %.4f %.4f %.4f\n" % (b_atom.Central_Atom.Element,b_atom.Central_Atom.Position[0],b_atom.Central_Atom.Position[1],b_atom.Central_Atom.Position[2]))
                try1 = b_atom.Central_Atom.Position - 1.08*self.Ring_List[Ring_Normal].Normal_Vector
                distance1 = (np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                try2 = b_atom.Central_Atom.Position + 1.08*self.Ring_List[Ring_Normal].Normal_Vector
                distance2 = (np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                if distance1 > distance2:
                    Normal_H = try1
                else:
                    Normal_H = try2
                Bonded_H_1 = (b_atom.Same_Ring_Bonded_Atom_List[0].Position - b_atom.Central_Atom.Position)/np.linalg.norm(b_atom.Same_Ring_Bonded_Atom_List[0].Position - b_atom.Central_Atom.Position)*1.08 + b_atom.Central_Atom.Position
                Bonded_H_2 = (b_atom.Same_Ring_Bonded_Atom_List[1].Position - b_atom.Central_Atom.Position)/np.linalg.norm(b_atom.Same_Ring_Bonded_Atom_List[1].Position - b_atom.Central_Atom.Position)*1.08 + b_atom.Central_Atom.Position
                f.write("H %.4f %.4f %.4f\n" % (Normal_H[0],Normal_H[1],Normal_H[2]))
                f.write("H %.4f %.4f %.4f\n" % (Bonded_H_1[0],Bonded_H_1[1],Bonded_H_1[2]))
                f.write("H %.4f %.4f %.4f\n" % (Bonded_H_2[0],Bonded_H_2[1],Bonded_H_2[2]))
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            if not b_atom.Is_Linked:
                f.write("H %.4f %.4f %.4f\n" % (b_atom.H_Atom.Position[0],b_atom.H_Atom.Position[1],b_atom.H_Atom.Position[2]))
        f.close()

    def Hydrogenated_Improper_Alternate(self,Ring_Index,Ring_Normal,XYZ_Filename):
        self.Update_All_Normal_Vectors()
        f = open(XYZ_Filename,'w')
        f.write("%d\n\n" % (len(self.Ring_List[Ring_Index].Atom_List) + 5))
        for atom in self.Ring_List[Ring_Index].Atom_List:
            f.write("%s %.4f %.4f %.4f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        for b_atom in self.Ring_List[Ring_Normal].Bonded_Atoms:
            if b_atom.Is_Linked:
                f.write("%s %.4f %.4f %.4f\n" % (b_atom.Central_Atom.Element,b_atom.Central_Atom.Position[0],b_atom.Central_Atom.Position[1],b_atom.Central_Atom.Position[2]))
                try1 = b_atom.Central_Atom.Position - 1.08*self.Ring_List[Ring_Normal].Normal_Vector
                distance1 = (np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try1 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                try2 = b_atom.Central_Atom.Position + 1.08*self.Ring_List[Ring_Normal].Normal_Vector
                distance2 = (np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Position) + np.linalg.norm(try2 - b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Position))/2
                if distance1 <= distance2:
                    Normal_H = try1
                else:
                    Normal_H = try2
                Bonded_H_1 = (b_atom.Same_Ring_Bonded_Atom_List[0].Position - b_atom.Central_Atom.Position)/np.linalg.norm(b_atom.Same_Ring_Bonded_Atom_List[0].Position - b_atom.Central_Atom.Position)*1.08 + b_atom.Central_Atom.Position
                Bonded_H_2 = (b_atom.Same_Ring_Bonded_Atom_List[1].Position - b_atom.Central_Atom.Position)/np.linalg.norm(b_atom.Same_Ring_Bonded_Atom_List[1].Position - b_atom.Central_Atom.Position)*1.08 + b_atom.Central_Atom.Position
                f.write("H %.4f %.4f %.4f\n" % (Normal_H[0],Normal_H[1],Normal_H[2]))
                f.write("H %.4f %.4f %.4f\n" % (Bonded_H_1[0],Bonded_H_1[1],Bonded_H_1[2]))
                f.write("H %.4f %.4f %.4f\n" % (Bonded_H_2[0],Bonded_H_2[1],Bonded_H_2[2]))
        for b_atom in self.Ring_List[Ring_Index].Bonded_Atoms:
            if not b_atom.Is_Linked:
                f.write("H %.4f %.4f %.4f\n" % (b_atom.H_Atom.Position[0],b_atom.H_Atom.Position[1],b_atom.H_Atom.Position[2]))
        f.close()

    def Zero_Bonded_Impropers(self):
        # replaces every interring bond impropers' Ki with 0
        Bonded_Atom_List = []
        for ring in self.Ring_List:
            for b_atom in ring.Bonded_Atoms:
                if b_atom.Is_Linked:
                    Bonded_Atom_List.append(b_atom)
        for bond in self.Interring_Bond_List:
            atom1_ID = bond.Bond_Master.Atom_ID
            atom2_ID = bond.Bond_Slave.Atom_ID
            for imp in self.Improper_List:
                if (imp.Improper_Master.Atom_ID == atom1_ID and imp.Improper_Slave3.Atom_ID == atom2_ID) or (imp.Improper_Master.Atom_ID == atom2_ID and imp.Improper_Slave3.Atom_ID == atom1_ID):
                    imp.Ki = 0.0

    def Replace_Interring_Dihedrals(self,index1,index2,coeffs_list,zero_coeffs = [0.0,0.0,0.0,0.0],different_dihedrals=1,double_dih=False):
        #replaces one interring bond dihedral with new coeffs (coeffs, standard list) and sets the rest to 0 energy. The specific dihedral to give the coeffs to is specified by index1 and index2, corresponding to the position in the Same_Ring_Bonded_Atom_List for each ring the end atom is
        Bonded_Atom_List = []
        for ring in self.Ring_List:
            for b_atom in ring.Bonded_Atoms:
                if b_atom.Is_Linked:
                    Bonded_Atom_List.append(b_atom)
        n=0
        for bond in self.Interring_Bond_List:
            atom1_ID = bond.Bond_Master.Atom_ID
            atom2_ID = bond.Bond_Slave.Atom_ID
            for dih in self.Dihedral_List:
                if (dih.Dihedral_Master1.Atom_ID == atom1_ID and dih.Dihedral_Master2.Atom_ID == atom2_ID):
                    for b_atom in Bonded_Atom_List:
                        if b_atom.Central_Atom.Atom_ID == atom1_ID:
                            bonded_atom_ID_1 = b_atom.Same_Ring_Bonded_Atom_List[index1[n]].Atom_ID
                            bonded_atom_ID_2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[index2[n]].Atom_ID
                            if double_dih:
                                if index1[n] == 0:
                                    alt_index1 = 1
                                else:
                                    alt_index1 = 0
                                if index2[n] == 0:
                                    alt_index2 = 1
                                else:
                                    alt_index2 = 0
                                alt_bonded_atom_ID_1 = b_atom.Same_Ring_Bonded_Atom_List[alt_index1].Atom_ID
                                alt_bonded_atom_ID_2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[alt_index2].Atom_ID
                            break
                    if double_dih and dih.Dihedral_Slave1.Atom_ID == alt_bonded_atom_ID_1 and dih.Dihedral_Slave2.Atom_ID == alt_bonded_atom_ID_2:
                        dih.Coeffs = coeffs_list[n]
                    if dih.Dihedral_Slave1.Atom_ID == bonded_atom_ID_1 and dih.Dihedral_Slave2.Atom_ID == bonded_atom_ID_2:
                        dih.Coeffs = coeffs_list[n]
                        n+=1
                        if n >= different_dihedrals:
                            n=0
                    else:
                        dih.Coeffs = zero_coeffs
                elif (dih.Dihedral_Master1.Atom_ID == atom2_ID and dih.Dihedral_Master2.Atom_ID == atom1_ID):
                    for b_atom in Bonded_Atom_List:
                        if b_atom.Central_Atom.Atom_ID == atom1_ID:
                            bonded_atom_ID_1 = b_atom.Same_Ring_Bonded_Atom_List[index1[n]].Atom_ID
                            bonded_atom_ID_2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[index2[n]].Atom_ID
                            if double_dih:
                                if index1[n] == 0:
                                    alt_index1 = 1
                                else:
                                    alt_index1 = 0
                                if index2[n] == 0:
                                    alt_index2 = 1
                                else:
                                    alt_index2 = 0
                                alt_bonded_atom_ID_1 = b_atom.Same_Ring_Bonded_Atom_List[alt_index1].Atom_ID
                                alt_bonded_atom_ID_2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[alt_index2].Atom_ID
                            break
                    if dih.Dihedral_Slave2.Atom_ID == bonded_atom_ID_1 and dih.Dihedral_Slave1.Atom_ID == bonded_atom_ID_2:
                        dih.Coeffs = coeffs_list[n]
                        n+=1
                        if n >= different_dihedrals:
                            n=0
                    elif double_dih and dih.Dihedral_Slave2.Atom_ID == alt_bonded_atom_ID_1 and dih.Dihedral_Slave1.Atom_ID == alt_bonded_atom_ID_2:
                        dih.Coeffs = coeffs_list[n]
                    else:
                        dih.Coeffs = zero_coeffs
            #Molecule.Assign_Lammps([self])

    def Calculate_Rg(self):
        Mass_Weighted_Sum = 0.0
        for Atom_Obj in self.Atom_List:
            Mass_Weighted_Sum += Atom_Obj.Mass*((Atom_Obj.Position[0] - self.COM[0])**2 + (Atom_Obj.Position[1] - self.COM[1])**2  + (Atom_Obj.Position[2] - self.COM[2])**2 )
        Rg2 = Mass_Weighted_Sum/self.MW
        self.Rg = np.sqrt(Rg2)
        print "The Radius of Gyration is", self.Rg
        return self.Rg

    def Calculate_PL(self,PL_Atom_ID,Polymer_Name,Make_Figure=False):
        Vectors = []
        Ps = []
        atoms = []
        for ring1,ring2 in zip(self.Ring_List[:-1],self.Ring_List[1:]):
            atoms.append(ring1.Atom_List[PL_Atom_ID].Atom_ID)
            Vectors.append((ring1.Atom_List[PL_Atom_ID].Position - ring2.Get_Atom(PL_Atom_ID).Position)/np.linalg.norm(ring1.Atom_List[PL_Atom_ID].Position - ring2.Get_Atom(PL_Atom_ID).Position))
        for i,v_1 in enumerate(Vectors):
            Cos_Thetas = []
            for v_2 in Vectors[i:]:
                Cos_Thetas.append(np.dot(v_1,v_2))
            Ls = range(len(Vectors) - i)
            for v_2 in Vectors[i::-1]:
                Cos_Thetas.append(np.dot(v_1,v_2))
            Ls = Ls + range(i + 1)
            Ps.append(scipy.optimize.curve_fit(Molecule.exponential,np.array(Ls),np.array(Cos_Thetas))[0][0])

            if Make_Figure:
                fig,ax = plt.subplots(1,1)
                plt.scatter(Ls,Cos_Thetas,marker = 's',c='k')
                x = np.linspace(0,15)
                fit = [Molecule.exponential(l,Ps[-1]) for l in x]
                plt.plot(x,fit,color='r')
                plt.tight_layout()
                fig.savefig('%s_%d_Persistence_Length' % (Polymer_Name,i))
                plt.close(fig)
                os.system("mkdir ./Figures/PL_Figures")
                os.system("scp %s_%d_Persistence_Length.png ./Figures/PL_Figures" % (Polymer_Name,i))
                os.system("rm -f %s_%d_Persistence_Length.png" % (Polymer_Name,i))
        return np.mean(np.array(Ps)),np.std(Ps)

    def CVFF(self,x,K,d,n):
        return K*(1+d*np.cos(np.pi/180*n*x))

    def Randomize_Torsions(self,Delocalization_Energies,Nonbonded_Energies,C_Probs=""):
        #Delocalization_Energies and Nonbonded_Energies need to be energies in kcal/mol for dihedral and OOP in 5 degree implements from 0 to 355 in first dimension and 0 to 60 in OOP
        OOP_1 = np.linspace(-15.0,15.0,len(Delocalization_Energies))
        OOP_2 = np.linspace(-15.0,15.0,len(Delocalization_Energies))
        Dih = np.linspace(0,360 - 360/len(Delocalization_Energies[0]),len(Delocalization_Energies[0]))

        Probabilities = np.zeros((len(self.Ring_List)-1,len(OOP_1),len(OOP_2),len(Dih)))
        Total_Energies = np.zeros((len(self.Ring_List)-1,len(OOP_1),len(OOP_2),len(Dih)))
        Cumulative_Probabilities = np.zeros((len(self.Ring_List)-1,len(OOP_1),len(OOP_2),len(Dih)))
        if C_Probs == "":
            for b_atom in self.Ring_List[0].Bonded_Atoms:
                if b_atom.Is_Linked:
                    b_atom_1 = b_atom
                    b_atom_2 = b_atom.Interring_Bond_Atom
                    for i,oop1 in enumerate(OOP_1):
                        for j,oop2 in enumerate(OOP_2):
                            for k,d in enumerate(Dih):
                                oop1_energy = self.CVFF(oop1,b_atom_1.K,b_atom_1.d,b_atom_1.n)
                                oop2_energy = self.CVFF(oop2,b_atom_2.K,b_atom_2.d,b_atom_2.n)
                                oop_index = (abs(oop1) + abs(oop2))/5
                                del_e = Delocalization_Energies[int(oop_index)][k]
                                non_e = Nonbonded_Energies[int(oop_index)][k]
                                Total_Energies[0][i][j][k] = oop1_energy + oop2_energy + del_e + non_e
                                Probabilities[0][i][j][k] = math.exp(-(oop1_energy + oop2_energy + del_e + non_e)/(.001986*300))
                    Probabilities[0] = Probabilities[0]/np.sum(Probabilities[0])
                    c_prob = 0
                    for i,oop1 in enumerate(OOP_1):
                        for j,oop2 in enumerate(OOP_2):
                            for k,d in enumerate(Dih):
                                c_prob += Probabilities[0][i][j][k]
                                Cumulative_Probabilities[0][i][j][k] = c_prob

            for q,ring in enumerate(self.Ring_List[1:-1]):
                for b_atom in ring.Bonded_Atoms:
                    if b_atom != b_atom_2:
                        b_atom_1 = b_atom
                        b_atom_2 = b_atom.Interring_Bond_Atom
                        for i,oop1 in enumerate(OOP_1):
                            for j,oop2 in enumerate(OOP_2):
                                for k,d in enumerate(Dih):
                                    oop1_energy = self.CVFF(oop1,b_atom_1.K,b_atom_1.d,b_atom_1.n)
                                    oop2_energy = self.CVFF(oop2,b_atom_2.K,b_atom_2.d,b_atom_2.n)
                                    oop_index = (abs(oop1) + abs(oop2))/5
                                    del_e = Delocalization_Energies[int(oop_index)][k]
                                    non_e = Nonbonded_Energies[int(oop_index)][k]
                                    Probabilities[q+1][i][j][k] = math.exp(-(oop1_energy + oop2_energy + del_e + non_e)/(.001986*300))
                        Probabilities[q+1] = Probabilities[q+1]/np.sum(Probabilities[q+1])
                        c_prob = 0
                        for i,oop1 in enumerate(OOP_1):
                            for j,oop2 in enumerate(OOP_2):
                                for k,d in enumerate(Dih):
                                    c_prob += Probabilities[q+1][i][j][k]
                                    Cumulative_Probabilities[q+1][i][j][k] = c_prob
                        break

        else:
            Cumulative_Probabilities = C_Probs

        for q,c_probs in enumerate(Cumulative_Probabilities):
            Chosen = False
            random.seed()
            r = random.random()
            for i,oop1 in enumerate(OOP_1):
                for j,oop2 in enumerate(OOP_2):
                    for k,d in enumerate(Dih):
                        if r < c_probs[i][j][k]:
                            self.Rotate_Ring("Dih",d,self.Ring_List[q],self.Ring_List[q+1])
                            self.Rotate_Ring("OOP",oop1,self.Ring_List[q],self.Ring_List[q+1])
                            self.Rotate_Ring("OOP",oop2,self.Ring_List[q+1],self.Ring_List[q])
                            Chosen = True
                            break
                    if Chosen:
                        break
                if Chosen:
                    break
        return Cumulative_Probabilities

    def Parameterize_Interring_Bonds(self):
        #Creates interring bonds, angles, and dihedrals after LigParGen read in
        for b_bond in self.Interring_Bond_List:
            self.Bond_List.append(Bond.Bond(b_bond.Bond_Master, b_bond.Bond_Slave, b_bond.req))
            self.Bond_List[-1].kb = 469.0
            #b_bond.req = 1.424
        #Create interring angles
        for ring in self.Ring_List:
            for b_atom in ring.Bonded_Atoms:
                if b_atom.Is_Linked:
                    if b_atom.Same_Ring_Bonded_Atom_List[0].Element == "S":
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[0],111.17))
                        self.Angle_List[-1].ka = 70.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                    else:
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[0], 107.300))
                        self.Angle_List[-1].ka = 70.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                    if b_atom.Same_Ring_Bonded_Atom_List[1].Element == "S":
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[1],111.17))
                        self.Angle_List[-1].ka = 70.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                    else:
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[1], 107.300))
                        self.Angle_List[-1].ka = 70.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                else:
                    self.Bond_List.append(Bond.Bond(b_atom.Central_Atom,b_atom.H_Atom,1.08))
                    self.Bond_List[-1].kb = 367.0
                    self.Bond_List[-1].Bond_ID = len(self.Bond_List)
                    b_atom.H_Atom.Sigma = 2.42
                    b_atom.H_Atom.Epsilon = 0.03
                    b_atom.H_Atom.Charge = .2
                    if b_atom.Same_Ring_Bonded_Atom_List[0].Element == "S":
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.H_Atom,b_atom.Same_Ring_Bonded_Atom_List[0],125.0))
                        self.Angle_List[-1].ka = 35.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                    else:
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.H_Atom,b_atom.Same_Ring_Bonded_Atom_List[0], 132.100))
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                        self.Angle_List[-1].ka = 35.0
                    if b_atom.Same_Ring_Bonded_Atom_List[1].Element == "S":
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.H_Atom,b_atom.Same_Ring_Bonded_Atom_List[1],125.0))
                        self.Angle_List[-1].ka = 35.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
                    else:
                        self.Angle_List.append(Angle.Angle(b_atom.Central_Atom,b_atom.H_Atom,b_atom.Same_Ring_Bonded_Atom_List[1], 132.100))
                        self.Angle_List[-1].ka = 35.0
                        self.Angle_List[-1].Angle_ID = len(self.Angle_List)
        #Create interring dihedrals
        for ring in self.Ring_List:
            b_atom = ring.Bonded_Atoms[1]
            if b_atom.Is_Linked:
                self.Dihedral_List.append(Dihedral.Dihedral(b_atom.Central_Atom, b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[0],b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0], 0.0))
                self.Dihedral_List[-1].Coeffs = [0.0,0.0,0.0,0.0]
                self.Dihedral_List[-1].Dihedral_ID = len(self.Dihedral_List)
                self.Dihedral_List.append(Dihedral.Dihedral(b_atom.Central_Atom, b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[0],b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1], 0.0))
                self.Dihedral_List[-1].Coeffs = [0.0, 0.0, 0.0, 0.0]
                self.Dihedral_List[-1].Dihedral_ID = len(self.Dihedral_List)
                self.Dihedral_List.append(Dihedral.Dihedral(b_atom.Central_Atom, b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[1], b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0], 0.0))
                self.Dihedral_List[-1].Coeffs = [0.0, 0.0, 0.0, 0.0]
                self.Dihedral_List[-1].Dihedral_ID = len(self.Dihedral_List)
                self.Dihedral_List.append(Dihedral.Dihedral(b_atom.Central_Atom, b_atom.Interring_Bond_Atom.Central_Atom,b_atom.Same_Ring_Bonded_Atom_List[1],b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1], 0.0))
                self.Dihedral_List[-1].Coeffs = [0.0, 0.0, 0.0, 0.0]
                self.Dihedral_List[-1].Dihedral_ID = len(self.Dihedral_List)

    def Read_From_Data_File(self,Input_File,No_Position_Update = False):
        Element_Dict = { 12.011:"C", 1.008:"H", 18.998:"F", 15.999:"O", 32.06:"S", 14.007:"N", 28.085:"Si",28.086:"Si", 35.453:"Cl"}
        Full = False
        j = -1
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        f = open(Input_File,'r')
        File_Lines = f.readlines()
        for Line in File_Lines:
            j += 1
            Line = Line.strip().split()
            if len(Line) > 1:
                if Line[1] == 'atoms':
                    self.N = int(Line[0])
                    #print self.N, "Atoms"
                if Line[1] == 'atom':
                    Atom_Types = int(Line[0])
                    Mass_List = np.empty(Atom_Types, dtype=float)
                    Pair_Coeffs = np.empty((Atom_Types,2), dtype=float)
                if Line[1] == 'bonds':
                    Num_Bonds = int(Line[0])
                if Line[1] == 'bond':
                    Bond_Types = int(Line[0])
                    Bond_Coeffs = np.empty((Bond_Types,2), dtype=float)
                if Line[1] == 'angles':
                    Num_Angles = int(Line[0])
                if Line[1] == 'angle':
                    Angle_Types = int(Line[0])
                    Angle_Coeffs = np.empty((Angle_Types,2),dtype=float)
                if Line[1] == 'dihedrals':
                    Num_Dihedrals = Line[0]
                if Line[1] == 'dihedral':
                    Dihedral_Types = int(Line[0])
                    #print('dihedral types')
                    #print (Dihedral_Types)
                    #Dihedral_Coeffs = [np.empty((Dihedral_Types,9),dtype=float)]
                    Dihedral_Coeffs = []
                    Dihedral_Styles = []
                if Line[1] == 'impropers':
                    Num_Impropers = Line[0]
                if Line[1] == 'improper':
                    Improper_Types = int(Line[0])
                    Improper_Coeffs = np.empty((Improper_Types,3),dtype=float)
                #self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
                #for Atom_Obj in self.Atom_List:
                    #print Atom_Obj.Atom_ID, Atom_Obj.Position
                if Line[0] == 'Pair':
                    #print "Getting Pair Coefficients"
                    for i in range(Atom_Types):
                        Pair_Coeffs[i,0] = float(File_Lines[i+j+2].split()[1])
                        Pair_Coeffs[i,1] = float(File_Lines[i+j+2].split()[2])
                if Line[0] == 'Bond':
                    #print "Getting Bond Coefficients"
                    for i in range(Bond_Types):
                        Bond_Coeffs[i,0] = float(File_Lines[i+j+2].split()[1])
                        Bond_Coeffs[i,1] = float(File_Lines[i+j+2].split()[2])
                if Line[0] == 'Angle':
                    #print "Getting Angle Coefficients"
                    for i in range(Angle_Types):
                        Angle_Coeffs[i,0] = float(File_Lines[i+j+2].split()[1])
                        Angle_Coeffs[i,1] = float(File_Lines[i+j+2].split()[2])
                if Line[0].strip() == 'Dihedral':
                    for i in range(Dihedral_Types):
                        dummy_coeffs = []
                        for k in range(1,len(File_Lines[i+j+2].split())):
                            dummy_coeffs.append(float(File_Lines[i+j+2].split()[k]))
                        Dihedral_Coeffs.append(np.array(dummy_coeffs))
                        Dihedral_Styles.append(File_Lines[i+j+2].split()[0])
                            #Dihedral_Coeffs[i,k-2] = float(File_Lines[i+j+2].split(' ')[k])
                        
                        """Dihedral_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                        Dihedral_Coeffs[i,2] = float(File_Lines[i+j+2].split(' ')[3])
                        Dihedral_Coeffs[i,3] = float(File_Lines[i+j+2].split(' ')[4])
                        Dihedral_Coeffs[i,4] = float(File_Lines[i+j+2].split(' ')[5])
                        Dihedral_Coeffs[i,5] = float(File_Lines[i+j+2].split(' ')[6])
                        Dihedral_Coeffs[i,6] = float(File_Lines[i+j+2].split(' ')[7])
                        Dihedral_Coeffs[i,7] = float(File_Lines[i+j+2].split(' ')[8])
                        Dihedral_Coeffs[i,8] = float(File_Lines[i+j+2].split(' ')[9])"""
                        
                if Line[0] == 'Improper':
                    #print "Getting Improper Coefficients"
                    for i in range(Improper_Types):
                        Improper_Coeffs[i,0] = float(File_Lines[i+j+2].split()[1])
                        Improper_Coeffs[i,1] = int(File_Lines[i+j+2].split()[2])
                        Improper_Coeffs[i,2] = int(File_Lines[i+j+2].split()[3])



            if len(Line) == 1 or (len(Line) > 1 and Line[1] == "#"):
                if Line == ['']:
                    continue
                if Line[0] == 'Atoms':
                    for i in range(self.N):
                        Type = int(File_Lines[i+j+2].split()[2])
                        Mass = Mass_List[Type-1]
                        Element = Element_Dict[Mass]
                        atom = self.Get_Atom(int(File_Lines[i+j+2].split()[0]))
                        """self.Atom_List[i].Atom_ID = int(File_Lines[i+j+2].split()[0])
                        self.Atom_List[i].Epsilon = Pair_Coeffs[Type-1,0]
                        self.Atom_List[i].Sigma = Pair_Coeffs[Type-1,1]
                        self.Atom_List[i].Charge = float(File_Lines[i+j+2].split()[3])
                        if not No_Position_Update:
                            self.Atom_List[i].Position = np.array([float(File_Lines[i+j+2].split()[4]),float(File_Lines[i+j+2].split()[5]),float(File_Lines[i+j+2].split()[6])])
                        if Full:
                            self.Atom_List[i].OPLS_Type = 2000 + Type
                            self.Atom_List[i].LAMMPS_Type = 2000 + Type
                        else:
                            self.Atom_List[i].OPLS_Type = Type
                            self.Atom_List[i].LAMMPS_Type = Type"""
                        atom.Atom_ID = int(File_Lines[i+j+2].split()[0])
                        atom.Epsilon = Pair_Coeffs[Type-1,0]
                        atom.Sigma = Pair_Coeffs[Type-1,1]
                        atom.Charge = float(File_Lines[i+j+2].split()[3])
                        if not No_Position_Update:
                            atom.Position = np.array([float(File_Lines[i+j+2].split()[4]),float(File_Lines[i+j+2].split()[5]),float(File_Lines[i+j+2].split()[6])])
                        if Full:
                            atom.OPLS_Type = 2000 + Type
                            atom.LAMMPS_Type = 2000 + Type
                        else:
                            atom.OPLS_Type = Type
                            atom.LAMMPS_Type = Type
                if Line == ['Masses']:
                    #print "Extracting Masses"
                    for i in range(Atom_Types):
                        Mass_List[i] = float(File_Lines[i+j+2].split()[1])
                if Line == ['Bonds']:
                    #print "Extracting Bonds"
                    for i in range(int(Num_Bonds)):
                        Bond_Info = File_Lines[i+j+2].strip().split()
                        #Master = self.Atom_List[int(Bond_Info[2])-1]
                        Master = self.Get_Atom(int(Bond_Info[2]))
                        #Slave = self.Atom_List[int(Bond_Info[3])-1]
                        Slave = self.Get_Atom(int(Bond_Info[3]))
                        Kb = Bond_Coeffs[int(Bond_Info[1])-1, 0]
                        Req = Bond_Coeffs[int(Bond_Info[1])-1,1]
                        Bond_ID = int(Bond_Info[0])
                        self.Bond_List.append(Bond.Bond(Master, Slave,Req))
                        self.Bond_List[i].kb = Kb
                        self.Bond_List[i].Bond_ID = Bond_ID
                        self.Bond_List[i].LAMMPS_Type = Bond_ID = int(Bond_Info[1])
                if Line == ['Angles']:
                    #print "Extracting Angles"
                    for i in range(int(Num_Angles)):
                        Angle_Info = File_Lines[i+j+2].strip('\n').split()
                        #Slave1 = self.Atom_List[int(Angle_Info[2])-1]
                        Slave1 = self.Get_Atom(int(Angle_Info[2]))
                        #Master = self.Atom_List[int(Angle_Info[3])-1]
                        Master = self.Get_Atom(int(Angle_Info[3]))
                        #Slave2 = self.Atom_List[int(Angle_Info[4])-1]
                        Slave2 = self.Get_Atom(int(Angle_Info[4]))
                        Ka = Angle_Coeffs[int(Angle_Info[1])-1,0]
                        Th0 = Angle_Coeffs[int(Angle_Info[1])-1,1]
                        Angle_ID = int(Angle_Info[0])
                        self.Angle_List.append(Angle.Angle(Master, Slave1, Slave2, Th0))
                        self.Angle_List[i].ka = Ka
                        self.Angle_List[i].Angle_ID = Angle_ID
                        self.Angle_List[i].LAMMPS_Type = int(Angle_Info[1])
                if Line == ['Dihedrals']:
                    #print "Extracting Dihedrals"
                    for i in range(int(Num_Dihedrals)):
                        Dihedral_Info = File_Lines[i+j+2].strip('\n').split()
                        #Slave1 = self.Atom_List[int(Dihedral_Info[2])-1]
                        Slave1 = self.Get_Atom(int(Dihedral_Info[2]))
                        #Master1 = self.Atom_List[int(Dihedral_Info[3])-1]
                        Master1 = self.Get_Atom(int(Dihedral_Info[3]))
                        #Master2 = self.Atom_List[int(Dihedral_Info[4])-1]
                        Master2 = self.Get_Atom(int(Dihedral_Info[4]))
                        #Slave2 = self.Atom_List[int(Dihedral_Info[5])-1]
                        Slave2 = self.Get_Atom(int(Dihedral_Info[5]))
                        Coeffs = Dihedral_Coeffs[int(Dihedral_Info[1])-1]
                        #Coeffs = []
                        #Style = Dihedral_Styles[int(Dihedral_Info[1])-1]
                        Style = int(Dihedral_Info[1])
                        Dihedral_ID = int(Dihedral_Info[0])
                        self.Dihedral_List.append(Dihedral.Dihedral(Master1, Master2, Slave1, Slave2, 0.0))
                        self.Dihedral_List[i].Coeffs = Coeffs
                        self.Dihedral_List[i].Dihedral_ID = Dihedral_ID
                        self.Dihedral_List[i].Style = Style
                        self.Dihedral_List[i].LAMMPS_Type = int(Dihedral_Info[1])

                if Line == ['Impropers']:
                    #print "Extracting Impropers"
                    for i in range(int(Num_Impropers)):
                        Improper_Info = File_Lines[i+j+2].strip('\n').split()
                        #Master = self.Atom_List[int(Improper_Info[2])-1]
                        Master = self.Get_Atom(int(Improper_Info[2]))
                        #Slave1 = self.Atom_List[int(Improper_Info[3])-1]
                        Slave1 = self.Get_Atom(int(Improper_Info[3]))
                        #Slave2 = self.Atom_List[int(Improper_Info[4])-1]
                        Slave2 = self.Get_Atom(int(Improper_Info[4]))
                        #Slave3 = self.Atom_List[int(Improper_Info[5])-1]
                        Slave3 = self.Get_Atom(int(Improper_Info[5]))
                        Coeff = Improper_Coeffs[int(Improper_Info[1])-1,0]
                        Improper_ID = int(Improper_Info[0])
                        self.Improper_List.append(Improper.Improper(Master, Slave1, Slave2, Slave3, Coeff, 180.0, Improper_ID))
                        self.Improper_List[i].LAMMPS_Type = int(Improper_Info[1])

        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
                        
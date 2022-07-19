# Open relevant modules

from logging import raiseExceptions
import Bond
import numpy as np
import copy
import Molecule
import os
import subprocess
import Write_Inputs
import Write_Submit_Script
import Bonded_Atom
import Cluster_IO
import Aux_Ring
import math
import Conjugated_Polymer
import time
import Configure

class Ring(Molecule.Molecule):
   
    def __init__(self,Supermolecule,Atom_Numbers,Core_Atom_List,Bond_Atoms,Bond_Atom_Vectors,Name,ID,Polymer_Name,Nickname,Symmetric = True,Cluster_Location=Configure.orca_dict["Base_Cluster_Location"]+'/Optimized_Monomers',Scheduler = "Torque",Aux_Ring_List = [],Plumed_Rings = []):
        self.Atom_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.Core_Atom_List = []
        self.Name = Name
        self.Ring_ID = ID
        self.Symmetric = Symmetric
        self.Polymer_Name = Polymer_Name
        self.Nickname = Nickname
        self.Aux_Ring_List = Aux_Ring_List
        self.Plumed_Rings = [[],[]]
        for atom_num in Atom_Numbers:
            self.Atom_List.append(copy.deepcopy(Supermolecule.Get_Atom(atom_num)))
            if atom_num in Core_Atom_List:
                self.Core_Atom_List.append(self.Atom_List[-1])
        for i,plumed_ring in enumerate(Plumed_Rings):
            if len(Plumed_Rings) == 2:
                for atom_num in plumed_ring:
                    for atom in self.Atom_List:
                        if atom_num == atom.Atom_ID:
                            self.Plumed_Rings[i].append(atom)
        for aux_ring in Aux_Ring_List:
            self.Aux_Ring_List.append(Aux_Ring.Aux_Ring_List(Supermolecule,aux_ring[0],aux_ring[1],aux_ring[2],aux_ring[3],self,aux_ring[4],Polymer_Name,aux_ring[5]))

        #Calculate Normal Vector using non-hydrogen atoms

        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2:]:
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)

        #Link hydrogens to Bond_Atoms
        self.Bonded_Atoms = np.empty(len(Bond_Atoms),dtype = object)
        max_atom_id = len(Supermolecule.Atom_List)
        for i,atoms,bond_atom_vec in zip(list(range(len(Bond_Atoms))),Bond_Atoms,Bond_Atom_Vectors):
            Same_Ring_Atom_List = atoms[2:]
            self.Bonded_Atoms[i] = Bonded_Atom.Bonded_Atom(self.Get_Atom(atoms[0]),bond_atom_vec,atoms[1],(max_atom_id + i + 1),self,Same_Ring_Atom_List)

        new_atom_id = 1
        for atom in self.Atom_List:
            atom.Atom_ID = new_atom_id 
            new_atom_id += 1
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Atom_ID = new_atom_id
            new_atom_id += 1

        #Optimize geometry
        self.Optimize_H_Positions(Cluster_Location,Shared_File_Location = Configure.orca_dict["Shared_File_Location"])
        self.Optimize_H_Positions_Hydrogenated(Cluster_Location,Shared_File_Location = Configure.orca_dict["Shared_File_Location"])


    def Add_Plumed_Rings(self,Plumed_Rings):
        self.Plumed_Rings = Plumed_Rings

    def Link_Plumed_Rings(self):
        for i,p_atom in enumerate(self.Plumed_Rings[0]):
            for atom in self.Atom_List:
                if p_atom.Atom_ID == atom.Atom_ID:
                    self.Plumed_Rings[0][i] = atom
        for j,p_atom_2 in enumerate(self.Plumed_Rings[1]):
            for atom in self.Atom_List:
                if p_atom_2.Atom_ID == atom.Atom_ID:
                    self.Plumed_Rings[1][j] = atom
        print("meh")

    def Find_Coul_Cutoff(self):
        coul_cutoff = 0
        for atom1 in self.Atom_List:
            for atom2 in self.Atom_List:
                if atom1 != atom2 and np.linalg.norm(atom1.Position - atom2.Position) > coul_cutoff:
                    coul_cutoff = np.linalg.norm(atom1.Position - atom2.Position)
        coul_cutoff += 2.5
        return coul_cutoff

    def Update_Positions_Data_File(self,Updated_Data_File):
        f = open(Updated_Data_File,'r')
        print(Updated_Data_File)
        lines = f.readlines()
        bond_atom_cutoff = (-1*len(self.Bonded_Atoms))
        regular_atom_lines = lines[2:bond_atom_cutoff]
        bonded_atom_lines = lines[bond_atom_cutoff:]
        if len(self.Atom_List) != len(regular_atom_lines):
            print((len(self.Atom_List)))
            print((len(regular_atom_lines)))
            raise Exception("Length mismatch in regular atoms and xyz file")
        for atom,line in zip(self.Atom_List,regular_atom_lines):
            atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
        if len(self.Bonded_Atoms) != len(bonded_atom_lines):
            print((len(self.Atom_List)))
            print((len(regular_atom_lines)))
            raise Exception("Length mismatch in bonded atoms and xyz file")
        for b_atom,line in zip(self.Bonded_Atoms,bonded_atom_lines):
            b_atom.H_Atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
            b_atom.H_Bond_Vector = b_atom.H_Atom.Position - b_atom.Central_Atom.Position

    def Update_Positions_Orca_Output(self,Orca_Output_File,Hydrogenated=False):
        f = open(Orca_Output_File,'r')
        lines = f.readlines()
        read_atoms = False
        off_count = 0
        for line in lines:
            if read_atoms and len(line.strip().split()) == 4:
                atom_list.append(line)
            elif read_atoms:
                off_count +=1
            if off_count >= 2:
                read_atoms = False
            if len(line.strip().split()) >= 3 and line.strip().split()[0] == "CARTESIAN" and line.strip().split()[1] == "COORDINATES" and line.strip().split()[2] == "(ANGSTROEM)":
                read_atoms = True
                off_count = 0
                atom_list = []
        if not Hydrogenated:
            bond_atom_cutoff = (-1*len(self.Bonded_Atoms))
            regular_atom_lines = atom_list[:bond_atom_cutoff]
            bonded_atom_lines = atom_list[bond_atom_cutoff:]
        else:
            regular_atom_lines = atom_list
            bonded_atom_lines = []
        if len(self.Atom_List) != len(regular_atom_lines) and not Hydrogenated:
            print((len(self.Atom_List)))
            print((len(regular_atom_lines)))
            raise Exception("Length mismatch in regular atoms and xyz file")
        for atom,line in zip(self.Atom_List,regular_atom_lines):
            atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
        if len(self.Bonded_Atoms) != len(bonded_atom_lines) and not Hydrogenated:
            raise Exception("Length mismatch in bonded atoms and xyz file")
        if not Hydrogenated:
            for b_atom,line in zip(self.Bonded_Atoms,bonded_atom_lines):
                b_atom.H_Atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
                b_atom.H_Bond_Vector = b_atom.H_Atom.Position - b_atom.Central_Atom.Position
        self.Update_Normal_Vector()

    def Read_In_LigParGen(self,index):
        temp_ring = copy.deepcopy(self)
        temp_polymer = Conjugated_Polymer.Conjugated_Polymer([copy.deepcopy(self)])
        temp_polymer.Read_From_Data_File("./LigParGen_Files/%s_%d.lmp" % (self.Name,index),No_Position_Update=True)
        self.Atom_List = copy.deepcopy(temp_polymer.Atom_List)
        self.Bond_List = copy.deepcopy(temp_polymer.Bond_List)
        self.Angle_List = copy.deepcopy(temp_polymer.Angle_List)
        self.Dihedral_List = copy.deepcopy(temp_polymer.Dihedral_List)
        self.Improper_List = copy.deepcopy(temp_polymer.Improper_List)
        Atom_List_Length = len(temp_ring.Atom_List)
        # remove interactions with removed hydrogens
        for atom in self.Atom_List[::-1]:
            if atom.Atom_ID > Atom_List_Length:
                self.Atom_List.remove(atom)
        for bond in self.Bond_List[::-1]:
            if bond.Bond_Main.Atom_ID > Atom_List_Length or bond.Bond_Node.Atom_ID > Atom_List_Length:
                self.Bond_List.remove(bond)
        for angle in self.Angle_List[::-1]:
            if angle.Angle_Main.Atom_ID > Atom_List_Length or angle.Angle_Node1.Atom_ID > Atom_List_Length or angle.Angle_Node2.Atom_ID > Atom_List_Length:
                self.Angle_List.remove(angle)
        for dih in self.Dihedral_List[::-1]:
            if dih.Dihedral_Main1.Atom_ID > Atom_List_Length or dih.Dihedral_Main2.Atom_ID > Atom_List_Length or dih.Dihedral_Node1.Atom_ID > Atom_List_Length or dih.Dihedral_Node2.Atom_ID > Atom_List_Length:
                self.Dihedral_List.remove(dih)
        for imp in self.Improper_List[::-1]:
            if imp.Improper_Main.Atom_ID > Atom_List_Length or imp.Improper_Node3.Atom_ID > Atom_List_Length or imp.Improper_Node1.Atom_ID > Atom_List_Length or imp.Improper_Node2.Atom_ID > Atom_List_Length:
                self.Improper_List.remove(imp)

    def Update_Positions_Orca_Output_Hydrogenated(self,Orca_Output_File):
        f = open(Orca_Output_File,'r')
        lines = f.readlines()
        read_atoms = False
        off_count = 0
        for line in lines:
            if read_atoms and len(line.strip().split()) == 4:
                atom_list.append(line.strip().split())
            elif read_atoms:
                off_count +=1
            if off_count >= 2:
                read_atoms = False
            if len(line.strip().split()) >= 3 and line.strip().split()[0] == "CARTESIAN" and line.strip().split()[1] == "COORDINATES" and line.strip().split()[2] == "(ANGSTROEM)":
                read_atoms = True
                off_count = 0
                atom_list = []

        for h_atom,c_atom in zip(self.Hydrogenated_H_List,self.Hydrogenated_C_List):
            h_atom_new_position = np.array([float(atom_list[h_atom.Atom_ID-1][1]),float(atom_list[h_atom.Atom_ID-1][2]),float(atom_list[h_atom.Atom_ID-1][3])])
            c_atom_new_position = np.array([float(atom_list[c_atom.Atom_ID-1][1]),float(atom_list[c_atom.Atom_ID-1][2]),float(atom_list[c_atom.Atom_ID-1][3])])
            h_atom.Position = c_atom.Position + (h_atom_new_position - c_atom_new_position)

    def Relink(self):
        for bond in self.Bond_List:
            bond.Bond_Main = self.Get_Atom(bond.Bond_Main.Atom_ID)
            bond.Bond_Node = self.Get_Atom(bond.Bond_Node.Atom_ID)
        for angle in self.Angle_List:
            angle.Angle_Main = self.Get_Atom(angle.Angle_Main.Atom_ID)
            angle.Angle_Node1 = self.Get_Atom(angle.Angle_Node1.Atom_ID)
            angle.Angle_Node2 = self.Get_Atom(angle.Angle_Node2.Atom_ID)
        for dih in self.Dihedral_List:
            dih.Dihedral_Main1 = self.Get_Atom(dih.Dihedral_Main1.Atom_ID)
            dih.Dihedral_Main2 = self.Get_Atom(dih.Dihedral_Main2.Atom_ID)
            dih.Dihedral_Node1 = self.Get_Atom(dih.Dihedral_Node1.Atom_ID)
            dih.Dihedral_Node2 = self.Get_Atom(dih.Dihedral_Node2.Atom_ID)
        for imp in self.Improper_List:
            imp.Improper_Main = self.Get_Atom(imp.Improper_Main.Atom_ID)
            imp.Improper_Node1 = self.Get_Atom(imp.Improper_Node1.Atom_ID)
            imp.Improper_Node2 = self.Get_Atom(imp.Improper_Node2.Atom_ID)
            imp.Improper_Node3 = self.Get_Atom(imp.Improper_Node3.Atom_ID)
        for b_atom in self.Bonded_Atoms:
            print((len(b_atom.Same_Ring_Bonded_Atom_List)))
            b_atom.Central_Atom = self.Get_Atom(b_atom.Central_Atom.Atom_ID)
            for i in range(len(b_atom.Same_Ring_Bonded_Atom_List)):
                b_atom.Same_Ring_Bonded_Atom_List[i] = self.Get_Atom(b_atom.Same_Ring_Bonded_Atom_List[i].Atom_ID)
                #print("Relinked b atoms")

    def Update_Positions_Data_File(self,Updated_Data_File):
        f = open(Updated_Data_File,'r')
        lines = f.readlines()
        bond_atom_cutoff = (-1*len(self.Bonded_Atoms))
        regular_atom_lines = lines[2:bond_atom_cutoff]
        bonded_atom_lines = lines[bond_atom_cutoff:]
        if len(self.Atom_List) != len(regular_atom_lines):
            raise Exception("Length mismatch in regular atoms and xyz file")
        for atom,line in zip(self.Atom_List,regular_atom_lines):
            atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
        if len(self.Bonded_Atoms) != len(bonded_atom_lines):
            raise Exception("Length mismatch in bonded atoms and xyz file")
        for b_atom,line in zip(self.Bonded_Atoms,bonded_atom_lines):
            b_atom.H_Atom.Position = np.array([float(line.strip().split()[1]),float(line.strip().split()[2]),float(line.strip().split()[3])])
            b_atom.H_Bond_Vector = b_atom.H_Atom.Position - b_atom.Central_Atom.Position
        self.Update_Normal_Vector()

    def Set_Up_Bonds(self,Orca_Output_File):
        f = open(Orca_Output_File,'r')
        lines = f.readlines()
        read_bonds = False
        for line in lines:
            if read_bonds and len(line.strip().split()) >= 1:
                for bond in line.strip().split("B("):
                    if len(bond.strip().split()) >= 3 and float(bond.strip().split()[-1]) > 0.5:
                        Main_ID = int(bond.strip().split()[0].split('-')[0]) + 1
                        if any(atom.Atom_ID == Main_ID for atom in self.Atom_List):
                            Main_Atom = self.Get_Atom(Main_ID)
                            Node_ID = int(bond.strip().split(',')[1].split()[0].split('-')[0].strip(")")) + 1
                            if any(atom.Atom_ID == Node_ID for atom in self.Atom_List):
                                Node_Atom = self.Get_Atom(Node_ID)
                                req = np.linalg.norm(Main_Atom.Position - Node_Atom.Position)
                                self.Bond_List.append(Bond.Bond(Main_Atom,Node_Atom,req))
            elif read_bonds:
                read_bonds = False
            if len(line.strip().split()) >= 3 and line.strip().split()[0] == "Mayer" and line.strip().split()[1] == "bond" and line.strip().split()[2] == "orders":
                read_bonds = True
                self.Bond_List = []


## Need to rewrite this function without the hardcoded variables, 
## with the ability to read an input file with the proper inputs
# Also need to rewrite as a qchem script
    # def xyz_write(f, el, x, y, z):
    #     f.write("%s\t%f\t%f\t%f\n" % (el,x,y,z))
    #     return

    # def rob_Optimize_H_Positions(self, cluster_data_file):
    #     f = open("%s_%s_%d.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID),'w')
    #     num_atoms = len(self.Atom_List) + len(self.Bonded_Atoms)
    #     f.write("%d\n\n" % (num_atoms))
    #     for atom in self.Atom_List:
    #         xyz_write(f, atom.Element,atom.Position[0],atom.Position[1],atom.Position[2])
        


    def Optimize_H_Positions(self,Cluster_Location,config_dict = {}, Shared_File_Location = ""):
        # config_dict should be a dictionary containing user-specific compute cluster information
        config_dict = {
            "Cluster_Login" : "0", # temp
            'Base_Cluster_Directory':"0",
            'Scheduler_Type':"0",
            'End_Condition':"0",
            "Executable_Location":"0",
            "OpenMPI_Location":"0",
        }

        #if not isinstance(config_dict, dict):
        #    raise ValueError('config_dict must be a dictionary')

        keys = [
            'Cluster_Login', 
            'Base_Cluster_Directory', 
            'Scheduler_Type', 
            'End_Condition', 
            'Executable_Location',
            'OpenMPI_Location'
            ]

        for key in keys:
            dat = config_dict[key]
            if len(dat) == 0 or not isinstance(dat, str):
                raise ValueError('config_dict must define a value for ' + key)
        
        config_dict['File_Name'] = "sub_%s_%s_%d_Optimize_Monomer" % (self.Polymer_Name,self.Name,self.Ring_ID)
        config_dict['In_File'] = "%s_%s_%d_Optimize_Monomer.inp" % (self.Polymer_Name,self.Name,self.Ring_ID)
        config_dict['End_File'] = "%s_%s_%d_Optimize_Monomer.out" % (self.Polymer_Name,self.Name,self.Ring_ID)
        config_dict['Job_Type'] = "Orca"
        config_dict['Folder_Name'] = "Optimized_Monomers"
        config_dict['Job_Name'] = "%s_%s_%d_Optimize_Monomer" % (self.Polymer_Name,self.Name,self.Ring_ID)
        config_dict['H_Position_File_Name'] = "./Optimized_Monomers/%s_%s_%d_H_Positions.txt" % (self.Polymer_Name,self.Name,self.Ring_ID)
        config_dict['Cluster_Location'] = config_dict['Base_Cluster_Directory'] + config_dict['Job_Name']
        
        File_Name = config_dict['File_Name']
        In_File = config_dict['In_File']
        End_File = config_dict['End_File']
        Job_Type = config_dict['Job_Type']
        Folder_Name = config_dict['Folder_Name']
        Job_Name = config_dict['Job_Name']

        Cluster_Login = Configure.orca_dict["Cluster_Login"]
        Base_Cluster_Location = Configure.orca_dict["Base_Cluster_Location"]
        Scheduler_Type = Configure.orca_dict["Scheduler_Type"]
        End_Condition = Configure.orca_dict["End_Condition"]
        #Executable_Location = Configure.qchem_dict["Executable_Location"]
        #OpenMP_Location = Configure.qchem_dict["OpenMP_Location"]

        print("I'm making %s_%s_%d.xyz!"% (self.Polymer_Name,self.Name,self.Ring_ID))#TODO: temporary, delete later
        f = open("%s_%s_%d.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID),'w')
        num_atoms = len(self.Atom_List) + len(self.Bonded_Atoms)
        
        f.write("%d\n\n" % (num_atoms))
        for atom in self.Atom_List:
            f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        for atom in self.Bonded_Atoms:
            f.write("%s\t%f\t%f\t%f\n" % ("H",atom.H_Atom.Position[0],atom.H_Atom.Position[1],atom.H_Atom.Position[2]))
        f.close()
        Monomer = Molecule.Molecule("%s_%s_%d.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID))
        if not os.path.isdir('Optimized_Monomers'):
            os.mkdir("Optimized_Monomers")


        Write_Inputs.Write_Orca_Optimize_Geometry(config_dict['In_File'],Monomer,H_Only = True)
        #Write_Submit_Script.Write_TORQUE(File_Name,In_File,Job_Name,1,Cluster_Location,Job_Type,Executable_Path = Executable_Location,OMP_Path = OpenMP_Location)
        print("Writing slurm; File_Name: %s, In_File: %s, Job_Name: %s" %(File_Name, In_File, Job_Name))#temp
        Write_Submit_Script.Write_SLURM(File_Name,In_File,Job_Name,1,Cluster_Location,Job_Type)
        Copy_File_List = [config_dict['File_Name'],config_dict['In_File']]
        if not os.path.exists("./Optimized_Monomers/%s" % End_File):
            print("started submit job, in ring.py")#TODO: temporary, delete later
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,File_Name,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            print("started return_info, in ring.py")#TODO: temporary, delete later
            Cluster_IO.Return_Info(End_File,End_File,Folder_Name,Job_Type,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Shared_File_Location = Shared_File_Location)
        print(End_File)
        #raise Exception("You've reached the end of the line :(, not implemented right now") 
        self.Update_Positions_Orca_Output("./Optimized_Monomers/%s" % End_File)
        self.Set_Up_Bonds("./Optimized_Monomers/%s" % End_File)


        
        #self.Update_H_Positions_Orca("./Optimized_Monomers/%s" % Log_File)
        return "%s_%s_%d.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID)

    def Optimize_H_Positions_Hydrogenated(self,Cluster_Location,Shared_File_Location = ""):
        """Aromatic_H_List = []
        Aromatic_C_List = []
        self.Hydrogenated_H_List = []
        self.Hydrogenated_C_List = []
        f = open("Position_Check_1.txt",'w')
        for atom in self.Atom_List:
            f.write("%d %s %.4f %.4f %.4f" % ())"""
        """for bond in self.Bond_List:
            if bond.Bond_Main.Element == "H" and bond.Bond_Node in self.Core_Atom_List:
                Aromatic_H_List.append(bond.Bond_Main)
                Aromatic_C_List.append(bond.Bond_Node)
            elif bond.Bond_Node.Element == "H" and bond.Bond_Main in self.Core_Atom_List:
                Aromatic_H_List.append(bond.Bond_Node)
                Aromatic_C_List.append(bond.Bond_Main)

        for i,h_atom,c_atom in zip(range(len(Aromatic_H_List)),Aromatic_H_List,Aromatic_C_List):
            self.Atom_List.remove(h_atom)
            h_atom_copy1 = copy.deepcopy(h_atom)
            h_atom_copy2 = copy.deepcopy(h_atom)
            h_atom_vector = h_atom.Position - c_atom.Position
            h_atom_unit_vector = h_atom_vector/np.linalg.norm(h_atom_vector)
            New_X = self.Normal_Vector
            New_Y = h_atom_unit_vector
            New_X = New_X - np.dot(New_X,New_Y) * New_Y
            New_Z = np.cross(New_X,New_Y)
            New_Z = New_Z/np.linalg.norm(New_Z)
            Rotation_Matrix = np.array([New_X,New_Y,New_Z])
            Inv_Rotation_Matrix = np.linalg.inv(Rotation_Matrix)
            Copy1_Rotation_Matrix = np.array([[math.cos(math.radians(55)),-math.sin(math.radians(55)),0],[math.sin(math.radians(55)),math.cos(math.radians(55)),0],[0,0,1]])
            Copy2_Rotation_Matrix = np.array([[math.cos(math.radians(-55)),-math.sin(math.radians(-55)),0],[math.sin(math.radians(-55)),math.cos(math.radians(-55)),0],[0,0,1]])
            h_atom_copy1.Position = c_atom.Position + np.dot(Inv_Rotation_Matrix,np.dot(Copy1_Rotation_Matrix,np.dot(Rotation_Matrix,h_atom_vector)))
            h_atom_copy2.Position = c_atom.Position + np.dot(Inv_Rotation_Matrix,np.dot(Copy2_Rotation_Matrix,np.dot(Rotation_Matrix,h_atom_vector)))
            self.Hydrogenated_H_List.append(h_atom_copy1)
            self.Hydrogenated_H_List.append(h_atom_copy2)
            self.Hydrogenated_C_List.append(c_atom)
            self.Hydrogenated_C_List.append(c_atom)
            self.Atom_List.append(h_atom_copy1)
            h_atom_copy2.Atom_ID = len(self.Atom_List) + 1
            self.Atom_List.append(h_atom_copy2)

        for b_atom in self.Bonded_Atoms:
            h_atom_copy1 = copy.deepcopy(b_atom.H_Atom)
            h_atom_copy2 = copy.deepcopy(b_atom.H_Atom)
            h_atom_vector = b_atom.H_Bond_Vector
            h_atom_unit_vector = h_atom_vector/np.linalg.norm(h_atom_vector)
            New_X = self.Normal_Vector
            New_Y = h_atom_unit_vector
            New_X = New_X - np.dot(New_X,New_Y) * New_Y
            New_Z = np.cross(New_X,New_Y)
            New_Z = New_Z/np.linalg.norm(New_Z)
            Rotation_Matrix = np.array([New_X,New_Y,New_Z])
            Inv_Rotation_Matrix = np.linalg.inv(Rotation_Matrix)
            Copy1_Rotation_Matrix = np.array([[math.cos(math.radians(55)),-math.sin(math.radians(55)),0],[math.sin(math.radians(55)),math.cos(math.radians(55)),0],[0,0,1]])
            Copy2_Rotation_Matrix = np.array([[math.cos(math.radians(-55)),-math.sin(math.radians(-55)),0],[math.sin(math.radians(-55)),math.cos(math.radians(-55)),0],[0,0,1]])
            h_atom_copy1.Position = b_atom.Central_Atom.Position + np.dot(Inv_Rotation_Matrix,np.dot(Copy1_Rotation_Matrix,np.dot(Rotation_Matrix,h_atom_vector)))
            h_atom_copy2.Position = b_atom.Central_Atom.Position + np.dot(Inv_Rotation_Matrix,np.dot(Copy2_Rotation_Matrix,np.dot(Rotation_Matrix,h_atom_vector)))
            self.Hydrogenated_H_List.append(h_atom_copy1)
            self.Hydrogenated_H_List.append(h_atom_copy2)
            self.Hydrogenated_C_List.append(b_atom.Central_Atom)
            self.Hydrogenated_C_List.append(b_atom.Central_Atom)
            h_atom_copy1.Atom_ID = len(self.Atom_List) + 1
            self.Atom_List.append(h_atom_copy1)
            h_atom_copy2.Atom_ID = len(self.Atom_List) + 1
            self.Atom_List.append(h_atom_copy2)

        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)

        f = open("%s_%s_%d_Hydrogenated.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID),'w')
        num_atoms = len(self.Atom_List)
        f.write("%d\n\n" % (num_atoms))
        for atom in self.Atom_List:
            f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        f.close()

        Monomer = Molecule.Molecule("%s_%s_%d_Hydrogenated.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID))
        os.system("mkdir ./Optimized_Monomers")
        File_Name = "sub_%s_%s_%d_Optimize_Monomer_Hydrogenated" % (self.Polymer_Name,self.Name,self.Ring_ID)
        In_File = "%s_%s_%d_Optimize_Monomer_Hydrogenated.inp" % (self.Polymer_Name,self.Name,self.Ring_ID)
        End_File = "%s_%s_%d_Optimize_Monomer_Hydrogenated.out" % (self.Polymer_Name,self.Name,self.Ring_ID)
        Job_Type = "Orca"
        Folder_Name = "Optimized_Monomers"
        Job_Name = "%s_%s_%d_Optimize_Monomer_Hydrogenated" % (self.Polymer_Name,self.Name,self.Ring_ID)
        Cluster_Login = "andrewk@tscc-login.sdsc.edu"
        Base_Cluster_Location = '/oasis/tscc/scratch/andrewk/'
        Scheduler_Type = "TORQUE"
        End_Condition = "Opt_Orca"
        Executable_Location = "/home/andrewk/orca_4_2_0_linux_x86-64_openmpi314"
        OpenMP_Location = "/home/andrewk/openmpi-3.1.4"
        H_Position_File_Name = "./Optimized_Monomers/%s_%s_%d_H_Positions.txt" % (self.Polymer_Name,self.Name,self.Ring_ID)

        Write_Inputs.Write_Orca_Optimize_Geometry(In_File,Monomer,H_Only = True)
        Write_Submit_Script.Write_TORQUE(File_Name,In_File,Job_Name,1,Cluster_Location,Job_Type,Executable_Path = Executable_Location,OMP_Path = OpenMP_Location)
        Copy_File_List = [File_Name,In_File]
        Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,File_Name,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
        Cluster_IO.Return_Info(End_File,End_File,Folder_Name,Job_Type,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Shared_File_Location = Shared_File_Location)
        self.Update_Positions_Orca_Output_Hydrogenated("./Optimized_Monomers/%s" % End_File)

        f = open("%s_%s_%d_Hydrogenated.xyz" % (self.Polymer_Name,self.Name,self.Ring_ID),'w')
        num_atoms = len(self.Atom_List)
        f.write("%d\n\n" % (num_atoms))
        for atom in self.Atom_List:
            f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
        f.close()

        for h_atom in self.Hydrogenated_H_List:
            self.Atom_List.remove(h_atom)
        for h_atom in Aromatic_H_List:
            self.Atom_List.append(h_atom)
        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
        self.Aromatic_H_List = Aromatic_H_List
        self.Aromatic_C_List = Aromatic_C_List"""

        self.Hydrogenated_H_List = []
        self.Hydrogenated_C_List = []
        for b_atom in self.Bonded_Atoms:
            h_atom_copy1 = copy.deepcopy(b_atom.H_Atom)
            h_atom_copy2 = copy.deepcopy(b_atom.H_Atom)
            h_atom_vector = b_atom.H_Bond_Vector
            h_atom_unit_vector = h_atom_vector/np.linalg.norm(h_atom_vector)
            New_X = self.Normal_Vector
            New_Y = h_atom_unit_vector
            New_X = New_X - np.dot(New_X,New_Y) * New_Y
            New_Z = np.cross(New_X,New_Y)
            New_Z = New_Z/np.linalg.norm(New_Z)
            Rotation_Matrix = np.array([New_X,New_Y,New_Z])
            Inv_Rotation_Matrix = np.linalg.inv(Rotation_Matrix)
            Copy1_Rotation_Matrix = np.array([[math.cos(math.radians(55)),-math.sin(math.radians(55)),0],[math.sin(math.radians(55)),math.cos(math.radians(55)),0],[0,0,1]])
            Copy2_Rotation_Matrix = np.array([[math.cos(math.radians(-55)),-math.sin(math.radians(-55)),0],[math.sin(math.radians(-55)),math.cos(math.radians(-55)),0],[0,0,1]])
            h_atom_copy1.Position = b_atom.Central_Atom.Position + np.dot(Inv_Rotation_Matrix,np.dot(Copy1_Rotation_Matrix,np.dot(Rotation_Matrix,h_atom_vector)))
            h_atom_copy2.Position = b_atom.Central_Atom.Position + np.dot(Inv_Rotation_Matrix,np.dot(Copy2_Rotation_Matrix,np.dot(Rotation_Matrix,h_atom_vector)))
            self.Hydrogenated_H_List.append(h_atom_copy1)
            self.Hydrogenated_H_List.append(h_atom_copy2)
            self.Hydrogenated_C_List.append(b_atom.Central_Atom)
            self.Hydrogenated_C_List.append(b_atom.Central_Atom)
            h_atom_copy1.Atom_ID = len(self.Atom_List) + 1
            h_atom_copy2.Atom_ID = len(self.Atom_List) + 1

        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)


    #Link Bond_Atoms to other rings and calculate angles with other rings
    def Link_Rings(self,Link_Ring):
        Linked = False
        for b_atom in Link_Ring.Bonded_Atoms:
            for b_atom2 in self.Bonded_Atoms:
                if not b_atom.Is_Linked and not b_atom2.Is_Linked and b_atom.Bonded_Vector[0] == -1*b_atom2.Bonded_Vector[0] and b_atom.Bonded_Vector[1] == -1*b_atom2.Bonded_Vector[1] and b_atom.Bonded_Vector[2] == -1*b_atom2.Bonded_Vector[2]:
                    b_atom.Add_Ring(self)
                    b_atom2.Add_Ring(Link_Ring)
                    center_position = b_atom2.Central_Atom.Position
                    bond_vector = b_atom.Bonded_Vector
                    Linked = True
                    break
                else:
                    continue
            break

        if not Linked:
            alignment = 0
            for b_atom in Link_Ring.Bonded_Atoms:
                for b_atom2 in self.Bonded_Atoms:
                    if not b_atom.Is_Linked and not b_atom2.Is_Linked and np.dot(b_atom.Bonded_Vector/np.linalg.norm(b_atom.Bonded_Vector),-1*b_atom2.Bonded_Vector/np.linalg.norm(b_atom2.Bonded_Vector)) > alignment:
                        alignment = np.dot(b_atom.Bonded_Vector,-1*b_atom2.Bonded_Vector)
                        temp_b_atom = b_atom
                        temp_b_atom2 = b_atom2
                        center_position = b_atom2.Central_Atom.Position
                        bond_vector = b_atom.Bonded_Vector
                
            temp_b_atom.Add_Ring(self)
            temp_b_atom2.Add_Ring(Link_Ring)
            b_atom = b_atom
       
            if alignment == 0:
                print("Ring 1:")
                for b_atom in Link_Ring.Bonded_Atoms:
                    print((b_atom.Bonded_Vector))
                    print((b_atom.Is_Linked))
                print("Ring 2:")
                for b_atom in self.Bonded_Atoms:
                    print((b_atom.Bonded_Vector))
                    print((b_atom.Is_Linked))
                raise Exception("Rings not properly linked")

        return center_position,bond_vector,b_atom

    def Link_Bonded_Atoms(self):
        for b_atom in self.Bonded_Atoms:
            if b_atom.Is_Linked:
                for b_atom_2 in b_atom.Bonded_Ring.Bonded_Atoms:
                    if b_atom_2.Is_Linked and b_atom_2.Bonded_Ring == self:
                        b_atom.Interring_Bond_Atom = b_atom_2

    #Update Normal Vector
    def Update_Normal_Vector(self):
        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2:]:
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    if np.linalg.norm(self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))) > np.linalg.norm(self.Normal_Vector - (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))):
                        self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
                    else:
                        self.Normal_Vector = self.Normal_Vector - (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)

    def Show_Normal_Vector(self):
        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2:]:
                    #print("Positions")
                    #print(atom1.Position)
                    #print(atom2.Position)
                    #print(atom3.Position)
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    #print("Vectors")
                    #print(vec1)
                    #print(vec2)
                    #print(np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
                    if np.linalg.norm(self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))) > np.linalg.norm(self.Normal_Vector - (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))):
                        self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
                    else:
                        self.Normal_Vector = self.Normal_Vector - (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
                    #print("New Normal:")
                    #print(self.Normal_Vector)
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)
        #print(self.Normal_Vector)
    
    def Rotate_Ring(self,Rotation_Matrix):
        for atom in self.Atom_List:
            atom.Position = np.matmul(Rotation_Matrix,atom.Position)
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Position = np.matmul(Rotation_Matrix,b_atom.H_Atom.Position)
            b_atom.H_Bond_Vector = np.matmul(Rotation_Matrix,b_atom.H_Bond_Vector)
        for h_atom in self.Hydrogenated_H_List:
            h_atom.Position = np.matmul(Rotation_Matrix,h_atom.Position)
        self.Update_Normal_Vector()

    def Rotate_Ring_Not_Interring(self,Rotation_Matrix):
        for atom in self.Atom_List:
            atom.Position = np.matmul(Rotation_Matrix,atom.Position)
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Position = np.matmul(Rotation_Matrix,b_atom.H_Atom.Position)
            b_atom.H_Bond_Vector = np.matmul(Rotation_Matrix,b_atom.H_Bond_Vector)
        self.Update_Normal_Vector()

    def Translate_Ring(self,Center_Position):
        for atom in self.Atom_List:
            atom.Position = atom.Position - Center_Position
        for b_atom in self.Bonded_Atoms:
            b_atom.H_Atom.Position = b_atom.H_Atom.Position - Center_Position
        for h_atom in self.Hydrogenated_H_List:
            h_atom.Position = h_atom.Position - Center_Position


    def Add_Bond_List(self):
        for atom in self.Atom_List:
            for bond in self.Bond_List:
                if bond.Bond_Main == atom:
                    atom.Bond_List.append(bond.Bond_Node)
                elif bond.Bond_Node == atom:
                    atom.Bond_List.append(bond.Bond_Main)
                if len(atom.Bond_List) == 4 or (atom.Element == "H" and len(atom.Bond_List) == 1):
                    break
        for b_atom in self.Bonded_Atoms:
            if b_atom.Is_Linked:
                b_atom.Central_Atom.Bond_List.append(b_atom.Interring_Bond_Atom.Central_Atom)
            else:
                b_atom.Central_Atom.Bond_List.append(b_atom.H_Atom)
                b_atom.H_Atom.Bond_List.append(b_atom.Central_Atom)

    def Improper_Bend_Control_Submit(self):
        In_File_List = []
        Copy_File_List = []
        End_File_List = []
        Initial_Position = self.Bonded_Atoms[0].Central_Atom.Position
        self.Translate_Ring(copy.deepcopy(Initial_Position))
        for i in range(10):
            #f = open("%s_XYZ_Improper_Bend_Phi_%d.xyz" % (self.Name,i*10),'w')
            f = open("%s_XYZ_Improper_Bend_Phi_%d.xyz" % (self.Name,i*5),'w')
            f.write("%d\n\n" % (len(self.Atom_List)+2))
            Align_X = self.Normal_Vector
            Align_Y = self.Bonded_Atoms[0].H_Bond_Vector/np.linalg.norm(self.Bonded_Atoms[0].H_Bond_Vector)
            Align_X = (Align_X - np.dot(Align_X,Align_Y)*Align_Y)/np.linalg.norm(Align_X - np.dot(Align_X,Align_Y)*Align_Y)
            Align_Z = np.cross(Align_X,Align_Y)/np.linalg.norm(np.cross(Align_X,Align_Y))
            Align_Matrix = [Align_X,Align_Y,Align_Z]

            for atom in self.Atom_List:
                f.write("%s\t%f\t%f\t%f\n" % (atom.Element,np.matmul(Align_Matrix,atom.Position)[0],np.matmul(Align_Matrix,atom.Position)[1],np.matmul(Align_Matrix,atom.Position)[2]))
            #Rotation_Angle = math.radians(i*10)
            Rotation_Angle = math.radians(i*5)
            Rotate_Matrix = [[math.cos(Rotation_Angle),-math.sin(Rotation_Angle),0],[math.sin(Rotation_Angle),math.cos(Rotation_Angle),0],[0,0,1]]
            f.write("%s\t%f\t%f\t%f\n" % (self.Bonded_Atoms[0].H_Atom.Element,np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,self.Bonded_Atoms[0].H_Atom.Position))[0],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,self.Bonded_Atoms[0].H_Atom.Position))[1],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,self.Bonded_Atoms[0].H_Atom.Position))[2]))
            f.write("%s\t%f\t%f\t%f\n" % (self.Bonded_Atoms[1].H_Atom.Element,np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[0],np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[1],np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[2]))
            f.close()
            """Bent_Molecule = Molecule.Molecule("%s_XYZ_Improper_Bend_Phi_%d.xyz" % (self.Name,i*10))
            os.system("scp %s_XYZ_Improper_Bend_Phi_%d.xyz ./XYZ_Files" % (self.Name,i*10))
            os.system("rm -f ./%s_XYZ_Improper_Bend_Phi_%d.xyz" % (self.Name,i*10))"""
            Bent_Molecule = Molecule.Molecule("%s_XYZ_Improper_Bend_Phi_%d.xyz" % (self.Name,i*5))
            os.system("scp %s_XYZ_Improper_Bend_Phi_%d.xyz ./XYZ_Files" % (self.Name,i*5))
            os.remove("%s_XYZ_Improper_Bend_Phi_%d.xyz" % (self.Name,i*5))

            Job_Type = "QChem"
            Folder_Name = "Improper_Bend_Test"
            #End_File = "%s_Improper_Bend_Phi_%d.out" % (self.Name,i*10)
            End_File = "%s_Improper_Bend_Phi_%d.out" % (self.Name,i*5)
            Cluster_Login = Configure.orca_dict["Cluster_Login"]
            Base_Cluster_Location = Configure.orca_dict["Base_Cluster_Location"]
            Cluster_Location=Base_Cluster_Location+"/Improper_Bend_Test"
            Scheduler_Type = Configure.orca_dict["Scheduler_Type"]
            End_Condition = Configure.orca_dict["End_Condition"]
            Shared_File_Location = Configure.orca_dict["Shared_File_Location"]
            """Job_Name = "%s_Improper_Bend_Phi_%d" % (self.Name,i*10)
            In_File = "%s_Improper_Bend_Phi_%d.qcin" % (self.Name,i*10)
            Sub_File = "sub_%s_Improper_Bend_Phi_%d" % (self.Name,i*10)"""
            Job_Name = "%s_Improper_Bend_Phi_%d" % (self.Name,i*5)
            In_File = "%s_Improper_Bend_Phi_%d.qcin" % (self.Name,i*5)
            Sub_File = "sub_%s_Improper_Bend_Phi_%d" % (self.Name,i*5)
            qos = "premium"
            Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
            if Return_File == "":
                Return_File = End_File
            Copy_File_List.append(In_File)
            End_File_List.append(Return_File)
            if not Finished:
                In_File_List.append(In_File)
                Write_Inputs.Write_QChem_SPE(In_File,Bent_Molecule)

        if len(In_File_List) != 0:
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            for file in Copy_File_List:
                os.system("scp %s ./Rotation_Run_Input_Copies" % file)
                os.remove(file)
        self.Translate_Ring(-1*Initial_Position)

    def Improper_Bend_Methyl_Control_Submit(self,Num_Rotations,Step):
        In_File_List = []
        Copy_File_List = []
        End_File_List = []
        Rotated_Flag = False
        Rotation_B_Atom = self.Bonded_Atoms[0]
        Initial_Position = Rotation_B_Atom.Central_Atom.Position
        self.Translate_Ring(copy.deepcopy(Initial_Position))
        for i in range(Num_Rotations):
            #f = open("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*10),'w')
            print("Making %s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz"% (self.Name,i*Step) ) #TODO: temporary?
            f = open("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*Step),'w')
            f.write("%d\n\n" % (len(self.Atom_List)+5))
            Align_X = self.Normal_Vector
            Align_Y = Rotation_B_Atom.Bonded_Vector/np.linalg.norm(Rotation_B_Atom.Bonded_Vector)
            Align_X = (Align_X - np.dot(Align_X,Align_Y)*Align_Y)/np.linalg.norm(Align_X - np.dot(Align_X,Align_Y)*Align_Y)
            Align_Z = np.cross(Align_X,Align_Y)/np.linalg.norm(np.cross(Align_X,Align_Y))
            Align_Matrix = [Align_X,Align_Y,Align_Z]

            for atom in self.Atom_List:
                f.write("%s\t%f\t%f\t%f\n" % (atom.Element,np.matmul(Align_Matrix,atom.Position)[0],np.matmul(Align_Matrix,atom.Position)[1],np.matmul(Align_Matrix,atom.Position)[2]))
            #Rotation_Angle = math.radians(i*10)
            Rotation_Angle = math.radians(i*Step)
            Rotate_Matrix = [[math.cos(Rotation_Angle),-math.sin(Rotation_Angle),0],[math.sin(Rotation_Angle),math.cos(Rotation_Angle),0],[0,0,1]]

            """Rotated_Flag = False
            for b_atom in self.Bonded_Atoms:
                if not Rotated_Flag and b_atom.Is_Linked:
                    f.write("%s\t%f\t%f\t%f\n" % ("C",np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,b_atom.Bonded_Vector))[0],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,b_atom.Bonded_Vector))[1],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,b_atom.Bonded_Vector))[2]))
                    Methyl_Bond = np.array([0.0,1.08,0.0])
                    rotations = [math.radians(0),math.radians(120),math.radians(240)]
                    oop_rot = math.radians(69)
                    OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
                    Methyl_Bond = np.matmul(OOP_Matrix,Methyl_Bond)
                    print("Hello")
                    for rot in rotations:
                        print("Made it")
                        Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                        New_H_Atom_Position = np.matmul(Align_Matrix,b_atom.Bonded_Vector) + np.matmul(Dih_Matrix,Methyl_Bond)
                        f.write("%s\t%f\t%f\t%f\n" % ("H",np.matmul(Rotate_Matrix,New_H_Atom_Position)[0],np.matmul(Rotate_Matrix,New_H_Atom_Position)[1],np.matmul(Rotate_Matrix,New_H_Atom_Position)[2]))
                    Rotated_Flag = True
                else:
                    f.write("%s\t%f\t%f\t%f\n" % (b_atom.H_Atom.Element,np.matmul(Align_Matrix,b_atom.H_Atom.Position)[0],np.matmul(Align_Matrix,b_atom.H_Atom.Position)[1],np.matmul(Align_Matrix,b_atom.H_Atom.Position)[2]))"""
            f.write("%s\t%f\t%f\t%f\n" % ("C",np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector))[0],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector))[1],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector))[2]))
            Methyl_Bond = np.array([0.0,1.08,0.0])
            rotations = [math.radians(0),math.radians(120),math.radians(240)]
            oop_rot = math.radians(69)
            OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
            Methyl_Bond = np.matmul(OOP_Matrix,Methyl_Bond)
            for rot in rotations:
                Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                New_H_Atom_Position = np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector) + np.matmul(Dih_Matrix,Methyl_Bond)
                f.write("%s\t%f\t%f\t%f\n" % ("H",np.matmul(Rotate_Matrix,New_H_Atom_Position)[0],np.matmul(Rotate_Matrix,New_H_Atom_Position)[1],np.matmul(Rotate_Matrix,New_H_Atom_Position)[2]))
            f.write("%s\t%f\t%f\t%f\n" % (self.Bonded_Atoms[1].H_Atom.Element,np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[0],np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[1],np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[2]))

            f.close()
            """Bent_Molecule = Molecule.Molecule("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*10))
            os.system("scp %s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz ./XYZ_Files" % (self.Name,i*10))
            os.system("rm -f ./%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*10))"""
            Bent_Molecule = Molecule.Molecule("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*Step))
            os.system("scp %s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz ./XYZ_Files" % (self.Name,i*Step))
            os.remove("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*Step))

            Job_Type = "QChem"
            Folder_Name = "Improper_Bend_Test"
            #End_File = "%s_Improper_Bend_Methyl_Phi_%d.out" % (self.Name,i*10)
            End_File = "%s_Improper_Bend_Methyl_Phi_%d.out" % (self.Name,i*Step)
            Cluster_Login = Configure.qchem_dict["Cluster_Login"]
            Base_Cluster_Location = Configure.qchem_dict["Base_Cluster_Location"]
            Cluster_Location=Base_Cluster_Location+"/" + Folder_Name
            Scheduler_Type = Configure.qchem_dict["Scheduler_Type"]
            End_Condition = Configure.qchem_dict["End_Condition"]
            Shared_File_Location = Configure.qchem_dict["Shared_File_Location"]
            """Job_Name = "%s_Improper_Bend_Methyl_Phi_%d" % (self.Name,i*10)
            In_File = "%s_Improper_Bend_Methyl_Phi_%d.qcin" % (self.Name,i*10)
            Sub_File = "sub_%s_Improper_Bend_Methyl_Phi_%d" % (self.Name,i*10)"""
            Job_Name = "%s_Improper_Bend_Methyl_Phi_%d" % (self.Name,i*Step)
            In_File = "%s_Improper_Bend_Methyl_Phi_%d.qcin" % (self.Name,i*Step)
            Sub_File = "sub_%s_Improper_Bend_Methyl_Phi_%d" % (self.Name,i*Step)
            qos = "overrun time-min=00:30:00"
            Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
            if Return_File == "":
                Return_File = End_File
            Copy_File_List.append(In_File)
            End_File_List.append(Return_File)
            if not Finished:
                In_File_List.append(In_File)
                Write_Inputs.Write_QChem_SPE(In_File,Bent_Molecule)

        #print(In_File_List)
        if len(In_File_List) != 0:
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File_List,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            for file in Copy_File_List:
                os.system("scp %s ./Rotation_Run_Input_Copies" % file)
                os.remove(file)

        self.Translate_Ring(-1*Initial_Position)

    def Improper_Bend_Methyl_Control_Fine_Submit(self):
        In_File_List = []
        Copy_File_List = []
        End_File_List = []
        Rotated_Flag = False
        Rotation_B_Atom = self.Bonded_Atoms[0]
        Initial_Position = Rotation_B_Atom.Central_Atom.Position
        self.Translate_Ring(copy.deepcopy(Initial_Position))
        for i in range(11):
            print("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*2))
            f = open("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*2),'w')
            f.write("%d\n\n" % (len(self.Atom_List)+5))
            Align_X = self.Normal_Vector
            Align_Y = Rotation_B_Atom.Bonded_Vector/np.linalg.norm(Rotation_B_Atom.Bonded_Vector)
            Align_X = (Align_X - np.dot(Align_X,Align_Y)*Align_Y)/np.linalg.norm(Align_X - np.dot(Align_X,Align_Y)*Align_Y)
            Align_Z = np.cross(Align_X,Align_Y)/np.linalg.norm(np.cross(Align_X,Align_Y))
            Align_Matrix = [Align_X,Align_Y,Align_Z]

            for atom in self.Atom_List:
                f.write("%s\t%f\t%f\t%f\n" % (atom.Element,np.matmul(Align_Matrix,atom.Position)[0],np.matmul(Align_Matrix,atom.Position)[1],np.matmul(Align_Matrix,atom.Position)[2]))
            Rotation_Angle = math.radians(i*2)
            Rotate_Matrix = [[math.cos(Rotation_Angle),-math.sin(Rotation_Angle),0],[math.sin(Rotation_Angle),math.cos(Rotation_Angle),0],[0,0,1]]

            """Rotated_Flag = False
            for b_atom in self.Bonded_Atoms:
                if not Rotated_Flag and b_atom.Is_Linked:
                    f.write("%s\t%f\t%f\t%f\n" % ("C",np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,b_atom.Bonded_Vector))[0],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,b_atom.Bonded_Vector))[1],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,b_atom.Bonded_Vector))[2]))
                    Methyl_Bond = np.array([0.0,1.08,0.0])
                    rotations = [math.radians(0),math.radians(120),math.radians(240)]
                    oop_rot = math.radians(69)
                    OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
                    Methyl_Bond = np.matmul(OOP_Matrix,Methyl_Bond)
                    print("Hello")
                    for rot in rotations:
                        print("Made it")
                        Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                        New_H_Atom_Position = np.matmul(Align_Matrix,b_atom.Bonded_Vector) + np.matmul(Dih_Matrix,Methyl_Bond)
                        f.write("%s\t%f\t%f\t%f\n" % ("H",np.matmul(Rotate_Matrix,New_H_Atom_Position)[0],np.matmul(Rotate_Matrix,New_H_Atom_Position)[1],np.matmul(Rotate_Matrix,New_H_Atom_Position)[2]))
                    Rotated_Flag = True
                else:
                    f.write("%s\t%f\t%f\t%f\n" % (b_atom.H_Atom.Element,np.matmul(Align_Matrix,b_atom.H_Atom.Position)[0],np.matmul(Align_Matrix,b_atom.H_Atom.Position)[1],np.matmul(Align_Matrix,b_atom.H_Atom.Position)[2]))"""
            f.write("%s\t%f\t%f\t%f\n" % ("C",np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector))[0],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector))[1],np.matmul(Rotate_Matrix,np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector))[2]))
            Methyl_Bond = np.array([0.0,1.08,0.0])
            rotations = [math.radians(0),math.radians(120),math.radians(240)]
            oop_rot = math.radians(69)
            OOP_Matrix = np.array([[math.cos(oop_rot),-math.sin(oop_rot),0],[math.sin(oop_rot),math.cos(oop_rot),0],[0,0,1]])
            Methyl_Bond = np.matmul(OOP_Matrix,Methyl_Bond)
            for rot in rotations:
                Dih_Matrix = [[math.cos(rot),0,-math.sin(rot)],[0,1,0],[math.sin(rot),0,math.cos(rot)]]
                New_H_Atom_Position = np.matmul(Align_Matrix,Rotation_B_Atom.Bonded_Vector) + np.matmul(Dih_Matrix,Methyl_Bond)
                f.write("%s\t%f\t%f\t%f\n" % ("H",np.matmul(Rotate_Matrix,New_H_Atom_Position)[0],np.matmul(Rotate_Matrix,New_H_Atom_Position)[1],np.matmul(Rotate_Matrix,New_H_Atom_Position)[2]))
            f.write("%s\t%f\t%f\t%f\n" % (self.Bonded_Atoms[1].H_Atom.Element,np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[0],np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[1],np.matmul(Align_Matrix,self.Bonded_Atoms[1].H_Atom.Position)[2]))

            f.close()
            Bent_Molecule = Molecule.Molecule("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*2))
            os.system("scp %s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz ./XYZ_Files" % (self.Name,i*2))
            os.remove("%s_XYZ_Improper_Bend_Methyl_Phi_%d.xyz" % (self.Name,i*2))

            Job_Type = "QChem"
            Folder_Name = "Improper_Bend_Test"
            End_File = "%s_Improper_Bend_Methyl_Phi_%d.out" % (self.Name,i*2)
            Cluster_Login = Configure.orca_dict["Cluster_Login"]
            Base_Cluster_Location = Configure.orca_dict["Base_Cluster_Location"]
            Cluster_Location=Base_Cluster_Location+"/Improper_Bend_Test"
            Scheduler_Type = Configure.orca_dict["Scheduler_Type"]
            End_Condition = Configure.orca_dict["End_Condition"]
            Shared_File_Location = Configure.orca_dict["Shared_File_Location"]
            Job_Name = "%s_Improper_Bend_Methyl_Phi_%d" % (self.Name,i*2)
            In_File = "%s_Improper_Bend_Methyl_Phi_%d.qcin" % (self.Name,i*2)
            Sub_File = "sub_%s_Improper_Bend_Methyl_Phi_%d" % (self.Name,i*2)
            qos = "premium"
            #print(End_File)
            Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
            if Return_File == "":
                Return_File = End_File
            #print(Finished)
            if not Finished:
                End_File_List.append(Return_File)
                Copy_File_List.append(In_File)
                In_File_List.append(In_File)
                Write_Inputs.Write_QChem_SPE(In_File,Bent_Molecule)

        if len(In_File_List) != 0:
            #print(Copy_File_List)
            End_File = Copy_File_List[-1]
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            for file in Copy_File_List:
                os.system("scp %s ./Rotation_Run_Input_Copies" % file)
                os.remove(file)

        self.Translate_Ring(-1*Initial_Position)


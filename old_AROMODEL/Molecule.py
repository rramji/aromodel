#! usr/bin/python


# Import relevant modules
import numpy as np
import os
import subprocess
import shutil
import time
import random
import glob
# import Class structure
import Atom
import Bond
import Angle
import Dihedral
import Improper
import Configure
import matplotlib.pyplot as plt
import copy

Element_Dict = { 12.011:"C", 1.008:"H", 18.998:"F", 15.999:"O", 32.06:"S", 14.007:"N", 28.086:"Si", 35.453:"Cl"}

class Molecule(object):
    """
    Class defining a molecule
    instance variables (self) :
        N = Number of Atoms
        Name = Molecule Name
        Atom_list[N] = List of Atom objects
        Bond_List[Num_Bonds] = list of bond objects
        Angle_List[Num_Angles] = list of angle objects
        Dihedral_List[Num_Dihedrals] = list of Dihedral objects
        Improper_List[Num_Impropers] = list of Improper objects
        Ring_List[Num_Rings] = list of Ring objects
        COM = Center of Mass
        MOL_ID = ID Number for molecule
        UnConverged = Flag for bypassing Orca convergence (Default = False)
    """

    def __init__(self, File_Name,Full = False):
        File = open(File_Name,'r') # File_Name is the name of an .xyz file outputted by Avogadro
        File_Lines = File.readlines()
        
        # Define Instance Variables
        self.Name = File_Name.split('/')[-1].split('.')[0]
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.Ring_List = []
        self.MW = 0.0
        self.COM = np.zeros(3,float)
        self.Mol_ID = 0
        self.Missing_Dihedrals = 0
        self.UnConverged = True # Unconverged Orca Optimization
        
        if File_Name.split('/')[-1].split('.')[-1] != 'data' and File_Name.split('/')[-1].split('.')[0] != 'data':
            self.N = int(File_Lines[0].strip('\n')) # Integer
            self.Atom_List = np.empty(self.N, dtype=object) # Numpy object array
            #print "Setting up molecule"
            #print "Molecule Name:", self.Name
            #print self.N, "Atoms in ", self.Name
            #print "----------------------------------"
            for i in range(self.N):
                Line = File_Lines[2+i]
                Line = Line.strip('\n').split()
                Element = Line[0]
                Position = np.array( [ float(Line[1]), float(Line[2]), float(Line[3]) ], dtype=float )
                self.Atom_List[i] = Atom.Atom(Position, Element, i+1) # Instantiate Atom_List with Atom objects
                self.MW += self.Atom_List[i].Mass
            
            #print "Initial XYZ Coordinates:\n"
            for Atom_Obj in self.Atom_List:
                print Atom_Obj.Element, Atom_Obj.Position
            #print "----------------------------------"
            
            # Compute center of Mass
            Mass_Weighted_Sum = np.zeros(3,dtype=float)
            for Atom_Obj in self.Atom_List:
                Mass_Weighted_Sum += Atom_Obj.Position*Atom_Obj.Mass
                self.MW += Atom_Obj.Mass
            
            #print "Molecular Weight is ", self.MW, "Grams/Mole"
            self.COM = Mass_Weighted_Sum/self.MW
            self.COM = np.asarray(self.COM)
            #print "COM is", self.COM
            
            Mass_Weighted_Sum = 0.0
            for Atom_Obj in self.Atom_List:
                Mass_Weighted_Sum += Atom_Obj.Mass*((Atom_Obj.Position[0] - self.COM[0])**2 + (Atom_Obj.Position[1] - self.COM[1])**2  + (Atom_Obj.Position[2] - self.COM[2])**2 )
            Rg2 = Mass_Weighted_Sum/self.MW
            self.Rg = np.sqrt(Rg2)
            #print "The Radius of Gyration is", self.Rg
                    
                
            
            # Zero COM
            for Atom_Obj in self.Atom_List:
                Atom_Obj.Position -= self.COM
            
            self.COM -= self.COM
            #print self.COM

        else:
            j = -1
            for Line in File_Lines:
                j += 1
                Line = Line.strip('\n').split(' ')
                #print(Line)
                if len(Line) > 1:
                    if len(Line) > 2:
                        if Line[2] == 'xlo':
                            #print(Line[1])
                            #self.Box_Size = float(Line[1].split('e')[0])*(10**(int(float(Line[1].split('+')[1]))))
                            self.Box_Size = float(Line[1])
                    if Line[1] == 'atoms':
                        self.N = int(Line[0])
                        #print self.N, "Atoms"
                    if Line[1] == 'atom':
                        Atom_Types = int(Line[0])
                        self.Atom_List = np.empty(self.N, dtype=object)
                        self.I_Flags_List = np.empty(self.N, dtype=object)
                        Mass_List = np.empty(Atom_Types, dtype=float)
                        Pair_Coeffs = np.empty((Atom_Types,2), dtype=float)
                    if Line[1] == 'bonds':
                        Num_Bonds = Line[0]
                    if Line[1] == 'bond':
                        Bond_Types = int(Line[0])
                        Bond_Coeffs = np.empty((Bond_Types,2), dtype=float)
                    if Line[1] == 'angles':
                        Num_Angles = Line[0]
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
                    if Line[0] == 'Atoms':
                        #print 'Getting Atoms'
                        for i in range(self.N):
                            Temp_Position = [float(File_Lines[i+j+2].split(' ')[4]), float(File_Lines[i+j+2].split(' ')[5]), float(File_Lines[i+j+2].split(' ')[6]) ]
                            try:
                                Temp_I_Flags = [ int(File_Lines[i+j+2].split(' ')[7]), int(File_Lines[i+j+2].split(' ')[8]), int(File_Lines[i+j+2].split(' ')[9])]
                                Temp_I_Flags = np.asarray(Temp_I_Flags, dtype = int)
                                Temp_Position[0] = Temp_Position[0] + Temp_I_Flags[0] * self.Box_Size
                                Temp_Position[1] = Temp_Position[1] + Temp_I_Flags[1] * self.Box_Size
                                Temp_Position[2] = Temp_Position[2] + Temp_I_Flags[2] * self.Box_Size
                                self.I_Flags_List[i] = Temp_I_Flags
                            except:
                                print("No Image Flags")
                            Temp_Position = np.asarray(Temp_Position, dtype=float)

                            #I_Flags = np.asarray(Temp_I_Flags, dtype=int)
                            #print "IMAGE FLAGS", I_Flags
                            #Temp_Position += I_Flags*Box_Length
                            Temp_Charge = float(File_Lines[i+j+2].split(' ')[3])
                            Type = int(File_Lines[i+j+2].split(' ')[2])
                            Mass = Mass_List[Type-1]
                            Element = Element_Dict[Mass]
                            Atom_Id = int(File_Lines[i+j+2].split(' ')[0])
                            #print Element, Mass, Atom_Id
                            self.Atom_List[i] = Atom.Atom(Temp_Position, Element, Atom_Id)
                            self.Atom_List[i].Charge = Temp_Charge
                            self.Atom_List[i].Epsilon = Pair_Coeffs[Type-1,0]
                            self.Atom_List[i].Sigma = Pair_Coeffs[Type-1,1]
                            if Full:
                                self.Atom_List[i].OPLS_Type = 2000 + Type
                                self.Atom_List[i].LAMMPS_Type = 2000 + Type
                            else:
                                self.Atom_List[i].OPLS_Type = Type
                                self.Atom_List[i].LAMMPS_Type = Type
                        self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
                        #for Atom_Obj in self.Atom_List:
                            #print Atom_Obj.Atom_ID, Atom_Obj.Position
                    try:
                        if Line[2] == 'xlo':
                            Box_Length = float(Line[1]) - float(Line[0])
                            #print Box_Length, "Box_Length"
                    except:
                        continue

                    if Line[0] == 'Pair':
                        #print "Getting Pair Coefficients"
                        for i in range(Atom_Types):
                            Pair_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                            Pair_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                    if Line[0] == 'Bond':
                        #print "Getting Bond Coefficients"
                        for i in range(Bond_Types):
                            Bond_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                            Bond_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                    if Line[0] == 'Angle':
                        #print "Getting Angle Coefficients"
                        for i in range(Angle_Types):
                            Angle_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                            Angle_Coeffs[i,1] = float(File_Lines[i+j+2].split(' ')[2])
                    if Line[0] == 'Dihedral':
                        #print "Getting Dihedral Coefficients"
                        #print(len(Dihedral_Coeffs))
                        for i in range(Dihedral_Types):
                            dummy_coeffs = []
                            for k in range(1,len(File_Lines[i+j+2].split(' '))):
                                dummy_coeffs.append(float(File_Lines[i+j+2].split(' ')[k]))
                            #print(dummy_coeffs)
                            Dihedral_Coeffs.append(np.array(dummy_coeffs))
                            Dihedral_Styles.append(File_Lines[i+j+2].split(' ')[0])
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
                            Improper_Coeffs[i,0] = float(File_Lines[i+j+2].split(' ')[1])
                            Improper_Coeffs[i,1] = int(File_Lines[i+j+2].split(' ')[2])
                            Improper_Coeffs[i,2] = int(File_Lines[i+j+2].split(' ')[3])



                if len(Line) == 1:
                    if Line == ['']:
                        continue
                    if Line == ['Masses']:
                        #print "Extracting Masses"
                        for i in range(Atom_Types):
                            Mass_List[i] = float(File_Lines[i+j+2].split(' ')[1])
                    if Line == ['Bonds']:
                        #print "Extracting Bonds"
                        for i in range(int(Num_Bonds)):
                            Bond_Info = File_Lines[i+j+2].strip('\n').split(' ')
                            Master = self.Atom_List[int(Bond_Info[2])-1]
                            Slave = self.Atom_List[int(Bond_Info[3])-1]
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
                            Angle_Info = File_Lines[i+j+2].strip('\n').split(' ')
                            Slave1 = self.Atom_List[int(Angle_Info[2])-1]
                            Master = self.Atom_List[int(Angle_Info[3])-1]
                            Slave2 = self.Atom_List[int(Angle_Info[4])-1]
                            Ka = Angle_Coeffs[int(Angle_Info[1])-1,0]
                            Th0 = Angle_Coeffs[int(Angle_Info[1])-1,1]
                            Angle_ID = int(Angle_Info[0])
                            self.Angle_List.append(Angle.Angle(Master, Slave1, Slave2, Th0))
                            self.Angle_List[i].ka = Ka
                            self.Angle_List[i].Angle_ID = Angle_ID
                            self.Angle_List[i].LAMMPS_Type = int(Angle_Info[1])
                    #print('Check')
                    if Line == ['Dihedrals']:
                        #print "Extracting Dihedrals"
                        for i in range(int(Num_Dihedrals)):
                            Dihedral_Info = File_Lines[i+j+2].strip('\n').split(' ')
                            Slave1 = self.Atom_List[int(Dihedral_Info[2])-1]
                            Master1 = self.Atom_List[int(Dihedral_Info[3])-1]
                            Master2 = self.Atom_List[int(Dihedral_Info[4])-1]
                            Slave2 = self.Atom_List[int(Dihedral_Info[5])-1]
                            Coeffs = Dihedral_Coeffs[int(Dihedral_Info[1])-1]
                            #Coeffs = []
                            #Style = Dihedral_Styles[int(Dihedral_Info[1])-1]
                            Style = int(Dihedral_Info[1])
                            Dihedral_ID = int(Dihedral_Info[0])
                            self.Dihedral_List.append(Dihedral.Dihedral(Master1, Master2, Slave1, Slave2, 0.0))
                            self.Dihedral_List[i].Coeffs = Coeffs
                            self.Dihedral_List[i].Dihedral_ID = Dihedral_ID
                            self.Dihedral_List[i].Style = Style
                            self.Dihedral_List[i].LAMMPS_Type = Dihedral_Info[1]

                    if Line == ['Impropers']:
                        #print "Extracting Impropers"
                        for i in range(int(Num_Impropers)):
                            Improper_Info = File_Lines[i+j+2].strip('\n').split(' ')
                            Master = self.Atom_List[int(Improper_Info[2])-1]
                            Slave1 = self.Atom_List[int(Improper_Info[3])-1]
                            Slave2 = self.Atom_List[int(Improper_Info[4])-1]
                            Slave3 = self.Atom_List[int(Improper_Info[5])-1]
                            Coeff = Improper_Coeffs[int(Improper_Info[1])-1,0]
                            Improper_ID = int(Improper_Info[0])
                            self.Improper_List.append(Improper.Improper(Master, Slave1, Slave2, Slave3, Coeff, 180.0, Improper_ID))
                            self.Improper_List[i].LAMMPS_Type = int(Improper_Info[1])
                            
                                
            # Compute center of Mass
            Mass_Weighted_Sum = np.zeros(3,dtype=float)
            #print len(self.Atom_List)
            for Atom_Obj in self.Atom_List:
                Mass_Weighted_Sum += Atom_Obj.Position*Atom_Obj.Mass
                self.MW += Atom_Obj.Mass
            
            #print "Molecular Weight is ", self.MW, "Grams/Mole"
            self.COM = Mass_Weighted_Sum/self.MW
            #print "COM is", self.COM
            
            # Zero COM
            for Atom_Obj in self.Atom_List:
                Atom_Obj.Position -= self.COM
            
            #self.COM -= self.COM
            #print self.COM


    def Adjust_COM(self,center=np.array([0.0,0.0,0.0]),Rotate=True):
        # This adjusts the center of mass and gives the molecule a random orientation
        x = random.random()*2*3.1415
        y = random.random()*2*3.1415
        z = random.random()*2*3.1415
        
        C1 = np.cos(x)
        S1 = np.sin(x)
        C2 = np.cos(y)
        S2 = np.sin(y)
        C3 = np.cos(z)
        S3 = np.sin(z)
        
        for Atom_Obj in self.Atom_List:
            if Rotate:
                # First rotation
                xt = Atom_Obj.Position[0]
                yt = Atom_Obj.Position[1]
                Atom_Obj.Position[0] = xt*C1 - yt*S1
                Atom_Obj.Position[1] = xt*S1 + yt*C1
                # Second rotation
                xt = Atom_Obj.Position[0]
                zt = Atom_Obj.Position[2]
                Atom_Obj.Position[0] = xt*C2 - zt*S2
                Atom_Obj.Position[2] = xt*S2 + zt*C2
                #Third rotation
                yt = Atom_Obj.Position[1]
                zt = Atom_Obj.Position[2]
                Atom_Obj.Position[1] = yt * C3 - zt * S3
                Atom_Obj.Position[2] = yt * S3 + zt * C3
            Atom_Obj.Position += self.COM + center
        return 

    def Write_Data_File(self,Filename,Exclude_All_Interring=False,Exclude_Interring_Bonds=False,Exclude_Interring_Angles=False,Exclude_Interring_Torsions=False,Exclude_Interring_Dihedrals=False,Interring_Angles_Only=False,Bond_Atoms=[],Soft_Potential=False):
        if Exclude_All_Interring or Exclude_Interring_Bonds:
            Changed_Bonds = []
            for bond in self.Bond_List:
                if bond.Bond_Master in Bond_Atoms or bond.Bond_Slave in Bond_Atoms:
                    Changed_Bonds.append(bond.kb)
                    bond.kb = 0.0

        if Interring_Angles_Only:
            Changed_Bonds = []
            for bond in self.Bond_List:
                Changed_Bonds.append(bond.kb)
                bond.kb = 0.0

        if Exclude_All_Interring or Exclude_Interring_Angles:
            Changed_Angles = []
            for ang in self.Angle_List:
                if (ang.Angle_Master in Bond_Atoms and ang.Angle_Slave1 in Bond_Atoms) or (ang.Angle_Master in Bond_Atoms and ang.Angle_Slave2 in Bond_Atoms):
                    Changed_Angles.append(ang.ka)
                    ang.ka = 0.0

        if Interring_Angles_Only:
            Changed_Angles = []
            for ang in self.Angle_List:
                if not ((ang.Angle_Master in Bond_Atoms and ang.Angle_Slave1 in Bond_Atoms) or (ang.Angle_Master in Bond_Atoms and ang.Angle_Slave2 in Bond_Atoms)):
                    Changed_Angles.append(ang.ka)
                    ang.ka = 0.0

        if Exclude_All_Interring or Exclude_Interring_Torsions or Exclude_Interring_Dihedrals:
            Changed_Dihs = []
            Lammps_Max = 1
            Zero_Dih = 0
            for dih in self.Dihedral_List:
                if dih.LAMMPS_Type > Lammps_Max:
                    Lammps_Max = dih.LAMMPS_Type
                print([0 for i in range(len(dih.Coeffs))])
                if list(dih.Coeffs) == [0 for i in range(len(dih.Coeffs))] or np.array_equal(dih.Coeffs,np.zeros(len(dih.Coeffs))):
                    Zero_Dih = dih.LAMMPS_Type
            for dih in self.Dihedral_List:
                match_count = 0
                if dih.Dihedral_Master1 in Bond_Atoms:
                    match_count += 1
                if dih.Dihedral_Master2 in Bond_Atoms:
                    match_count += 1
                if dih.Dihedral_Slave1 in Bond_Atoms:
                    match_count += 1
                if dih.Dihedral_Slave2 in Bond_Atoms:
                    match_count += 1
                if match_count >= 2:
                    Changed_Dihs.append(copy.deepcopy(dih.Coeffs))
                    dih.Coeffs = np.zeros(len(dih.Coeffs))
                    dih.LAMMPS_Type = Zero_Dih

        if Exclude_All_Interring or Exclude_Interring_Torsions:
            Changed_Imps = []
            for imp in self.Improper_List:
                if (imp.Improper_Master in Bond_Atoms and imp.Improper_Slave3 in Bond_Atoms):
                    Changed_Imps.append(imp.Ki)
                    imp.Ki = 0.0

        if Interring_Angles_Only:
            Changed_Dihs = []
            for dih in self.Dihedral_List:
                Changed_Dihs.append(copy.deepcopy(dih.Coeffs))
                dih.Coeffs = np.zeros(len(dih.Coeffs))

            Changed_Imps = []
            for imp in self.Improper_List:
                Changed_Imps.append(imp.Ki)
                imp.Ki = 0.0

        f = open(Filename,'w')
        f.write("LAMMPS data file\n\n")
        f.write("%d atoms\n" % len(self.Atom_List))
        type_count = 0 
        type_list = []
        Atom_Info = []
        for atom in self.Atom_List:
            if atom.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(atom.LAMMPS_Type)
                Atom_Info.append((atom.LAMMPS_Type,atom.Mass,atom.Epsilon,atom.Sigma))
        Atom_Info = sorted(Atom_Info, key=lambda atom_type: atom_type[0])

        f.write("%d atom types\n" % type_count)

        f.write("%d bonds\n" % len(self.Bond_List))
        type_count = 0 
        type_list = []
        Bond_Info = []
        for bond in self.Bond_List:
            if bond.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(bond.LAMMPS_Type)
                Bond_Info.append((bond.LAMMPS_Type,bond.kb,bond.req))
        Bond_Info = sorted(Bond_Info, key=lambda bond_type: bond_type[0])
        f.write("%d bond types\n" % type_count)
        f.write("%d angles\n" % len(self.Angle_List))
        type_count = 0 
        type_list = []
        Angle_Info = []
        for ang in self.Angle_List:
            if ang.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(ang.LAMMPS_Type)
                Angle_Info.append((ang.LAMMPS_Type,ang.ka,ang.Angle_Eq))
        Angle_Info = sorted(Angle_Info, key=lambda angle_type: angle_type[0])
        f.write("%d angle types\n" % type_count)
        f.write("%d dihedrals\n" % len(self.Dihedral_List))
        type_count = 0 
        type_list = []
        Dihedral_Info = []
        for dih in self.Dihedral_List:
            if dih.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(dih.LAMMPS_Type)
                Dihedral_Info.append((int(dih.LAMMPS_Type),dih.Coeffs))
        Dihedral_Info = sorted(Dihedral_Info, key=lambda dih_type: dih_type[0])
        if Exclude_All_Interring or Exclude_Interring_Torsions:
            dih_type = 1
            for dih in Dihedral_Info:
                if dih[0] == dih_type:
                    dih_type += 1
                else:
                    while dih[0] > dih_type:
                        Dihedral_Info.append((dih_type,np.zeros(len(dih[1]))))
                        dih_type += 1
                        type_count += 1
                        print(str(dih_type) + " " + str(type_count) + "\n")
                    dih_type += 1
        Dihedral_Info = sorted(Dihedral_Info, key=lambda dih_type: dih_type[0])
        f.write("%d dihedral types\n" % Dihedral_Info[-1][0])
        f.write("%d impropers\n" % len(self.Improper_List))
        type_count = 0 
        type_list = []
        Improper_Info = []
        for imp in self.Improper_List:
            if imp.LAMMPS_Type not in type_list:
                type_count += 1
                type_list.append(imp.LAMMPS_Type)
                Improper_Info.append((int(imp.LAMMPS_Type),imp.Ki,imp.d,imp.n))
        Improper_Info = sorted(Improper_Info, key=lambda imp_type: imp_type[0])
        f.write("%d improper types\n\n" % type_count)
        f.write("\n\n\n0.0000 300.0000 xlo xhi\n0.0000 300.0000 ylo yhi\n0.0000 300.0000 zlo zhi\n\n")
        f.write("Masses\n\n")
        for atom in Atom_Info:
            f.write("%d %.3f\n" % (atom[0],atom[1]))
        f.write("\nPair Coeffs\n\n")
        if not Soft_Potential:
            for atom in Atom_Info:
                f.write("%d %f %f\n" % (atom[0],atom[2],atom[3]))
        else:
            f.write("* * 0.0\n")
        f.write("\nBond Coeffs\n\n")
        for bond in Bond_Info:
            f.write("%d %f %f\n" % (bond[0],bond[1],bond[2]))
        f.write("\nAngle Coeffs\n\n")
        for ang in Angle_Info:
            f.write("%d %f %f\n" % (ang[0],ang[1],ang[2]))
        f.write("\nDihedral Coeffs\n\n")
        for dih in Dihedral_Info:
            f.write("%s" % dih[0])
            for c in dih[1]:
                f.write(" %s" % c)
            f.write("\n")
        f.write("\nImproper Coeffs\n\n")
        for imp in Improper_Info:
            f.write("%d %f %d %d\n" % (imp[0],imp[1],imp[2],imp[3]))
        f.write("\nAtoms\n\n")
        for atom in self.Atom_List:
            f.write("%d 1 %d %.16f %.16f %.16f %.16f %d %d %d\n" % (atom.Atom_ID,atom.LAMMPS_Type,atom.Charge,atom.Position[0],atom.Position[1],atom.Position[2],atom.Image_Flags[0],atom.Image_Flags[1],atom.Image_Flags[2]))
        f.write("\nBonds\n\n")
        for bond in self.Bond_List:
                f.write("%d %d %d %d\n" % (bond.Bond_ID,bond.LAMMPS_Type,bond.Bond_Master.Atom_ID,bond.Bond_Slave.Atom_ID))
        f.write("\nAngles\n\n")
        for ang in self.Angle_List:
                f.write("%d %d %d %d %d\n" % (ang.Angle_ID,ang.LAMMPS_Type,ang.Angle_Slave1.Atom_ID,ang.Angle_Master.Atom_ID,ang.Angle_Slave2.Atom_ID))
        f.write("\nDihedrals\n\n")
        for dih in self.Dihedral_List:
                f.write("%d %d %d %d %d %d\n" % (dih.Dihedral_ID,int(dih.LAMMPS_Type),dih.Dihedral_Slave1.Atom_ID,dih.Dihedral_Master1.Atom_ID,dih.Dihedral_Master2.Atom_ID,dih.Dihedral_Slave2.Atom_ID))
        f.write("\nImpropers\n\n")
        for imp in self.Improper_List:
                f.write("%d %d %d %d %d %d\n" % (imp.Improper_ID,int(imp.LAMMPS_Type),imp.Improper_Master.Atom_ID,imp.Improper_Slave1.Atom_ID,imp.Improper_Slave2.Atom_ID,imp.Improper_Slave3.Atom_ID))
        f.close()

        if Exclude_All_Interring or Exclude_Interring_Bonds:
            i = 0
            for bond in self.Bond_List:
                if bond.Bond_Master in Bond_Atoms or bond.Bond_Slave in Bond_Atoms:
                    bond.kb = Changed_Bonds[i]
                    i+=1

        if Exclude_All_Interring or Exclude_Interring_Angles:
            i = 0
            for ang in self.Angle_List:
                if (ang.Angle_Master in Bond_Atoms and ang.Angle_Slave1 in Bond_Atoms) or (ang.Angle_Master in Bond_Atoms and ang.Angle_Slave2 in Bond_Atoms):
                    ang.ka = Changed_Angles[i]
                    i+=1

        if Exclude_All_Interring or Exclude_Interring_Torsions:
            i = 0
            for dih in self.Dihedral_List:
                match_count = 0
                if dih.Dihedral_Master1 in Bond_Atoms:
                    match_count += 1
                if dih.Dihedral_Master2 in Bond_Atoms:
                    match_count += 1
                if dih.Dihedral_Slave1 in Bond_Atoms:
                    match_count += 1
                if dih.Dihedral_Slave2 in Bond_Atoms:
                    match_count += 1
                if match_count >= 2:
                    dih.Coeffs = Changed_Dihs[i]
                    i+=1

            i = 0
            for imp in self.Improper_List:
                if (imp.Improper_Master in Bond_Atoms and imp.Improper_Slave3 in Bond_Atoms):
                    imp.Ki = Changed_Imps[i]
                    i += 1

    def Set_Up_FF(self, run_orca=True, local= True):
        if run_orca:
            print "Setting up Orca input script"
        
            # Write Orca Input File
            File_Name = self.Name + ".inp"
            File = open(File_Name, 'w')
            #File.write('! RKS B3LYP 6-31+G** NormalSCF Opt NOSOSCF CHELPG PAL8\n\n')
            File.write('! RKS B3LYP 6-31+G** NormalSCF Opt NOSOSCF CHELPG\n\n')
            File.write('*xyz 0 1\n')
            for Atom_Obj in self.Atom_List:
                File.write('%s %.5f %.5f %.5f\n' % ( Atom_Obj.Element, Atom_Obj.Position[0], Atom_Obj.Position[1],Atom_Obj.Position[2]))
            File.write('*')
            File.close()
            Finished = False
            
            File_Out = self.Name + ".out"
            if local:
                #Run subprocess on local machine
                
                
                os.system('mkdir Orca')
                os.system('mv %s ./Orca' % File_Name)
                os.chdir('./Orca')
                File_Out = self.Name + ".out"
                try:
                    File = open(File_Out,'r')
                except:
                    os.system('/Users/andrewkleinschmidt/Library/Orca/orca %s > %s' %(File_Name, File_Out)) # Run Orca Job
            
            
            else:
                print "Running Orca Geometry Optimization on Comet"
                cmd = "mkdir " + Configure.Comet_Path % self.Name
                subprocess.call(["ssh", Configure.Comet_Login, cmd])
                subtemp = Configure.Template_Path + "sub_orca_temp"
                submit = "submit_orca"
                Path = Configure.Comet_Path % self.Name
                # Write submit script
                with open(subtemp) as f:
                    template = f.read()
                s = template.format(Comet_Path=Path, Orca_Path = Configure.Orca_Path, name = self.Name )
                with open(submit ,"w") as f:
                    f.write(s)
            
                # Copy Files over  to Comet
                os.system( Configure.c2c % (submit, self.Name))
                os.system( Configure.c2c % (File_Name, self.Name))
                # Run job
                os.system( Configure.c2l % (self.Name, File_Out))
                try:
                    File = open(File_Out,'r')
                except:
                    subprocess.call(["ssh",Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])


                i = 0
                while not Finished:
                    os.system( Configure.c2l % (self.Name, File_Out))
                    try:
                        File = open(File_Out,'r')
                        File_Lines = File.readlines()
                        print File_Lines[-1]
                        if File_Lines[-1].split(' ')[0] == "TOTAL" or self.UnConverged:
                            Finished = True
                        else:
                            print "Not Finished"
                            i += 10
                            print "Sleeping process", i, "Minutes"
                            time.sleep(600)
                            self.Unconverged = True
                    except:
                        print "Sleeping process", i, "miniutes"
                        time.sleep(600)
                        i += 10
    
    
            
    
    
        else:
            os.chdir('./Orca')
        
        File_Out = self.Name + ".out"
        
        # Extract info from Orca output file
        Orca_File = open(File_Out, 'r')
        File_Lines = Orca_File.readlines()
        print "Extracting Redundant Coordinates..."
        for i in range(len(File_Lines)):
            Line = File_Lines[i].strip('\n').split()
            #print Line
            Found = False
            try:
                if Line[0] == "Redundant" and self.UnConverged:
                    Redundant = True
                    j = i + 6
                    k = 0
                    while Redundant:
                        R_Line = File_Lines[j]
                        R_Line = R_Line.split()
                        print R_Line
                        try:
                            if int(R_Line[0].strip('.')) >= 1:
                                j = j + 1
                                Type = R_Line[1].split('(')[0]
                                k += 1
                                # Extract Bonds
                                if Type == "B":
                                    Master_ID = int(R_Line[2].split(',')[0])
                                    Slave_ID = int(R_Line[3].split(')')[0])
                                    req = float(R_Line[-1])
                                    self.Atom_List[Master_ID].Bond_List.append(Slave_ID + 1)
                                    self.Atom_List[Slave_ID].Bond_List.append(Master_ID + 1)
                                    self.Bond_List.append( Bond.Bond(self.Atom_List[Master_ID], self.Atom_List[Slave_ID], req))
                                # Extract Angles
                                if Type == "A":
                                    Master_ID = int(R_Line[3].split(',')[0])
                                    Slave1_ID = int(R_Line[2].split(',')[0])
                                    Slave2_ID = int(R_Line[4].split(')')[0])
                                    Angle_Eq = float(R_Line[-1])
                                    self.Angle_List.append( Angle.Angle(self.Atom_List[Master_ID], self.Atom_List[Slave1_ID], self.Atom_List[Slave2_ID], Angle_Eq))
                                if Type == "D":
                                    Master1_ID = int(R_Line[3].split(',')[0])
                                    Master2_ID = int(R_Line[4].split(',')[0])
                                    Slave1_ID = int(R_Line[2].split(',')[0])
                                    Slave2_ID = int(R_Line[5].split(')')[0])
                                    Dihedral_Eq = float(R_Line[-1])
                                    self.Dihedral_List.append(Dihedral.Dihedral(self.Atom_List[Master1_ID],self.Atom_List[Master2_ID], self.Atom_List[Slave1_ID], self.Atom_List[Slave2_ID], Dihedral_Eq))
                        except:
                            Redundant = False
                            Found = True
                            print k
            
                if Found:
                    break
                
                elif Line[0] == "Redundant" and (File_Lines[i+2].strip('\n').split()[1] == "Optimized"):
                    print "Found redundant internal coordinates"
                    Redundant = True
                    j = i + 7
                    while Redundant:
                        R_Line = File_Lines[j]
                        R_Line = R_Line.split()
                        try:
                            if int(R_Line[0].strip('.')) >= 1:
                                j = j + 1
                                Type = R_Line[1].split('(')[0]
                                # Extract Bonds
                                if Type == "B":
                                    Master_ID = int(R_Line[2].split(',')[0])
                                    Slave_ID = int(R_Line[3].split(')')[0])
                                    req = float(R_Line[-1])
                                    self.Atom_List[Master_ID].Bond_List.append(Slave_ID + 1)
                                    self.Atom_List[Slave_ID].Bond_List.append(Master_ID + 1)
                                    self.Bond_List.append( Bond.Bond(self.Atom_List[Master_ID], self.Atom_List[Slave_ID], req))
                                # Extract Angles
                                if Type == "A":
                                    Master_ID = int(R_Line[3].split(',')[0])
                                    Slave1_ID = int(R_Line[2].split(',')[0])
                                    Slave2_ID = int(R_Line[4].split(')')[0])
                                    Angle_Eq = float(R_Line[-1])
                                    self.Angle_List.append( Angle.Angle(self.Atom_List[Master_ID], self.Atom_List[Slave1_ID], self.Atom_List[Slave2_ID], Angle_Eq))
                                if Type == "D":
                                    Master1_ID = int(R_Line[3].split(',')[0])
                                    Master2_ID = int(R_Line[4].split(',')[0])
                                    Slave1_ID = int(R_Line[2].split(',')[0])
                                    Slave2_ID = int(R_Line[5].split(')')[0])
                                    Dihedral_Eq = float(R_Line[-1])
                                    self.Dihedral_List.append(Dihedral.Dihedral(self.Atom_List[Master1_ID],self.Atom_List[Master2_ID], self.Atom_List[Slave1_ID], self.Atom_List[Slave2_ID], Dihedral_Eq))
                        except:
                            Redundant = False
    
                if Line[0] == "CHELPG" and len(Line) == 2 and not self.UnConverged:
                    for j in range(self.N):
                        Chelp_Line = File_Lines[i+2+j].split()
                        index = int(Chelp_Line[0])
                        charge = float(Chelp_Line[3])
                        self.Atom_List[index].Charge = charge
                        print("NOW YOU SEE ME")
                
                if Line[0] == "CARTESIAN" and Line[2] == "(ANGSTROEM)":
                    for j in range(self.N):
                        XYZ_Line = File_Lines[i+2+j].split()
                        Position_Temp = np.array( [float(XYZ_Line[1]) ,float(XYZ_Line[2]) ,float(XYZ_Line[3])], dtype = float)
                        self.Atom_List[j].Position = Position_Temp
            except:
                continue
        print "Redundant Internal Coordinates and ChelpG partial charges Extracted"
        print "Bond_List = ", len(self.Bond_List)
        
        if not run_orca:
            os.chdir('..')
        
        if local:
            os.chdir('..')

        for atom in self.Atom_List:
            b_list = []
            for b_atom in atom.Bond_List:
                b_list.append(self.Get_Atom(b_atom))
            atom.Bond_List = b_list


        
        print "Optimized XYZ Coordinates:\n"
        for Atom_Obj in self.Atom_List:
            # Finds OPLS Types and Classes
            print Atom_Obj.Atom_ID, Atom_Obj.Element, Atom_Obj.Position, sorted([ self.Get_Atom(Atomobj).Element for Atomobj in Atom_Obj.Bond_List ])
            Atom_Obj.Assign_OPLS_ID()
        print "----------------------------------"

        
        # Find all the ring units in the molecule and add them to the ring list
#self.Ring_List = Ring.create_rings(self.Atom_List)
        
        return

    def Run_Dihedral_Scan(self, File_List):
        Energy = []

        for file in File_List:
            file = open(file, 'r')
            file = file.readlines()
            i = 0
            for line in file:
                try:
                    line = line.split()
                    coords = np.array([float(line[1]), float(line[2]), float(line[3])])
                    self.Atom_List[i].Position = coords
                    i += 1
                except:
                    continue
            D_Sys = System.System([self], [1], 100.0, "Dihedral")
            D_Sys.Gen_Rand()
            D_Sys.Write_LAMMPS_Data( Dihedral=True)
            os.system("cp %sin.Dihedral_Energy ./" % Configure.Template_Path)
            lmp = lammps()
            lmp.file("in.Dihedral_Energy")
            Output_File = open("log.lammps", 'r')
            Output_File = Output_File.readlines()
            for Line in Output_File:
                try:
                    if len(Line.split()) == 2 and float(Line.split()[0]) == float(Line.split()[1]):
                        Energy.append(float(Line.split()[0]))
                except:
                    continue

        Energy = np.asarray(Energy)

        return Energy- Energy[0]

    def Scan_Dihedrals(self, Dihedral_Scan_List, local = False):
        i = 0
        for Dihedral_Angle in Dihedral_Scan_List:
            i += 1
            print Dihedral_Angle
            print self.Dihedral_List
            j = 0
            for Dihedral_Obj in self.Dihedral_List:
                Compare = [ Dihedral_Obj.Dihedral_Slave1.Atom_ID, Dihedral_Obj.Dihedral_Master1.Atom_ID, Dihedral_Obj.Dihedral_Master2.Atom_ID, Dihedral_Obj.Dihedral_Slave2.Atom_ID]
                print Compare
                if sorted(Compare) == sorted(Dihedral_Angle):
                    print "Found Dihedral, Dihedral ID = ", j
                    Dihedral_Eq = Dihedral_Obj.Dihedral_Eq
                    Dihedral_Index = j
                j += 1
        
            print Dihedral_Index
            # Write Orca Input File
            File_Name = self.Name + "_Dih_%d.inp" % i
            File = open(File_Name, 'w')
            #File.write('! RKS B3LYP 6-31+G** NormalSCF Opt NOSOSCF CHELPG PAL8\n\n')
            File.write('! RKS B3LYP 6-31+G** NormalSCF Opt NOSOSCF CHELPG\n\n')
            File.write('%scf MaxIter 500 end\n')
            File.write('%geom Scan\n')
            if abs(Dihedral_Eq) <= 180.0 and abs(Dihedral_Eq) >=170:
                File.write('\tD %d %d %d %d = %.1f, %.1f, 37\n' % ( Dihedral_Angle[0]-1, Dihedral_Angle[1]-1, Dihedral_Angle[2]-1, Dihedral_Angle[3]-1, Dihedral_Eq, -Dihedral_Eq))
                Dihedral_Angles = np.arange(-180,190,10)
                print Dihedral_Angles
            elif abs(Dihedral_Eq) >= 0.0 and abs(Dihedral_Eq) <= 10.0:
                File.write('\tD %d %d %d %d = %.1f, %.1f, 37\n' % ( Dihedral_Angle[0]-1, Dihedral_Angle[1]-1, Dihedral_Angle[2]-1, Dihedral_Angle[3]-1, Dihedral_Eq, 360.0))
                Dihedral_Angles = np.arange(0,370,10)
                print Dihedral_Angles
            else:
                print "Error: Edit Molecule.py/Scan_Dihedrals()"
                time.sleep(600)
            File.write('end end\n\n')
            File.write('*xyz 0 1\n')
            for Atom_Obj in self.Atom_List:
                File.write('%s %.5f %.5f %.5f\n' % ( Atom_Obj.Element, Atom_Obj.Position[0], Atom_Obj.Position[1],Atom_Obj.Position[2]))
            File.write('*')
            
            File.close()
            Finished = False
            File_Out = self.Name + "_Dih_%d.out" % i
            Dihedral_Name = self.Name + "_Dih_%d" % i

            if local:
                #Run subprocess on local machine
                print "Runnning Orca on local machine"
                os.system('pwd')
                os.system('cd Orca')
                os.system('mkdir ./%s' % Dihedral_Name)
                os.system('mv %s ./%s' % (File_Name, Dihedral_Name))
                os.system('cp %s ./%s' % (File_Out, Dihedral_Name))
                os.chdir('./%s' % Dihedral_Name )
                try:
                    File = open(File_Out,'r')
                except:
                    os.system('/Users/andrewkleinschmidt/Library/Orca/orca %s > %s' %(File_Name, File_Out)) # Run Orca Job

            else:
                print "Running Orca Dihedral Scan %d on Comet" % i
                cmd = "mkdir " + Configure.Comet_Path % self.Name + "/Dihedral_%d" % i
                subprocess.call(["ssh", Configure.Comet_Login, cmd])
                subtemp = Configure.Template_Path + "sub_orca_temp"
                submit = "submit_orca"
                Path = Configure.Comet_Path % self.Name + "/Dihedral_%d" % i
                # Write submit script
                with open(subtemp) as f:
                    template = f.read()
                    s = template.format(Comet_Path=Path, Orca_Path = Configure.Orca_Path, name = Dihedral_Name )
                with open(submit ,"w") as f:
                    f.write(s)
            
                # Copy Files over  to Comet
                Dihedral_Path = self.Name + "/Dihedral_%d" % i
            
                os.system( Configure.c2c % (submit, Dihedral_Path))
                os.system( Configure.c2c % (File_Name, Dihedral_Path))
                # Run job
                os.system( Configure.c2l % (Dihedral_Path, File_Out))
                try:
                    File = open(File_Out,'r')
                except:
                    subprocess.call(["ssh",Configure.Comet_Login, Configure.SBATCH % (Dihedral_Path, submit)])
            
            
                # Continuously check to see if  is finished
                j = 0
                while not Finished:
                    os.system( Configure.c2l % (Dihedral_Path, File_Out))
                    try:
                        File = open(File_Out,'r')
                        File_Lines = File.readlines()
                        print File_Lines[-1]
                        if File_Lines[-1].split(' ')[0] == "TOTAL" or self.UnConverged:
                            Finished = True
                        else:
                            print "Not Finished"
                            j += 10
                            print "Sleeping process", j, "Minutes"
                            time.sleep(600)
                    except:
                        print "Sleeping process", j, "miniutes"
                        time.sleep(600)
                        j  += 10

                try:
                    os.mkdir("Dihedral_%d" % i)
                    os.chdir("Dihedral_%d" % i)
                except:
                    os.chdir("Dihedral_%d" % i)
            
                # Copy XYZ trajectory to working directory
                XYZ_File = self.Name +  "_Dih_%d" % i + ".*.xyz"
            

                os.system( Configure.c2l % (Dihedral_Path, XYZ_File))
            
            XYZ_File = self.Name +  "_Dih_%d" % i + ".*.xyz"
            File_List = glob.glob(XYZ_File)
            Traj_File = open('Traj.xyz', 'w')
            # Concatenate files
            for File in File_List:
                File_Obj = open(File, 'r')
                File_Lines = File_Obj.readlines()
                for Line in File_Lines:
                    Traj_File.write(Line)
                print
                index = int(File.split('.')[1].strip('0'))
                print File_Lines[-1]
                print index
                if index < 37 and File_Lines[-1] != ">\n":
                    Traj_File.write('>\n')
                        
            Traj_File.close()
            MP2_Name = "MP2"
            
            if local == True:
                print "Running MP2 on local machine"
                os.system('pwd')
                os.system('cp %s ./' % (Configure.Template_Path + "MP2.inp"))
                File_Name = "MP2.inp"
                File_Out = "MP2.out"
                try:
                    File = open(File_Out,'r')
                    print "Found File"
                except:
                    print "submitted job"
                    os.system('/Users/andrewkleinschmidt/Library/Orca/orca %s > %s' %(File_Name, File_Out)) # Run Orca Job
                    
            else:
                # Prepare new submit script
                with open(subtemp) as f:
                    template = f.read()
                    s = template.format(Comet_Path=Path, Orca_Path = Configure.Orca_Path, name= MP2_Name)
                with open(submit, "w") as f:
                    f.write(s)

                # Copy over to comet
                os.system( Configure.c2c % (submit, Dihedral_Path))
                os.system( Configure.c2c % ("Traj.xyz", Dihedral_Path))
                os.system( Configure.c2c % (Configure.Template_Path + "MP2.inp", Dihedral_Path))
                File_Out = "MP2.out"
                os.system( Configure.c2l % (Dihedral_Path, File_Out))

                # Run Job
                try:
                    File = open(File_Out,'r')
                except:
                    subprocess.call(["ssh",Configure.Comet_Login, Configure.SBATCH % (Dihedral_Path, submit)])

                # Continuously check to see if job is finished
                j = 0
                Finished = False
                while not Finished:
                    os.system( Configure.c2l % (Dihedral_Path, File_Out))
                    try:
                        File = open(File_Out,'r')
                        File_Lines = File.readlines()
                        print File_Lines[-1]
                        if File_Lines[-1].split(' ')[0] == "TOTAL":
                            print "Found Last Line"
                            Finished = True
                        else:
                            print "Not Finished"
                            j += 10
                            print "Sleeping process", j, "Minutes"
                            time.sleep(600)
                    except:
                        print "Sleeping process", j, "miniutes"
                        time.sleep(600)
                        j  += 10
            
            #Copy back over output and store the results for the MP2 relaxed energy scan
            
            File = open("MP2.out", 'r')
            Found = False
            Dihedral_MP2_Energy = np.zeros(37, dtype=float)
            i = 0
            for line in File:
                print line.split()
                if Found and i < 37:
                    Dihedral_MP2_Energy[i] = float(line.split()[1])
                    print Dihedral_MP2_Energy
                    i += 1
            
                
                if line.split() == ['The', 'Calculated', 'Surface', 'using', 'the', 'MP2', 'energy']:
                    Found = True
                

            print Dihedral_MP2_Energy
            Dihedral_MP2_Energy = (Dihedral_MP2_Energy - Dihedral_MP2_Energy[0])*630
            OPLS_Energy = self.Run_Dihedral_Scan(self, File_List)
            
            """
            plt.xlim((Dihedral_Angles[0], Dihedral_Angles[-1]))
            plt.plot(Dihedral_Angles, Dihedral_MP2_Energy, label = "MP2")
            
            plt.ylim((min(Dihedral_MP2_Energy), max(Dihedral_MP2_Energy)))
            plt.ylabel('Relative Energy (kcal/mol)')
            plt.xlabel('Dihedral Angle (Degrees)')
            plt.legend()
            plt.show()
            
            plt.plot(Dihedral_Angles, OPLS_Energy, label = "OPLS")
            plt.ylabel('Relative Energy (kcal/mol)')
            plt.xlabel('Dihedral Angle (Degrees)')
            plt.xlim((Dihedral_Angles[0], Dihedral_Angles[-1]))
            plt.legend()
            plt.show()
            """
            Dihedral_Corrected = Dihedral_MP2_Energy - OPLS_Energy
            Dihedral_Radians = Dihedral_Angles*(3.14/180.0)
            self.Dihedral_List[Dihedral_Index].Fit_Parameters( Dihedral_Corrected, Dihedral_Radians)
          
            """
            plt.plot(Dihedral_Angles, OPLS_Energy, label = "OPLS")
            plt.ylabel('Relative Energy (kcal/mol)')
            plt.xlabel('Dihedral Angle (Degrees)')
            plt.xlim((Dihedral_Angles[0], Dihedral_Angles[-1]))
            plt.legend()
            plt.show()
            """

            os.chdir("..")
                

        return


    def Get_Atom(self,Atom_ID):
        ID = -1
        for i,atom in enumerate(self.Atom_List):
            if atom.Atom_ID == Atom_ID:
                ID = i
                break
        return self.Atom_List[i]

    def Get_Bond(self,Bond_ID):
        ID = -1
        for i,bond in enumerate(self.Bond_List):
            if bond.Bond_ID == Bond_ID:
                ID = i
                break
        return self.Bond_List[i]

    def Unwrap_Coordinates(self,Cutoff,Box_Size):
    #moves atoms to be adjacent if separated by periodic boundary conditions. If atoms are farther apart than Cutoff in any dimension, they are moved by Box_Size
        for atom1 in self.Atom_List:
            for atom2 in self.Atom_List:
                if abs(atom1.Position[0] - atom2.Position[0]) > Cutoff:
                    if atom1.Position[0] > atom2.Position[0]:
                        while abs(atom1.Position[0] - atom2.Position[0]) > Cutoff:
                            atom1.Position[0] += -1*Box_Size
                    else:
                        while abs(atom1.Position[0] - atom2.Position[0]) > Cutoff:
                            atom2.Position[0] += -1*Box_Size
                if abs(atom1.Position[1] - atom2.Position[1]) > Cutoff:
                    if atom1.Position[1] > atom2.Position[1]:
                        while abs(atom1.Position[1] - atom2.Position[1]) > Cutoff:
                            atom1.Position[1] += -1*Box_Size
                    else:
                        while abs(atom1.Position[1] - atom2.Position[1]) > Cutoff:
                            atom2.Position[1] += -1*Box_Size
                if abs(atom1.Position[2] - atom2.Position[2]) > Cutoff:
                    if atom1.Position[2] > atom2.Position[2]:
                        while abs(atom1.Position[2] - atom2.Position[2]) > Cutoff:
                            atom1.Position[2] += -1*Box_Size
                    else:
                        while abs(atom1.Position[2] - atom2.Position[2]) > Cutoff:
                            atom2.Position[2] += -1*Box_Size
# Functions Operating on sets of Molecule objects

def Assign_Lammps(Moltemp_List):
    """
        Function that inputs a list of molecule templates. It searches through all the atoms, bonds, angles etc. to find the unique types of interactions
        present in an arbitrary system object made  of molecule templates.
        
        returns a list of unique params for writing to a LAMMPS data file. these are lists defined such that the i-1 element corresponds to the ith LAMMPS_Type
    """
    
    print "Finding unique Atoms"
    Unique_Atoms = []
    Atom_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Atom_Obj in Moltemp_Obj.Atom_List:
            i = 1
            for Type in Unique_Atoms:
                if Atom_Obj.OPLS_Type != Type:
                    i += 1
                if Atom_Obj.OPLS_Type == Type:
                    Atom_Obj.LAMMPS_Type = i
            if i > len(Unique_Atoms):
                Unique_Atoms.append(Atom_Obj.OPLS_Type)
                Atom_Params.append([Atom_Obj.Mass, Atom_Obj.Sigma, Atom_Obj.Epsilon])
                Atom_Obj.LAMMPS_Type = i
                print "Found new atom type:", i, "OPLS_ID =", Atom_Obj.OPLS_Type, Atom_Params[i-1]
                print "Mass is", Atom_Obj.Mass, "Element is", Atom_Obj.Element

    print "Finding unique Bonds"
    Bond_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Bond_Obj in Moltemp_Obj.Bond_List:
            i = 1
            for Params in Bond_Params:
                if Params[0] != Bond_Obj.kb or Params[1] != Bond_Obj.req:
                    i += 1
                if Params[0] == Bond_Obj.kb and Params[1] == Bond_Obj.req:
                    Bond_Obj.LAMMPS_Type = i
            if i > len(Bond_Params):
                Bond_Params.append([Bond_Obj.kb, Bond_Obj.req])
                Bond_Obj.LAMMPS_Type = i
                print "Found new bond type:", i, Bond_Obj.kb, Bond_Obj.req

    print "Finding unique Angles"
    Angle_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Angle_Obj in Moltemp_Obj.Angle_List:
            i = 1
            for Params in Angle_Params:
                if Params[0] != Angle_Obj.ka or Params[1] != Angle_Obj.Angle_Eq:
                    i += 1
                if Params[0] == Angle_Obj.ka and Params[1] == Angle_Obj.Angle_Eq:
                    Angle_Obj.LAMMPS_Type = i
            if i > len(Angle_Params):
                Angle_Params.append([Angle_Obj.ka, Angle_Obj.Angle_Eq])
                Angle_Obj.LAMMPS_Type = i
                print "Found new angle type:", i, Angle_Params[i-1]


    print "Finding unique Dihedrals"
    Dihedral_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Dihedral_Obj in Moltemp_Obj.Dihedral_List:
            i = 1
            for Params in Dihedral_Params:
                if len(Dihedral_Obj.Coeffs) != len(Params[0]):
                    i += 1
                else:
                    Match_Check = 0
                    for j in range(len(Dihedral_Obj.Coeffs)):
                        if Dihedral_Obj.Coeffs[j] != Params[0][j]:
                            Match_Check = 1
                            break
                    if Match_Check == 1:
                        i+=1
                    else:
                        Dihedral_Obj.LAMMPS_Type = i
                """if Dihedral_Obj.Coeffs[0] != Params[0] or Dihedral_Obj.Coeffs[1] != Params[1] or Dihedral_Obj.Coeffs[2] != Params[2] or Dihedral_Obj.Coeffs[3] != Params[3] or Dihedral_Obj.Coeffs[4] != Params[4] or Dihedral_Obj.Coeffs[5] != Params[5] or Dihedral_Obj.Coeffs[6] != Params[6] or Dihedral_Obj.Coeffs[7] != Params[7] or Dihedral_Obj.Coeffs[8] != Params[8]:
                    i += 1
                if Dihedral_Obj.Coeffs[0] == Params[0] and Dihedral_Obj.Coeffs[1] == Params[1] and Dihedral_Obj.Coeffs[2] == Params[2] and Dihedral_Obj.Coeffs[3] == Params[3] and Dihedral_Obj.Coeffs[4] == Params[4] and Dihedral_Obj.Coeffs[5] == Params[5] and Dihedral_Obj.Coeffs[6] == Params[6] and Dihedral_Obj.Coeffs[7] == Params[7] and Dihedral_Obj.Coeffs[8] == Params[8]:
                    Dihedral_Obj.LAMMPS_Type = i"""
            if i > len(Dihedral_Params):
                Dihedral_Params.append((Dihedral_Obj.Coeffs,Dihedral_Obj.Style))
                Dihedral_Obj.LAMMPS_Type = i
                print "Found new dihedral type:", i, Dihedral_Params[i-1]


    print "Finding unique Impropers"
    Improper_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Improper_Obj in Moltemp_Obj.Improper_List:
            i = 1
            for Params in Improper_Params:
                if Improper_Obj.Ki != Params[0] or Improper_Obj.Improper_Eq != Params[1]:
                    i += 1
                else:
                    Improper_Obj.LAMMPS_Type = i
            try:
                if i > len(Improper_Params):
                    Improper_Params.append([Improper_Obj.Ki, Improper_Obj.Improper_Eq])
                    Improper_Obj.LAMMPS_Type = i
                    print "Found improper type:", i, Improper_Obj.Ki, Improper_Obj.Improper_Eq
            except:
                continue
        


    return Atom_Params, Bond_Params, Angle_Params, Dihedral_Params, Improper_Params

def exponential(L,P):
    return np.exp(L/P)









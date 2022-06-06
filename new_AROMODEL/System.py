#! usr/bin/python

import Molecule
import random
import numpy as np
from copy import deepcopy
import os
import subprocess
import time
import Configure
import Parallel
import pickle
import Conjugated_Polymer
import math
# Class defining an MD system for simulation with LAMMPS


class System(object):
    """
    Class defining an MD system for simulation with LAMMPS
    instance variables:
        Molecule_List = List of Molecule objects
        Composition_List = List of integers defining the number of each molecule type in the system
        Box_Size = float (Angstrom)
        Atom_Params = [[Mass, Sigma, Epsilon]] List of floats
        Bond_Params = [[KB, RO]] List of floats
        Angle_Params = [[KA, A0]] List of floats
        Dihedral_Params = [[V1,V2,V3,V4]] List of floats
        Improper_Params = [[KI, AI]] List of floats
        Num_Atoms = float
        
    """

    def __init__(self, Moltemp_List, Composition_List, Box_Size, Name):
        self.Name = Name
    
        self.Moltemp_List = Moltemp_List
        self.Atom_Params, self.Bond_Params, self.Angle_Params, self.Dihedral_Params, self.Improper_Params = Molecule.Assign_Lammps(Moltemp_List)
        self.Composition_List = Composition_List
        self.Box_Size = Box_Size
        Num_Mol = 0
        for Comp in self.Composition_List:
            Num_Mol += Comp
        print "Number of Molecules = ", Num_Mol
        self.Molecule_List = []
        self.Current_Restart = ""
        self.Temperature = 800
        self.PQEq = False
        self.Hybrid_Dih = False
        

        return

    def Run_Packmol(self,Name="Test"):
        for Molecule_Obj in self.Moltemp_List:
            FileName_P = '%s_Packmol.xyz' % (Molecule_Obj.Name)
            if len(FileName_P) > 30:
                FileName_P = FileName_P[:26]
                FileName_P = FileName_P + ('.xyz')
            File = open(FileName_P, 'w')
            File.write('%d\n' % len(Molecule_Obj.Atom_List))
            for Atom_Obj in Molecule_Obj.Atom_List:
                File.write('\n%s\t%.4f\t%.4f\t%.4f' % (Atom_Obj.Element,Atom_Obj.Position[0],Atom_Obj.Position[1],Atom_Obj.Position[2]))
            File.close()
        Output_Name = ''
        for i in range(len(self.Moltemp_List)):
            Output_Name = Output_Name + ('%s_%d' % (self.Moltemp_List[i].Name,self.Composition_List[i]))
        if len(Output_Name) > 30:
            Output_Name = Output_Name[:25]
        File = open(Output_Name + '.inp', 'w')
        File.write('tolerance 4.0\noutput %s.xyz\nfiletype xyz\n' % Output_Name)
        for z,Molecule_Obj in enumerate(self.Moltemp_List):
            if Name != "Test":
                Molecule_Obj.Write_XYZ(Name = Name + "_" + str(z) + '_Packmol')
            if z == 0:
                if Name == "Test":
                    File.write('structure %s\n\tnumber %d\n\tcenter\nfixed %d. %d. %d. 0. 0. 0.\nend structure\n' % (Molecule_Obj.Name + '_Packmol.xyz',self.Composition_List[z],int(self.Box_Size/2),int(self.Box_Size/2),int(self.Box_Size/2)))
                else:
                    File.write('structure %s\n\tnumber %d\n\tcenter\nfixed %d. %d. %d. 0. 0. 0.\nend structure\n' % (Name + '_' + str(z) + '_Packmol.xyz',self.Composition_List[z],int(self.Box_Size / 2),int(self.Box_Size / 2), int(self.Box_Size / 2)))
            else:
                if Name == "Test":
                    File.write('structure %s\n\tnumber %d\n\tinside cube 0. 0. 0. %d.0\nend structure\n' % (Molecule_Obj.Name + '_Packmol.xyz',self.Composition_List[z],int(self.Box_Size)-2))
                else:
                    File.write('structure %s\n\tnumber %d\n\tinside cube 0. 0. 0. %d.0\nend structure\n' % (Name + '_' + str(z) + '_Packmol.xyz', self.Composition_List[z], int(self.Box_Size) - 2))
        File.close()
        os.system('packmol < %s.inp' % Output_Name)
        Mol_ID = 1
        File = open('%s.xyz' % Output_Name, 'r')
        Coords = File.readlines()
        Atom_Counter = 2
        for i,Molecule_Obj in enumerate(self.Moltemp_List):
            for j in range(self.Composition_List[i]):
                Temp_Mol = deepcopy(Molecule_Obj)
                print("Depositing %d" % Mol_ID)
                Temp_Mol.Mol_ID = Mol_ID
                for k in range(len(Molecule_Obj.Atom_List)):
                    Temp_Mol.Atom_List[k].Position[0] = Coords[Atom_Counter].split()[1]
                    Temp_Mol.Atom_List[k].Position[1] = Coords[Atom_Counter].split()[2]
                    Temp_Mol.Atom_List[k].Position[2] = Coords[Atom_Counter].split()[3]
                    Atom_Counter += 1
                Mol_ID += 1
                self.Molecule_List.append(Temp_Mol)
        #equilibrate charge
        Q = np.longdouble(0.0)
        Num_Atoms = 0
        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                Q += Atom_Obj.Charge
                if Atom_Obj.Charge != 0.0:
                    Num_Atoms += 1.0
        print "The total charge of the system is ", Q

        dQ = np.longdouble(Q / Num_Atoms)
        print "dQ =", dQ

        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                if Atom_Obj.Charge != 0.0:
                    Atom_Obj.Charge -= dQ

    def Assign_PQEq(self):
        self.PQEq_Params = []
        for Atom_Type in Atom_Params:
            self.PQEq_Params.append(Atom.Find_PQEq_Params(Atom_Type[0]))
        self.PQEq = True

    def Gen_Rand(self):
        i = 0
        k = 1
        for Molecule_Obj in self.Moltemp_List:
            print "Genenerated", self.Composition_List[i], Molecule_Obj.Name, "Molecules"
            for j in range(self.Composition_List[i]):
                Temp_Mol = deepcopy(Molecule_Obj)
                #Temp_Mol = Molecule_Obj
                
               
                Temp_Mol.COM = np.asarray([random.random()*self.Box_Size, random.random()*self.Box_Size, random.random()*self.Box_Size], dtype = float)
                Temp_Mol.Mol_ID = k
                k += 1
                Temp_Mol.Adjust_COM()
                self.Molecule_List.append(Temp_Mol)
                
                """Deposited = False
                while not Deposited:
                    
                    Temp_Mol.COM = np.asarray([random.random()*(self.Box_Size-5)+5, random.random()*(self.Box_Size-5)+5, random.random()*(self.Box_Size-5)+5], dtype = float)
                    Temp_Mol.Adjust_COM()
                    Temp_Mol.Mol_ID = k
                    
                    
                    
                    if len(self.Molecule_List) == 0:
                        Temp_Mol.COM = np.asarray([self.Box_Size/2.0, self.Box_Size/2.0, self.Box_Size/2.0], dtype = float)
                        Temp_Mol.Adjust_COM()
                        print "No Overlap"
                        self.Molecule_List.append(Temp_Mol)
                        Deposited = True
                        k+= 1
                    
                
                    e = 0
                    for Mol_Obj in self.Molecule_List:
                        Distance = np.linalg.norm(Temp_Mol.COM - Mol_Obj.COM)
                        if Distance > 3:
                            e += 1

                        else:
                            print "Overlap"
                            print Distance
                            break

        
                    if e == len(self.Molecule_List):
                        Temp_Mol.Adjust_COM()
                        self.Molecule_List.append(Temp_Mol)
                        print "Depositing", j
                        Deposited = True
                        k += 1
                    else:
                        print "Didn't Deposit
                    #Deposited = True
                    for mol in self.Molecule_List:
                        if np.linalg.norm(Temp_Mol.COM - mol.COM < 100):
                            for atom1 in Temp_Mol.Atom_List:
                                for atom2 in mol.Atom_List:
                                    if np.linalg.norm(atom1.Position - atom2.Position) < 4:
                                        Deposited = False
                                        print(atom1.Position)
                                        print(atom2.Position)
                                    if not Deposited:
                                        break
                                if not Deposited:
                                    break
                        if not Deposited:
                            break
                    if Deposited:
                        self.Molecule_List.append(deepcopy(Temp_Mol))
                        k += 1
                    else:
                        print("Overlap")
                
        
            i += 1"""
        
    
        Q = np.longdouble(0.0)
        Num_Atoms = 0
        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                Q += Atom_Obj.Charge
                if Atom_Obj.Charge != 0.0:
                    Num_Atoms += 1.0
        print "The total charge of the system is ", Q
        
        dQ = np.longdouble(abs(Q/Num_Atoms))
        print "dQ =", dQ

        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                if Atom_Obj.Charge < 0.0:
                    Atom_Obj.Charge -= dQ
                if Atom_Obj.Charge > 0.0:
                    Atom_Obj.Charge -= dQ

        return

    def Write_LAMMPS_Data(self, Dihedral = False, Data_File = ""):
        """
        Function for writing LAMMPS data file 
        """
        if Data_File == "":
            Data_File = self.Name + ".data"
        if Dihedral:
            self.Data_File = "Dihedral.data"
        else:
            self.Data_File = Data_File

        File = open(self.Data_File, 'w')
        File.write('LAMMPS data file via System.Write_LAMMPS_Data()\n\n')
        
        # Find number of atoms, bonds, dihedrals, impropers in the system
        self.Num_Atoms = 0
        self.Num_Bonds = 0
        self.Num_Angles = 0
        self.Num_Dihedrals = 0
        self.Num_Impropers = 0
        i = 0
        for Moltemp_Obj in self.Moltemp_List:
            self.Num_Atoms += len(Moltemp_Obj.Atom_List)*self.Composition_List[i]
            self.Num_Bonds += len(Moltemp_Obj.Bond_List)*self.Composition_List[i]
            self.Num_Angles += len(Moltemp_Obj.Angle_List)*self.Composition_List[i]
            self.Num_Dihedrals += len(Moltemp_Obj.Dihedral_List)*self.Composition_List[i]
            self.Num_Impropers += len(Moltemp_Obj.Improper_List)*self.Composition_List[i]
            i += 1
        
        
        File.write('%d atoms\n' % self.Num_Atoms)
        File.write('%d atom types\n' % len(self.Atom_Params))
        File.write('%d bonds\n' % self.Num_Bonds)
        File.write('%d bond types\n' % len(self.Bond_Params))
        File.write('%d angles\n' % self.Num_Angles)
        File.write('%d angle types\n' % len(self.Angle_Params))
        if self.Num_Dihedrals > 0:
            File.write('%d dihedrals\n' % self.Num_Dihedrals)
            File.write('%d dihedral types\n' % len(self.Dihedral_Params))
    
        if self.Num_Impropers > 0:
            File.write('%d impropers\n' % self.Num_Impropers)
            File.write('%d improper types\n' % len(self.Improper_Params))

        File.write('\n\n0.0000 %.4f xlo xhi\n' % self.Box_Size)
        File.write('0.0000 %.4f ylo yhi\n' % self.Box_Size)
        File.write('0.0000 %.4f zlo zhi\n' % self.Box_Size)

        File.write('\n\nMasses\n\n')
        i = 1
        for Params in self.Atom_Params:
            File.write('%d %.3f\n' % ( i, Params[0]))
            i += 1

        File.write('\n\nPair Coeffs # lj/cut/coul/long\n\n')
        i = 1
        if self.PQEq:
            File.write('* * coul/pqeqgauss 0.0000 0.0000 0.000 0 0.00 0.000 0.0 0 #dummy')
            for Atom_Param,PQEq_Param in self.Atom_Params,self.PQEq_Params:
                File.write('%d lj/cut/coul/long %.3f %.3f\n' % (i,Atom_Param[2],Atom_Param[3]))
                #File.write('%d %d coul/pqeqgauss %.4f %.4f %.3f %d %.2f %.3f %.1f %d' % )
        for Params in self.Atom_Params:

            File.write('%d %.3f %.3f\n' % (i, Params[2], Params[1]))
            i += 1

        File.write('\n\nBond Coeffs # harmonic\n\n')
        i = 1
        for Params in self.Bond_Params:
            File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
            i += 1

        File.write('\n\nAngle Coeffs # harmonic\n\n')
        i = 1
        for Params in self.Angle_Params:
            File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
            i += 1

        if self.Num_Dihedrals > 0:
            for Param in self.Dihedral_Params:
                if Param[1] != self.Dihedral_Params[0][1]:
                    self.Hybrid_Dih = True
            File.write('\n\nDihedral Coeffs # opls\n\n')
            i = 1
            for Params in self.Dihedral_Params:
                """if self.Hybrid_Dih:
                    File.write('%d %s' % (i, Params[1]))
                else:"""
                File.write("%d" % i)
                for j in range(len(Params[0])):
                    File.write(' %.4f' % Params[0][j])
                File.write('\n')
                i += 1

        if self.Num_Impropers > 0:
            File.write('\n\nImproper Coeffs # harmonic\n\n')
            i = 1
            for Params in self.Improper_Params:
                File.write('%d %.4f -1 2\n' % (i, Params[0]))
                i += 1

        File.write('\n\nAtoms # full\n\n')
        i = 1
        j = 1
        for Molecule_Obj in self.Molecule_List:
            for Atom_Obj in Molecule_Obj.Atom_List:
                Atom_Obj.System_ID = i
                File.write('%d %d %d %.12f %.4f %.4f %.4f\n' % (Atom_Obj.System_ID, Molecule_Obj.Mol_ID, Atom_Obj.LAMMPS_Type, Atom_Obj.Charge, Atom_Obj.Position[0], Atom_Obj.Position[1], Atom_Obj.Position[2]))
                i += 1
            j += 1
        

        File.write('\n\nBonds\n\n')
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Bond_Obj in Molecule_Obj.Bond_List:
                Bond_Obj.System_ID = i
                File.write('%d %d %d %d\n' % ( Bond_Obj.System_ID, Bond_Obj.LAMMPS_Type, Bond_Obj.Bond_Master.System_ID, Bond_Obj.Bond_Slave.System_ID))
                i += 1


        File.write('\n\nAngles\n\n')
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Angle_Obj in Molecule_Obj.Angle_List:
                Angle_Obj.System_ID = i
                File.write('%d %d %d %d %d\n' % (Angle_Obj.System_ID, Angle_Obj.LAMMPS_Type, Angle_Obj.Angle_Slave1.System_ID, Angle_Obj.Angle_Master.System_ID,  Angle_Obj.Angle_Slave2.System_ID))
                i += 1
        if self.Num_Dihedrals > 0:
            File.write('\n\nDihedrals\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Dihedral_Obj in Molecule_Obj.Dihedral_List:
                    Dihedral_Obj.System_ID = i
                    File.write('%d %d %d %d %d %d\n' % (Dihedral_Obj.System_ID, Dihedral_Obj.LAMMPS_Type, Dihedral_Obj.Dihedral_Slave1.System_ID, Dihedral_Obj.Dihedral_Master1.System_ID, Dihedral_Obj.Dihedral_Master2.System_ID, Dihedral_Obj.Dihedral_Slave2.System_ID))
                    i += 1

        if self.Num_Impropers > 0:
            File.write('\n\nImpropers\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Improper_Obj in Molecule_Obj.Improper_List:
                    Improper_Obj.System_ID = i
                    File.write('%d %d %d %d %d %d\n' % (Improper_Obj.System_ID, Improper_Obj.LAMMPS_Type, Improper_Obj.Improper_Master.System_ID, Improper_Obj.Improper_Slave1.System_ID, Improper_Obj.Improper_Slave2.System_ID, Improper_Obj.Improper_Slave3.System_ID))
                    i += 1

    def Write_LAMMPS_Data_Imp_Only(self, Dihedral = False, Fold = -1):
        """
        Function for writing LAMMPS data file 
        """
        if Dihedral:
            self.Data_File = "Dihedral.data"
        else:
            self.Data_File = self.Name + ".data"

        File = open(self.Data_File, 'w')
        File.write('LAMMPS data file via System.Write_LAMMPS_Data()\n\n')
        
        # Find number of atoms, bonds, dihedrals, impropers in the system
        self.Num_Atoms = 0
        self.Num_Bonds = 0
        self.Num_Angles = 0
        self.Num_Dihedrals = 0
        self.Num_Impropers = 0
        i = 0
        for Moltemp_Obj in self.Moltemp_List:
            self.Num_Atoms += len(Moltemp_Obj.Atom_List)*self.Composition_List[i]
            self.Num_Bonds += len(Moltemp_Obj.Bond_List)*self.Composition_List[i]
            self.Num_Angles += len(Moltemp_Obj.Angle_List)*self.Composition_List[i]
            self.Num_Dihedrals += len(Moltemp_Obj.Dihedral_List)*self.Composition_List[i]
            self.Num_Impropers += len(Moltemp_Obj.Improper_List)*self.Composition_List[i]
            i += 1
        
        
        File.write('%d atoms\n' % self.Num_Atoms)
        File.write('%d atom types\n' % len(self.Atom_Params))
        File.write('%d bonds\n' % self.Num_Bonds)
        File.write('%d bond types\n' % len(self.Bond_Params))
        if Fold == -1:
            File.write('0 angles\n')
            File.write('0 angle types\n')
        else:
            File.write('1 angles\n')
            File.write('1 angle types\n')
        if self.Num_Dihedrals > 0:
            File.write('0 dihedrals\n')
            File.write('0 dihedral types\n')
    
        if self.Num_Impropers > 0:
            File.write('%d impropers\n' % self.Num_Impropers)
            File.write('%d improper types\n' % len(self.Improper_Params))

        File.write('\n\n0.0000 %.4f xlo xhi\n' % self.Box_Size)
        File.write('0.0000 %.4f ylo yhi\n' % self.Box_Size)
        File.write('0.0000 %.4f zlo zhi\n' % self.Box_Size)

        File.write('\n\nMasses\n\n')
        i = 1
        for Params in self.Atom_Params:
            File.write('%d %.3f\n' % ( i, Params[0]))
            i += 1

        """File.write('\n\nPair Coeffs # lj/cut/coul/long\n\n')
        i = 1
        for Params in self.Atom_Params:
            File.write('%d %.3f %.3f\n' % (i, Params[2], Params[1]))
            i += 1"""

        File.write('\n\nBond Coeffs # harmonic\n\n')
        i = 1
        for Params in self.Bond_Params:
            File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
            i += 1

        if Fold != -1:
            File.write('\n\nAngle Coeffs # harmonic\n\n')
            i = 1
            File.write('%d %.4f %.4f\n' % (i, 12000.0000, 0.0000)) # Its possible that im missing a factor of 2 here
            i += 1

        if self.Num_Impropers > 0:
            File.write('\n\nImproper Coeffs # harmonic\n\n')
            i = 1
            for Params in self.Improper_Params:
                File.write('%d %.4f -1 2\n' % (i, Params[0]))
                i += 1

        File.write('\n\nAtoms # full\n\n')
        i = 1
        j = 1
        for Molecule_Obj in self.Molecule_List:
            for Atom_Obj in Molecule_Obj.Atom_List:
                Atom_Obj.System_ID = i
                File.write('%d %d %d %.8f %.4f %.4f %.4f\n' % (Atom_Obj.System_ID, Molecule_Obj.Mol_ID, Atom_Obj.LAMMPS_Type, Atom_Obj.Charge, Atom_Obj.Position[0], Atom_Obj.Position[1], Atom_Obj.Position[2]))
                i += 1
            j += 1
        

        File.write('\n\nBonds\n\n')
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Bond_Obj in Molecule_Obj.Bond_List:
                Bond_Obj.System_ID = i
                File.write('%d %d %d %d\n' % ( Bond_Obj.System_ID, Bond_Obj.LAMMPS_Type, Bond_Obj.Bond_Master.System_ID, Bond_Obj.Bond_Slave.System_ID))
                i += 1

        if Fold != -1:
            File.write('\n\nAngles\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                File.write('%d %d %d %d %d\n' % (1, 1, Molecule_Obj.Atom_List[0].Atom_ID, Fold,  Molecule_Obj.Atom_List[-1].Atom_ID))
                i += 1

        if self.Num_Impropers > 0:
            File.write('\n\nImpropers\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Improper_Obj in Molecule_Obj.Improper_List:
                    Improper_Obj.System_ID = i
                    File.write('%d %d %d %d %d %d\n' % (Improper_Obj.System_ID, Improper_Obj.LAMMPS_Type, Improper_Obj.Improper_Master.System_ID, Improper_Obj.Improper_Slave1.System_ID, Improper_Obj.Improper_Slave2.System_ID, Improper_Obj.Improper_Slave3.System_ID))
                    i += 1

    def Run_Lammps_Soft(self, Nodes = 1, GPU = False):
        cmd = "mkdir " + Configure.Comet_Path % self.Name
        subprocess.call(["ssh", Configure.Comet_Login, cmd])

        In_Temp = Configure.Template_Path + "in.init_temp_soft"
        In_File = "in.soft_%s" % self.Name
        Sim_Name = "soft_%s" % self.Name
        Pair_Num = 1
        Pair_Coeffs = ''
        for Params in self.Atom_Params:
            Pair_Coeffs = Pair_Coeffs + ('pair_coeff %d %d %.3f %.3f\n' % (Pair_Num, Pair_Num, Params[2], Params[1]))
            Pair_Num += 1
        with open(In_Temp) as f:
            template = f.read()
        s = template.format(System_Name = self.Name,Pair_Params = Pair_Coeffs)
        with open(In_File,'w') as f:
            f.write(s)
        """
        Current Rule of thumb for Parallel Job Submission:
        MPI: 1 processor per 1000 atoms
        MPI + GPU: 4 GPU and 24 Processors --> Only use for >100K particles
        Need to do more benchmarks with current system --> good task for CJ
        """
        
        if not GPU:
            sub_temp = Configure.Template_Path + "sub_Lammps"
        else:
            sub_temp = Configure.Template_Path + "GPU_Sub"
        NProcs = 24
        submit = "sub_%s" % self.Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)

        File_Out1 = 'log.%s' % Sim_Name
        File_Out = 'restart.%s_Soft' % self.Name

        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (In_File, self.Name))
        os.system( Configure.c2c % (self.Data_File, self.Name))
        os.system( Configure.c2l % (self.Name, File_Out1))

        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])

        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10
        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        return

    def Run_Lammps_Init(self, Nodes = 1, GPU = False, Generate_Models = False, Minimize_Only = False, Compress = False, Condense = False, Shared = False, Implicit = False, PQEQ = False):
        cmd = "mkdir " + Configure.Comet_Path % self.Name
        subprocess.call(["ssh", Configure.Comet_Login, cmd])
        
        sections = 1

        # Set up input file
        if Generate_Models:
            sections = 3
            if Compress:
                In_Temp = Configure.Template_Path + "in.generate_models_compressed_initial"
            elif Implicit:
                In_Temp = Configure.Template_Path + "in.generate_models_implicit"
            else:
                In_Temp = Configure.Template_Path + "in.generate_models_initial"
        elif Minimize_Only:
            In_Temp = Configure.Template_Path + "in.minimize_fold"
        elif Compress:
            In_Temp = Configure.Template_Path + "in.init_temp_compress"
        elif Condense:
            In_Temp = Configure.Template_Path + "in.init_temp_condense"
        else:
            In_Temp = Configure.Template_Path + "in.init_temp"
        In_File = "in.init_%s" % self.Name
        Sim_Name = "init_%s" % self.Name
        for section in range(0,sections):
            if Generate_Models:
                if section == 1:
                    In_Temp = Configure.Template_Path + "in.generate_models_middle"
                if section == 2:
                    In_Temp = Configure.Template_Path + "in.generate_models_final"
            Dih_Coeff = ''
            Dih_Num = 1
            if self.Hybrid_Dih:
                for Params in self.Dihedral_Params:
                    Dih_Coeff = Dih_Coeff + ('dihedral_coeff %d %s' % (Dih_Num, Params[1]))
                    for j in range(len(Params[0])):
                        Dih_Coeff = Dih_Coeff + (' %.4f' % Params[0][j])
                    Dih_Coeff = Dih_Coeff + '\n'
                    Dih_Num += 1
            seed = int(random.random() * 10000)
            File_Out1 = 'log.%s' % Sim_Name
            if Condense:
                File_Out = 'restart.%s_Condensed_%d' % (self.Name, section)
            else:
                File_Out = 'restart.%s_800_%d_1' % (self.Name, section)
            with open(In_Temp) as f:
                template = f.read()
            s = template.format(System_Name = self.Name,Dihedral_Coeffs = Dih_Coeff, Seed = seed,Restart_Out = File_Out, Restart_In = self.Current_Restart)
            with open(In_File,'w') as f:
                f.write(s)
            """
            Current Rule of thumb for Parallel Job Submission:
            MPI: 1 processor per 1000 atoms
            MPI + GPU: 4 GPU and 24 Processors --> Only use for >100K particles
            Need to do more benchmarks with current system --> good task for CJ
            """
            

            if GPU:
                sub_temp = Configure.Template_Path + "GPU_Sub"
            elif Shared:
                sub_temp = Configure.Template_Path + "sub_Shared"
            else:
                sub_temp = Configure.Template_Path + "sub_Lammps"
            if Shared:
                NProcs = 8
            else:
                NProcs = 24
            submit = "sub_%s" % self.Name

            with open(sub_temp) as f:
                template = f.read()
                s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
            with open(submit,'w') as f:
                f.write(s)


            # Copy over to Comet
            os.system( Configure.c2c % (submit, self.Name))
            os.system( Configure.c2c % (In_File, self.Name))
            os.system( Configure.c2c % (self.Data_File, self.Name))
            os.system( Configure.c2l % (self.Name, File_Out1))

            try:
                File = open(File_Out1,'r')
            except:
                subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])

            Finished = False
            i = 0
            while not Finished:
                os.system( Configure.c2l % (self.Name, File_Out))
                try:
                    File = open(File_Out,'r')
                    Finished = True
                except:
                    print "Sleeping process", i, "minutes"
                    if not Minimize_Only:
                        time.sleep(600)
                    else:
                        time.sleep(30)
                    i += 10
            if Minimize_Only:
                os.system( Configure.c2l % (self.Name, "Minimized_Structure.data"))
            os.system( 'rm %s' % File_Out)
            self.Current_Restart = File_Out
        return

    def Run_Lammps_NPT(self, GPU = False, Implicit_Solv = False, Generate_Models = False, Compress = False, Temp_Out = 0.0, time_steps = 1000000, Nodes = 1):
        """
            Function for running NPT dynamics for the system in lammps
            
        """
        Temp_In = self.Temperature
        try:
            count = int(self.Current_Restart.split('_')[-1])
        except:
            count = 1
        count += 1
        sections = 1
        if Implicit_Solv:
            NPT_Temp = Configure.Template_Path + "in.NPT_Temp_Solv"    
        elif Generate_Models:
            if Compress:
                NPT_Temp = Configure.Template_Path + "in.generate_models_compressed"
                sections = 3
            else:
                NPT_Temp = Configure.Template_Path + "in.generate_models"
                sections = 3
        else:
            NPT_Temp = Configure.Template_Path + "in.NPT_Temp"
        for section in range(0,sections):
            if sections > 1:
                NPT_In = "in.NPT_%s_%d_%d_%d" % (self.Name, count, Temp_In, sections)
            else:
                NPT_In = "in.NPT_%s_%d_%d" % (self.Name, count, Temp_In)
            if Generate_Models:
                if section == 0:
                    NPT_Temp = Configure.Template_Path + 'in.generate_models_initial'
                elif section == 1:
                    NPT_Temp = Configure.Template_Path + 'in.generate_models_middle'
                elif section == 2:
                    NPT_Temp = Configure.Template_Path + 'in.generate_models_final'
            Sim_Name = NPT_In.split('.')[-1]
            if Temp_Out != 0.0:
                NPT_In += "_%d" % Temp_Out
                Sim_Name = NPT_In.split('.')[-1]
                self.Temperature = Temp_Out
            if Temp_Out == 0.0:
                Temp_Out = Temp_In
            if sections == 1:
                New_Restart = 'restart.%s_%d_%d' % (self.Name, Temp_Out, count)
            else:
                New_Restart = 'restart.%s_%d_%d' % (self.Name, Temp_Out, section)
            Dih_Coeff = ''
            Dih_Num = 1
            if self.Hybrid_Dih:
                for Params in self.Dihedral_Params:
                    Dih_Coeff = Dih_Coeff + ('dihedral_coeff %d %s' % (Dih_Num, Params[1]))
                    for j in range(len(Params[0])):
                        Dih_Coeff = Dih_Coeff + (' %.4f' % Params[0][j])
                    Dih_Coeff = Dih_Coeff + '\n'
                    Dih_Num += 1

            with open(NPT_Temp) as f:
                template = f.read()
            s = template.format(Name = self.Name, count = count, Coeffs = Dih_Coeff, Temp_In = Temp_In, Temp_Out = Temp_Out, Restart_In = self.Current_Restart, Restart_Out = New_Restart, Steps = time_steps)
            with open(NPT_In,'w') as f:
                f.write(s)
                #f.insert(33, 'fix def1 all print 100 "${temperature} ${volume} ${dens} ${Enthalpy}" append Thermo_{Temp_In}_{Temp_Out}.txt screen no')
            
            
            if not GPU:
                sub_temp = Configure.Template_Path + "sub_Lammps"
            else:
                sub_temp = Configure.Template_Path + "GPU_Sub"
            NProcs = 4
            submit = "sub_%s" % self.Name
            with open(sub_temp) as f:
                template = f.read()
                s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
            with open(submit,'w') as f:
                f.write(s)

            File_Out1 = 'log.%s' % Sim_Name
            File_Out = New_Restart
            Traj_File = self.Name + "_%d.lammpstrj" % count

            # Copy over to Comet
            os.system( Configure.c2c % (submit, self.Name))
            os.system( Configure.c2c % (NPT_In, self.Name))
            os.system( Configure.c2l % (self.Name, File_Out1))
            
            try:
                File = open(File_Out1,'r')
            except:
                subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])
            
            Finished = False
            i = 0
            while not Finished:
                os.system( Configure.c2l % (self.Name, File_Out))
                try:
                    File = open(File_Out,'r')
                    Finished = True
                except:
                    print "Sleeping process", i, "minutes"
                    time.sleep(600)
                    i += 10

            os.system( 'rm %s' % File_Out)
            self.Current_Restart = File_Out
            os.system( Configure.c2l % (self.Name,  Traj_File))
        return

    def Run_Lammps_NVT(self, Old_Restart, GPU = False, Temp_Out = 0.0, time_steps = 1000000, Nodes = 1):
        """
            Function for running NPT dynamics for the system in lammps
            
        """
        cmd = "mkdir " + Configure.Comet_Path % self.Name
        subprocess.call(["ssh", Configure.Comet_Login, cmd])
        Temp_In = self.Temperature
        self.Current_Restart = Old_Restart
        count = int(self.Current_Restart.split('_')[-1])
        count += 1
        NVT_Temp = Configure.Template_Path + "in.NVT_Temp"
        NVT_In = "in.NVT_%s_%d_%d" % (self.Name, count, Temp_In)
        Sim_Name = NVT_In.split('.')[-1]
        if Temp_Out != 0.0:
            NVT_In += "_%d" % Temp_Out
            Sim_Name = NVT_In.split('.')[-1]
            self.Temperature = Temp_Out
        if Temp_Out == 0.0:
            Temp_Out = Temp_In
        New_Restart = 'restart.%s_%d_%d' % (self.Name, Temp_Out, count)
    

        with open(NVT_Temp) as f:
            template = f.read()
        s = template.format(Name = self.Name, count = count, Temp_In = Temp_In, Temp_Out = Temp_Out, Restart_In = self.Current_Restart, Restart_Out = New_Restart, Steps = time_steps)
        with open(NVT_In,'w') as f:
            f.write(s)
            #f.insert(33, 'fix def1 all print 100 "${temperature} ${volume} ${dens} ${Enthalpy}" append Thermo_{Temp_In}_{Temp_Out}.txt screen no')
        
        
        if not GPU:
            sub_temp = Configure.Template_Path + "sub_Lammps"
        else:
            sub_temp = Configure.Template_Path + "GPU_Sub"
        NProcs = 24
        submit = "sub_%s" % self.Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)

        File_Out1 = 'log.%s' % Sim_Name
        File_Out = New_Restart
        Traj_File = self.Name + "_%d.lammpstrj" % count

        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (NVT_In, self.Name))
        os.system( Configure.c2c % (Old_Restart, self.Name))
        os.system( Configure.c2l % (self.Name, File_Out1))
        
        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])
        
        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10

        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        os.system( Configure.c2l % (self.Name,  Traj_File))
        return

    def Run_Lammps_Strain(self, Restart_In, Restart_Out, Index, Nodes = 1):
        " Function for running uniaxial straining simulations on the system in LAMMPS"

        Strain_Temp = Configure.Template_Path + "in.strain_Temp"
        Strain_In = "in.%s_strain_%s" % (self.Name, Index)
        Sim_Name = Strain_In.split('.')[-1]
        
        # Prepare input script
        with open(Strain_Temp) as f:
            Template = f.read()
        s = Template.format( index = Index, Restart_In= Restart_In, Restart_Out = Restart_Out)
        with open(Strain_In, 'w') as f:
            f.write(s)

        # Prepare Submission Script
        sub_temp = Configure.Template_Path + "sub_Lammps"
        NProcs = 24
        submit = "sub_%s_strain" % self.Name
        with open(sub_temp) as f:
            template = f.read()
        s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)

        File_Out1 = 'log.%s' % Sim_Name
        File_Out = Restart_Out
        Traj_File1 = "DAStrain_%d.lammpstrj" % Index
        Strain_File1 = "SS_%d.txt" %  Index

        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (Strain_In, self.Name))
        os.system( Configure.c2l % (self.Name, File_Out1))

        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])

        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10

        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        os.system( Configure.c2l % (self.Name,  Traj_File1))
        os.system( Configure.c2l % (self.Name,  Strain_File1))
        
        
        
        return

    def Change_Ring_Positions(self):
        #for Conjugated_Polymer class only
        for mol in self.Molecule_List:
            for ring in mol.Ring_List:
                for r_atom in ring.Atom_List:
                    for atom in mol.Atom_List:
                        if atom.Atom_ID == r_atom.Atom_ID:
                            r_atom.Position = atom.Position
                            break

    def Update_Atom_Positions(self,File_Name,Image_Flags=False):
        #Updates Atom positions based on input file
        f = open(File_Name,'r')
        lines = f.readlines()
        Read_Flag = False
        for line in lines:
            if len(line.split()) > 0 and line.split()[0].strip() == 'Atoms':
                Read_Flag = True
            elif Read_Flag and len(line.split()) == 7 and not Image_Flags:
                #temp_atom = self.Get_Atom_System(int(line.split()[0].strip()))
                temp_atom = self.Molecule_List[(int(line.split()[0].strip())-1)//len(self.Molecule_List[0].Atom_List)].Atom_List[(int(line.split()[0].strip()) - 1) % len(self.Molecule_List[0].Atom_List)] #Gets the number atom in the proper molecule for a system with all identical molecules. // is floor division
                temp_atom.Position[0] = float(line.split()[4].strip())
                temp_atom.Position[1] = float(line.split()[5].strip())
                temp_atom.Position[2] = float(line.split()[6].strip())
            elif Read_Flag and len(line.split()) == 10 and Image_Flags:
                #temp_atom = self.Get_Atom_System(int(line.split()[0].strip()))
                temp_atom = self.Molecule_List[(int(line.split()[0].strip())-1)//len(self.Molecule_List[0].Atom_List)].Atom_List[(int(line.split()[0].strip()) - 1) % len(self.Molecule_List[0].Atom_List)] #Gets the number atom in the proper molecule for a system with all identical molecules. // is floor division
                temp_atom.Position[0] = float(line.split()[4].strip())
                temp_atom.Position[1] = float(line.split()[5].strip())
                temp_atom.Position[2] = float(line.split()[6].strip())
                temp_atom.Image_Flags[0] = float(line.split()[7].strip())
                temp_atom.Image_Flags[1] = float(line.split()[8].strip())
                temp_atom.Image_Flags[2] = float(line.split()[9].strip())
            elif len(line.split()) > 0 and line.split()[0] == 'Bonds':
                Read_Flag = False

    def Get_Atom_System(self,Atom_ID):
        mol_atom_ID = -1
        mol_ID = -1
        found = False
        for j,mol in enumerate(self.Molecule_List):
            for i,atom in enumerate(mol.Atom_List):
                if atom.System_ID == Atom_ID:
                    mol_atom_ID = i
                    mol_ID = j
                    found = True
                    break
            if found:
                break
        return self.Molecule_List[mol_ID].Atom_List[mol_atom_ID]

    def Assign_System_ID(self):
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Atom_Obj in Molecule_Obj.Atom_List:
                Atom_Obj.System_ID = i
                i += 1
                for Nested_Mol in self.Molecule_List:
                    found = False
                    for Nested_Ring in Nested_Mol.Ring_List:
                        for Nested_Atom in Nested_Ring.Plumed_Rings[0]:
                            if Nested_Atom.Position[0] == Atom_Obj.Position[0] and Nested_Atom.Position[1] == Atom_Obj.Position[1] and Nested_Atom.Position[2] == Atom_Obj.Position[2]:
                                Nested_Atom.System_ID = Atom_Obj.System_ID
                                found = True
                                break
                        for Nested_Atom in Nested_Ring.Plumed_Rings[1]:
                            if found:
                                break
                            if Nested_Atom.Position[0] == Atom_Obj.Position[0] and Nested_Atom.Position[1] == Atom_Obj.Position[1] and Nested_Atom.Position[2] == Atom_Obj.Position[2]:
                                Nested_Atom.System_ID = Atom_Obj.System_ID
                                found = True
                                break

                        if found:
                            break
                    if found:
                        found = False
                        break

    def Unwrap_Image_Flags(self):
        for mol in self.Molecule_List:
            for atom in mol.Atom_List:
                atom.Position[0] += self.Box_Size * atom.Image_Flags[0]
                atom.Position[1] += self.Box_Size * atom.Image_Flags[1]
                atom.Position[2] += self.Box_Size * atom.Image_Flags[2]

    def Write_Plumed_Files(self,Write_Flags,Base_File_Name,Fit_Energies,Force_y,Force_x,Delocalization_Modifier=1.0,Nonbonded_Modifier=1.0,Non_Interacting=False,Reverse=False,Folder_Name="./Full_Polymer_Plumed",Bias=True):
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

        for flag,mol in zip(Write_Flags,self.Molecule_List):
            plumed_index = 0
            if flag:
                oop_blocks = np.linspace(0,math.pi,200)
                dih_blocks = np.linspace(-math.pi,math.pi,200)

                if Non_Interacting == True:
                    for atom in mol.Atom_List:
                        atom.Charge = 0.000000
                        atom.Sigma = 0.0
                        atom.Epsilon = 0.0




                torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (index*2,mol.Ring_List[0].Plumed_Rings[plumed_index][0].System_ID))
                for atom in mol.Ring_List[0].Core_Atom_List[1:]:
                    torsion_file.write(",%d" % atom.System_ID)
                torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c1,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (index*2,mol.Ring_List[0].Plumed_Rings[plumed_index][0].System_ID,mol.Ring_List[0].Plumed_Rings[plumed_index][1].System_ID))

                Parameter_Index = 0
                for i in range(0,len(mol.Ring_List)-1):
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
                        f.write("#! FIELDS DIH_%d OOP_%d ext_%d.bias der_DIH_%d der_OOP_%d\n#! SET min_DIH_%d -pi\n#! SET max_DIH_%d pi\n#! SET nbins_DIH_%d 199\n#! SET periodic_DIH_%d true\n#! SET min_OOP_%d 0\n#! SET max_OOP_%d pi\n#! SET nbins_OOP_%d 199\n#! SET periodic_OOP_%d false\n" % (index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1,index-1))
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
                        os.system("scp %s.dat ./%s" % (Bias_File_Name,Folder_Name))
                        os.system("rm -f %s.dat" % (Bias_File_Name))
                        if Non_Interacting:
                            os.system("scp %s_Nonbonded.dat ./%s" % (Bias_File_Name,Folder_Name))
                            os.system("rm -f %s_Nonbonded.dat" % (Bias_File_Name))
                    torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (index*2-3, mol.Ring_List[i].Plumed_Rings[0][0].System_ID)) #Writes the first atom in the first plumed ring
                    for atom in mol.Ring_List[i].Plumed_Rings[0][1:]:
                        torsion_file.write(",%d" % atom.System_ID)
                    torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (index*2-2, mol.Ring_List[i+1].Plumed_Rings[1][0].System_ID)) #Writes the first atom in the second plumed ring
                    for atom in mol.Ring_List[i+1].Plumed_Rings[1][1:]:
                        torsion_file.write(",%d" % atom.System_ID)
                    for b_atom in mol.Ring_List[i].Bonded_Atoms:
                        if b_atom.Is_Linked and b_atom.Interring_Bond_Atom.Self_Ring.Ring_ID == i+2:
                            atom1 = b_atom.Interring_Bond_Atom.Central_Atom.System_ID
                            atom2 = b_atom.Central_Atom.System_ID
                            ring1_batom1 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].System_ID
                            ring1_batom2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].System_ID
                            ring2_batom1 = b_atom.Same_Ring_Bonded_Atom_List[0].System_ID
                            ring2_batom2 = b_atom.Same_Ring_Bonded_Atom_List[1].System_ID
                    torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c%d,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (index*2-3,index*2-3,mol.Ring_List[i].Plumed_Rings[0][0].System_ID,mol.Ring_List[i].Plumed_Rings[0][1].System_ID))
                    torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c%d,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (index*2-2,index*2-2,mol.Ring_List[i+1].Plumed_Rings[1][0].System_ID,mol.Ring_List[i+1].Plumed_Rings[1][1].System_ID))
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
                    torsion_file.write("PRINT ARG=DIH_%d,OOP_%d FILE=/scratch/andrewk/job_SLURM_JOBID/colvar_%s%d.txt STRIDE=100\n\n" % (index-1,index-1,Conventional_Add,index-1))
                    if Bias:
                        torsion_file.write("EXTERNAL ARG=DIH_%d,OOP_%d FILE=/scratch/andrewk/job_SLURM_JOBID/%s.dat LABEL=ext_%d\n\n" % (index-1,index-1,Bias_File_Name,index-1))
                    if Non_Interacting:
                        torsion_file.write("EXTERNAL ARG=DIH_%d,OOP_%d FILE=/scratch/andrewk/job_SLURM_JOBID/%s_Nonbonded.dat LABEL=ext_%d_nb\n\n" % (i,i,Bias_File_Name,i))
                    Parameter_Index += 1
                    if Parameter_Index == len(Fit_Energies)/2:
                        Parameter_Index = 0

                os.system("mkdir ./%s" % (Folder_Name))

                index += 1


            torsion_file.write("WHOLEMOLECULES")
            for q in range(len(self.Molecule_List)):
                torsion_file.write(" ENTITY%d=%d-%d" % (q,q*len(self.Molecule_List[q].Atom_List)+1,q*len(self.Molecule_List[q].Atom_List) + len(self.Molecule_List[q].Atom_List)))
                for i in range(q*len(self.Molecule_List[q].Ring_List)*2+1, q*len(self.Molecule_List[q].Ring_List)*2+2*len(self.Molecule_List[q].Ring_List)-1):
                    torsion_file.write(",c%d,c%d_normal" % (i, i))
            Num_Atoms += len(mol.Atom_List)
            Num_Rings += len(mol.Ring_List)

        torsion_file.close()
        os.system("scp %s.dat ./%s" % (Base_Calculate_Torsions_String, Folder_Name))
        os.system("rm -f %s.dat" % Base_Calculate_Torsions_String)
        for i in range(1,Num_Rings):
            if Non_Interacting:
                os.system("scp %s_Bias_File_%d_Nonbonded.dat ./%s" % (Base_File_Name,i,Folder_Name))
                os.system("rm -f %s_Bias_File_%d_Nonbonded.dat" % (Base_File_Name,i))


def Run_Glass_Transition(system, Interval, Ramp_Steps = 50000, Equil_Steps = 50000, T_End = 100, Nodes = 1):
    T_init = system.Temperature
    Range = abs(T_init - T_End)
    Steps = Range/Interval
    T_out = T_init
    for i in range(Steps):
        T_in = T_out
        T_out = T_out + Interval
        system.Run_Lammps_NPT(Temp_Out= T_out, time_steps = Ramp_Steps,Nodes=Nodes)
        system.Run_Lammps_NPT( time_steps = Equil_Steps, Nodes=Nodes)
        File_Out = "Thermo_%d_%d" % (T_out, T_out) + ".txt"
        os.system( Configure.c2l % (system.Name, File_Out))
    return


def Uniaxial_Strain(system, Restart_In, steps = 20):
    Restart_In_Temp = Restart_In
    for i in range(steps):
        index = i+1
        Restart_Out = "restart." + system.Name + "_Strain_%d" % (i+1)
        system.Run_Lammps_Strain(Restart_In_Temp, Restart_Out, index )
        Restart_In_Temp = Restart_Out

    return










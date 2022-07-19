import numpy as np
from sympy import Mod
import Atom
import Molecule
import Ring
import OPLS as op
import System
import Configure
import Conjugated_Polymer
import Cluster_IO
import Write_Inputs
import Write_Submit_Script
import math
import copy
import scipy
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib as mpl
from symfit import parameters, variables, sin, cos, Fit
import re
import time
import dill

def Rotate_Ring(Planar_Polymer,Rotation_List):
    #A function where rotates one of the rings in the multimer through an improper or dihedral rotation
    return

def Run_LAMMPS_No_Inter(Rotated_Polymer):
    #Runs lammps locally with no interring interactions (e.g. LJ forces)
    return

def Run_Quantum(Rotated_Polymer):
    #Runs on the supercomputer quantum mechanical simulations
    return

def Run_Nonbonded(Rotated_Polymer):
    #Runs on the supercomputer with one monomer methylated to block conjugation effects
    return

def CVFF(x,K,d,n):
    #A CVFF function for evaluating or fitting an improper twist
    return K*(1+d*np.cos(np.pi/180*n*x))

def Return_LAMMPS_No_Inter(Rotated_Polymer):
    #Get back information from a local Lammps run with no interring forces, e.g. LJ
    return

def Return_Quantum(Rotated_Polymer):
    #Get back information from a quantum mechanical run about different types of interactions
    return

def OPLS(x,a,b,c,d,e):
    #Returns an OPLS measurement or fitting, potentially offset in y-direction, of how the rings move
    return .5*a*(1+np.cos(np.radians(x)))+.5*b*(1-np.cos(np.radians(2*x)))+.5*c*(1+np.cos(np.radians(3*x)))+.5*d*(1-np.cos(np.radians(4*x)))+e

def OPLS_LAMMPS(x,a,b,c,d):
    #Returns an OPLS measurement or fitting with *no* y-direction offset
    return .5*a*(1+np.cos(np.radians(x)))+.5*b*(1-np.cos(np.radians(2*x)))+.5*c*(1+np.cos(np.radians(3*x)))+.5*d*(1-np.cos(np.radians(4*x)))

def OPLS_No_Offset(x,a,b,c,d,e,f,g,h):
    #Returns OPLS function without y-offset *but* with x-offset
    return .5*a*(1+np.cos(np.radians(x)-e))+.5*b*(1-np.cos(np.radians(2*x)-f))+.5*c*(1+np.cos(np.radians(3*x)-g))+.5*d*(1-np.cos(np.radians(4*x)-h))

def OPLS_Derivative(x,a,b,c,d):
    #Returns the first derivative of the OPLS function
    return .5*a*(-np.sin(np.radians(x))) + b*(np.sin(np.radians(2*x))) + 1.5*c*(-np.sin(np.radians(3*x))) + 2*d*(np.sin(np.radians(4*x)))

def Fourier(x,a,c,e,g):
    #Returns an offset OPLS series but without a 4th order part, with y-offset
    return (a*(1+np.cos(np.radians(x))) + c*(1+np.cos(np.radians(2*x))) + e*(1+np.cos(np.radians(3*x)))) + g

def Harmonic(x,a,b,c,d,e,f,g):
    return (a*np.cos(x)+b*np.cos(x)**2+c*np.cos(x)**3+d*np.cos(x)**4+e*np.cos(x)**5+f*np.cos(x)**6+g)

def Cosine_Harmonic(x,a,b,c):
    return a*(1+b*np.cos(np.radians(x))) + c

def Sine_Harmonic(x,a,b,c):
    return a*(1+b*np.sin(np.radians(x))) + c

def Gaussian_And_Sigmoid(x,a,b,c,d,e,f,g,h,i):
    return a*math.e**(-(x-b)**2/c) + d*(1-np.power((x-e)*f,g))/(1-np.power((x-e)*f,h)) + i

def Multi_Gaussian(x,a,b,c,d,e,f,g,h,i,j):
    return a*math.e**(-(x-b)**2/c) + d*math.e**(-(x-e)**2/f) + g*math.e**(-(x-h)**2/i) + j

def Log(x,a,b,c,d,e):
    return a*(np.log(b*x+c)/np.log(d))+e

def Exp(x,a,b,c,d,e):
    return a*np.power(b,(c*x+d))+e

def Power(x,a,b,c,d):
    return a*np.power(x-b,int(c))+d

def Sigmoid(x,a,b,c,d,e,f):
    return a*(1-np.power((x-b)*c,d))/(1-np.power((x-b)*c,e)) + f

def Run_Parameters(Folder_Name):
    #Returns regularly the current lab parameters for the supercomputer (e.g. location)
    Cluster_Login = Configure.orca_dict["Cluster_Login"]
    Base_Cluster_Location = Configure.orca_dict["Base_Cluster_Location"]
    Cluster_Location = Base_Cluster_Location + "/" + Folder_Name
    Scheduler_Type = Configure.orca_dict["Scheduler_Type"]
    End_Condition = Configure.orca_dict["End_Condition"]
    Shared_File_Location = Configure.orca_dict["Shared_File_Location"]
    qos = "overrun --time-min=00:30:00"

    return Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos

def Run_Parameters_LAMMPS(Folder_Name):
    #What supercomputer (if any) to run LAMMPS on
    Cluster_Login = Configure.lammps_dict["Cluster_Login"]
    Base_Cluster_Location = Configure.lammps_dict["Base_Cluster_Location"]
    Cluster_Location = Base_Cluster_Location + "/" + Folder_Name
    Scheduler_Type = Configure.lammps_dict["Scheduler_Type"]
    End_Condition = Configure.lammps_dict["End_Condition"]
    Shared_File_Location = Configure.lammps_dict["Shared_File_Location"]
    qos = "overrun --time-min=00:30:00"

    return Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos

def Return_Filenames(Base_Name,Reversed_Name=""):
    #Returns a list of filenames for the supercomputer for a given base name and possibly reversed name
    End_File = Base_Name + ".out"
    Reversed_End_File = Reversed_Name + ".out"
    Job_Name = Base_Name
    Reversed_Job_Name = Reversed_Name
    In_File = Base_Name + ".qcin"
    Reversed_In_File = Reversed_Name + ".qcin"
    Sub_File = "sub_" + Base_Name
    Reversed_Sub_File = "sub_" + Reversed_Name

    return End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File

def fourier_series(x, f, n=0):
    """
    Returns a symbolic fourier series of order `n`.

    :param n: Order of the fourier series.
    :param x: Independent variable
    :param f: Frequency of the fourier series
    """
    # Make the parameter objects for all the terms
    a0, cos_a = parameters(','.join(['a{}'.format(i) for i in range(0, n + 1)]))[0],parameters(','.join(['a{}'.format(i) for i in range(0, n + 1)]))[1:]
    sin_b = parameters(','.join(['b{}'.format(i) for i in range(1, n + 1)]))
    # Construct the series
    series = a0 + sum(ai * cos(i * f * x) + bi * sin(i * f * x)
                     for i, (ai, bi) in enumerate(zip(cos_a, sin_b), start=1))
    return series

def fourier_nonfit(x, a_params, b_params, a0, n=0):
    """
    Returns a symbolic fourier series of order `n`.

    :param n: Order of the fourier series.
    :param x: Independent variable
    :param f: Frequency of the fourier series
    """
    # Make the parameter objects for all the terms
    series = a0
    for a,b,i in zip(a_params,b_params,range(1,n+1)):
        series += a * cos(i * x) + b * sin(i *  x)
    return series

def fourier_series_derivative(x, a_params, b_params, n=0):
    """
    Returns a symbolic fourier series of order `n`.

    :param n: Order of the fourier series.
    :param x: Independent variable
    :param f: Frequency of the fourier series
    """
    # Make the parameter objects for all the terms

    # Construct the series
    series = 0
    for a,b,i in zip(a_params,b_params,range(1,n+1)):
        series += -1 * a * i * sin(i * x) + b * i * cos(i * x)
    return series

def Make_Surface_Plot(x_axis,y_axis,Surface_Data,File_Name,Title="",xlabel="",ylabel="",Tight_Layout = True,vmin=0,vmax=10):
#This function takes in x and y-information nad surface data and prints a pretty plot from it (3-D data)
    fig,ax = plt.subplots(1,1)
    x,y = np.meshgrid(x_axis,y_axis)
    c = ax.pcolor(x,y,Surface_Data,cmap = 'seismic',vmin=vmin,vmax=vmax)
    ax.set_title(Title,fontdict = {'fontsize':24})
    plt.xlabel(xlabel,size = 24)
    plt.ylabel(ylabel,size = 24)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.tick_params(length=4,width=4)
    cbar = fig.colorbar(c,ax=ax)
    cbar.ax.tick_params(labelsize = 20)
    if Tight_Layout:
        plt.tight_layout()
    fig.savefig(File_Name)

    plt.close(fig)

def Read_Input(Input_File,XYZ_File,Polymer_Name):
#This function reads in the parameter file, assigns the atoms to rings, adds available LJ, coulombic, bonded, angular, dihedral, and improper potentials, and tells the program whether it needs to parameterize missing bond potentials or partial charges. Returns Ring_List: NumPy array of Ring objects categorizing all available atoms into separate rings; Paramaterize_Bond: Boolean that equals "True" if bond parameters for interring bonds have not been specified; Paramaterize_Charges: Boolean that equals "True" if partial charges for atoms have not been specified
    Supermolecule = Molecule.Molecule(XYZ_File)
    #Open file and read lines
    f = open(Input_File)
    lines = f.readlines()
    Ring_List = []
    Ring_Read = False
    Bonded_Read = False
    Aux_Read = False
    for line in lines:
        if len(line.strip().split()) > 0 and line.strip().split()[0] == "***":
            Bonded_Read = False
            for b_atoms in Bonded_Atom_List:
                Bonded_Atom_Vectors.append(Supermolecule.Get_Atom(b_atoms[1]).Position - Supermolecule.Get_Atom(b_atoms[0]).Position)
            New_Ring = Ring.Ring(Supermolecule,Ring_Atom_List,Core_Atom_List,Bonded_Atom_List,Bonded_Atom_Vectors,Ring_Name,Ring_ID,Polymer_Name,Ring_Nickname,Symmetric = Symmetry,Aux_Ring_List=Aux_Rings)
            print("Bond Length Check")
            print(len(New_Ring.Bond_List))
            Ring_List.append(New_Ring)
        if Bonded_Read:
            String_Bonded_Atom_List = line.strip().split()
            Bonded_Atom_List.append([int(i) for i in String_Bonded_Atom_List])
        if Aux_Read:
            if Aux_Count <= Aux_Ring_Num:
                Aux_Rings.append([[],[],int(line.strip().split()[1]),int(line.strip().split()[2]),line.strip().split()[0],line.strip().split()[3]])
                Aux_Count += 1
            else:
                Aux_Read = False
                Ring_Read = True
        if Ring_Read:
            if len(line.strip().split()) == 1 or line.strip().split()[-1] == "CORE":
                Ring_Atom_List.append(int(line.strip().split()[0]))
                if line.strip().split()[-1] == "CORE":
                    Core_Atom_List.append(int(line.strip().split()[0]))
            elif line.strip().split()[-1][:3] == "AUX":
                Aux_Rings[int(line.strip().split()[-1][-1])][0].append(int(line.strip().split()[0]))
                if len(line.strip().split()[-1]) > 5:
                    Aux_Rings[int(line.strip().split()[-1][-1])][1].append(int(line.strip().split()[0]))
            elif line.strip().split()[0] == "Bonded":
                Ring_Read = False
                Bonded_Read = True
        if len(line.strip().split()) > 0 and line.strip().split()[0] == "Ring":
            Ring_ID = int(line.strip().split()[1])
            if line.strip().split()[2].strip() == "True":
                Symmetry = True
            else:
                Symmetry = False
            Ring_Name = line.strip().split()[3].strip()
            Ring_Nickname = line.strip().split()[4].strip()
            Aux_Ring_Num = int(line.strip().split()[5].strip())
            if Aux_Ring_Num != 0:
                Aux_Read = True
                Aux_Count = 1
            else:
                Ring_Read = True
            Ring_Atom_List = []
            Bonded_Atom_List = []
            Bonded_Atom_Vectors = []
            Core_Atom_List = []
            Aux_Rings = []
            #if Aux_Ring_Num != 0:

    return Ring_List,False,False

    #Determine if interring bond potentials are specified in the input file and set Parameterize_Bonds

    #Determine if partial charges are specified in the input file and set Parameterize_Charges

    #Create Ring objects and append them to the list

def Set_Up_Nonbonded(Rotated_Polymer,Solvent):
    Density = 1.49
    MW = 119.38
    Solvent = Molecule.Molecule(Solvent)
    Coul_Cutoff = Rotated_Polymer.Find_Coul_Cutoff()
    Coul_Cutoff = round(Coul_Cutoff)
    Coul_Cutoff += 2
    Num_Solv = int(round((2*Coul_Cutoff+3)**3/10*Density/MW*6.022))
    Solvent.Set_Up_FF()
    op.Assign_OPLS(Solvent, ChelpG = False)
    #Rotated_Polymer.Zero_Total_Charge()
    Solvent_System = System.System([Rotated_Polymer,Solvent], [1,Num_Solv], (2*Coul_Cutoff+1), "%s_Nonbonded" % (Rotated_Polymer.Name))
    Solvent_System.Run_Packmol()
    Solvent_System.Write_LAMMPS_Data(Data_File = Rotated_Polymer.Name + "_Nonbonded.data")

    Rotated_Polymer.Write_XYZ()
    return Rotated_Polymer.Name + "_Nonbonded.data"

def Run_SPE_Methyl_Impropers(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name):
    Improper_Run_List = []
    """for ring in Ring_List:
        if ring.Name not in Improper_Run_List:
            Improper_Run_List.append(ring.Name)
            ring.Improper_Bend_Control_Submit()"""

    Phi_Rotation = np.linspace(0,Max_OOP,Rotated_Shape[0])[1]

    for ring in Ring_List:
        if ring.Name not in Improper_Run_List:
            Improper_Run_List.append(ring.Name)
            print(ring.Name)
            ring.Improper_Bend_Methyl_Control_Submit(Rotated_Shape[0],Phi_Rotation)
    
def Run_SPE_Dimers(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name):
    Phi_Rotation = np.linspace(0,Max_OOP,Rotated_Shape[0])[1]
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    Full_Nontorsional_Energy_List = []
    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))

    Run_List = []
    k = 0
    Ring_By_Ring_End_File_Matrices = []
    Ring_By_Ring_Improper_File_Matrices = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Reversed_Dimer = Conjugated_Polymer.Conjugated_Polymer([ring2,ring1])
        XYZ_Filename = Dimer.Write_XYZ() 
        if not os.path.isdir("XYZ_Files"):
            os.mkdir('XYZ_Files')
        os.system("cp %s ./XYZ_Files/" % XYZ_Filename)
        os.remove(XYZ_Filename)
        Reversed_End_File_Matrix = []
        End_File_Matrix = []
        Improper_File_Matrix = []
        Reversed_Improper_File_Matrix = []
        Nontorsional_Energy_Matrix = []
        Reversed_Nontorsional_Energy_Matrix = []
        for j in range(Rotated_Shape[0]):
            In_File_List = []
            Reversed_In_File_List = []
            Copy_File_List = []
            Reversed_Copy_File_List = []
            End_File_List = []
            Improper_File_List = []
            Reversed_Improper_File_List = []
            Reversed_End_File_List = []
            Nontorsional_Energy_List = []
            Reversed_Nontorsional_Energy_List = []
            for i in range(36): #TODO: figure out what 36 is
                time.sleep(2)#TODO:rapid repeated calls to expanse server seems to cause issues, this slows it down. Should look for a better solution in the future
                XYZ_Filename = Dimer.Write_XYZ()
                os.system("cp %s ./XYZ_Files/" % XYZ_Filename)
                os.remove(XYZ_Filename)
                XYZ_Filename = Reversed_Dimer.Write_XYZ()
                os.system("cp %s ./XYZ_Files/" % XYZ_Filename)
                os.remove(XYZ_Filename)
                print()
                End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s" % Dimer.Name,Reversed_Name = "%s" % Reversed_Dimer.Name)
                Job_Type = "QChem"
                Folder_Name = "Rotation_Test"
                #Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)
                #TODO: Confirm this should actually use qchem
                Cluster_Login = Configure.qchem_dict["Cluster_Login"]
                Base_Cluster_Location = Configure.qchem_dict["Base_Cluster_Location"]
                Cluster_Location = Configure.qchem_dict["Cluster_Location"]
                Scheduler_Type = Configure.qchem_dict["Scheduler_Type"]
                End_Condition = Configure.qchem_dict["End_Condition"]
                Shared_File_Location = Configure.qchem_dict["Shared_File_Location"]
                N = 2
                Percent_Change = 0.0
                """if ring1.Symmetric and ring2.Symmetric:
                    Symmetry_End_File = Dimer.Return_Symmetry_Name() + ".out"
                    Check_List = [End_File,Symmetry_End_File]
                else:
                    Symmetry_End_File = ""
                    Check_List = End_File"""
                """Symmetry_End_File = ""
                Check_List = End_File
                if End_File in Run_List:
                    End_File_List.append(End_File_List)
                elif Symmetry_End_File in Run_List:
                    End_File_List.append(Symmetry_End_File)
                if End_File not in Run_List and Symmetry_End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(Check_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = Check_List)
                    if Return_File == "":
                        Return_File = End_File
                    End_File_List.append(Return_File)
                    if not Finished:
                        Copy_File_List.append(In_File)
                        In_File_List.append(In_File)
                        Run_List.append(End_File)
                        if Symmetry_End_File != "":
                            Run_List.append(Symmetry_End_File)
                        Write_Inputs.Write_QChem_SPE(In_File,Dimer)"""
                if End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
                    if Return_File == "":
                        Return_File = End_File
                    End_File_List.append(Return_File)
                    Improper_File_List.append(ring2.Name + "_Improper_Bend_Methyl_Phi_%d.out" % (Max_OOP/(Rotated_Shape[0]-1)*j))
                    Reversed_Improper_File_List.append(ring1.Name + "_Improper_Bend_Methyl_Phi_%d.out" % (Max_OOP/(Rotated_Shape[0]-1)*j))
                    if not Finished:
                        Copy_File_List.append(In_File)
                        In_File_List.append(In_File)
                        Run_List.append(End_File)
                        Write_Inputs.Write_QChem_SPE(In_File,Dimer)

                Dimer.Read_From_Data_File("./LigParGen_Files/%s_%s.lmp" % (ring1.Name,ring2.Name),No_Position_Update = True)
                #Replaced "~/lammps-11Aug17/src/lmp_serial" with Configure.local_dict["Lammps_Path"]
                Nontorsional_Energy_List.append(Dimer.Calculate_Internal_Energy(Polymer_Name,Configure.local_dict["Lammps_Path"],Interring_Angles_Only=True,dielectric=4.9))

                if Reversed_End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(Reversed_End_File,Folder_Name,Reversed_Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = Reversed_End_File)
                    if Return_File == "":
                        Return_File = Reversed_End_File
                    Reversed_End_File_List.append(Return_File)
                    if not Finished:
                        Reversed_Copy_File_List.append(Reversed_In_File)
                        Reversed_In_File_List.append(Reversed_In_File)
                        Run_List.append(Reversed_End_File)
                        Write_Inputs.Write_QChem_SPE(Reversed_In_File,Reversed_Dimer)

                Reversed_Dimer.Read_From_Data_File("./LigParGen_Files/%s_%s.lmp" % (ring2.Name,ring1.Name),No_Position_Update = True)
                #Reversed_Nontorsional_Energy_List.append(Reversed_Dimer.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/whi",Interring_Angles_Only=True,dielectric=4.9))
                Reversed_Nontorsional_Energy_List.append(Reversed_Dimer.Calculate_Internal_Energy(Polymer_Name,Configure.local_dict["Lammps_Path"],Interring_Angles_Only=True,dielectric=4.9))
                """if Symmetry_End_File == "":
                    Dimer.Read_From_Data_File("./LigParGen_Files/%s_%s.lmp" % (ring1.Name,ring2.Name))
                    Nontorsional_Energy_List.append(Dimer.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/lmp_serial",Exclude_Interring_Torsions=True,dielectric=4.9))
                else:
                    if Return_File != Symmetry_End_File or End_File == Symmetry_End_File:
                        Dimer.Read_From_Data_File("./LigParGen_Files/%s_%s.lmp" % (ring1.Name,ring2.Name))
                    Nontorsional_Energy_List.append(Dimer.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/lmp_serial",Exclude_Interring_Torsions=True,dielectric=4.9,Symmetric=True))"""
                Dimer.Rotate_Ring("Dih",10,Dimer.Ring_List[0],Dimer.Ring_List[1])
                Reversed_Dimer.Rotate_Ring("Dih",10,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            if len(In_File_List) != 0:
                End_File = In_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
                Copy_File_List.append(Sub_File)
                Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Copy_File_List:
                    os.system("cp %s ./Rotation_Run_Input_Copies/" % file)
                    os.remove(file)

            if len(Reversed_In_File_List) != 0:
                Reversed_End_File = Reversed_In_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Reversed_Sub_File,Reversed_In_File_List,Reversed_Job_Name,Cluster_Location,Job_Type)
                Reversed_Copy_File_List.append(Reversed_Sub_File)
                Cluster_IO.Submit_Job(Reversed_Copy_File_List,Folder_Name,Reversed_Sub_File,Reversed_End_File,Reversed_Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = Reversed_End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Reversed_Copy_File_List:
                    os.system("cp %s ./Rotation_Run_Input_Copies/" % file)
                    os.remove(file)
            Dimer.Rotate_Ring("OOP",Phi_Rotation,Dimer.Ring_List[0],Dimer.Ring_List[1])
            #Dimer.Rotate_Ring("OOP",5,Dimer.Ring_List[0],Dimer.Ring_List[1])
            Reversed_Dimer.Rotate_Ring("OOP",Phi_Rotation,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            #Reversed_Dimer.Rotate_Ring("OOP",5,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            End_File_Matrix.append(End_File_List)
            Improper_File_Matrix.append(Improper_File_List)
            Reversed_Improper_File_Matrix.append(Reversed_Improper_File_List)
            Reversed_End_File_Matrix.append(Reversed_End_File_List)
            Nontorsional_Energy_Matrix.append(np.array(Nontorsional_Energy_List))
            Reversed_Nontorsional_Energy_Matrix.append(np.array(Reversed_Nontorsional_Energy_List))
            print("taking a quick 20 sec nap to let the server cool down")
            #time.sleep(20)
            #print("nap time over") #TODO: temporary, yeah.

        Ring_By_Ring_End_File_Matrices.append(End_File_Matrix)
        Ring_By_Ring_End_File_Matrices.append(Reversed_End_File_Matrix)
        Ring_By_Ring_Improper_File_Matrices.append(Improper_File_Matrix)
        Ring_By_Ring_Improper_File_Matrices.append(Reversed_Improper_File_Matrix)
        Full_Nontorsional_Energy_List.append(np.array(Nontorsional_Energy_Matrix) - np.amin(np.array(Nontorsional_Energy_Matrix)[0]))
        Full_Nontorsional_Energy_List.append(np.array(Reversed_Nontorsional_Energy_Matrix) - np.amin(np.array(Reversed_Nontorsional_Energy_Matrix)[0]))
        

    """for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Copy_File_List = []
        In_File_List = []
        for i in range(11):
            Job_Type = "QChem"
            Folder_Name = "Rotation_Test"
            End_File = "%s.out" % Dimer.Name
            Cluster_Login = "andrewk@cori.nersc.gov"
            Base_Cluster_Location = '/global/cscratch1/sd/andrewk'
            Cluster_Location='/global/cscratch1/sd/andrewk/Rotation_Test'
            Scheduler_Type = "SLURM"
            End_Condition = "SPE_QChem"
            Shared_File_Location = "/Users/andrewkleinschmidt/Shared_File_Location"
            Job_Name = "%s" % Dimer.Name
            In_File = "%s.qcin" % Dimer.Name
            Sub_File = "sub_%s" % Dimer.Name
            N = 2
            qos = "premium"
            if End_File not in Run_List:
                Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
                if Return_File == "":
                    Return_File = End_File
                if not Finished:
                    Run_List.append(Symmetry_End_File)
                    Write_Inputs.Write_QChem_SPE(In_File,Dimer)
                    In_File_List.append(In_File)
                    Copy_File_List.append(In_File)
            Dimer.Rotate_Ring("OOP",2,Dimer.Ring_List[0],Dimer.Ring_List[1])
        if len(In_File_List) != 0:
            End_File = In_File_List[-1]
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,32,Cluster_Location,Job_Type,walltime = 5,queue = qos,proc_per_node = 32,constraint = 'haswell')
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            for file in Copy_File_List:
                os.system("scp %s ./Rotation_Run_Input_Copies" % file)
                os.system("rm -f %s" % file)"""


    return Ring_By_Ring_End_File_Matrices,Full_Nontorsional_Energy_List,Ring_By_Ring_Improper_File_Matrices

def Run_SPE_Dimers_Hydrogenated(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Alternate = False):
    Phi_Rotation = np.linspace(0,Max_OOP,Rotated_Shape[0])[1]
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    Full_Nontorsional_Energy_List = []
    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))

    Run_List = []
    k = 0
    Ring_By_Ring_End_File_Matrices = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Reversed_Dimer = Conjugated_Polymer.Conjugated_Polymer([ring2,ring1])
        End_File_Matrix = []
        Reversed_End_File_Matrix = []
        for j in range(Rotated_Shape[0]):
            In_File_List = []
            Reversed_In_File_List = []
            Copy_File_List = []
            Reversed_Copy_File_List = []
            End_File_List = []
            Reversed_End_File_List = []
            Nontorsional_Energy_List = []
            Reversed_Nontorsional_Energy_List = []
            for i in range(36):
                if not Alternate:
                    Hydrogenated_Dimer,_ = Dimer.Create_Hydrogenated_Copy(0,1)
                    XYZ_Filename = Hydrogenated_Dimer.Write_XYZ()
                else:
                    Hydrogenated_Dimer = Dimer.Create_Hydrogenated_Copy_Alternate(0,1)
                    XYZ_Filename = Hydrogenated_Dimer.Write_XYZ()

                if not Alternate:
                    Reversed_Hydrogenated_Dimer,_ = Reversed_Dimer.Create_Hydrogenated_Copy(0,1)
                    Reversed_XYZ_Filename = Reversed_Hydrogenated_Dimer.Write_XYZ()
                else:
                    Reversed_Hydrogenated_Dimer = Reversed_Dimer.Create_Hydrogenated_Copy_Alternate(0,1)
                    Reversed_XYZ_Filename = Reversed_Hydrogenated_Dimer.Write_XYZ()

                os.system("cp %s ./Hydrogenated_XYZ_Files/" % XYZ_Filename)
                os.remove(XYZ_Filename)
                if XYZ_Filename != Reversed_XYZ_Filename:
                    os.system("cp %s ./Hydrogenated_XYZ_Files/" % Reversed_XYZ_Filename)
                    os.remove(Reversed_XYZ_Filename)
                if not Alternate:
                    End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated" % Hydrogenated_Dimer.Name,"%s_Hydrogenated" % Reversed_Hydrogenated_Dimer.Name)
                else:
                    End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated_Alternate" % Hydrogenated_Dimer.Name,"%s_Hydrogenated_Alternate" % Reversed_Hydrogenated_Dimer.Name)
                Job_Type = "QChem"
                Folder_Name = "Hydrogenated_Rotation_Test"
                Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

                """if ring1.Symmetric and ring2.Symmetric:
                    Symmetry_End_File = Hydrogenated_Dimer.Return_Symmetry_Name() + "_Hydrogenated.out"
                    Check_List = [End_File,Symmetry_End_File]
                else:
                    Symmetry_End_File = ""
                    Check_List = End_File"""
                if End_File in Run_List:
                    End_File_List.append(End_File)
                if Reversed_End_File in Run_List:
                    Reversed_End_File_List.append(Reversed_End_File)
                """elif Symmetry_End_File in Run_List:
                    End_File_List.append(Symmetry_End_File)"""
                if End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
                    if Return_File == "":
                        Return_File = End_File
                    End_File_List.append(Return_File)
                    if not Finished:
                        Copy_File_List.append(In_File)
                        In_File_List.append(In_File)
                        Run_List.append(End_File)
                        Write_Inputs.Write_QChem_SPE(In_File,Hydrogenated_Dimer)

                if Reversed_End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(Reversed_End_File,Folder_Name,Reversed_Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = Reversed_End_File)
                    if Return_File == "":
                        Return_File = Reversed_End_File
                    Reversed_End_File_List.append(Return_File)
                    if not Finished:
                        Reversed_Copy_File_List.append(Reversed_In_File)
                        Reversed_In_File_List.append(Reversed_In_File)
                        Run_List.append(Reversed_End_File)
                        Write_Inputs.Write_QChem_SPE(Reversed_In_File,Reversed_Hydrogenated_Dimer)
                Dimer.Rotate_Ring("Dih",10,Dimer.Ring_List[0],Dimer.Ring_List[1])
                Reversed_Dimer.Rotate_Ring("Dih",10,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            if len(In_File_List) != 0:
                End_File = In_File_List[-1]
                #Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,32,Cluster_Location,Job_Type,walltime = 5,queue = qos,proc_per_node = 32,constraint = 'haswell')
                Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
                Copy_File_List.append(Sub_File)
                Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Copy_File_List:
                    os.remove(file)

            if len(Reversed_In_File_List) != 0:
                Reversed_End_File = Reversed_In_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Reversed_Sub_File,Reversed_In_File_List,Reversed_Job_Name,Cluster_Location,Job_Type)
                Reversed_Copy_File_List.append(Reversed_Sub_File)
                Cluster_IO.Submit_Job(Reversed_Copy_File_List,Folder_Name,Reversed_Sub_File,Reversed_End_File,Reversed_Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = Reversed_End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Reversed_Copy_File_List:
                    os.remove(file)
            Dimer.Rotate_Ring("OOP",Phi_Rotation,Dimer.Ring_List[0],Dimer.Ring_List[1])
            Reversed_Dimer.Rotate_Ring("OOP",Phi_Rotation,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            End_File_Matrix.append(End_File_List)
            Reversed_End_File_Matrix.append(Reversed_End_File_List)
        Ring_By_Ring_End_File_Matrices.append(End_File_Matrix)
        Ring_By_Ring_End_File_Matrices.append(Reversed_End_File_Matrix)

    """for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Copy_File_List = []
        In_File_List = []
        for i in range(11):
            Hydrogenated_Dimer = Dimer.Create_Hydrogenated_Copy(0,1)
            XYZ_Filename = Hydrogenated_Dimer.Write_XYZ()
            os.system("scp %s ./Hydrogenated_XYZ_Files" % XYZ_Filename)
            os.system("rm -f ./%s" % XYZ_Filename)
            Job_Type = "QChem"
            Folder_Name = "Hydrogenated_Rotation_Test"
            End_File = "%s_Hydrogenated.out" % Hydrogenated_Dimer.Name
            Cluster_Login = "andrewk@cori.nersc.gov"
            Base_Cluster_Location = '/global/cscratch1/sd/andrewk'
            Cluster_Location='/global/cscratch1/sd/andrewk/Rotation_Test'
            Scheduler_Type = "SLURM"
            End_Condition = "SPE_QChem"
            Shared_File_Location = "/Users/andrewkleinschmidt/Shared_Files_Dihedral_Parameterization"
            Job_Name = "%s_Hydrogenated" % Hydrogenated_Dimer.Name
            In_File = "%s_Hydrogenated.qcin" % Hydrogenated_Dimer.Name
            Sub_File = "sub_%s_Hydrogenated" % Hydrogenated_Dimer.Name
            N = 2
            qos = "premium"
            if End_File not in Run_List:
                Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
                if Return_File == "":
                    Return_File = End_File
                if not Finished:
                    Run_List.append(Symmetry_End_File)
                    Write_Inputs.Write_QChem_SPE(In_File,Hydrogenated_Dimer)
                    In_File_List.append(In_File)
                    Copy_File_List.append(In_File)
            Dimer.Rotate_Ring("OOP",2,Dimer.Ring_List[0],Dimer.Ring_List[1])
        if len(In_File_List) != 0:
            End_File = In_File_List[-1]
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,32,Cluster_Location,Job_Type,walltime = 5,queue = qos,proc_per_node = 32,constraint = 'haswell')
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            for file in Copy_File_List:
                os.system("scp %s ./Rotation_Run_Input_Copies" % file)
                os.system("rm -f %s" % file)"""


    return Ring_By_Ring_End_File_Matrices,Full_Nontorsional_Energy_List

def Run_SPE_Impropers_Hydrogenated(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Alternate = False):
    Phi_Rotation = np.linspace(0,Max_OOP,Rotated_Shape[0])[1]
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    Full_Nontorsional_Energy_List = []
    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))

    Run_List = []
    k = 0
    Ring_By_Ring_End_File_Matrices = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Reversed_Dimer = Conjugated_Polymer.Conjugated_Polymer([ring2,ring1])
        End_File_Matrix = []
        Reversed_End_File_Matrix = []
        for j in range(Rotated_Shape[0]):
            In_File_List = []
            Reversed_In_File_List = []
            Copy_File_List = []
            Reversed_Copy_File_List = []
            End_File_List = []
            Reversed_End_File_List = []
            Nontorsional_Energy_List = []
            Reversed_Nontorsional_Energy_List = []
            for i in range(36):
                if not Alternate:
                    XYZ_Filename = Dimer.Name + "_Hydrogenated_Improper.xyz"
                    Dimer.Hydrogenated_Improper(1,0,XYZ_Filename)
                else:
                    XYZ_Filename = Dimer.Name + "_Hydrogenated_Improper_Alternate.xyz"
                    Dimer.Hydrogenated_Improper_Alternate(1,0,XYZ_Filename)
                Improper_Molecule = Molecule.Molecule(XYZ_Filename)
                os.system("cp %s ./Hydrogenated_Improper_XYZ_Files/" % XYZ_Filename)
                os.remove(XYZ_Filename)

                if not Alternate:
                    Reversed_XYZ_Filename = Reversed_Dimer.Name + "_Hydrogenated_Improper.xyz"
                    Reversed_Dimer.Hydrogenated_Improper(1,0,Reversed_XYZ_Filename)
                else:
                    Reversed_XYZ_Filename = Reversed_Dimer.Name + "_Hydrogenated_Improper_Alternate.xyz"
                    Reversed_Dimer.Hydrogenated_Improper_Alternate(1,0,Reversed_XYZ_Filename)
                Reversed_Improper_Molecule = Molecule.Molecule(Reversed_XYZ_Filename)
                os.system("cp %s ./Hydrogenated_Improper_XYZ_Files/" % Reversed_XYZ_Filename)
                os.remove(Reversed_XYZ_Filename)
                if not Alternate:
                    End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated_Improper" % Dimer.Name,"%s_Hydrogenated_Improper" % Reversed_Dimer.Name)
                else:
                    End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated_Improper_Alternate" % Dimer.Name,"%s_Hydrogenated_Improper_Alternate" % Reversed_Dimer.Name)
                Job_Type = "QChem"
                Folder_Name = "Hydrogenated_Impropers_Test"

                Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

                """if ring1.Symmetric and ring2.Symmetric:
                    Symmetry_End_File = Dimer.Return_Symmetry_Name() + "_Hydrogenated_Improper.out"
                    Check_List = [End_File,Symmetry_End_File]
                else:
                    Symmetry_End_File = ""
                    Check_List = End_File"""
                if End_File in Run_List:
                    End_File_List.append(End_File_List)
                """elif Symmetry_End_File in Run_List:
                    End_File_List.append(Symmetry_End_File)"""
                if End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
                    if Return_File == "":
                        Return_File = End_File
                    End_File_List.append(Return_File)
                    if not Finished:
                        Copy_File_List.append(In_File)
                        In_File_List.append(In_File)
                        Run_List.append(End_File)
                        Write_Inputs.Write_QChem_SPE(In_File,Improper_Molecule)

                if Reversed_End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(Reversed_End_File,Folder_Name,Reversed_Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = Reversed_End_File)
                    if Return_File == "":
                        Return_File = Reversed_End_File
                    Reversed_End_File_List.append(Return_File)
                    if not Finished:
                        Reversed_Copy_File_List.append(Reversed_In_File)
                        Reversed_In_File_List.append(Reversed_In_File)
                        Run_List.append(Reversed_End_File)
                        Write_Inputs.Write_QChem_SPE(Reversed_In_File,Reversed_Improper_Molecule)

                Dimer.Rotate_Ring("Dih",10,Dimer.Ring_List[0],Dimer.Ring_List[1])
                Reversed_Dimer.Rotate_Ring("Dih",10,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])

            if len(In_File_List) != 0:
                End_File = In_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
                Copy_File_List.append(Sub_File)
                Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Copy_File_List:
                    os.remove(file)

            if len(Reversed_In_File_List) != 0:
                Reversed_End_File = Reversed_In_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Reversed_Sub_File,Reversed_In_File_List,Reversed_Job_Name,Cluster_Location,Job_Type)
                Reversed_Copy_File_List.append(Reversed_Sub_File)
                Cluster_IO.Submit_Job(Reversed_Copy_File_List,Folder_Name,Reversed_Sub_File,Reversed_End_File,Reversed_Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = Reversed_End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Reversed_Copy_File_List:
                    os.remove(file)
            #Dimer.Rotate_Ring("OOP",10,Dimer.Ring_List[0],Dimer.Ring_List[1])
            Dimer.Rotate_Ring("OOP",Phi_Rotation,Dimer.Ring_List[0],Dimer.Ring_List[1])
            #Reversed_Dimer.Rotate_Ring("OOP",10,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            Reversed_Dimer.Rotate_Ring("OOP",Phi_Rotation,Reversed_Dimer.Ring_List[0],Reversed_Dimer.Ring_List[1])
            End_File_Matrix.append(End_File_List)
            Reversed_End_File_Matrix.append(Reversed_End_File_List)
        Ring_By_Ring_End_File_Matrices.append(End_File_Matrix)
        Ring_By_Ring_End_File_Matrices.append(Reversed_End_File_Matrix)
        print(End_File)

    """for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Copy_File_List = []
        In_File_List = []
        for i in range(11):
            Hydrogenated_Dimer = Dimer.Create_Hydrogenated_Copy(0,1)
            XYZ_Filename = Hydrogenated_Dimer.Write_XYZ()
            os.system("scp %s ./Hydrogenated_XYZ_Files" % XYZ_Filename)
            os.system("rm -f ./%s" % XYZ_Filename)
            Job_Type = "QChem"
            Folder_Name = "Hydrogenated_Rotation_Test"
            End_File = "%s_Hydrogenated.out" % Hydrogenated_Dimer.Name
            Cluster_Login = "andrewk@cori.nersc.gov"
            Base_Cluster_Location = '/global/cscratch1/sd/andrewk'
            Cluster_Location='/global/cscratch1/sd/andrewk/Rotation_Test'
            Scheduler_Type = "SLURM"
            End_Condition = "SPE_QChem"
            Shared_File_Location = "/Users/andrewkleinschmidt/Shared_Files_Dihedral_Parameterization"
            Job_Name = "%s_Hydrogenated" % Hydrogenated_Dimer.Name
            In_File = "%s_Hydrogenated.qcin" % Hydrogenated_Dimer.Name
            Sub_File = "sub_%s_Hydrogenated" % Hydrogenated_Dimer.Name
            N = 2
            qos = "premium"
            if End_File not in Run_List:
                Finished,Return_File = Cluster_IO.Check_Finished(End_File,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = End_File)
                if Return_File == "":
                    Return_File = End_File
                if not Finished:
                    Run_List.append(Symmetry_End_File)
                    Write_Inputs.Write_QChem_SPE(In_File,Hydrogenated_Dimer)
                    In_File_List.append(In_File)
                    Copy_File_List.append(In_File)
            Dimer.Rotate_Ring("OOP",2,Dimer.Ring_List[0],Dimer.Ring_List[1])
        if len(In_File_List) != 0:
            End_File = In_File_List[-1]
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,32,Cluster_Location,Job_Type,walltime = 5,queue = qos,proc_per_node = 32,constraint = 'haswell')
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            for file in Copy_File_List:
                os.system("scp %s ./Rotation_Run_Input_Copies" % file)
                os.system("rm -f %s" % file)"""


    return Ring_By_Ring_End_File_Matrices,Full_Nontorsional_Energy_List

def Return_SPE_Dimers(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_End_File_Matrices,Ring_By_Ring_Nontorsional_Energy,Ring_By_Ring_Hydrogenated_Energy):
    Full_Raw_Energies = []
    Full_Corrected_Energies = []
    Full_Raw_Energies_Lookup = {}
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])

    Job_Type = "QChem"
    Folder_Name = "Rotation_Test"
    Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

    print("Dimer Matrix size:")
    print(np.asarray(Ring_By_Ring_End_File_Matrices).shape)
    print("Dimer Nonbonded Matrix size:")
    print(np.asarray(Ring_By_Ring_Nontorsional_Energy).shape)

    Ring1_List = []
    Ring2_List = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring2_List.append(ring2)
        Ring1_List.append(ring2)
        Ring2_List.append(ring1)


    for ring1,ring2,End_File_Matrix,Nontorsional_Energy_Matrix,Hydrogenated_Energy in zip(Ring1_List,Ring2_List,Ring_By_Ring_End_File_Matrices,Ring_By_Ring_Nontorsional_Energy,Ring_By_Ring_Hydrogenated_Energy):
        Raw_Energies = []
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s" % Dimer.Name)
        End_File = "%s.out" % Dimer.Name
        Job_Name = "%s" % Dimer.Name
        for End_File_List in End_File_Matrix:
            Raw_Energies.append(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Shared_File_Location = Shared_File_Location,Return_Energy_QChem = True)[0])
            Nontorsional_Energy_List = []
        Raw_Energies = np.asarray(Raw_Energies)
        Raw_Energies = Raw_Energies - np.amin(Raw_Energies)
        Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
        OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)
        Corrected_Energy_Matrix = Raw_Energies - Nontorsional_Energy_Matrix
        Corrected_Energy_Matrix = Corrected_Energy_Matrix - np.amin(Corrected_Energy_Matrix[0])
        Sub_Hydrog_Energy_Matrix = Raw_Energies - Hydrogenated_Energy
        Sub_Hydrog_Energy_Matrix = Sub_Hydrog_Energy_Matrix - np.amin(Sub_Hydrog_Energy_Matrix[0])
        Normalized_Energy_List = []
        for energy_list in Sub_Hydrog_Energy_Matrix:
            norm_energies = energy_list - np.amin(energy_list)
            Normalized_Energy_List.append(norm_energies)
        Full_Raw_Energies.append(Raw_Energies)
        Full_Raw_Energies_Lookup[(ring1.Name,ring2.Name)] = Raw_Energies
        Full_Corrected_Energies.append(Corrected_Energy_Matrix)
        Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Raw_Energies,'%s_%s_Raw_Torsional_Energies' % (ring1.Name,ring2.Name),Title='RI-MP2 Energies (Raw)',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')
        Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Nontorsional_Energy_Matrix - np.amin(Nontorsional_Energy_Matrix),'%s_%s_Nonbonded_Energies' % (ring1.Name,ring2.Name),Title='Nonbonded Energies',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')
        Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Corrected_Energy_Matrix,'%s_%s_Modified_Torsional_Energies' % (ring1.Name,ring2.Name),Title='RI-MP2 Energies (Modified)',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')
        Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Sub_Hydrog_Energy_Matrix,'%s_%s_Sub_Hydrog_Torsional_Energies' % (ring1.Name,ring2.Name),Title='RI-MP2 Energies (Sub Hydrog)',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')
        Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Normalized_Energy_List,'%s_%s_Sub_Hydrog_Torsional_Energies_Normalized' % (ring1.Name,ring2.Name),Title='RI-MP2 Energies (Sub Hydrog)',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')

    return Full_Raw_Energies,Full_Raw_Energies_Lookup,Full_Corrected_Energies

def Return_SPE_Dimers_Hydrogenated(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_End_File_Matrices,Alternate = False):
    Full_Raw_Energies = []
    Full_Corrected_Energies = []
    Full_Raw_Energies_Lookup = {}
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])

    Job_Type = "QChem"
    Folder_Name = "Hydrogenated_Rotation_Test"
    Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

    Ring1_List = []
    Ring2_List = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring2_List.append(ring2)
        Ring1_List.append(ring2)
        Ring2_List.append(ring1)

    for ring1,ring2,End_File_Matrix in zip(Ring1_List,Ring2_List,Ring_By_Ring_End_File_Matrices):
        Raw_Energies = []
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        if not Alternate:
            End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated" % Dimer.Name)
        else:
            End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated_Alternate" % Dimer.Name)
        for End_File_List in End_File_Matrix:
            Raw_Energies.append(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Shared_File_Location = Shared_File_Location,Return_Energy_QChem = True)[0])
        Raw_Energies = np.asarray(Raw_Energies)
        Raw_Energies = Raw_Energies - np.amin(Raw_Energies)
        Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
        OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)

        if not Alternate:
            Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Raw_Energies,'%s_%s_Raw_Torsional_Energies_Hydrogenated' % (ring1.Name,ring2.Name),Title='Hydrogenated RI-MP2 Energies (Raw)',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')
        else:
            Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Raw_Energies,'%s_%s_Raw_Torsional_Energies_Hydrogenated_Alternate' % (ring1.Name,ring2.Name),Title='Hydrogenated RI-MP2 Energies (Raw)',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')

        Full_Raw_Energies.append(Raw_Energies)

    return Full_Raw_Energies

def Return_SPE_Dimers_Impropers_Hydrogenated(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_End_File_Matrices,Alternate = False):
    Full_Raw_Energies = []
    Full_Corrected_Energies = []
    Full_Raw_Energies_Lookup = {}
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])

    Job_Type = "QChem"
    Folder_Name = "Hydrogenated_Impropers_Test"
    Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

    print("Dimer Matrix size:")
    print("Hydrogenated")

    Ring1_List = []
    Ring2_List = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring2_List.append(ring2)
        Ring1_List.append(ring2)
        Ring2_List.append(ring1)

    for ring1,ring2,End_File_Matrix in zip(Ring1_List,Ring2_List,Ring_By_Ring_End_File_Matrices):
        Raw_Energies = []
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        if not Alternate:
            End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated_Improper" % Dimer.Name)
        else:
            End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Hydrogenated_Improper_Alternate" % Dimer.Name)
        for End_File_List in End_File_Matrix:
            Raw_Energies.append(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Shared_File_Location = Shared_File_Location,Return_Energy_QChem = True)[0])
        Raw_Energies = np.asarray(Raw_Energies)
        Raw_Energies = Raw_Energies - np.amin(Raw_Energies)
        print(Raw_Energies)
        print(Raw_Energies.shape)
        Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
        OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)

        if not Alternate:
            Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Raw_Energies,'%s_%s_Raw_Torsional_Energies_Hydrogenated_Impropers' % (ring1.Name,ring2.Name),Title='Hydrogenated Improper RI-MP2 Energies',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')
        else:
            Make_Surface_Plot(Dih_Rotations_Degrees,OOP_Rotations_Degrees,Raw_Energies,'%s_%s_Raw_Torsional_Energies_Hydrogenated_Impropers_Alternate' % (ring1.Name,ring2.Name),Title='Hydrogenated Improper RI-MP2 Energies',xlabel='Dihedral Angle (degrees)',ylabel='Out of Plane Angle (degrees)')

        Full_Raw_Energies.append(Raw_Energies)

    return Full_Raw_Energies

def Run_SPE_Trimers_Dih(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name):
    Ring_By_Ring_End_File_Matrices = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []
    Full_Trimer_Nontorsional_Energy_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(copy.deepcopy(ring))

    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Dual_Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Dual_Offset_Ring_List.append(copy.deepcopy(Ring_List[1]))

    k = 0
    Run_List = []

    for ring1,ring2,ring3 in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List):
        Trimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2,ring3])
        Trimer.Create_Dihydrogenated_Monomer(1)
        End_File_Matrix = []
        Nontorsional_Energy_Matrix = []
        for j in range(Rotated_Shape[1]):
            In_File_List = []
            Copy_File_List = []
            End_File_List = []
            Nontorsional_Energy_List = []
            for i in range(Rotated_Shape[1]):
                Trimer.Rotate_Ring("Dih",0,Trimer.Ring_List[0],Trimer.Ring_List[1])
                XYZ_Filename = Trimer.Write_XYZ()
                os.system("cp %s ./XYZ_Files/" % XYZ_Filename)
                os.remove(XYZ_Filename)
                Job_Type = "QChem"
                Folder_Name = "Multi_Ring_Rotation_Test"
                End_File = "%s.out" % Trimer.Name
                Cluster_Login = Configure.qchem_dict["Cluster_Login"]
                Base_Cluster_Location = Configure.qchem_dict["Base_Cluster_Location"]
                Cluster_Location=Base_Cluster_Location + "/Multi_Ring_Rotation_Test"
                Scheduler_Type = Configure.qchem_dict["Scheduler_Type"]
                End_Condition = Configure.qchem_dict["End_Condition"]
                Shared_File_Location = Configure.qchem_dict["Shared_File_Location"]
                Job_Name = "%s" % Trimer.Name
                In_File = "%s.qcin" % Trimer.Name
                Sub_File = "sub_%s" % Trimer.Name
                N = 2
                if k < 5:
                    qos = "regular"
                else:
                    qos = "regular"
                Percent_Change = 0.0
                if ring1.Symmetric and ring2.Symmetric:
                    Symmetry_End_File = Trimer.Return_Symmetry_Name() + ".out"
                    Check_List = [End_File,Symmetry_End_File]
                else:
                    Symmetry_End_File = ""
                    Check_List = End_File
                if End_File in Run_List:
                    End_File_List.append(End_File)
                elif Symmetry_End_File in Run_List:
                    End_File_List.append(Symmetry_End_File)
                elif End_File not in Run_List and Symmetry_End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(Check_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = Check_List)
                    if Return_File == "":
                        Return_File = End_File
                    End_File_List.append(Return_File)
                    if not Finished:
                        Copy_File_List.append(In_File)
                        In_File_List.append(In_File)
                        Run_List.append(End_File)
                        if Symmetry_End_File != "":
                            Run_List.append(Symmetry_End_File)
                        Write_Inputs.Write_QChem_SPE(In_File,Trimer)
                if Symmetry_End_File == "":
                    Trimer.Read_From_Data_File("./LigParGen_Files/%s_%s_%s.lmp" % (ring1.Name,ring2.Name,ring3.Name),No_Position_Update = True)
                    Nontorsional_Energy_List.append(Trimer.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/lmp_serial",Exclude_Interring_Torsions=True,dielectric=4.9))
                else:
                    if Return_File != Symmetry_End_File or End_File == Symmetry_End_File:
                        Trimer.Read_From_Data_File("./LigParGen_Files/%s_%s_%s.lmp" % (ring1.Name,ring2.Name,ring3.Name),No_Position_Update = True)
                    Nontorsional_Energy_List.append(Trimer.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/lmp_serial",Exclude_Interring_Torsions=True,dielectric=4.9,Symmetric=True))

                Trimer.Rotate_Ring("Dih",10,Trimer.Ring_List[0],Trimer.Ring_List[1])
            if len(In_File_List) != 0:
                End_File = Copy_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
                Copy_File_List.append(Sub_File)
                Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Copy_File_List:
                    os.system("cp %s ./Rotation_Run_Input_Copies/" % file)
                    os.remove(file)
            Trimer.Rotate_Ring("Dih",10,Trimer.Ring_List[2],Trimer.Ring_List[1])
            End_File_Matrix.append(End_File_List)
            Nontorsional_Energy_Matrix.append(np.array(Nontorsional_Energy_List))
        Ring_By_Ring_End_File_Matrices.append(End_File_Matrix)
        print(End_File)
        Full_Trimer_Nontorsional_Energy_List.append(np.array(Nontorsional_Energy_Matrix) - np.amin(np.array(Nontorsional_Energy_Matrix)[0]))

    print("Matrix size:")
    print(np.asarray(Ring_By_Ring_End_File_Matrices).shape)

    return Ring_By_Ring_End_File_Matrices,Full_Trimer_Nontorsional_Energy_List

def Run_Paired_Hydrogenation_Energy(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name):
    Ring_By_Ring_Dual_Hydrogenated_End_File_Matrices = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []

    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(copy.deepcopy(ring))

    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Dual_Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Dual_Offset_Ring_List.append(copy.deepcopy(Ring_List[1]))

    k = 0
    Run_List = []

    for ring1,ring2,ring3 in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List):
        Trimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2,ring3])
        if Trimer.Name == "Furan_Phi_0_Theta_0_Dimethyl_Diketopyrrolopyrrole_Phi_0_Theta_0_Furan":
            Trimer.Rotate_Ring("Dih",180,Trimer.Ring_List[2],Trimer.Ring_List[1])
            Trimer.Dih_Rotation = np.zeros(len(Trimer.Ring_List)-1)
            Trimer.Refresh_Name()
        End_File_Matrix = []
        Nontorsional_Energy_Matrix = []
        Nontorsional_Energy_List = []
        Name_List = Trimer.Create_Dihydrogenated_Monomer(1)

        Check_List = [file.split(".")[0] + ".out" for file in Name_List]
        End_File_List = [file.split(".")[0] + ".out" for file in Name_List]
        Run_List = [file.split(".")[0] + ".out" for file in Name_List]
        In_File_List = [file.split(".")[0] + ".qcin" for file in Name_List]
        Copy_File_List = [file.split(".")[0] + ".qcin" for file in Name_List]

        print(End_File_List)


        for file,in_file in zip(Name_List,In_File_List):
            Dihydrogenated_Ring = Molecule.Molecule(file)
            Write_Inputs.Write_QChem_SPE(in_file,Dihydrogenated_Ring)

        Job_Type = "QChem"
        Folder_Name = "Dual_Hydrogenated_Test"
        End_File = End_File_List[-1]
        Cluster_Login = Configure.qchem_dict["Cluster_Login"]
        Base_Cluster_Location = Configure.qchem_dict["Base_Cluster_Location"]
        Cluster_Location = Configure.qchem_dict["Cluster_Location"]
        Scheduler_Type = Configure.qchem_dict["Scheduler_Type"]
        End_Condition = Configure.qchem_dict["End_Condition"]
        Shared_File_Location = Configure.qchem_dict["Shared_File_Location"]
        Job_Name = "%s_Dual_Hydrogenated" % Trimer.Ring_List[1].Name
        In_File = "%s_Dual_Hydrogenated.qcin" % Trimer.Ring_List[1].Name
        Sub_File = "sub_%s_Dual_Hydrogenated" % Trimer.Ring_List[1].Name
        #qos = "debug"
        #Symmetry_End_File = ""

        if len(In_File_List) != 0:
            Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
            Copy_File_List.append(Sub_File)
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
            if not os.path.isdir('Rotation_Run_Input_Copies'):
                os.mkdir('Rotation_Run_Input_Copies')
            for file in Copy_File_List:
                os.system("cp %s ./Rotation_Run_Input_Copies/" % file)
                os.remove(file)

        Ring_By_Ring_Dual_Hydrogenated_End_File_Matrices.append(End_File_List)

    return Ring_By_Ring_Dual_Hydrogenated_End_File_Matrices

def Return_SPE_Methyl_Impropers(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_Improper_File_Matrices):
    #Check which ring1 vs ring2 correspond
    Offset_Ring_List = []

    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))

    Phi_Rotation = np.linspace(0,Max_OOP,Rotated_Shape[0])[1]

    Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
    OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)

    Dimer_Names = []

    Ring1_List = []
    Ring2_List = []
    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring2_List.append(ring2)
        Ring1_List.append(ring2)
        Ring2_List.append(ring1)

    Ring_By_Ring_Methyl_Impropers = []
    Ring_By_Ring_Methyl_Improper_Lists = []
    for ring1,ring2,Improper_File_Matrix in zip(Ring1_List,Ring2_List,Ring_By_Ring_Improper_File_Matrices):
        Improper_Energies = []
        End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Improper_Bend_Methyl" % ring1.Name)
        Folder_Name = "Improper_Bend_Test"
        #Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)
        #TODO: Confirm if this should actually use qchem
        Cluster_Login = Configure.qchem_dict["Cluster_Login"]
        Base_Cluster_Location = Configure.qchem_dict["Base_Cluster_Location"]
        Cluster_Location = Configure.qchem_dict["Cluster_Location"]
        Scheduler_Type = Configure.qchem_dict["Scheduler_Type"]
        End_Condition = Configure.qchem_dict["End_Condition"]
        Shared_File_Location = Configure.qchem_dict["Shared_File_Location"]

        for End_File_List in Improper_File_Matrix:
            Improper_Energies.append(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Shared_File_Location = Shared_File_Location,Return_Energy_QChem = True)[0])
        Improper_Energies = np.asarray(Improper_Energies)
        Improper_Energies = Improper_Energies - Improper_Energies[0][0]
        Ring_By_Ring_Methyl_Impropers.append(Improper_Energies)
        """if ring1.Name + ring2.Name not in Dimer_Names:
            Dimer_Names.append(ring1.Name + ring2.Name)
            Dimer_Names.append(ring2.Name + ring1.Name)
            End_File_List = []
            Reversed_End_File_List = []

            for i in range(Rotated_Shape[0]):
                Job_Type = "QChem"
                Folder_Name = "Improper_Bend_Test"
                End_File,Reversed_End_File,Job_Name,Reversed_Job_Name,In_File,Reversed_In_File,Sub_File,Reversed_Sub_File = Return_Filenames("%s_Improper_Bend_Methyl_Phi_%d" % (ring2.Name,i*Phi_Rotation),Reversed_Name = "%s_Improper_Bend_Methyl_Phi_%d" % (ring1.Name,i*Phi_Rotation))
                End_File_List.append("%s_Improper_Bend_Methyl_Phi_%d.out" % (ring2.Name,i*Phi_Rotation))
                Reversed_End_File_List.append("%s_Improper_Bend_Methyl_Phi_%d.out" % (ring1.Name,i*Phi_Rotation))
                Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

            Methyl_Impropers = (np.array(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Return_Energy_QChem = True,Shared_File_Location = Shared_File_Location)[0]))
            Reversed_Methyl_Impropers = (np.array(Cluster_IO.Return_Info_Batch(Reversed_End_File_List,Reversed_End_File_List,Folder_Name,Reversed_Job_Name,Cluster_Login,Cluster_Location,Return_Energy_QChem = True,Shared_File_Location = Shared_File_Location)[0]))
            Methyl_Impropers = Methyl_Impropers - np.amin(Methyl_Impropers)
            Methyl_Improper_Matrix = []
            for meth_imp in Methyl_Impropers:
                Methyl_Improper_Matrix.append(np.ones(Rotated_Shape[1])*meth_imp)
            Methyl_Improper_Matrix = np.asarray(Methyl_Improper_Matrix)

            Reversed_Methyl_Impropers = Reversed_Methyl_Impropers - np.amin(Reversed_Methyl_Impropers)
            Reversed_Methyl_Improper_Matrix = []
            for rev_meth_imp in Methyl_Impropers:
                Reversed_Methyl_Improper_Matrix.append(np.ones(Rotated_Shape[1])*rev_meth_imp)
            Reversed_Methyl_Improper_Matrix = np.asarray(Reversed_Methyl_Improper_Matrix)

            Ring_By_Ring_Methyl_Impropers.append(Methyl_Improper_Matrix)
            Ring_By_Ring_Methyl_Improper_Lists.append(Methyl_Impropers)
            Ring_By_Ring_Methyl_Impropers.append(Reversed_Methyl_Improper_Matrix)
        Ring_By_Ring_Methyl_Improper_Lists.append(Reversed_Methyl_Impropers)"""
        fig,ax = plt.subplots(1,1)
        x,y = np.meshgrid(Dih_Rotations_Degrees,OOP_Rotations_Degrees)
        c = ax.pcolor(x,y,Improper_Energies,cmap = 'seismic',vmin=0,vmax=10)
        ax.set_title('Conjugated Energies (Raw)',fontdict = {'fontsize':24})
        plt.xlabel('Dihedral Angle (degrees)',size = 24)
        plt.ylabel('OOP Angle (degrees)',size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        fig.savefig('%s_Conjugated_Energies' % (ring2.Name))
        plt.close(fig)

        """fig,ax = plt.subplots(1,1)
        plt.scatter(np.linspace(0,90,Rotated_Shape[0]),Methyl_Impropers)
        ax.set_title('Conjugated Energies (Raw)',fontdict = {'fontsize':24})
        plt.xlabel('Dihedral Angle (degrees)',size = 24)
        plt.ylabel('OOP Angle (degrees)',size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        fig.savefig('%s_Conjugated_Energies_XY' % (ring2.Name))
        plt.close(fig)

        fig,ax = plt.subplots(1,1)
        x,y = np.meshgrid(Dih_Rotations_Degrees,OOP_Rotations_Degrees)
        c = ax.pcolor(x,y,Reversed_Methyl_Improper_Matrix,cmap = 'seismic',vmin=0,vmax=10)
        ax.set_title('Conjugated Energies (Raw)',fontdict = {'fontsize':24})
        plt.xlabel('Dihedral Angle (degrees)',size = 24)
        plt.ylabel('OOP Angle (degrees)',size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        fig.savefig('%s_Conjugated_Energies' % (ring1.Name))
        plt.close(fig)"""

    color_list = ['k','b','r','c','y']
    fig,ax = plt.subplots(1,1)
    for meth_imp,color in zip(Ring_By_Ring_Methyl_Improper_Lists,color_list):
        plt.scatter(np.linspace(0,90,Rotated_Shape[0]),meth_imp,marker = "s",c = color)

    print(np.array(Ring_By_Ring_Methyl_Impropers).shape)

    plt.ylabel('Energy (kcal/mol)',size = 24)
    plt.xlabel('Out Of Plane Angle ($^\circ$)',size = 24)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.tick_params(length=4,width=4)
    plt.tight_layout()
    fig.savefig('%s_Conjugated_Energies_Comparison' % (Polymer_Name))
    plt.close(fig)

    return Ring_By_Ring_Methyl_Impropers,Ring_By_Ring_Methyl_Improper_Lists,Dimer_Names

def Run_SPE_Trimers_Hydrogenated_Dih(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name):
    Ring_By_Ring_End_File_Matrices = []
    Ring_By_Ring_Syn_Anti_Matrices = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []

    Full_Trimer_Nontorsional_Energy_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(copy.deepcopy(ring))

    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Dual_Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Dual_Offset_Ring_List.append(copy.deepcopy(Ring_List[1]))

    k = 0
    Run_List = []

    for ring1,ring2,ring3 in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List):
        Trimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2,ring3])
        if Trimer.Name == "Furan_Phi_0_Theta_0_Dimethyl_Diketopyrrolopyrrole_Phi_0_Theta_0_Furan":
            Trimer.Rotate_Ring("Dih",180,Trimer.Ring_List[2],Trimer.Ring_List[1])
            Trimer.Dih_Rotation = np.zeros(len(Trimer.Ring_List)-1)
            Trimer.Refresh_Name()
        End_File_Matrix = []
        Nontorsional_Energy_Matrix = []
        Syn_Anti_Matrix = []
        for j in range(Rotated_Shape[1]):
            In_File_List = []
            Copy_File_List = []
            End_File_List = []
            Nontorsional_Energy_List = []
            Syn_Anti_List = []
            for i in range(Rotated_Shape[1]):
                Hydrogenated_Trimer,Syn = Trimer.Create_Hydrogenated_Copy(1,0)
                Syn_Anti_List.append(Syn)
                XYZ_Filename = Hydrogenated_Trimer.Write_XYZ()
                os.system("cp %s ./Hydrogenated_XYZ_Files/" % XYZ_Filename)
                os.remove(XYZ_Filename)
                Job_Type = "QChem"
                Folder_Name = "Multi_Ring_Hydrogenated_Rotation_Test"
                End_File = "%s_Hydrogenated.out" % Hydrogenated_Trimer.Name
                Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)
                Job_Name = "%s_Hydrogenated" % Hydrogenated_Trimer.Name
                In_File = "%s_Hydrogenated.qcin" % Hydrogenated_Trimer.Name
                Sub_File = "sub_%s_Hydrogenated" % Hydrogenated_Trimer.Name
                N = 2
                if k < 5:
                    qos = "regular"
                else:
                    qos = "regular"
                Percent_Change = 0.0
                if ring1.Symmetric and ring2.Symmetric:
                    Symmetry_End_File = Hydrogenated_Trimer.Return_Symmetry_Name() + "_Hydrogenated.out"
                    Check_List = [End_File,Symmetry_End_File]
                else:
                    Symmetry_End_File = ""
                    Check_List = End_File
                if End_File in Run_List:
                    End_File_List.append(End_File)
                elif Symmetry_End_File in Run_List:
                    End_File_List.append(Symmetry_End_File)
                elif End_File not in Run_List and Symmetry_End_File not in Run_List:
                    Finished,Return_File = Cluster_IO.Check_Finished(Check_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Analyze_File = Check_List)
                    if Return_File == "":
                        Return_File = End_File
                    End_File_List.append(Return_File)
                    if not Finished:
                        Copy_File_List.append(In_File)
                        In_File_List.append(In_File)
                        Run_List.append(End_File)
                        if Symmetry_End_File != "":
                            Run_List.append(Symmetry_End_File)
                        Write_Inputs.Write_QChem_SPE(In_File,Hydrogenated_Trimer)

                Trimer.Rotate_Ring("Dih",10,Trimer.Ring_List[0],Trimer.Ring_List[1])
            if len(In_File_List) != 0:
                End_File = Copy_File_List[-1]
                Write_Submit_Script.Write_SLURM_Batch(Sub_File,In_File_List,Job_Name,Cluster_Location,Job_Type)
                Copy_File_List.append(Sub_File)
                Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Shared_File_Location = Shared_File_Location)
                k+=1
                for file in Copy_File_List:
                    os.system("cp %s ./Rotation_Run_Input_Copies/"%file)
                    os.remove(file)
            Trimer.Rotate_Ring("Dih",10,Trimer.Ring_List[2],Trimer.Ring_List[1])
            End_File_Matrix.append(End_File_List)
            Syn_Anti_Matrix.append(Syn_Anti_List)
        Ring_By_Ring_End_File_Matrices.append(End_File_Matrix)
        Ring_By_Ring_Syn_Anti_Matrices.append(Syn_Anti_Matrix)
        print(End_File)

    return Ring_By_Ring_End_File_Matrices,Ring_By_Ring_Syn_Anti_Matrices

def Return_SPE_Trimers_Hydrogenated_Dih(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_End_File_Matrices,Dimer_Energies,Ring_By_Ring_Trimer_Nontorsional_Energy):
    Full_Raw_Energies = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []

    Folder_Name = "Multi_Ring_Hydrogenated_Rotation_Test"
    Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)

    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[1])


    for ring1,ring2,ring3,End_File_Matrix,trimer_nontorsional_matrix in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List,Ring_By_Ring_End_File_Matrices,Ring_By_Ring_Trimer_Nontorsional_Energy):
        Raw_Energies = []
        Trimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2,ring3])
        Job_Name = "%s" % Trimer.Name
        for End_File_List in End_File_Matrix:
            print(End_File_List[1])
            Raw_Energies.append(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Return_Energy_QChem = True,Shared_File_Location = Shared_File_Location)[0])
        Raw_Energies = np.asarray(Raw_Energies)
        Raw_Energies = Raw_Energies - np.amin(Raw_Energies)
        fig,ax = plt.subplots(1,1)
        Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
        OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Raw_Energies,'%s_%s_%s_Hydrogenated_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Hydrogenated)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Full_Raw_Energies.append(Raw_Energies)

    return Full_Raw_Energies

def Return_SPE_Trimers_Dih(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_End_File_Matrices,Dimer_Energies,Ring_By_Ring_Trimer_Nontorsional_Energy):
    Full_Raw_Energies = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []
    #TODO: Figure out what Folder name this is supposed to get
    Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name) 
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[1])


    for ring1,ring2,ring3,End_File_Matrix,trimer_nontorsional_matrix in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List,Ring_By_Ring_End_File_Matrices,Ring_By_Ring_Trimer_Nontorsional_Energy):
        Raw_Energies = []
        Trimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2,ring3])
        Job_Name = "%s" % Trimer.Name
        for End_File_List in End_File_Matrix:
            print(End_File_List[1])
            Raw_Energies.append(Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Return_Energy_QChem = True,Shared_File_Location = Shared_File_Location)[0])
        Raw_Energies = np.asarray(Raw_Energies)
        Raw_Energies = Raw_Energies - np.amin(Raw_Energies)

        fig,ax = plt.subplots(1,1)
        Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
        OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Raw_Energies,'%s_%s_%s_Raw_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Hydrogenated)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Full_Raw_Energies.append(Raw_Energies)

        Corrected_Energies = Raw_Energies - trimer_nontorsional_matrix
        Corrected_Energies = Corrected_Energies - np.amin(Corrected_Energies)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Corrected_Energies,'%s_%s_%s_Corrected_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Corrected)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Additive_Energy_List = []

        for energy_list,file_list in zip(Raw_Energies,End_File_Matrix):
            add_energies = []
            for e,f in zip(energy_list,file_list):
                Ring1_Name = f.split("_Phi")[0]
                theta_pattern = re.compile("_Theta_\d+_")
                Ring2_Name = re.split(theta_pattern,f)[1].split("_Phi")[0]
                Ring3_Name = re.split(theta_pattern,f)[-1].split(".")[0]
                Torsion1_Index = int(f.split("Theta_")[1].split("_")[0])/10
                Torsion2_Index = int(f.split("Theta_")[-1].split("_")[0])/10
                add_energies.append(e-Dimer_Energies[(Ring1_Name,Ring2_Name)][0][Torsion1_Index]-Dimer_Energies[(Ring2_Name,Ring3_Name)][0][Torsion2_Index])
            Additive_Energy_List.append(add_energies)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Additive_Energy_List,'%s_%s_%s_Additive_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Relative to Additive)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))


        Relative_Energy_List = []

        for energy_list in Raw_Energies[::-1]:
            rel_energies = energy_list - np.amin(energy_list)
            rel_energies = rel_energies - Raw_Energies[0]
            Relative_Energy_List.append(rel_energies)

        Relative_Energy_List = Relative_Energy_List[::-1]

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Relative_Energy_List,'%s_%s_%s_Relative_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Relative to Planar)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Relative_Energy_List_Percent = np.zeros((len(Relative_Energy_List),len(Relative_Energy_List[0])))

        for energy_list,r_energy_list,i in zip(Raw_Energies,Relative_Energy_List,range(len(Raw_Energies))):
            for e,r_e,j in zip(energy_list,r_energy_list,range(len(energy_list))):
                if e < .01:
                    Relative_Energy_List_Percent[i][j] = 0
                else:
                    Relative_Energy_List_Percent[i][j] = r_e/e

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Relative_Energy_List_Percent,'%s_%s_%s_Relative_Dihedral_Energies_Percent' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Relative to Planar, %)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Normalized_Energy_List = []

        for energy_list in Raw_Energies:
            norm_energies = energy_list - np.amin(energy_list)
            Normalized_Energy_List.append(norm_energies)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Normalized_Energy_List,'%s_%s_%s_Normalized_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Per-line normalized)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        for i in range(len(Normalized_Energy_List)):
            fig,ax = plt.subplots(1,1)
            for energy_list in Normalized_Energy_List[:i]:
                plt.scatter(np.linspace(0,350,Rotated_Shape[1]),energy_list,alpha=0.2,marker = 's',c='k')

            plt.scatter(np.linspace(0,350,Rotated_Shape[1]),Normalized_Energy_List[i],marker = 's',c='k')
            plt.xlabel('Dihedral Angle ($^\circ$)',size = 24)
            plt.ylabel('Energy (kcal/mol)',size = 24)
            ax.tick_params(axis="x", labelsize=18)
            ax.tick_params(axis="y", labelsize=18)
            ax.tick_params(length=4,width=4)
            plt.ylim([-1,12])
            plt.tight_layout()
            fig.savefig('%s_%s_%s_Cooperative_Energies_Row_By_Row_%d' % (ring1.Name,ring2.Name,ring3.Name,i))
            plt.close(fig)
            os.system("cp %s_%s_%s_Cooperative_Energies_Row_By_Row_%d.png ./Figures/"%(ring1.Name,ring2.Name,ring3.Name,i))
            os.remove("%s_%s_%s_Cooperative_Energies_Row_By_Row_%d.png"%(ring1.Name,ring2.Name,ring3.Name,i))

        Normalized_Energy_List = []

        for energy_list in np.matrix.transpose(Raw_Energies):
            norm_energies = energy_list - np.amin(energy_list)
            Normalized_Energy_List.append(norm_energies)


        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Normalized_Energy_List,'%s_%s_%s_Normalized_Dihedral_Energies' % (ring3.Name,ring2.Name,ring1.Name),Title='RI-MP2 Energies (Per-line normalized)',xlabel='%s %s Dihedral (degrees)' % (ring3.Nickname,ring2.Nickname),ylabel='%s %s Dihedral (degrees)' % (ring2.Nickname,ring1.Nickname))

        for i in range(len(Normalized_Energy_List)):
            fig,ax = plt.subplots(1,1)
            for energy_list in Normalized_Energy_List[:i]:
                plt.scatter(np.linspace(0,350,Rotated_Shape[1]),energy_list,alpha=0.2,marker = 's',c='k')

            plt.scatter(np.linspace(0,350,Rotated_Shape[1]),Normalized_Energy_List[i],marker = 's',c='k')
            plt.xlabel('Dihedral Angle ($^\circ$)',size = 24)
            plt.ylabel('Energy (kcal/mol)',size = 24)
            ax.tick_params(axis="x", labelsize=18)
            ax.tick_params(axis="y", labelsize=18)
            ax.tick_params(length=4,width=4)
            plt.ylim([-1,12])
            plt.tight_layout()
            fig.savefig('%s_%s_%s_Cooperative_Energies_Row_By_Row_%d' % (ring3.Name,ring2.Name,ring1.Name,i))
            plt.close(fig)
            os.system("scp %s_%s_%s_Cooperative_Energies_Row_By_Row_%d.png ./Figures/" % (ring3.Name,ring2.Name,ring1.Name,i))
            os.remove("%s_%s_%s_Cooperative_Energies_Row_By_Row_%d.png" % (ring3.Name,ring2.Name,ring1.Name,i))

        Normalized_Energy_List = []

        for energy_list in Corrected_Energies:
            norm_energies = energy_list - np.amin(energy_list)
            Normalized_Energy_List.append(norm_energies)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Normalized_Energy_List,'%s_%s_%s_Normalized_Corrected_Dihedral_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies (Relative to Planar, %)',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))


        f = open('%s_%s_%s_Deflection_Angles' % (ring1.Name,ring2.Name,ring3.Name),'w')
        Current_External_Index = 0
        for i,ne_energy_list in enumerate(Normalized_Energy_List):
            Current_Index = 0
            Reverse_Flag = False
            Weights = np.zeros(10)
            for e in ne_energy_list:
                Weights[Current_Index] += (math.exp(-1*e/.593))
                if Current_Index == 0:
                    Reverse_Flag = False
                if Current_Index == 9:
                    Reverse_Flag = True
                if Reverse_Flag:
                    Current_Index += -1
                else:
                    Current_Index += 1
            Weights = Weights/sum(Weights)
            Averaged_Deflection = 0.0
            for j,w in zip(range(10),Weights):
                Averaged_Deflection += j*10*w
            f.write("%d %.2f\n" % (i*10,Averaged_Deflection))
        f.close()


        Normalized_Energy_List = []

        for energy_list in np.matrix.transpose(Corrected_Energies):
            norm_energies = energy_list - np.amin(energy_list)
            Normalized_Energy_List.append(norm_energies)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Normalized_Energy_List,'%s_%s_%s_Normalized_Corrected_Dihedral_Energies' % (ring3.Name,ring2.Name,ring1.Name),Title='RI-MP2 Energies (Per-line normalized)',xlabel='%s %s Dihedral (degrees)' % (ring3.Nickname,ring2.Nickname),ylabel='%s %s Dihedral (degrees)' % (ring2.Nickname,ring1.Nickname))

        f = open('%s_%s_%s_Deflection_Angles' % (ring3.Name,ring2.Name,ring1.Name),'w')
        for i,ne_energy_list in enumerate(Normalized_Energy_List):
            Current_Index = 0
            Reverse_Flag = False
            Weights = np.zeros(10)
            for e in ne_energy_list:
                Weights[Current_Index] += (math.exp(-1*e/.593))
                if Current_Index == 0:
                    Reverse_Flag = False
                if Current_Index == 9:
                    Reverse_Flag = True
                if Reverse_Flag:
                    Current_Index += -1
                else:
                    Current_Index += 1
            Weights = Weights/sum(Weights)
            Averaged_Deflection = 0.0
            for j,w in zip(range(10),Weights):
                Averaged_Deflection += j*10*w
            f.write("%d %.2f\n" % (i*10,Averaged_Deflection))
        f.close()

    Full_Raw_Energies = np.asarray(Full_Raw_Energies)

    return Full_Raw_Energies

def Merge_Hydrogenated_Energies(Ring_By_Ring_Hydrogenated_Energies,Ring_By_Ring_Hydrogenated_Alternate_Energies):
    Combined_Ring_By_Ring_Hydrogenated_Energies = []
    Combined_Ring_By_Ring_Hydrogenated_Energies_Tracker = []
    for reg_matrix,alt_matrix in zip(Ring_By_Ring_Hydrogenated_Energies,Ring_By_Ring_Hydrogenated_Alternate_Energies):
        Combined_Matrix = []
        Matrix_Tracker = []
        for reg_line,alt_line in zip(reg_matrix,alt_matrix):
            Combined_Line = []
            Line_Tracker = []
            for reg,alt in zip(reg_line,alt_line):
                print("Regular:")
                print(reg)
                print("Alternate:")
                print(alt)
                if reg < alt:
                    Combined_Line.append(reg)
                    Line_Tracker.append(0)
                else:
                    Combined_Line.append(alt)
                    Line_Tracker.append(10)
            Combined_Matrix.append(Combined_Line)
            Matrix_Tracker.append(Line_Tracker)
        Combined_Ring_By_Ring_Hydrogenated_Energies.append(Combined_Matrix)
        Combined_Ring_By_Ring_Hydrogenated_Energies_Tracker.append(Matrix_Tracker)

    Combined_Ring_By_Ring_Hydrogenated_Energies = np.asarray(Combined_Ring_By_Ring_Hydrogenated_Energies)

    return Combined_Ring_By_Ring_Hydrogenated_Energies,Combined_Ring_By_Ring_Hydrogenated_Energies_Tracker

def Merge_Hydrogenated_Improper_Energies(Ring_By_Ring_Hydrogenated_Improper_Energies,Ring_By_Ring_Hydrogenated_Alternate_Improper_Energies,Combined_Ring_By_Ring_Hydrogenated_Energies_Tracker):
    Combined_Ring_By_Ring_Hydrogenated_Improper_Energies = []
    Combined_Ring_By_Ring_Hydrogenated_Improper_Energies_Tracker = []
    for reg_matrix,alt_matrix,track_matrix in zip(Ring_By_Ring_Hydrogenated_Improper_Energies,Ring_By_Ring_Hydrogenated_Alternate_Improper_Energies,Combined_Ring_By_Ring_Hydrogenated_Energies_Tracker):
        Combined_Matrix = []
        Matrix_Tracker = []
        for reg_line,alt_line,track_line in zip(reg_matrix,alt_matrix,track_matrix):
            Combined_Line = []
            Line_Tracker = []
            for reg,alt,track in zip(reg_line,alt_line,track_line):
                if track == 0:
                    Combined_Line.append(reg)
                    Line_Tracker.append(0)
                else:
                    Combined_Line.append(alt)
                    Line_Tracker.append(10)
            Combined_Matrix.append(Combined_Line)
            Matrix_Tracker.append(Line_Tracker)
        Combined_Ring_By_Ring_Hydrogenated_Improper_Energies.append(Combined_Matrix)
        Combined_Ring_By_Ring_Hydrogenated_Improper_Energies_Tracker.append(Matrix_Tracker)

    Combined_Ring_By_Ring_Hydrogenated_Improper_Energies = np.asarray(Combined_Ring_By_Ring_Hydrogenated_Improper_Energies)

    return Combined_Ring_By_Ring_Hydrogenated_Improper_Energies,Combined_Ring_By_Ring_Hydrogenated_Improper_Energies_Tracker

def Merge_Coarse_Fine(Rotated_Shape,Max_OOP,Max_Dih,Fine_Rotated_Shape,Fine_Max_OOP,Fine_Max_Dih,E_delocalization,E_delocalization_Fine):
    OOP_Rotations_Degrees = np.linspace(0,Max_OOP,Rotated_Shape[0])
    OOP_Rotations_Degrees = OOP_Rotations_Degrees.tolist()
    OOP_Rotations_Degrees_Fine = np.linspace(0,Fine_Max_OOP,Fine_Rotated_Shape[0])
    OOP_Rotations_Degrees_Fine = OOP_Rotations_Degrees_Fine.tolist()

    i = 0
    j = 0
    Merged_OOP_Rotations_Degrees = OOP_Rotations_Degrees + OOP_Rotations_Degrees_Fine
    Merged_OOP_Rotations_Degrees = list(dict.fromkeys(Merged_OOP_Rotations_Degrees))
    Merged_OOP_Rotations_Degrees = sorted(Merged_OOP_Rotations_Degrees)

    print(Merged_OOP_Rotations_Degrees)

    Ring_By_Ring_Merged_Delocalization_Energies = []
    for e_del,e_del_f in zip(E_delocalization,E_delocalization_Fine):
        Merged_Delocalization_Energies = []
        for tors in Merged_OOP_Rotations_Degrees:
            if tors in OOP_Rotations_Degrees:
                Merged_Delocalization_Energies.append(e_del[OOP_Rotations_Degrees.index(tors)])
            else:
                Merged_Delocalization_Energies.append(e_del_f[OOP_Rotations_Degrees_Fine.index(tors)])
        Merged_Delocalization_Energies = np.array(Merged_Delocalization_Energies)
        Ring_By_Ring_Merged_Delocalization_Energies.append(Merged_Delocalization_Energies)

    Ring_By_Ring_Merged_Delocalization_Energies = np.array(Ring_By_Ring_Merged_Delocalization_Energies)

    Merged_OOP_Rotations_Degrees = np.array(Merged_OOP_Rotations_Degrees)

    print(Ring_By_Ring_Merged_Delocalization_Energies.shape)

    return Merged_OOP_Rotations_Degrees,Ring_By_Ring_Merged_Delocalization_Energies

def Return_Dually_Hydrogenated(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_Dual_Hydrogenated_End_File_Matrices,Syn_Anti_Matrices):
    Full_Raw_Energies = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []
    Dually_Hydrogenated_Energies = []

    Folder_Name = "Dual_Hydrogenated_Test"
    Cluster_Login,Base_Cluster_Location,Cluster_Location,Scheduler_Type,End_Condition,Shared_File_Location,qos = Run_Parameters(Folder_Name)
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[1])
    Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
    OOP_Rotations_Degrees = np.linspace(0,Max_OOP + Max_OOP/Rotated_Shape[0]-1,Rotated_Shape[0]+1)

    for ring1,ring2,ring3,End_File_List,Syn_Anti_Matrix in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List,Ring_By_Ring_Dual_Hydrogenated_End_File_Matrices,Syn_Anti_Matrices):
        Raw_Energies = []
        Trimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2,ring3])
        Job_Name = "%s" % Trimer.Name
        Syn_Anti_Energies = (Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Return_Energy_QChem = True,Shared_File_Location = Shared_File_Location)[0])
        Syn_Anti_Energies = Syn_Anti_Energies - np.amin(Syn_Anti_Energies)
        Dually_Hydrogenated_Matrix = []
        for direction_list in Syn_Anti_Matrix:
            Dually_Hydrogenated_List = []
            for direction in direction_list:

                Dually_Hydrogenated_List.append(Syn_Anti_Energies[direction])
            Dually_Hydrogenated_Matrix.append(Dually_Hydrogenated_List)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Dually_Hydrogenated_Matrix,'%s_%s_%s_Dually_Hydrogenated_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Dually_Hydrogenated_Energies.append(Dually_Hydrogenated_Matrix)

    return Dually_Hydrogenated_Energies

def Calculate_Trimer_Long_Range_Energy(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,Ring_By_Ring_Dimer_Delocalization_Energy,Ring_By_Ring_Hydrogenated_Improper_Energies,Ring_By_Ring_Total_Trimer_Energy,Ring_By_Ring_Hydrogenated_Trimer_Energy,Ring_By_Ring_End_File_Matrices,Dually_Hydrogenated_Energies,Ring_By_Ring_Nonbonded_Energies):
    Dih_Rotations_Degrees = np.linspace(0,Max_Dih + Max_Dih/Rotated_Shape[1]-1,Rotated_Shape[1]+1)
    Trimer_Delocalization_Energy = []

    Full_Raw_Energies = []
    Offset_Ring_List = []
    Dual_Offset_Ring_List = []

    Offset_Delocalization_Energy = []

    Dual_Offset_Hydrogenated_Improper_Ring_List = []

    for del_energy in Ring_By_Ring_Dimer_Delocalization_Energy[1:]:
        Offset_Delocalization_Energy.append(del_energy)
    Offset_Delocalization_Energy.append(Ring_By_Ring_Dimer_Delocalization_Energy[0])

    for hyd_ring in Ring_By_Ring_Hydrogenated_Improper_Energies[2:]:
        Dual_Offset_Hydrogenated_Improper_Ring_List.append(hyd_ring)

    Dual_Offset_Hydrogenated_Improper_Ring_List.append(Ring_By_Ring_Hydrogenated_Improper_Energies[0])
    Dual_Offset_Hydrogenated_Improper_Ring_List.append(Ring_By_Ring_Hydrogenated_Improper_Energies[1])

    Offset_Nonbonded_Energies = Ring_By_Ring_Nonbonded_Energies[1:]
    np.concatenate((np.array([Ring_By_Ring_Nonbonded_Energies[0]]),Offset_Nonbonded_Energies))

    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    for ring in Ring_List[2:]:
        Dual_Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[0])
    Dual_Offset_Ring_List.append(Ring_List[1])

    for ring1,ring2,ring3,hyd_ring1,hyd_ring3,End_File_Matrix,trimer_energy_matrix,trimer_nontorsional_matrix,del_energy12,del_energy23,dual_hyd_matrix,nonbonded_dimer_energies,o_nonbonded_dimer_energies,i in zip(Ring_List,Offset_Ring_List,Dual_Offset_Ring_List,Ring_By_Ring_Hydrogenated_Improper_Energies,Dual_Offset_Hydrogenated_Improper_Ring_List,Ring_By_Ring_End_File_Matrices,Ring_By_Ring_Total_Trimer_Energy,Ring_By_Ring_Hydrogenated_Trimer_Energy,Ring_By_Ring_Dimer_Delocalization_Energy,Offset_Delocalization_Energy,Dually_Hydrogenated_Energies,Ring_By_Ring_Nonbonded_Energies,Offset_Nonbonded_Energies,range(len(Ring_By_Ring_Total_Trimer_Energy))):

        Hydrogenated_Matrix_1 = []
        for angle_val in hyd_ring1[0]:
            Hydrogenated_Matrix_1.append(np.ones(len(hyd_ring1[0]))*angle_val)
        Hydrogenated_Matrix_1 = Hydrogenated_Matrix_1 - np.amin(Hydrogenated_Matrix_1)
        Hydrogenated_Matrix_3 = []
        for angle_val in hyd_ring3[0]:
            Hydrogenated_Matrix_3.append(np.ones(len(hyd_ring3[0]))*angle_val)
        Hydrogenated_Matrix_3 = Hydrogenated_Matrix_3 - np.amin(Hydrogenated_Matrix_3)

        Nonbonded_Energy_12 = []
        for angle_val in nonbonded_dimer_energies[0]:
            Nonbonded_Energy_12.append(np.ones(len(nonbonded_dimer_energies[0]))*angle_val)
        Nonbonded_Energy_12 = Nonbonded_Energy_12 - np.amin(Nonbonded_Energy_12)

        Nonbonded_Energy_23 = []
        for angle_val in o_nonbonded_dimer_energies[0]:
            Nonbonded_Energy_23.append(np.ones(len(o_nonbonded_dimer_energies[0]))*angle_val)
        Nonbonded_Energy_23 = Nonbonded_Energy_23 - np.amin(Nonbonded_Energy_23)

        Delocalization_Matrix_12 = []
        for angle_val in del_energy12[0]:
            Delocalization_Matrix_12.append(np.ones(len(del_energy12[0]))*angle_val)

        Delocalization_Matrix_23 = []
        for angle_val in del_energy23[0]:
            Delocalization_Matrix_23.append(np.ones(len(del_energy23[0]))*angle_val)

        if End_File_Matrix[0][0].split("Phi")[0].strip("_") == ring1.Name:
            Hydrogenated_Matrix_1 = np.transpose(Hydrogenated_Matrix_1)
            Delocalization_Matrix_12 = np.transpose(Delocalization_Matrix_12)
            Nonbonded_Energy_12 = np.transpose(Nonbonded_Energy_12)
        elif End_File_Matrix[0][0].split("Phi")[0].strip("_") == ring3.Name:
            Hydrogenated_Matrix_3 = np.transpose(Hydrogenated_Matrix_3)
            Delocalization_Matrix_23 = np.transpose(Delocalization_Matrix_23)
            Nonbonded_Energy_23 = np.transpose(Nonbonded_Energy_23)
        else:
            raise Exception("Name Does Not Match")

        Nonbonded_Energy = trimer_nontorsional_matrix - Hydrogenated_Matrix_1 - Hydrogenated_Matrix_3 - dual_hyd_matrix + Nonbonded_Energy_12 + Nonbonded_Energy_23

        Delocalization_Energy = trimer_energy_matrix - Nonbonded_Energy - Delocalization_Matrix_12 - Delocalization_Matrix_23

        Delocalization_Energy = Delocalization_Energy - np.amin(Delocalization_Energy)

        print(dual_hyd_matrix)

        Ring_By_Ring_Dual_Hydrogenated_End_File_Matrices = Trimer_Delocalization_Energy.append(Delocalization_Energy)

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Delocalization_Energy,'%s_%s_%s_Trimer_Delocalization_Energies' % (ring1.Name,ring2.Name,ring3.Name),Title='RI-MP2 Energies',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

        Make_Surface_Plot(Dih_Rotations_Degrees,Dih_Rotations_Degrees,Nonbonded_Energy,'%s_%s_%s_Trimer_Nonbonded_Energies' % (ring3.Name,ring2.Name,ring1.Name),Title='RI-MP2 Energies',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s Dihedral (degrees)' % (ring2.Name,ring3.Name))

def Weight_Impropers(Raw_Dimer_Energies,Temp):
    Dimer_Weights = []

    for dimer_matrix in Raw_Dimer_Energies:
        Ring_By_Ring_Dimer_Weights = []
        for dimer_row in dimer_matrix:
            weight = 0
            for num in dimer_row:
                weight += np.exp(num/(-Temp*0.001985875))
            Ring_By_Ring_Dimer_Weights.append(weight)
        Dimer_Weights.append(np.array(Ring_By_Ring_Dimer_Weights)/sum(Ring_By_Ring_Dimer_Weights))
    Dimer_Weights = np.array(Dimer_Weights)

    return Dimer_Weights

def Fit_Improper_Energies(Dimer_Delocalization_Energies,Ring_List,Rotated_Shape,Merged_OOP_Rotations_Degrees,Dimer_Weights,Coarse_Grid_Size=(13,72),Nonbonded=False):
    color = np.linspace(0.0,1.0,len(Merged_OOP_Rotations_Degrees))
    trimer_color = np.linspace(0.0,1.0,len(Merged_OOP_Rotations_Degrees))
    cmap = mpl.cm.get_cmap("coolwarm")
    Offset_Ring_List = []

    a_params_list = []
    b_params_list = []
    a0_params_list = []

    Merged_OOP_Rotations_Radians = np.array([np.radians(x) for x in Merged_OOP_Rotations_Degrees])

    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)

    Offset_Ring_List.append(Ring_List[0])

    Fit_Energies_list = []
    Fit_Coarse_Energies_list = []
    Force_x_list = []
    Force_y_list = []

    Ring1_List = []
    Ring2_List = []

    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring1_List.append(ring2)
        Ring2_List.append(ring2)
        Ring2_List.append(ring1)

    for deloc_matrix,ring1,ring2,weights,z in zip(Dimer_Delocalization_Energies,Ring1_List,Ring2_List,Dimer_Weights,range(len(Ring1_List))):
        deloc_matrix = deloc_matrix - np.amin(deloc_matrix)
        Parameter_Lists = []
        fig,ax = plt.subplots(1,1)
        for i,energy_list in enumerate(deloc_matrix):
            if np.mod(i,7) == 0:
                alpha = 1.0
            else:
                alpha = 0.2
            ax.scatter(np.linspace(0,350,Rotated_Shape[1]),energy_list,alpha=alpha,marker = 's',color=cmap(color[i]))
            opt,cov = scipy.optimize.curve_fit(OPLS,np.linspace(0,350,Rotated_Shape[1]),energy_list)
            Parameter_Lists.append(opt)
            x = np.linspace(0,360,10000)
            fit = []
            for j in x:
                fit.append(OPLS(j,opt[0],opt[1],opt[2],opt[3],opt[4]))
            plt.plot(x,fit,color=cmap(color[i]),alpha=alpha)

        plt.xlabel('Dihedral Angle ($^\circ$)',size = 24)
        plt.ylabel('Energy (kcal/mol)',size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        #plt.ylim([-1,15])
        plt.tight_layout()
        if Nonbonded:
            fig.savefig('%s_%s_Conjugated_Energies_Nonbonded' % (ring1.Name,ring2.Name))
            plt.close(fig)
            os.system("scp %s_%s_Conjugated_Energies_Nonbonded.png ./Figures/" % (ring1.Name,ring2.Name))
            os.remove("%s_%s_Conjugated_Energies_Nonbonded.png" % (ring1.Name,ring2.Name))
        else:
            fig.savefig('%s_%s_Conjugated_Energies' % (ring1.Name,ring2.Name))
            plt.close(fig)
            os.system("scp %s_%s_Conjugated_Energies.png ./Figures/" % (ring1.Name,ring2.Name))
            os.remove("%s_%s_Conjugated_Energies.png" % (ring1.Name,ring2.Name))

        opls_params = []
        a_params = []
        b_params = []
        a0_params = []
        fourier_order = 6
        for j in range(len(opt)):
            x, y = variables('x, y')
            w, = parameters('w')
            #model_dict = {y: fourier_series(x, f=w, n=6)}
            model_dict = {y: fourier_series(x, f=1, n=fourier_order)}
            fig,ax = plt.subplots(1,1)
            plt.scatter(Merged_OOP_Rotations_Radians,[q[j] for q in Parameter_Lists],color=cmap(color[i]),alpha=alpha)
            plt.scatter(-1*Merged_OOP_Rotations_Radians[:0:-1],[q[j] for q in Parameter_Lists][:0:-1],color=cmap(color[i]),alpha=alpha)
            Mirror_OOP_Radians = np.concatenate((-1*Merged_OOP_Rotations_Radians[:0:-1],Merged_OOP_Rotations_Radians))
            Mirror_y_coords = [q[j] for q in Parameter_Lists][:0:-1] + [p[j] for p in Parameter_Lists]
            Mirror_weights = [1/q for q in (weights[:0:-1])] + [1/p for p in weights]
            #Mirror_weights = [q*0+1 for q in (weights[:0:-1])] + [q*0+1 for p in weights] #TODO: Figure out whats the deal with this q*0 thing
            fit_obj = Fit(model_dict, x=Mirror_OOP_Radians, y=Mirror_y_coords,sigma_y=Mirror_weights,absolute_sigma=False)
            fit_result = fit_obj.execute()
            q = np.linspace(0,math.pi,1000)
            plt.plot(q, fit_obj.model(x=q, **fit_result.params).y, color='green', ls=':')
            opls_params.append(np.array(fit_obj.model(x=q, **fit_result.params).y))
            plt.xlabel('OPLS Parameter %d' % j,size = 24)
            plt.ylabel('OOP Rotation (Degree)',size = 24)
            plt.ylim([np.amin([q[j] for q in Parameter_Lists]) - .15*(np.amax([q[j] for q in Parameter_Lists]) - np.amin([q[j] for q in Parameter_Lists])),np.amax([q[j] for q in Parameter_Lists]) + .15*(np.amax([q[j] for q in Parameter_Lists]) - np.amin([q[j] for q in Parameter_Lists]))])
            ax.tick_params(axis="x", labelsize=18)
            ax.tick_params(axis="y", labelsize=18)
            ax.tick_params(length=4,width=4)
            plt.tight_layout()
            fig.savefig('%s_%s_Improper_OPLS_Parameters_%d_Fourier_Fit' % (ring1.Name,ring2.Name,j))
            plt.close(fig)
            os.system("scp %s_%s_Improper_OPLS_Parameters_%d_Fourier_Fit.png ./Figures/" % (ring1.Name,ring2.Name,j))
            os.remove("%s_%s_Improper_OPLS_Parameters_%d_Fourier_Fit.png" % (ring1.Name,ring2.Name,j))
            a_params.append([fit_result.params['a%d' % z] for z in range(1,fourier_order + 1)])
            b_params.append([fit_result.params['b%d' % z] for z in range(1,fourier_order + 1)])
            a0_params.append(fit_result.params['a0'])
        # print(a_params)
        # print(b_params)
        # print(a0_params)
        # print(Mirror_weights)

        Fit_Matrix = []
        p = np.linspace(0,360,1000)
        for a,b,c,d,e in zip(opls_params[0][:300],opls_params[1][:300],opls_params[2][:300],opls_params[3][:300],opls_params[4][:300]):
            Fit_List = []
            for dih in p:
                Fit_List.append(OPLS(dih,a,b,c,d,e))
            Fit_Matrix.append(Fit_List)

        print("0 OOP Dihedral Twist energy:")
        print(OPLS(90,fourier_nonfit(0,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0,a_params[4],b_params[4],a0_params[4],n=fourier_order)) - OPLS(0,fourier_nonfit(0,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0,a_params[4],b_params[4],a0_params[4],n=fourier_order)))
        print("30 degree OOP Dihedral Twist energy:")
        print(OPLS(90,fourier_nonfit(0.523599,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0.523599,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0.523599,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0.523599,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0.523599,a_params[4],b_params[4],a0_params[4],n=fourier_order)) - OPLS(0,fourier_nonfit(0.523599,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0.523599,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0.523599,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0.523599,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0.523599,a_params[4],b_params[4],a0_params[4],n=fourier_order)))
        print("30 degree vs 0 degree OOP energy:")
        print(OPLS(0,fourier_nonfit(0.523599,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0.523599,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0.523599,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0.523599,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0.523599,a_params[4],b_params[4],a0_params[4],n=fourier_order)) - OPLS(0,fourier_nonfit(0,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0,a_params[4],b_params[4],a0_params[4],n=fourier_order)))
        Make_Surface_Plot(p,np.linspace(0,math.pi,1000)[:300],Fit_Matrix,'%s_%s_Fit_Energies_%d' % (ring1.Name,ring2.Name,z),Title='RI-MP2 Energies',xlabel='%s %s Dihedral (degrees)' % (ring1.Name,ring2.Name),ylabel='%s %s OOP (degrees)' % (ring1.Name,ring2.Name),Tight_Layout = False)

        p = np.linspace(0,360-360/Coarse_Grid_Size[1],Coarse_Grid_Size[1])
        q = np.linspace(0,math.pi/3-math.pi/(3*Coarse_Grid_Size[0]),Coarse_Grid_Size[0])
        U = np.zeros(Coarse_Grid_Size)
        V = np.zeros(Coarse_Grid_Size)
        Coarse_Energy = np.zeros(Coarse_Grid_Size)
        for oop_num,oop in enumerate(q):
            for dih_num,dih in enumerate(p):
                U[oop_num][dih_num] = -1 * OPLS_Derivative(dih,fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order))
                print(oop)
                print(a_params[4])
                print(b_params[4])
                print("%.4f %.4f %.4f %.4f" % (fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order)))
                V[oop_num][dih_num] = -1 * OPLS(dih,fourier_series_derivative(oop,a_params[0],b_params[0],n=fourier_order),fourier_series_derivative(oop,a_params[1],b_params[1],n=fourier_order),fourier_series_derivative(oop,a_params[2],b_params[2],n=fourier_order),fourier_series_derivative(oop,a_params[3],b_params[3],n=fourier_order),fourier_series_derivative(oop,a_params[4],b_params[4],n=fourier_order))
                Coarse_Energy[oop_num][dih_num] = 4.184 * OPLS(dih,fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(oop,a_params[4],b_params[4],a0_params[4],n=fourier_order))



        #print(U)
        #print(V)
        fig,ax = plt.subplots(1,1)
        plt.quiver(np.linspace(0,360-360/Coarse_Grid_Size[1],Coarse_Grid_Size[1]),np.linspace(0,math.pi/3-math.pi/(3*Coarse_Grid_Size[0]),Coarse_Grid_Size[0]),U,V)
        ax.set_title("Forces",fontdict = {'fontsize':24})
        plt.xlabel('Dih',size = 24)
        plt.ylabel('OOP',size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        fig.savefig("Test_Quiver")
        fig.tight_layout()

        plt.close(fig)

        p = np.linspace(-180,180,200)
        q = np.linspace(0,math.pi/2,200)
        Force_x = np.zeros((200,200))
        Force_y = np.zeros((200,200))
        Energy = np.zeros((200,200))
        counter = 1#TODO: temporary, delete later
        print("should run for %d loop"%(len(q)*len(p)))
        for oop_num,oop in enumerate(q):
            for dih_num,dih in enumerate(p):
                Force_x[oop_num][dih_num] = -1 * 4.184 * OPLS_Derivative(dih,fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order))
                Force_y[oop_num][dih_num] = -1 * 4.184 * OPLS(dih,fourier_series_derivative(oop,a_params[0],b_params[0],n=fourier_order),fourier_series_derivative(oop,a_params[1],b_params[1],n=fourier_order),fourier_series_derivative(oop,a_params[2],b_params[2],n=fourier_order),fourier_series_derivative(oop,a_params[3],b_params[3],n=fourier_order),fourier_series_derivative(oop,a_params[4],b_params[4],n=fourier_order))
                Energy[oop_num][dih_num] = 4.184 * OPLS(dih,fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(oop,a_params[4],b_params[4],a0_params[4],n=fourier_order))
                print("working in loop # %d"%counter)
                counter += 1

        a_params_list.append(a_params)
        b_params_list.append(b_params)
        a0_params_list.append(a0_params)
        Fit_Energies_list.append(Energy)
        Fit_Coarse_Energies_list.append(Coarse_Energy)
        Force_x_list.append(Force_x)
        Force_y_list.append(Force_y)

    return Force_x_list,Force_y_list,Fit_Energies_list,Fit_Coarse_Energies_list,a_params_list,b_params_list,a0_params_list

def Fit_Conjugated_Energies(Conjugated_Energies,OOP_Angles):
    current_params = [0,0,0]
    Conjugated_Parameters = []
    for i,energy_matrix in enumerate(Conjugated_Energies):
        best_fit = 100000000000
        for n in [0,1,2,3,4,6]:
            for d in [1,-1]:
                optimized = np.array([e[0] for e in energy_matrix])
                K,_ = scipy.optimize.curve_fit(lambda x,K: CVFF(x,K,d,n),OOP_Angles,optimized)
                K = K[0]
                if sum([(CVFF(x,K,d,n) - energy)**2 for energy,x in zip([e[0] for e in energy_matrix],OOP_Angles)]) < best_fit:
                    current_params = [K,d,n]
                    best_fit = sum([(CVFF(x,K,d,n) - energy)**2 for energy,x in zip([e[0] for e in energy_matrix],OOP_Angles)])
        Conjugated_Parameters.append(current_params)
        fig,ax = plt.subplots(1,1)
        plt.scatter(OOP_Angles,[e[0] for e in energy_matrix])
        plt.plot(np.linspace(OOP_Angles[0],OOP_Angles[-1],1000),[CVFF(x,current_params[0],current_params[1],current_params[2]) for x in np.linspace(OOP_Angles[0],OOP_Angles[-1],1000)])
        plt.savefig("Conjugated_Energies_Fit_%d" % i)
        plt.close(fig)
        print(current_params)
    return Conjugated_Parameters

def Check_New_Torsional_Energy(Ring_List,Polymer_Name,Rotated_Shape,a_params_list,b_params_list,a0_params_list,OPLS_Fit_list,Conjugated_Parameters,Nonbonded_a_params_list,Nonbonded_b_params_list,Nonbonded_a0_params_list):
    if Polymer_Name == "P3HT_Input":
        Sub_Directory_List = ["./Dimer_Rigid/","./Dimer_Rigid/"]
        Colvar_File = "colvar_1_No_Torsion.txt"
        Planar_File_List = ["P3HT_Dimer_800_0_Rigid.data","P3HT_Dimer_800_0_Rigid.data"]
        Quantum_Planar_File_List = ["P3HT_Dimer_Verification_Rigid0.out","P3HT_Dimer_Verification_Rigid0.out"]
        OPLS_Dih_List_list = [[3,5,12,13],[3,5,12,13]]
        File_Name_Fragment_1_List = ["P3HT_Dimer_800_","P3HT_Dimer_800_"]
        File_Name_Fragment_2_List = ["_Rigid.data","_Rigid.data"]
        Improper_1_List_list = [[12,13,17,5],[12,13,17,5]] #P3HT
        Improper_2_List_list = [[5,6,3,12],[5,6,3,12]] #P3HT
        Out_File_Fragment_1_List = ["P3HT_Dimer_Verification_Rigid","P3HT_Dimer_Verification_Rigid"]
        Out_File_Fragment_2_List = [".out",".out"]
        Job_Name = "P3HT_Dimer_Verification"
        LAMMPS_Name = "P3HT"

    elif Polymer_Name == "PTB7_Input":
        Sub_Directory_List = ["./Dimer_Plumed_PTB7/F_Outside/","./Dimer_Plumed_PTB7/F_Inside/"]
        Colvar_File = "colvar_1.txt"
        Planar_File_List = ["PTB7_Dimer_F_Outside_800_0_No_Torsion.data","PTB7_Dimer_F_Inside_800_0_No_Torsion.data"]
        Quantum_Planar_File_List = ["PTB7_Dimer_Verification_Rigid_0.out","PTB7_Dimer_Verification_Rigid_F_Inside_0.out"]
        OPLS_Dih_List_list = [[3,5,12,13],[3,5,12,13]]
        File_Name_Fragment_1_List = ["PTB7_Dimer_F_Outside_800_","PTB7_Dimer_F_Inside_800_"]
        File_Name_Fragment_2_List = ["_No_Torsion.data","_No_Torsion.data"]
        Improper_1_List_list = [[23,22,21,26],[3,2,5,18]] #PTB7 may not be accurate
        Improper_2_List_list = [[26,27,29,23],[18,19,20,3]]
        Out_File_Fragment_1_List = ["PTB7_Dimer_Verification_Rigid_","PTB7_Dimer_Verification_Rigid_F_Inside_"]
        Out_File_Fragment_2_List = [".out",".out",".out",".out"]
        Job_Name = "PTB7_Dimer_Verification"
        LAMMPS_Name = "PTB7"

    elif Polymer_Name == "N2200_Input":
        Sub_Directory_List = ["./Dimer_Plumed_N2200/Naph_Thio/","./Dimer_Plumed_N2200/Naph_Thio/","./Dimer_Plumed_N2200/Naph_Thio/","./Dimer_Plumed_N2200/Naph_Thio/","./Dimer_Plumed_N2200/Naph_Thio/","./Dimer_Plumed_N2200/Naph_Thio/"]
        Colvar_File = "colvar_1.txt"
        Planar_File_List = ["N2200_Dimer_Naph_Thio_800_0_No_Torsion_Alternate_Plus_Charge.data","N2200_Dimer_Naph_Thio_800_0_No_Torsion_Alternate_Plus_Charge.data","N2200_Dimer_Naph_Thio_800_0_No_Torsion_Alternate_Plus_Charge.data","N2200_Dimer_Naph_Thio_800_0_No_Torsion_Alternate_Plus_Charge.data","N2200_Dimer_Naph_Thio_800_0_No_Torsion_Alternate_Plus_Charge.data","N2200_Dimer_Naph_Thio_800_0_No_Torsion_Alternate_Plus_Charge.data"]
        Quantum_Planar_File_List = ["N2200_Dimer_Verification_Rigid_Naph_Thio_0.out","N2200_Dimer_Verification_Rigid_Naph_Thio_0.out","N2200_Dimer_Verification_Rigid_Naph_Thio_0.out","N2200_Dimer_Verification_Rigid_Naph_Thio_0.out","N2200_Dimer_Verification_Rigid_Naph_Thio_0.out","N2200_Dimer_Verification_Rigid_Naph_Thio_0.out"]
        OPLS_Dih_List_list = [[10,11,32,33],[3,5,11,13],[3,5,9,10],[3,5,9,10],[3,5,11,13],[10,11,32,33]]
        File_Name_Fragment_1_List = ["N2200_Dimer_Naph_Thio_800_","N2200_Dimer_Naph_Thio_800_","N2200_Dimer_Naph_Thio_800_","N2200_Dimer_Naph_Thio_800_","N2200_Dimer_Naph_Thio_800_","N2200_Dimer_Naph_Thio_800_"]
        File_Name_Fragment_2_List = ["_No_Torsion_Alternate_Plus_Charge.data","_No_Torsion_Alternate_Plus_Charge.data","_No_Torsion_Alternate_Plus_Charge.data","_No_Torsion_Alternate_Plus_Charge.data","_No_Torsion_Alternate_Plus_Charge.data","_No_Torsion_Alternate_Plus_Charge.data"]
        Improper_1_List_list = [[32,33,38,11],[11,8,10,32],[32,33,38,11],[32,33,38,11],[11,8,10,32],[32,33,38,11]] #N2200
        Improper_2_List_list = [[11,8,10,32],[32,33,38,11],[32,33,38,11],[32,33,38,11],[32,33,38,11],[11,8,10,32]] #N2200
        Out_File_Fragment_1_List = ["N2200_Dimer_Verification_Rigid_Naph_Thio_","N2200_Dimer_Verification_Rigid_Naph_Thio_","N2200_Dimer_Verification_Rigid_Naph_Thio_","N2200_Dimer_Verification_Rigid_Naph_Thio_","N2200_Dimer_Verification_Rigid_Naph_Thio_","N2200_Dimer_Verification_Rigid_Naph_Thio_"]
        Out_File_Fragment_2_List = [".out",".out",".out",".out",".out",".out"]
        Job_Name = "N2200_Dimer_Verification"
        LAMMPS_Name = "N2200"

    else:
        raise Exception("Undefined Polymer Name")

    Improper_Params_1_list = Conjugated_Parameters #N2200
    Improper_Params_2_list = []
    for i in range(0,len(Conjugated_Parameters),2):
        Improper_Params_2_list.append(Conjugated_Parameters[i+1])
        Improper_Params_2_list.append(Conjugated_Parameters[i])

    LAMMPS_Start_Timestep = 0
    LAMMPS_End_Timestep = 9900000
    LAMMPS_Step = 100000
    Colvar_Start_Timestep = 0.0
    Colvar_End_Timestep = 9900.0
    Colvar_Step =  100.0
    Offset_Ring_List = []
    fourier_order = 6
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))
    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))

    Ring1_List = []
    Ring2_List = []

    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring1_List.append(ring2)
        Ring2_List.append(ring2)
        Ring2_List.append(ring1)

    for File_Name_Fragment_1,File_Name_Fragment_2,Out_File_Fragment_1,Out_File_Fragment_2,OPLS_Dih_List,OPLS_Fit,a_params,b_params,a0_params,ring1,ring2,Sub_Directory,Planar_File,Quantum_Planar_File,Improper_1_List,Improper_2_List,Improper_Params_1,Improper_Params_2,Nonbonded_a_params,Nonbonded_b_params,Nonbonded_a0_params,q in zip(File_Name_Fragment_1_List,File_Name_Fragment_2_List,Out_File_Fragment_1_List,Out_File_Fragment_2_List,OPLS_Dih_List_list,OPLS_Fit_list,a_params_list,b_params_list,a0_params_list,Ring1_List,Ring2_List,Sub_Directory_List,Planar_File_List,Quantum_Planar_File_List,Improper_1_List_list,Improper_2_List_list,Improper_Params_1_list,Improper_Params_2_list,Nonbonded_a_params_list,Nonbonded_b_params_list,Nonbonded_a0_params_list,range(len(File_Name_Fragment_1_List))):
        File_Name_List = []
        Out_File_List = []
        Colvar_OOP_List = []
        Colvar_DIH_List = []

        LAMMPS_Current_Timestep = LAMMPS_Start_Timestep
        while LAMMPS_Current_Timestep <= LAMMPS_End_Timestep:
            File_Name = (File_Name_Fragment_1 + "%d" + File_Name_Fragment_2) % (LAMMPS_Current_Timestep)
            File_Name_List.append(File_Name)
            Out_File_Name = (Out_File_Fragment_1 + "%d" + Out_File_Fragment_2) % (LAMMPS_Current_Timestep)
            Out_File_List.append(Out_File_Name)
            LAMMPS_Current_Timestep += LAMMPS_Step

        Colvar_Current_Timestep = Colvar_Start_Timestep
        f = open(Sub_Directory + Colvar_File,'r')
        Colvars = f.readlines()
        while Colvar_Current_Timestep <= Colvar_End_Timestep:
            for line in Colvars[3:]:
                if len(line.strip().split()) == 3 and all((ch.isdigit() or ch == ".") for ch in line.strip().split()[0].strip()) and float(line.strip().split()[0].strip()) == Colvar_Current_Timestep:
                    Colvar_DIH_List.append(float(line.strip().split()[1].strip()))
                    Colvar_OOP_List.append(float(line.strip().split()[2].strip()))
                    Colvar_Current_Timestep += Colvar_Step
                    break

        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        #Dimer.Read_From_Data_File("./LigParGen_Files/%s_%s.lmp" % (ring1.Name,ring2.Name))
        Dimer.Read_From_Data_File(Sub_Directory + Planar_File)
        Dimer.Unwrap_Coordinates(50.0,300.0)
        Nonbonded_Planar_Energy = Dimer.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/lmp_serial",dielectric=4.9,Exclude_Interring_Torsions=True)
        Nonbonded_Fit_Planar_Energy = OPLS(0,fourier_nonfit(0,Nonbonded_a_params[0],Nonbonded_b_params[0],Nonbonded_a0_params[0],n=fourier_order),fourier_nonfit(0,Nonbonded_a_params[1],Nonbonded_b_params[1],Nonbonded_a0_params[1],n=fourier_order),fourier_nonfit(0,Nonbonded_a_params[2],Nonbonded_b_params[2],Nonbonded_a0_params[2],n=fourier_order),fourier_nonfit(0,Nonbonded_a_params[3],Nonbonded_b_params[3],Nonbonded_a0_params[3],n=fourier_order),fourier_nonfit(0,Nonbonded_a_params[4],Nonbonded_b_params[4],Nonbonded_a0_params[4],n=fourier_order))
        Colvar_Planar_Energy = Nonbonded_Planar_Energy + OPLS(0,fourier_nonfit(0,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0,a_params[4],b_params[4],a0_params[4],n=fourier_order))
        Nonbonded_Fit_Colvar_Planar_Energy = Nonbonded_Fit_Planar_Energy + OPLS(0,fourier_nonfit(0,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(0,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(0,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(0,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(0,a_params[4],b_params[4],a0_params[4],n=fourier_order))
        Control_Planar_Energy = Nonbonded_Planar_Energy + OPLS_LAMMPS(0,OPLS_Fit[0],OPLS_Fit[1],OPLS_Fit[2],OPLS_Fit[3])

        Colvar_Energy_List = []
        Control_Energy_List = []
        Nonbonded_Fit_Colvar_Energy_List = []

        Quantum_Planar_Energy = Cluster_IO.Return_Info(Quantum_Planar_File,Quantum_Planar_File,Sub_Directory,Job_Name,"","",Shared_File_Location = Sub_Directory,Return_Energy_QChem = True)[0]
        Quantum_Energies = Cluster_IO.Return_Info_Batch(Out_File_List,Out_File_List,Sub_Directory,Job_Name,"","",Shared_File_Location = Sub_Directory,Return_Energy_QChem = True)[0]

        Nontorsional_Energies = []
        Conv_Dih_List = []
        Improper_1_Angles = []
        Improper_2_Angles = []
        for File_Name,Out_File_Name,oop,dih_rad in zip(File_Name_List,Out_File_List,Colvar_OOP_List,Colvar_DIH_List):
            print(File_Name)
            dih = dih_rad * 180 / math.pi

            Dimer.Read_From_Data_File(Sub_Directory + File_Name)
            Dimer.Unwrap_Coordinates(50.0,300.0)
            vec_1 = Dimer.Get_Atom(OPLS_Dih_List[0]).Position - Dimer.Get_Atom(OPLS_Dih_List[1]).Position
            vec_1 = vec_1/np.linalg.norm(vec_1)
            vec_2 = Dimer.Get_Atom(OPLS_Dih_List[1]).Position - Dimer.Get_Atom(OPLS_Dih_List[2]).Position
            vec_2 = vec_2/np.linalg.norm(vec_2)
            vec_3 = Dimer.Get_Atom(OPLS_Dih_List[2]).Position - Dimer.Get_Atom(OPLS_Dih_List[3]).Position
            vec_3 = vec_3/np.linalg.norm(vec_3)

            plane_1 = np.cross(vec_1,vec_2)/np.linalg.norm(np.cross(vec_1,vec_2))
            plane_2 = np.cross(vec_2,vec_3)/np.linalg.norm(np.cross(vec_2,vec_3))
            control_dih = np.arccos(np.dot(plane_1,plane_2)) * 180 / math.pi
            Conv_Dih_List.append(control_dih)

            vec_1 = Dimer.Get_Atom(Improper_1_List[0]).Position - Dimer.Get_Atom(Improper_1_List[1]).Position
            vec_1 = vec_1/np.linalg.norm(vec_1)
            vec_2 = Dimer.Get_Atom(Improper_1_List[1]).Position - Dimer.Get_Atom(Improper_1_List[2]).Position
            vec_2 = vec_2/np.linalg.norm(vec_2)
            vec_3 = Dimer.Get_Atom(Improper_1_List[2]).Position - Dimer.Get_Atom(Improper_1_List[3]).Position
            vec_3 = vec_3/np.linalg.norm(vec_3)

            plane_1 = np.cross(vec_1,vec_2)/np.linalg.norm(np.cross(vec_1,vec_2))
            plane_2 = np.cross(vec_2,vec_3)/np.linalg.norm(np.cross(vec_2,vec_3))
            imp_1 = np.arccos(np.dot(plane_1,plane_2)) * 180 / math.pi
            Improper_1_Angles.append(imp_1)

            vec_1 = Dimer.Get_Atom(Improper_2_List[0]).Position - Dimer.Get_Atom(Improper_2_List[1]).Position
            vec_1 = vec_1/np.linalg.norm(vec_1)
            vec_2 = Dimer.Get_Atom(Improper_2_List[1]).Position - Dimer.Get_Atom(Improper_2_List[2]).Position
            vec_2 = vec_2/np.linalg.norm(vec_2)
            vec_3 = Dimer.Get_Atom(Improper_2_List[2]).Position - Dimer.Get_Atom(Improper_2_List[3]).Position
            vec_3 = vec_3/np.linalg.norm(vec_3)

            plane_1 = np.cross(vec_1,vec_2)/np.linalg.norm(np.cross(vec_1,vec_2))
            plane_2 = np.cross(vec_2,vec_3)/np.linalg.norm(np.cross(vec_2,vec_3))
            imp_2 = np.arccos(np.dot(plane_1,plane_2)) * 180 / math.pi
            Improper_2_Angles.append(imp_2)

            Dimer.Name = File_Name.split('.data')[0]
            print([l.Position for l in Dimer.Atom_List])
            print(Dimer.Name)
            Nontorsional_Energy = Dimer.Calculate_Internal_Energy(LAMMPS_Name,"~/lammps-11Aug17/src/lmp_serial",dielectric=4.9,Exclude_Interring_Torsions=True)
            print(Nontorsional_Energy)
            Nontorsional_Energies.append(Nontorsional_Energy)
            Nontorsional_Fit_Energy = OPLS(dih,fourier_nonfit(oop,Nonbonded_a_params[0],Nonbonded_b_params[0],Nonbonded_a0_params[0],n=fourier_order),fourier_nonfit(oop,Nonbonded_a_params[1],Nonbonded_b_params[1],Nonbonded_a0_params[1],n=fourier_order),fourier_nonfit(oop,Nonbonded_a_params[2],Nonbonded_b_params[2],Nonbonded_a0_params[2],n=fourier_order),fourier_nonfit(oop,Nonbonded_a_params[3],Nonbonded_b_params[3],Nonbonded_a0_params[3],n=fourier_order),fourier_nonfit(oop,Nonbonded_a_params[4],Nonbonded_b_params[4],Nonbonded_a0_params[4],n=fourier_order))
            Colvar_Total_Energy = Nontorsional_Energy + OPLS(dih,fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(oop,a_params[4],b_params[4],a0_params[4],n=fourier_order)) + CVFF(Improper_1_Angles[-1],Improper_Params_1[0],Improper_Params_1[1],Improper_Params_1[2]) + CVFF(Improper_2_Angles[-1],Improper_Params_2[0],Improper_Params_2[1],Improper_Params_2[2])

            Nonbonded_Fit_Colvar_Energy = Nontorsional_Fit_Energy + OPLS(dih,fourier_nonfit(oop,a_params[0],b_params[0],a0_params[0],n=fourier_order),fourier_nonfit(oop,a_params[1],b_params[1],a0_params[1],n=fourier_order),fourier_nonfit(oop,a_params[2],b_params[2],a0_params[2],n=fourier_order),fourier_nonfit(oop,a_params[3],b_params[3],a0_params[3],n=fourier_order),fourier_nonfit(oop,a_params[4],b_params[4],a0_params[4],n=fourier_order)) + CVFF(Improper_1_Angles[-1],Improper_Params_1[0],Improper_Params_1[1],Improper_Params_1[2]) + CVFF(Improper_2_Angles[-1],Improper_Params_2[0],Improper_Params_2[1],Improper_Params_2[2])

            Control_Total_Energy = Nontorsional_Energy + OPLS_LAMMPS(control_dih,OPLS_Fit[0],OPLS_Fit[1],OPLS_Fit[2],OPLS_Fit[3])

            Control_Energy_List.append(Control_Total_Energy)
            Nonbonded_Fit_Colvar_Energy_List.append(Nonbonded_Fit_Colvar_Energy)
            Colvar_Energy_List.append(Colvar_Total_Energy)

        Colvar_Energy_List = np.array(Colvar_Energy_List) - Colvar_Planar_Energy
        Control_Energy_List = np.array(Control_Energy_List) - Control_Planar_Energy
        Nonbonded_Fit_Colvar_Energy_List = np.array(Nonbonded_Fit_Colvar_Energy_List) - Nonbonded_Fit_Colvar_Planar_Energy
        Quantum_Energies = np.array(Quantum_Energies) - Quantum_Planar_Energy

        Energy_Diff = Colvar_Energy_List - Control_Energy_List
        Energy_Diff = (Energy_Diff + 2.0)/4.0
        Energy_Diff = [float(e) for e in Energy_Diff]

        Nonbonded_Energy_Diff = Nonbonded_Fit_Colvar_Energy_List - Colvar_Energy_List
        Nonbonded_Energy_Diff = [float(e) for e in Nonbonded_Energy_Diff]

        Relative_Accuracies = []
        Nonbonded_Fit_Relative_Accuracies = []
        dih_degrees = np.array(Colvar_DIH_List) * 180 / math.pi
        oop_degrees = np.array(Colvar_OOP_List) * 180 / math.pi
        OOP_Diffs = []
        Dih_Diffs = []
        Histogram_Error_Lists = np.empty((6,6),dtype=object)
        for hist_row in Histogram_Error_Lists:
            for hist_cell in hist_row:
                hist_cell = []

        Histogram_Nonbonded_Fit_Error_Lists = copy.deepcopy(Histogram_Error_Lists)

        Out_File1_List = []
        Out_File2_List = []
        Colvar_OOP_List1 = []
        Colvar_OOP_List2 = []
        Colvar_Dih_List1 = []
        Colvar_Dih_List2 = []
        Conv_Dih_List1 = []
        Conv_Dih_List2 = []

        Quantum_Relative_Energies = []
        Control_Relative_Energies = []
        Colvar_Relative_Energies = []
        Nonbonded_Fit_Colvar_Relative_Energies = []

        for quantum1,colvar1,control1,nf_colvar1,z,dih1,oop1,conv_dih1,out_file1 in zip(Quantum_Energies[:-1],Colvar_Energy_List[:-1],Control_Energy_List[:-1],Nonbonded_Fit_Colvar_Energy_List[:-1],range(len(Quantum_Energies)),dih_degrees[:-1],oop_degrees[:-1],Conv_Dih_List[:-1],Out_File_List[:-1]):
            for quantum2,colvar2,control2,nf_colvar2,dih2,oop2,conv_dih2,out_file2 in zip(Quantum_Energies[z+1:],Colvar_Energy_List[z+1:],Control_Energy_List[z+1:],Nonbonded_Fit_Colvar_Energy_List[z+1:],dih_degrees[z+1:],oop_degrees[z+1:],Conv_Dih_List[z+1:],Out_File_List[z+1:]):
                quantum_relative_energy = quantum1 - quantum2
                Quantum_Relative_Energies.append(quantum_relative_energy)
                colvar_relative_energy = colvar1 - colvar2
                Colvar_Relative_Energies.append(colvar_relative_energy)
                control_relative_energy = control1 - control2
                Control_Relative_Energies.append(control_relative_energy)
                nf_colvar_relative_energy = nf_colvar1 - nf_colvar2
                Nonbonded_Fit_Colvar_Relative_Energies.append(nf_colvar_relative_energy)
                colvar_error = abs(quantum_relative_energy - colvar_relative_energy)
                control_error = abs(quantum_relative_energy - control_relative_energy)
                nf_colvar_error = abs(quantum_relative_energy - nf_colvar_relative_energy)
                Relative_Accuracies.append(colvar_error - control_error) #negative is good
                Nonbonded_Fit_Relative_Accuracies.append(nf_colvar_error - colvar_error)

                Out_File1_List.append(out_file1)
                Out_File2_List.append(out_file2)
                Colvar_OOP_List1.append(oop1)
                Colvar_OOP_List2.append(oop2)
                Colvar_Dih_List1.append(dih1)
                Colvar_Dih_List2.append(dih2)
                Conv_Dih_List1.append(conv_dih1)
                Conv_Dih_List2.append(conv_dih2)
                if abs(oop1-oop2) < 180:
                    OOP_Diffs.append(abs(oop1-oop2))
                else:
                    OOP_Diffs.append(360-abs(oop1-oop2))
                if abs(dih1-dih2) < 180:
                    Dih_Diffs.append(abs(dih1-dih2))
                else:
                    Dih_Diffs.append(360-abs(dih1-dih2))
                print(oop1-oop2)
                print(OOP_Diffs[-1])
                print(Dih_Diffs[-1])
                try:
                    Histogram_Error_Lists[int(math.floor(OOP_Diffs[-1]/5.0))][int(math.floor(Dih_Diffs[-1]/30.0))].append(colvar_error - control_error)
                except:
                    if int(math.floor(Dih_Diffs[-1]/30.0)) < Histogram_Error_Lists.shape[1] and int(math.floor(OOP_Diffs[-1]/5.0)) < Histogram_Error_Lists.shape[0]:
                        Histogram_Error_Lists[int(math.floor(OOP_Diffs[-1]/5.0))][int(math.floor(Dih_Diffs[-1]/30.0))] = [colvar_error - control_error]
                    else:
                        continue

                try:
                    Histogram_Nonbonded_Fit_Error_Lists[int(math.floor(OOP_Diffs[-1]/5.0))][int(math.floor(Dih_Diffs[-1]/30.0))].append(nf_colvar_error - colvar_error)
                except:
                    if int(math.floor(Dih_Diffs[-1]/30.0)) < Histogram_Nonbonded_Fit_Error_Lists.shape[1] and int(math.floor(OOP_Diffs[-1]/5.0)) < Histogram_Nonbonded_Fit_Error_Lists.shape[0]:
                        Histogram_Nonbonded_Fit_Error_Lists[int(math.floor(OOP_Diffs[-1]/5.0))][int(math.floor(Dih_Diffs[-1]/30.0))] = [nf_colvar_error - colvar_error]
                    else:
                        continue

        Histogram_Errors = np.zeros(Histogram_Error_Lists.shape,dtype=float)
        Histogram_Nonbonded_Fit_Errors = np.zeros(Histogram_Nonbonded_Fit_Error_Lists.shape,dtype=float)
        Histogram_Frequency = np.zeros(Histogram_Error_Lists.shape,dtype=float)
        for ave_row,list_row,nf_list_row,i in zip(Histogram_Errors,Histogram_Error_Lists,Histogram_Nonbonded_Fit_Error_Lists,range(Histogram_Errors.shape[0])):
            for ave,l,nf_l,j in zip(ave_row,list_row,nf_list_row,range(Histogram_Errors.shape[0])):
                if l != None:
                    Histogram_Errors[i][j] = np.mean(np.array(l))
                    Histogram_Frequency[i][j] = len(l)
                if nf_l != None:
                    Histogram_Nonbonded_Fit_Errors[i][j] = np.mean(np.array(nf_l))
        Histogram_Frequency = Histogram_Frequency/np.sum(Histogram_Frequency)

        Files = zip(Relative_Accuracies,Out_File1_List,Out_File2_List,Colvar_OOP_List1,Colvar_OOP_List2,Colvar_Dih_List1,Colvar_Dih_List2,Conv_Dih_List1,Conv_Dih_List2,Quantum_Relative_Energies,Colvar_Relative_Energies,Control_Relative_Energies)
        Files = sorted(Files)
        f = open("Relative_Accuracies_Best_%s_%d.txt" % (Polymer_Name,q),'w')
        for combo in Files[:100]:
            f.write("%.4f %s %s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n" % (combo[0],combo[1],combo[2],combo[3],combo[4],combo[5],combo[6],combo[7],combo[8],combo[9],combo[10],combo[11]))
        f.close()

        f = open("Relative_Accuracies_Worst_%s_%d.txt" % (Polymer_Name,q),'w')
        for combo in Files[-100:]:
            f.write("%.4f %s %s %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n" % (combo[0],combo[1],combo[2],combo[3],combo[4],combo[5],combo[6],combo[7],combo[8],combo[9],combo[10],combo[11]))
        f.close()

        Error_Colors = [float(err+2.0)/4.0 for err in Relative_Accuracies]

        f = open("Relative_Accuracies_Average_%s_%d.txt" % (Polymer_Name,q),'w')
        f.write("Relative Accuracy:\n")
        f.write("%.2f" % (sum(Relative_Accuracies)/len(Relative_Accuracies)))
        f.write("Std_Dev:\n%.2f" % (np.std(np.array(Relative_Accuracies,dtype=np.float64))))
        f.close()

        #Color_Energy_Diff = mpl.colors.Normalize(vmin=-1.5,vmax=1.5)
        for norm in Energy_Diff:
            if norm > 1.0:
                norm = 1.0
            elif norm < 0.0:
                norm = 0.0
        cmap = mpl.cm.get_cmap("seismic")
        colors_array = [cmap(c) for c in Energy_Diff]

        for norm in Error_Colors:
            if norm > 1.0:
                norm = 1.0
            elif norm < 0.0:
                norm = 0.0
        colors_array_error = [cmap(c) for c in Error_Colors]
        #scalar_map = matplotlib.cm.ScalarMappable(cmap=cmap,norm=Color_Energy_Diff)
        #colors_array = scalar_map.to_rgba(Color_Energy_Diff)

        fig,ax = plt.subplots(1,1)
        plt.scatter(np.array(Colvar_DIH_List) * 180 / math.pi,np.array(Colvar_OOP_List) * 180 / math.pi,color=colors_array,cmap=cmap,marker='s',edgecolors='k')
        plt.ylabel("Improper Angle",size = 24)
        plt.xlabel("Dihedral Angle",size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        plt.tight_layout()
        #cbar = fig.colorbar(c,ax=ax)
        #cbar.ax.tick_params(labelsize = 20)
        fig.savefig(Sub_Directory + "%s_Dimer_Check_Rigid_800_%d" % (Polymer_Name,q))


        fig,ax = plt.subplots(1,1)
        plt.scatter(Dih_Diffs,OOP_Diffs,color=colors_array_error,cmap=cmap,marker='s',s=3.0)
        plt.ylabel("Improper Angle Difference",size = 24)
        plt.xlabel("Dihedral Angle Difference",size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        plt.tight_layout()
        #cbar = fig.colorbar(c,ax=ax)
        #cbar.ax.tick_params(labelsize = 20)
        fig.savefig(Sub_Directory + "%s_Dimer_Check_Rigid_800_Errors_%d" % (Polymer_Name,q))

        Make_Surface_Plot(np.linspace(0.0,180.0+180.0/len(Histogram_Error_Lists),len(Histogram_Error_Lists)+1),np.linspace(0.0,30.0+30.0/len(Histogram_Error_Lists[0]),len(Histogram_Error_Lists[0])+1),Histogram_Errors,'%s_Error_Surface_%d' % (Polymer_Name,q),Title='Errors',xlabel='Dihedral Difference',ylabel='OOP Difference (degrees)',vmin=-6,vmax=6)
        Make_Surface_Plot(np.linspace(0.0,180.0+180.0/len(Histogram_Nonbonded_Fit_Error_Lists),len(Histogram_Nonbonded_Fit_Error_Lists)+1),np.linspace(0.0,30.0+30.0/len(Histogram_Nonbonded_Fit_Error_Lists[0]),len(Histogram_Nonbonded_Fit_Error_Lists[0])+1),Histogram_Nonbonded_Fit_Errors,'%s_Nonbonded_Fit_Error_Surface_%d' % (Polymer_Name,q),Title='NF Errors',xlabel='Dihedral Difference',ylabel='OOP Difference (degrees)',vmin=-6,vmax=6)
        Make_Surface_Plot(np.linspace(0.0,180.0+180.0/len(Histogram_Error_Lists),len(Histogram_Error_Lists)+1),np.linspace(0.0,30.0+30.0/len(Histogram_Error_Lists[0]),len(Histogram_Error_Lists[0])+1),Histogram_Frequency,'%s_Frequency_%d' % (Polymer_Name,q),Title='Frequencies',xlabel='Dihedral Difference',ylabel='OOP Difference (degrees)',vmin=0.0,vmax=0.15)


    plt.close(fig)

def Control_Fit_OPLS(Dimer_Energies,Rotated_Shape):
    OPLS_Fit = []
    for j,energies in enumerate(Dimer_Energies):
        OPLS_Fit_Params,cov = scipy.optimize.curve_fit(OPLS_LAMMPS,np.linspace(0,350,Rotated_Shape[1]),energies[0])
        translated_x = np.linspace(0,350,Rotated_Shape[1]) - 180
        translated_x = [x+360 if x < 0 else x for x in translated_x]
        OPLS_Fit_Params_Translated,cov = scipy.optimize.curve_fit(OPLS_LAMMPS,translated_x,energies[0])
        fig,ax = plt.subplots(1,1)
        #print(sum([(OPLS_LAMMPS(dih,OPLS_Fit_Params[0],OPLS_Fit_Params[1],OPLS_Fit_Params[2],OPLS_Fit_Params[3])-e)**2 for dih,e in zip(np.linspace(0,350,Rotated_Shape[1]),energies[0])]))
        #print(sum([(OPLS_LAMMPS(dih,OPLS_Fit_Params_Translated[0],OPLS_Fit_Params_Translated[1],OPLS_Fit_Params_Translated[2],OPLS_Fit_Params_Translated[3])-e)**2 for dih,e in zip(translated_x,energies[0])]))
        if sum([(OPLS_LAMMPS(dih,OPLS_Fit_Params[0],OPLS_Fit_Params[1],OPLS_Fit_Params[2],OPLS_Fit_Params[3])-e)**2 for dih,e in zip(np.linspace(0,350,Rotated_Shape[1]),energies[0])]) < sum([(OPLS_LAMMPS(dih,OPLS_Fit_Params_Translated[0],OPLS_Fit_Params_Translated[1],OPLS_Fit_Params_Translated[2],OPLS_Fit_Params_Translated[3])-e)**2 for dih,e in zip(translated_x,energies[0])]):
            OPLS_Fit.append(OPLS_Fit_Params)
            plt.scatter(np.linspace(0,350,Rotated_Shape[1]), energies[0], color='k')
            q = np.linspace(0,360,1000)
            y_values = [OPLS_LAMMPS(dih,OPLS_Fit_Params[0],OPLS_Fit_Params[1],OPLS_Fit_Params[2],OPLS_Fit_Params[3]) for dih in q]
        else:
            OPLS_Fit.append(OPLS_Fit_Params_Translated)
            plt.scatter(translated_x, energies[0], color='k')
            q = np.concatenate((np.linspace(180,360,500),np.linspace(0,180,500)))
            y_values = [OPLS_LAMMPS(dih,OPLS_Fit_Params_Translated[0],OPLS_Fit_Params_Translated[1],OPLS_Fit_Params_Translated[2],OPLS_Fit_Params_Translated[3]) for dih in q]
            #print(y_values)
        plt.plot(q, y_values, color='black')
        plt.xlabel('Dihedral Rotation (Degree)',size = 24)
        plt.ylabel('OOP Rotation (Degree)',size = 24)
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        plt.tight_layout()
        fig.savefig('Conventional_OPLS_Fit_%d' % j)
        plt.close()

    return OPLS_Fit

def Calculate_Conventional_Energies(Ring_List,Rotated_Shape,Max_Dih,Max_OOP,Polymer_Name,OPLS_Fit):
    #Gives back first
    OOP_Length = 10
    Dih_Length = 36
    Phi_Rotations = np.linspace(0,30,OOP_Length)
    Theta_Rotations = np.linspace(0,360,Dih_Length)

    #OPLS_Dih_List_list is a list of single dihedral angles considered "definitive" for each set of polymers. Note that they are offset by the additional hydrogens added by the program!
    if Polymer_Name == "P3HT_Input": 
        OPLS_Dih_List_list = [[3,5,12,13],[3,5,12,13]] #P3HT
    elif Polymer_Name == "PTB7_Input":
        OPLS_Dih_List_list = [[21,23,26,27],[2,3,18,19],[2,3,18,19],[21,23,26,27]] #PTB7
    elif Polymer_Name == "N2200_Input":
        OPLS_Dih_List_list = [[10,11,32,33],[3,5,11,13],[3,5,9,10],[3,5,9,10],[3,5,11,13],[10,11,32,33]] #N2200
    elif Polymer_Name == "PNDI_T_Input":
        OPLS_Dih_List_list = [[9,11,32,33],[9,11,32,33]]

    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    Full_Nontorsional_Energy_List = []

    Ring1_List = []
    Ring2_List = []

    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Ring1_List.append(ring1)
        Ring1_List.append(ring2)
        Ring2_List.append(ring2)
        Ring2_List.append(ring1)

    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))
    Run_List = []
    k = 0
    OPLS_Dih_Angle_Matrices = []
    OPLS_Energy_Matrices = []

    for ring1,ring2,OPLS_Dih_List,OPLS_Params in zip(Ring1_List,Ring2_List,OPLS_Dih_List_list,OPLS_Fit):
        OPLS_Dih_Angle_Matrix = []
        OPLS_Energy_Matrix = []
        Dimer = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        for i in range(OOP_Length):
            OPLS_Angle_List = []
            OPLS_Energy_List = []
            for j in range(Dih_Length):
                vec_1 = Dimer.Get_Atom(OPLS_Dih_List[0]).Position - Dimer.Get_Atom(OPLS_Dih_List[1]).Position
                vec_1 = vec_1/np.linalg.norm(vec_1)
                vec_2 = Dimer.Get_Atom(OPLS_Dih_List[1]).Position - Dimer.Get_Atom(OPLS_Dih_List[2]).Position
                vec_2 = vec_2/np.linalg.norm(vec_2)
                vec_3 = Dimer.Get_Atom(OPLS_Dih_List[2]).Position - Dimer.Get_Atom(OPLS_Dih_List[3]).Position
                vec_3 = vec_3/np.linalg.norm(vec_3)
                plane_1 = np.cross(vec_1,vec_2)/np.linalg.norm(np.cross(vec_1,vec_2))
                plane_2 = np.cross(vec_2,vec_3)/np.linalg.norm(np.cross(vec_2,vec_3))
                control_dih = np.arccos(np.dot(plane_1,plane_2)) * 180 / math.pi
                OPLS_Energy_List.append(OPLS_LAMMPS(control_dih,OPLS_Params[0],OPLS_Params[1],OPLS_Params[2],OPLS_Params[3]))
                OPLS_Angle_List.append(control_dih)
                Dimer.Rotate_Ring("Dih",360/Dih_Length,Dimer.Ring_List[0],Dimer.Ring_List[1])
                if i == 0:
                    if j == 0:
                        Planar_Energy = OPLS_Energy_List[-1]
                    if j == 10:
                        No_OOP_Twisted_Energy = OPLS_Energy_List[-1]
                if i == (OOP_Length -1):
                    if j == 0:
                        OOP_No_Dih_Energy = OPLS_Energy_List[-1]
                    if j == 10:
                        OOP_Twisted_Energy = OPLS_Energy_List[-1]
            print(Dimer.Name)
            print([a.Position for a in Dimer.Atom_List])
            OPLS_Energy_Matrix.append(OPLS_Energy_List)
            Dimer.Rotate_Ring("OOP",30/OOP_Length,Dimer.Ring_List[0],Dimer.Ring_List[1])
            OPLS_Dih_Angle_Matrix.append(OPLS_Angle_List)
        OPLS_Energy_Matrices.append(OPLS_Energy_Matrix)
        OPLS_Dih_Angle_Matrices.append(OPLS_Dih_Angle_Matrix)

    print(OPLS_Fit)
    print("0 OOP:")
    print(No_OOP_Twisted_Energy - Planar_Energy)
    print("30 OOP:")
    print(OOP_Twisted_Energy - OOP_No_Dih_Energy)
    print("OOP Diff:")
    print(OOP_No_Dih_Energy - Planar_Energy)
    for OPLS_Energy_Matrix,ring1,ring2,i in zip(OPLS_Energy_Matrices,Ring_List,Offset_Ring_List,range(len(OPLS_Energy_Matrices))):
        Make_Surface_Plot(Theta_Rotations,Phi_Rotations,OPLS_Energy_Matrix,'%s_%s_Conventional_Energies' % (ring1.Name,ring2.Name),Title='Dihedral-torsion',xlabel='Dihedral Angle',ylabel='OOP Angle (radians)')

def Calculate_Statistics(Input_File_Sidechains,XYZ_File_Sidechains,Polymer_Name,Deg_Polym,PL_Atom_ID):
    Ring_List,Parameterize_Bond,Parameterize_Charges = Read_Input(Input_File_Sidechains,XYZ_File_Sidechains,Polymer_Name)
    Sub_Directory = "./Full_Polymer_Plumed/Double_Deloc_Double_Nonbonded/"
    File_Name_Fragment_1_List = ["P3HT_Polymer_Double_Deloc_Double_Nonbonded_300_"]
    File_Name_Fragment_2_List = ["_Production.data"]
    LAMMPS_Start_Timestep = 1000100
    LAMMPS_End_Timestep = 10100100
    LAMMPS_Step = 100000
    Offset_Ring_List = []
    Rg_List = []

    Full_Polymer_Ring_List = []

    for i in range(Deg_Polym):
        Full_Polymer_Ring_List = Full_Polymer_Ring_List + copy.deepcopy(Ring_List)

    Polymer = Conjugated_Polymer.Conjugated_Polymer(Full_Polymer_Ring_List)

    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))

    for File_Name_Fragment_1,File_Name_Fragment_2,ring1,ring2 in zip(File_Name_Fragment_1_List,File_Name_Fragment_2_List,Ring_List,Offset_Ring_List):
        File_Name_List = []
        Colvar_OOP_List = []
        Colvar_DIH_List = []


        LAMMPS_Current_Timestep = LAMMPS_Start_Timestep
        while LAMMPS_Current_Timestep <= LAMMPS_End_Timestep:
            File_Name = (File_Name_Fragment_1 + "%d" + File_Name_Fragment_2) % (LAMMPS_Current_Timestep)
            File_Name_List.append(File_Name)
            LAMMPS_Current_Timestep += LAMMPS_Step

        Rgs = []
        PLs = []
        for File_Name in File_Name_List:
            Polymer.Read_From_Data_File(Sub_Directory + File_Name)
            Rgs.append(Polymer.Calculate_Rg())
            PLs.append(Polymer.Calculate_PL(0,Polymer_Name,Make_Figure=True)[0])
    f = open("%s_Polymer_Stats.txt" % Polymer_Name,'w')
    print(np.array(Rgs))
    print(np.array(PLs))
    f.write("Radius of gyration: %.2f" % np.mean(np.array(Rgs)))
    f.close()

    f = open("Check.txt",'w')
    f.write("Hi")
    f.close()

def Make_Solvated_Box(Polymer_Ring_List,Base_Name,Solvent_Name,Folder_Name = "Plumed_Scripts"):
    #Box_Size in Angstroms
    Polymer = Conjugated_Polymer.Conjugated_Polymer(Polymer_Ring_List)
    Polymer.Read_From_Data_File("P3HT_Polymer_300.0_4000000.data")
    Solvent = Molecule.Molecule("%s.xyz" % Solvent_Name)

    Box_Size = np.linalg.norm(Polymer.Get_Atom(1).Position - Polymer.Get_Atom(len(Polymer.Atom_List)).Position) * 4/3
    print(Box_Size)
    Box_Size = 60
    Density = 1.11 # g/mL
    MW = 112.56
    Box_Volume = (Box_Size**3) * (10**-24) # in mL
    Num_Solv = math.floor(Box_Volume * Density * (6.022*(10**23)))/MW
    print("Number of Solvents:")
    print(Num_Solv)
    Solvent.Set_Up_FF()
    op.Assign_OPLS(Solvent, ChelpG = False)
    Polymer.Zero_Total_Charge()
    Polymer.Name = "P3HT"
    Solvent_System = System.System([Polymer,Solvent], [1,int(Num_Solv)], Box_Size, "%s_Solvated" % (Base_Name))
    Solvent_System.Run_Packmol()
    Solvent_System.Write_LAMMPS_Data(Data_File = "%s_Solvated.data" % (Base_Name))

def Find_Charges(Ring_List,Polymer_Name):
    Extended_Ring_List = copy.deepcopy(Ring_List)
    for i in range(2):
        for ring in Ring_List:
            Extended_Ring_List.append(copy.deepcopy(ring))
    Multimer = Conjugated_Polymer.Conjugated_Polymer(Extended_Ring_List)
    #Multimer.Write_XYZ()
    Job_Name = "%s_Find_Charges" % Polymer_Name
    In_File = "%s_Find_Charges.qcin" % Polymer_Name
    Sub_File = "sub_%s_Find_Charges" % Polymer_Name
    Job_Type = "QChem"
    Folder_Name = "Find_Charges"
    End_File = "%s_Find_Charges.out" % Polymer_Name
    Cluster_Login = Configure.orca_dict["Cluster_Login"]
    Base_Cluster_Location = Configure.orca_dict["Base_Cluster_Location"]
    Cluster_Location=Base_Cluster_Location + "/Find_Charges"
    Scheduler_Type = Configure.orca_dict["Scheduler_Type"]
    End_Condition = Configure.orca_dict["End_Condition"]
    Copy_File_List = [In_File,Sub_File]
    #Write_Inputs.Write_QChem_SPE(In_File,Multimer,Exchange_Method = "B88",Correlation_Method = "P86",Basis = "def2-SVP",Implicit_Solvent_Dielectric = 0.0,ChelpG = True)
    Write_Inputs.Write_QChem_SPE(In_File,Multimer,Implicit_Solvent_Dielectric = 0.0,ChelpG = True)
    Write_Submit_Script.Write_SLURM(Sub_File,In_File,Job_Name,32,Cluster_Location,Job_Type,walltime = 4,queue = "regular",proc_per_node = 32,constraint = 'haswell')
    Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File)
    Charges = Cluster_IO.Return_Info(End_File,End_File,Folder_Name,Job_Type,Cluster_Login,Cluster_Location,End_Condition = End_Condition,Return_ChelpG_Charges_QChem = True)
    Charges = Charges[0]
    Ring_List_Atoms = []
    for ring in Ring_List:
        Ring_List_Atoms = Ring_List_Atoms + ring.Atom_List
    if len(Charges[len(Ring_List_Atoms):2*len(Ring_List_Atoms)]) != len(Ring_List_Atoms):
        print(len(Charges))
        print(len(Charges[len(Ring_List_Atoms):2*len(Ring_List_Atoms)]))
        print(len(Ring_List_Atoms))
        raise Exception("Mismatched Charge List Size")
    for q,atom in zip(Charges[len(Ring_List_Atoms)+1:2*len(Ring_List_Atoms)+1],Ring_List_Atoms):
        atom.Charge = q
    for b_atom in Ring_List[0].Bonded_Atoms:
        b_atom.H_Atom.Charge = Charges[len(Ring_List[0].Atom_List)]
    for b_atom in Ring_List[-1].Bonded_Atoms:
        b_atom.H_Atom.Charge = Charges[-1]

    Multimer.Write_XYZ()

def Strech_Bond(Ring_List,Polymer_Name):
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(ring)
    Offset_Ring_List.append(Ring_List[0])
    Job_Type = Configure.qchem_dict["Job_Type"]
    Folder_Name = Configure.qchem_dict["Folder_Name"]
    Cluster_Login = Configure.qchem_dict["Cluster_Login"]
    Base_Cluster_Location = Configure.qchem_dict["Base_Cluster_Location"]
    Cluster_Location= Configure.qchem_dict["Cluster_Location"]
    Scheduler_Type = Configure.qchem_dict["Scheduler_Type"]
    End_Condition = Configure.qchem_dict["End_Condition"]
    Executable_Location = Configure.qchem_dict["Executable_Location"]
    OpenMP_Location = Configure.qchem_dict["OpenMP_Location"]
    Run_List = []

    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Stretch = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        Stretch.Write_XYZ()
        Bonded_Atom_List = []
        for bond in Stretch.Interring_Bond_List:
            Bonded_Atom_List.append(bond.Bond_Main)
            Bonded_Atom_List.append(bond.Bond_Node)
        Bond_Eq = np.linalg.norm(Bonded_Atom_List[0].Position - Bonded_Atom_List[1].Position)
        Job_Name = "%s_%s_Interring_Bond" % (ring1.Name,ring2.Name)
        In_File = "%s_%s_Interring_Bond.inp" % (ring1.Name,ring2.Name)
        Sub_File = "sub_%s_%s_Interring_Bond" % (ring1.Name,ring2.Name)
        End_File = "%s_%s_Interring_Bond.out" % (ring1.Name,ring2.Name)
        for atom in Stretch.Atom_List:
            print("Atom ID:")
            print(atom.Atom_ID)
        """if ring1.Symmetric and ring2.Symmetric:
            Symmetry_End_File = "%s_%s_Interring_Bond.out" % (ring2.Name,ring1.Name)
        else:
            Symmetry_End_File = "" """
        Symmetry_End_File = ""
        Copy_File_List = [In_File,Sub_File]
        Write_Inputs.Write_Orca_Optimize_Geometry(In_File,Stretch,Bond_Eq = Bond_Eq,Bond_Atoms = [Bonded_Atom_List[0].Atom_ID,Bonded_Atom_List[1].Atom_ID])
        Write_Submit_Script.Write_TORQUE(Sub_File,In_File,Job_Name,4,Cluster_Location,Job_Type,walltime = 8,Executable_Path = Executable_Location,OMP_Path = OpenMP_Location)
        if End_File not in Run_List:
            Cluster_IO.Submit_Job(Copy_File_List,Folder_Name,Sub_File,End_File,Job_Name,Cluster_Login,Cluster_Location,Base_Cluster_Location,Scheduler_Type,End_Condition = End_Condition,Analyze_File = End_File,Symmetry_End_File = Symmetry_End_File,Symmetry_Analyze_File = Symmetry_End_File)
            Run_List.append(End_File)
            if Symmetry_End_File != "":
                Run_List.append(Symmetry_End_File)

    for ring1,ring2 in zip(Ring_List,Offset_Ring_List):
        Job_Name = "%s_%s_Interring_Bond" % (ring1.Name,ring2.Name)
        In_File = "%s_%s_Interring_Bond.inp" % (ring1.Name,ring2.Name)
        Sub_File = "sub_%s_%s_Interring_Bond" % (ring1.Name,ring2.Name)
        End_File = "%s_%s_Interring_Bond.out" % (ring1.Name,ring2.Name)
        """if ring1.Symmetric and ring2.Symmetric:
            Symmetry_End_File = "%s_%s_Interring_Bond.out" % (ring2.Name,ring1.Name)
        else:
            Symmetry_End_File = "" """
        Symmetry_End_File = ""
        Stretch = Conjugated_Polymer.Conjugated_Polymer([ring1,ring2])
        End_File_List = [End_File]
        for i in range(1,13):
            if i < 10:
                End_File_List.append(Job_Name + ".00%d.xyz" % i)
            else:
                End_File_List.append(Job_Name + ".0%d.xyz" % i)
        Cluster_IO.Return_Info_Batch(End_File_List,End_File_List,Folder_Name,Job_Name,Cluster_Login,Cluster_Location,Symmetry_Analyze_File = Symmetry_End_File,Symmetry_End_File = Symmetry_End_File)

        Molecule.Assign_Lammps([Stretch])
        Bond_Lengths = []
        Raw_Energies = []
        f = open("./%s/%s" % (Folder_Name,End_File),'r')
        print("./%s/%s" % (Folder_Name,End_File))
        lines = f.readlines()
        Energy_Flag = False
        for line in lines:
            if Energy_Flag:
                if len(line.split()) >= 2:
                    Bond_Lengths.append(float(line.strip().split()[0]))
                    Raw_Energies.append(float(line.strip().split()[1]) * 627.5)
                else:
                    Energy_Flag = False
            if len(line.split()) > 2 and line.split()[-1] == "Energy'" and line.split()[-2] == "'Actual":
                Energy_Flag = True

        Bond_Lengths = np.asarray(Bond_Lengths)
        Raw_Energies = np.asarray(Raw_Energies)
        Raw_Energies = Raw_Energies - np.amin(Raw_Energies)
        fig = plt.figure()
        plt.scatter(Bond_Lengths,Raw_Energies)
        plt.xlabel('Bond Lengths (A)')
        plt.ylabel('Raw Energies (kcal/mol)')
        fig.savefig('%s_%s_Raw_Bond_Energies' % (ring1.Name,ring2.Name))
        plt.close(fig)


        Correction_Energy = []
        for file in End_File_List[-12:]:
            Stretch.Update_From_XYZ("./%s/%s" % (Folder_Name,file))
            Correction_Energy.append(Stretch.Calculate_Internal_Energy(Polymer_Name,"~/lammps-11Aug17/src/lmp_serial",Exclude_Interring_Bonds=True))
        Correction_Energy = np.asarray(Correction_Energy)
        Final_Energies = Raw_Energies - Correction_Energy
        Final_Energies = Final_Energies - np.amin(Final_Energies)
        fig = plt.figure()
        plt.scatter(Bond_Lengths,Final_Energies)
        plt.xlabel('Bond Lengths (A)')
        plt.ylabel('Final Energies (kcal/mol)')
        fig.savefig('%s_%s_Final_Bond_Energies' % (ring1.Name,ring2.Name))
        plt.close(fig)

def Make_Example_Files(Ring_List):
    Planar = Conjugated_Polymer.Conjugated_Polymer(Ring_List)
    Planar.Write_XYZ()
    Planar.Rotate_Ring("Dih",10,Planar.Ring_List[0],Planar.Ring_List[1])
    Planar.Write_XYZ()
    Planar1 = Conjugated_Polymer.Conjugated_Polymer(Ring_List)
    Planar2 = Conjugated_Polymer.Conjugated_Polymer(Ring_List)
    Planar2.Translate_Polymer([4.0,0.0,0.0])
    f = open("Ordered_Pi_Stack.xyz",'w')
    f.write("%d\n\n" % (len(Planar1.Atom_List)*2))
    for atom in Planar1.Atom_List:
        f.write("%s %f %f %f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
    for atom in Planar2.Atom_List:
        f.write("%s %f %f %f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
    f.close()
    Planar1.Rotate_Ring("OOP",30,Planar1.Ring_List[0],Planar1.Ring_List[1])
    Planar1.Rotate_Ring("Dih",50,Planar1.Ring_List[0],Planar1.Ring_List[1])
    #Planar1.Rotate_Ring("Dih",330,Planar1.Ring_List[1],Planar1.Ring_List[0])
    f = open("Disordered_Pi_Stack.xyz",'w')
    f.write("%d\n\n" % (len(Planar1.Atom_List)*2))
    for atom in Planar1.Atom_List:
        f.write("%s %f %f %f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
    for atom in Planar2.Atom_List:
        f.write("%s %f %f %f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
    f.close()

def Set_Up_Folders():
    os.system("mkdir ./Bond_Parameters")
    os.system("mkdir ./Figures")
    os.system("mkdir ./Rotation_Run_Input_Copies")
    os.system("mkdir ./XYZ_Files")
    os.system("mkdir ./Hydrogenated_XYZ_Files")
    os.system("mkdir ./Alternate_Hydrogenated_XYZ_Files")
    os.system("mkdir ./Nontorsional_Outputs")
    os.system("mkdir ./Nontorsional_Inputs")
    os.system("mkdir ./Hydrogenated_Improper_XYZ_Files")
    os.system("mkdir ./Alternate_Hydrogenated_Improper_XYZ_Files")
    os.system("mkdir ./Multi_Ring_Hydrogenated_Rotation_Test")

def Write_All_Dimers_Plumed(Ring_List,Fit_Energies_list,Force_y_list,Force_x_list,Polymer_Name):
    
    Offset_Ring_List = []
    for ring in Ring_List[1:]:
        Offset_Ring_List.append(copy.deepcopy(ring))
    Offset_Ring_List.append(copy.deepcopy(Ring_List[0]))

    for ring1,ring2,Fit_Energies,Force_y,Force_x,i in zip(Ring_List,Offset_Ring_List,Fit_Energies_list,Force_y_list,Force_x_list,range(len(Ring_List))):
        Write_Example_Plumed_Script([ring1,ring2],[Fit_Energies],[Force_y],[Force_x],"%s_Dimers_%d" % (Polymer_Name,i),Folder_Name="Plumed_Scripts")

def Write_Example_Plumed_Script(Polymer_Ring_List,Fit_Energies,Force_y,Force_x,Base_Name,Folder_Name="Plumed_Scripts",Non_Interacting=False,Nonbonded_Fit_Energies=[],Nonbonded_Force_y=[],Nonbonded_Force_x=[],Delocalization_Modifier=1,Nonbonded_Modifier=1):
    oop_blocks = np.linspace(0,2*math.pi,200)
    dih_blocks = np.linspace(-math.pi,math.pi,200)

    #TODO: temp, remove later
    import imp
    imp.reload(Conjugated_Polymer)
    imp.reload(Molecule)
    Polymer = Conjugated_Polymer.Conjugated_Polymer(Polymer_Ring_List)

    Molecule.Assign_Lammps([Polymer]) #where's the output? -Leon
    Bond_Atoms = []
    for bond in Polymer.Interring_Bond_List:
        if bond.Bond_Main not in Bond_Atoms:
            Bond_Atoms.append(bond.Bond_Main)
        if bond.Bond_Node not in Bond_Atoms:
            Bond_Atoms.append(bond.Bond_Node)
    if Non_Interacting == True:
        for atom in Polymer.Atom_List:
            atom.Charge = 0.000000
            atom.Sigma = 0.0
            atom.Epsilon = 0.0

    Polymer.Adjust_COM(center=np.array([150.0,150.0,150.0]))
    Polymer.Write_Data_File("Run_%s.data" % (Base_Name),Exclude_Interring_Dihedrals=True,Bond_Atoms=Bond_Atoms)
    Base_Calculate_Torsions_String = "%s_Calculate_Torsions" % (Base_Name)
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
    torsion_file = open(Base_Calculate_Torsions_String + ".dat", 'w')


    torsion_file.write("\n\nc1: CENTER ATOMS=%d" % Polymer.Ring_List[0].Core_Atom_List[0].Atom_ID)
    for atom in Polymer.Ring_List[0].Core_Atom_List[1:]:
        torsion_file.write(",%d" % atom.Atom_ID)
    torsion_file.write("\n\nc1_normal: GHOST ATOMS=c1,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (Polymer.Ring_List[0].Core_Atom_List[0].Atom_ID,Polymer.Ring_List[0].Core_Atom_List[1].Atom_ID))

    Parameter_Index = 0
    for i in range(1,len(Polymer.Ring_List)):
        Bias_File_Name = "%s_Bias_File_%d" % (Base_Name,i)
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
        f = open(Bias_File_Name+".dat", 'w')
        f.write("#! FIELDS DIH_%d OOP_%d ext_%d.bias der_DIH_%d der_OOP_%d\n#! SET min_DIH_%d -pi\n#! SET max_DIH_%d pi\n#! SET nbins_DIH_%d 199\n#! SET periodic_DIH_%d true\n#! SET min_OOP_%d 0\n#! SET max_OOP_%d 2*pi\n#! SET nbins_OOP_%d 199\n#! SET periodic_OOP_%d false\n" % (i,i,i,i,i,i,i,i,i,i,i,i,i))
        for q,energy_list in enumerate(Fit_Energies[Parameter_Index][:-1]):
            for j,e in enumerate(energy_list[:-1]):
                f.write("%.6f %.6f %.6f %.6f %.6f\n" % (dih_blocks[j],oop_blocks[q],e*Delocalization_Modifier,Force_x[Parameter_Index][q][j]*Delocalization_Modifier,Force_y[Parameter_Index][q][j]*Delocalization_Modifier))
            f.write("\n")
        f.close()

        if Non_Interacting:
            f = open(Bias_File_Name + "_Nonbonded.dat", 'w')
            f.write("#! FIELDS DIH_%d OOP_%d ext_%d_nb.bias der_DIH_%d der_OOP_%d\n#! SET min_DIH_%d -pi\n#! SET max_DIH_%d pi\n#! SET nbins_DIH_%d 199\n#! SET periodic_DIH_%d true\n#! SET min_OOP_%d 0\n#! SET max_OOP_%d 2*pi\n#! SET nbins_OOP_%d 199\n#! SET periodic_OOP_%d false\n" % (i,i,i,i,i,i,i,i,i,i,i,i,i))
            for q,energy_list in enumerate(Nonbonded_Fit_Energies[Parameter_Index][:-1]):
                for j,e in enumerate(energy_list[:-1]):
                    f.write("%.6f %.6f %.6f %.6f %.6f\n" % (dih_blocks[j],oop_blocks[q],e*Nonbonded_Modifier,Nonbonded_Force_x[Parameter_Index][q][j]*Nonbonded_Modifier,Nonbonded_Force_y[Parameter_Index][q][j]*Nonbonded_Modifier))
                f.write("\n")
            f.close()
        os.system("scp %s.dat ./%s/" % (Bias_File_Name,Folder_Name))
        os.remove("%s.dat" % (Bias_File_Name))
        if Non_Interacting:
            os.system("scp %s_Nonbonded.dat ./%s/" % (Bias_File_Name,Folder_Name))
            os.remove("%s_Nonbonded.dat" % (Bias_File_Name))

        torsion_file.write("\n\nc%d: CENTER ATOMS=%d" % (i+1,Polymer.Ring_List[i].Core_Atom_List[0].Atom_ID))
        for atom in Polymer.Ring_List[i].Core_Atom_List[1:]:
            torsion_file.write(",%d" % atom.Atom_ID)
        for b_atom in Polymer.Ring_List[i].Bonded_Atoms:
            if b_atom.Is_Linked and b_atom.Interring_Bond_Atom.Self_Ring.Ring_ID == i:
                atom1 = b_atom.Interring_Bond_Atom.Central_Atom.Atom_ID
                atom2 = b_atom.Central_Atom.Atom_ID
                ring1_batom1 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[0].Atom_ID
                ring1_batom2 = b_atom.Interring_Bond_Atom.Same_Ring_Bonded_Atom_List[1].Atom_ID
                ring2_batom1 = b_atom.Same_Ring_Bonded_Atom_List[0].Atom_ID
                ring2_batom2 = b_atom.Same_Ring_Bonded_Atom_List[1].Atom_ID
        torsion_file.write("\n\nc%d_normal: GHOST ATOMS=c%d,%d,%d COORDINATES=0.0,1.0,0.0\n\n" % (i+1,i+1,Polymer.Ring_List[i].Core_Atom_List[0].Atom_ID,Polymer.Ring_List[i].Core_Atom_List[1].Atom_ID))
        torsion_file.write("DIH_%d: TORSION VECTOR1=c%d,c%d_normal AXIS=%d,%d VECTOR2=c%d,c%d_normal\n\n" % (i,i,i,atom1,atom2,i+1,i+1))
        torsion_file.write("OOP_%d_1: TORSION ATOMS=%d,%d,%d,%d\n\n" % (i,atom1,ring1_batom1,ring1_batom2,atom2))
        torsion_file.write("OOP_%d_2: TORSION ATOMS=%d,%d,%d,%d\n\n" % (i,atom2,ring2_batom1,ring2_batom2,atom1))
        torsion_file.write("OOP_%d: MATHEVAL ARG=OOP_%d_1,OOP_%d_2 FUNC=abs(x)+abs(y) PERIODIC=NO\n\n" % (i,i,i))
        torsion_file.write("PRINT ARG=DIH_%d,OOP_%d FILE=colvar_%d.txt STRIDE=100\n\n" % (i,i,i))
        torsion_file.write("EXTERNAL ARG=DIH_%d,OOP_%d FILE=%s.dat LABEL=ext_%d\n\n" % (i,i,Bias_File_Name,i))
        if Non_Interacting:
            torsion_file.write("EXTERNAL ARG=DIH_%d,OOP_%d FILE=%s_Nonbonded.dat LABEL=ext_%d_nb\n\n" % (i,i,Bias_File_Name,i))
        Parameter_Index += 1
        if Parameter_Index == len(Fit_Energies):
            Parameter_Index = 0

    torsion_file.write("WHOLEMOLECULES ENTITY0=1-%d" % len(Polymer.Atom_List))
    for i in range(1,len(Polymer.Ring_List)):
        torsion_file.write(",c%d,c%d_normal" % (i,i))

    torsion_file.close()
    if not os.path.isdir(Folder_Name):
        os.mkdir(Folder_Name)

    os.system("scp %s.dat ./%s/" % (Base_Calculate_Torsions_String,Folder_Name))
    os.system("scp Run_%s.data ./%s" % (Base_Name,Folder_Name))
    os.remove("%s.dat" % Base_Calculate_Torsions_String)
    os.remove("Run_%s.data" % Base_Name)

    for i in range(1,len(Polymer.Ring_List)):
        if Non_Interacting:
            os.system("scp %s_Bias_File_%d_Nonbonded.dat ./%s/" % (Base_Name,i,Folder_Name))
            os.system("rm -f %s_Bias_File_%d_Nonbonded.dat" % (Base_Name,i))

def Write_Conventional_Data_File(Polymer_Name,Input_File_Sidechains,XYZ_File_Sidechains,Deg_Polym,OPLS_Fit,Folder_Name="Conventional_Scripts",Non_Interacting=False):

    Ring_List,Parameterize_Bond,Parameterize_Charges = Read_Input(Input_File_Sidechains,XYZ_File_Sidechains,Polymer_Name)

    Full_Polymer_Ring_List = []

    for i in range(Deg_Polym):
        Full_Polymer_Ring_List = Full_Polymer_Ring_List + copy.deepcopy(Ring_List)

    Polymer = Conjugated_Polymer.Conjugated_Polymer(Full_Polymer_Ring_List)
    Molecule.Assign_Lammps([Polymer])

    if Non_Interacting == True:
        for atom in Polymer.Atom_List:
            atom.Charge = 0.000000
            atom.Sigma = 0.0
            atom.Epsilon = 0.0

    print(Polymer.Interring_Bond_List)
    Bond_Atoms = []
    for bond in Polymer.Interring_Bond_List:
        if bond.Bond_Main not in Bond_Atoms:
            Bond_Atoms.append(bond.Bond_Main)
        if bond.Bond_Node not in Bond_Atoms:
            Bond_Atoms.append(bond.Bond_Node)

    Polymer.Adjust_COM(center=np.array([150.0,150.0,150.0]))
    Polymer.Replace_Interring_Dihedrals([0],[0],[o for o in OPLS_Fit[0]]) # Had to change some inputs to list([0]), to fit with the function definition better, not sure which one is actually the problem - Leon
    Polymer.Write_Data_File("Run_%s_Conventional.data" % Polymer_Name,Bond_Atoms=Bond_Atoms)
    
    os.system("mkdir ./%s" % (Folder_Name))

    os.system("scp Run_%s_Conventional.data ./%s/" % (Polymer_Name,Folder_Name))
    os.remove("Run_%s_Conventional.data" % (Polymer_Name))

def Compare_Rigidity_Types(Merged_Conjugation_Energies,Merged_Delocalization_Energies,Merged_Nonbonded_Energies,Merged_OOP_Rotations_Degrees):
    for E_conjugation,E_delocalization,E_nonbonded,i in zip(Merged_Conjugation_Energies,Merged_Delocalization_Energies,Merged_Nonbonded_Energies,range(len(Merged_Conjugation_Energies))):
        fig,ax = plt.subplots(1,1)
        Dihedral_Degrees = np.linspace(0,180,19)
        Total_Dihedral_Energy = E_delocalization[0][:19] + E_nonbonded[0][:19]
        Total_Dihedral_Energy = Total_Dihedral_Energy - np.amin(Total_Dihedral_Energy)
        #plt.plot(Dihedral_Degrees,E_conjugation[0][:10] - np.amin(E_conjugation[0][:10]),color='k',marker='s',markersize=19)
        ax.plot(Dihedral_Degrees,Total_Dihedral_Energy,color='k',marker='s',markersize=10,label="Total Energy")
        ax.plot(Dihedral_Degrees,E_delocalization[0][:19] - np.amin(E_delocalization[0][:19]),color='r',marker='s',markersize=10,label="Delocalization Energy")
        ax.plot(Dihedral_Degrees,E_nonbonded[0][:19] - np.amin(E_nonbonded[0][:19]),color='b',marker='s',markersize=10,label="Nonbonded Energy")
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        plt.ylim((-1,10))
        plt.xlabel("Dihedral Angle (Degrees)",fontsize=24)
        plt.ylabel("Energy (kcal/mol)",fontsize=24)
        ax.legend()
        plt.tight_layout()
        fig.savefig("Dihedral_Energy_Types_%d" %i)
        plt.close(fig)

        fig,ax = plt.subplots(1,1)
        Total_Improper_Energy = np.array([e[0] for e in E_conjugation] - np.amin([e[0] for e in E_conjugation])) + np.array([e[0] for e in E_delocalization]) - np.amin([e[0] for e in E_delocalization]) + np.array([e[0] for e in E_nonbonded]) - np.amin([e[0] for e in E_nonbonded])
        ax.plot(Merged_OOP_Rotations_Degrees,np.array([e[0] for e in E_conjugation] - np.amin([e[0] for e in E_conjugation])),color='c',marker='s',markersize=10,label="Conjugation Energy")
        ax.plot(Merged_OOP_Rotations_Degrees,np.array([e[0] for e in E_delocalization]) - np.amin([e[0] for e in E_delocalization]),color='r',marker='s',markersize=10,label="Delocalization Energy")
        ax.plot(Merged_OOP_Rotations_Degrees,np.array([e[0] for e in E_nonbonded]) - np.amin([e[0] for e in E_nonbonded]),color='b',marker='s',markersize=10,label="Nonbonded Energy")
        ax.plot(Merged_OOP_Rotations_Degrees,Total_Improper_Energy,color='k',marker='s',markersize='10',label="Total Energy")
        ax.tick_params(axis="x", labelsize=18)
        ax.tick_params(axis="y", labelsize=18)
        ax.tick_params(length=4,width=4)
        plt.ylim((-1,10))
        plt.xlabel("Improper Angle (Degrees)",fontsize=24)
        plt.ylabel("Energy (kcal/mol)",fontsize=24)
        ax.legend()
        plt.tight_layout()
        fig.savefig("Improper_Energy_Types_%d" %i,pad_inches=2.0)
        plt.close(fig)


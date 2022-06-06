import numpy as np
import Atom
import Molecule
import Ring
import OPLS as op
import System
import Conjugated_Polymer
import Cluster_IO
import Write_Inputs
import Write_Submit_Script
import math
import copy

Ring_List,Parameterize_Bond,Parameterize_Charges = Read_Input(Input_File,XYZ_File,Polymer_Name)
#Write XYZ files for LigParGen
for i,pring in enumerate(Polymer_Rings):
    pring.Write_XYZ("%s_%d.xyz" % (pring.Name,i))
#Manually run LigParGen, then you're ready for:
for i,pring in enumerate(Polymer_Rings):
    pring.Read_In_LigParGen(i)
#Now create a list of all the rings you want to include
Polymer_Ring_List = []
for i in range(15):
    for ring in Polymer_Rings:
        Polymer_Ring_List.append(copy.deepcopy(ring))
#Connect the rings together
for i in range(14):
    if (not Polymer_Ring_List[i].Bonded_Atoms[-1].Is_Linked) and not (Polymer_Ring_List[i+1].Bonded_Atoms[0].Is_Linked):
        Polymer_Ring_List[i].Bonded_Atoms[-1].Add_Ring(Polymer_Ring_List[i+1])
        Polymer_Ring_List[i+1].Bonded_Atoms[0].Add_Ring(Polymer_Ring_List[i])
        Polymer_Ring_List[i].Bonded_Atoms[-1].Add_Bonded_Atom(Polymer_Ring_List[i+1].Bonded_Atoms[0])
        Polymer_Ring_List[i+1].Bonded_Atoms[0].Add_Bonded_Atom(Polymer_Ring_List[i+1].Bonded_Atoms[-1])
#Create a system of 30 of those polymers and print a LAMMPS data file
Polymer = Conjugated_Polymer.Conjugated_Polymer(Polymer_Ring_List)
System_List = []
Comp_List = np.ones(30,dtype=int)
for poly in range(30):
    System_List.append(copy.deepcopy(Conventional_Polymer))
Film = System.System(System_List,Comp_List,160.0,"PTB7_Film_30_30mers_Dihedral_torsion")
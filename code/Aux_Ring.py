#! usr/bin/python

# Open relevant modules

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
import Ring

class Aux_Ring(Molecule.Molecule):
   
    def __init__(self,Supermolecule,Atom_Numbers,Core_Atom_List,Self_Bond_Atom,Other_Bond_Atom,Main_Ring,Name,Polymer_Name,Nickname,Symmetric = True):
        self.Atom_List = []
        self.Core_Atom_List = []
        self.Name = Name
        self.Symmetric = Symmetric
        self.Polymer_Name = Polymer_Name
        self.Nickname = Nickname
        self.Main_Ring = Main_Ring
        self.Self_Bond_Atom = self.Main_Ring.Get_Atom(Self_Bond_Atom)
        self.Other_Bond_Atom = self.Main_Ring.Get_Atom(Other_Bond_Atom)
        for atom_num in Atom_Numbers:
            self.Atom_List.append(self.Main_Ring.Get_Atom(atom_num))
            if atom_num in Core_Atom_List:
                self.Core_Atom_List.append(self.Main_Ring.Get_Atom(atom_num))

        #Calculate Normal Vector using non-hydrogen atoms

        self.Normal_Vector = np.array([0.0,0.0,0.0])
        for cut1,atom1 in enumerate(self.Core_Atom_List):
            for cut2,atom2 in enumerate(self.Core_Atom_List[cut1+1:]):
                for atom3 in self.Core_Atom_List[cut1+cut2+2:]:
                    vec1 = atom1.Position - atom2.Position
                    vec2 = atom2.Position - atom3.Position
                    self.Normal_Vector = self.Normal_Vector + (np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2)))
        self.Normal_Vector = self.Normal_Vector/np.linalg.norm(self.Normal_Vector)

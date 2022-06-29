import Molecule
import Configure
import sys

def main():
    DA_polymer = Molecule.Molecule(Configure.local_dict["User_Path"] + "/DPP-HD.xyz")
    DA_polymer.Set_Up_FF()

if __name__=='__main__': main()

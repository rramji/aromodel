#! usr/bin/python

import os
import subprocess
import time

def Write_Orca_SPE(File_Name,Test_Molecule,Method = "RI-MP2",Basis = "cc-pVTZ",Aux_Basis = "cc-pVTZ/C",convergence = "NORMALSCF",nproc = 4,Implicit_Solvent = "CPCM(Chloroform)",Memory = 3000):
	f = open(File_Name,'w')
	if nproc == 1:
		Parallel = ""
	else:
		Parallel = "PAL%d" % nproc
	f.write("! %s %s %s %s %s %s\n%%maxcore %d\n\n*xyz 0 1\n" % (Method,Basis,Aux_Basis,convergence,Parallel,Implicit_Solvent,Memory))
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("*")
	f.close()

def Write_Orca_Optimize_Geometry(File_Name,Test_Molecule,Method = "RI BP86",Basis = "def2-SVP",Aux_Basis = "def2/J",Dispersion_Correction = "D3BJ",convergence = "TIGHTSCF",nproc = 1,Memory = 3000,H_Only = False,Planarize = False,Planarize_Atoms = [],Aux_Torsions1 = [],Aux_Torsions2 = [],Bond_Eq = 0.0,Bond_Atoms = []):
	f = open(File_Name,'w')
	if nproc == 1:
		Parallel = ""
	else:
		Parallel = "PAL%d" % nproc
	if H_Only:
		H_Opt = "%geom\noptimizehydrogens true\nend\n\n"
	else:
		H_Opt = ""
	if Planarize:
		Plane = "%geom\nConstraints\n"
		for atoms in Planarize_Atoms:
			Add_Line = "{D %d %d %d %d 0.0 C}\n" % (atoms[0]-1,atoms[1]-1,atoms[2]-1,atoms[3]-1)
			Plane = Plane + Add_Line
		for atoms in Aux_Torsions1:
			Add_Line = "{D %d %d %d %d 150.0 C}\n" % (atoms[0]-1,atoms[1]-1,atoms[2]-1,atoms[3]-1)
			Plane = Plane + Add_Line
		for atoms in Aux_Torsions2:
			Add_Line = "{D %d %d %d %d 150.0 C}\n" % (atoms[0]-1,atoms[1]-1,atoms[2]-1,atoms[3]-1)
			Plane = Plane + Add_Line
		Plane = Plane + "end\nend\n"
	else:
		Plane = ""
	if Bond_Eq == 0.0:
		Bond_Stretch = ""
	else:
		Bond_Stretch = "%%geom Scan\nB %d %d = %.2f,%.2f,12\nend\nend\n" % (Bond_Atoms[0]-1,Bond_Atoms[1]-1,Bond_Eq-.05,Bond_Eq+.05)
	f.write("! %s %s %s %s %s %s Opt Grid3 FinalGrid5\n%%maxcore %d\n%s%s%s\n*xyz 0 1\n" % (Method,Basis,Aux_Basis,Dispersion_Correction,convergence,Parallel,Memory,H_Opt,Plane,Bond_Stretch))
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("*")
	f.close()

def Write_Orca_ChelpG(File_Name,Test_Molecule,Method = "RI BP86",Basis = "def2-SVP",Aux_Basis = "def2/J",Dispersion_Correction = "D3BJ",convergence = "TIGHTSCF",nproc = 4,Memory = 3000,H_Only = False,Planarize = False,Planarize_Atoms = []):
	f = open(File_Name,'w')
	if nproc == 1:
		Parallel = ""
	else:
		Parallel = "PAL%d" % nproc
	f.write("! RKS B3LYP 6-31+G** NormalSCF NOSOSCF CHELPG\n%%maxcore %d\n\n*xyz 0 1\n" % (Memory))
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("*")
	f.close()

def Write_QChem_SPE(File_Name,Test_Molecule,Exchange_Method = "HF",Correlation_Method = "pRIMP2",Basis = "cc-pvtz",Method = "rimp2",Aux_Basis = "rimp2-cc-pvtz",Memory = 110000,convergence = 6,Implicit_Solvent_Method = "PCM",Implicit_Solvent_Dielectric = 4.9,ChelpG=False):
	f = open(File_Name,'w')
	f.write("$molecule\n\t0 1\n")
	if ChelpG:
		ChelpG_Line = "\tCHELPG\t\tTRUE\n"
	else:
		ChelpG_Line = ''
	if Implicit_Solvent_Dielectric != 0.0:
		Solvent_Line = "\n\n$solvent dielectric 4.9 $end"
	else:
		Solvent_Line = ''
	for atom in Test_Molecule.Atom_List:
		f.write("%s\t%f\t%f\t%f\n" % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("$end\n\n$rem\n\tJOBTYPE\t\tSP\n\tEXCHANGE\t%s\n\tCORRELATION\t%s\n\tBASIS\t\t%s\n\tMETHOD\t\t%s\n\tAUX_BASIS\t%s\n\tSOLVENT_METHOD\t%s\n\tPURECART\t11111\n\tSYMMETRY\tfalse\n\tMEM_TOTAL\t%d\n\tSCF_CONVERGENCE = %d\n\tTHRESH=%d\n%s$end\n\n%s" % (Exchange_Method,Correlation_Method,Basis,Method,Aux_Basis,Implicit_Solvent_Method,Memory,convergence,convergence + 4,ChelpG_Line,Solvent_Line))
	f.close()

def Write_QChem_Optimize_Geometry(File_Name,Test_Molecule,Basis = "def2-SVP",Method = "BP86",Approximations = "\n\tRI-J\t\tTRUE",Aux_Basis = "def2-SVP-J",Memory = 12000,convergence = 6,Implicit_Solvent_Method = "PCM",Implicit_Solvent_Dielectric = 4.9):
	f = open(File_Name,'w')
	f.write("$molecule\n\t0 1\n")
	for atom in Test_Molecule.Atom_List:
		f.write('%s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("$end\n\n$rem\n\tJOBTYPE\t\tOPT\n\tBASIS\t\t%s\n\tMETHOD\t\t%s\n\tDFT_D\t\tD3_BJ\n\tAUX_BASIS\t%s\n\tPURECART\t11111\n\tSYMMETRY\tfalse\n\tMEM_TOTAL\t%d\n\tSCF_CONVERGENCE = %d\n\tTHRESH=%d\n$end" % (Basis,Method,Aux_Basis,Memory,convergence,convergence + 4))
	if Implicit_Solvent_Dielectric != "":
		f.write("\n\n$solvent\nDielectric %.2f\n$end" % Implicit_Solvent_Dielectric)
	f.close()


def Write_NWChem_SPE(File_Name,Test_Molecule,Name,Method = "rimp2",charge = 0,Basis = "cc-pVTZ",Aux_Basis = "cc-pVTZ",Implicit_Solvent_Dielectric = 4.9,ChelpG = False,DFT_Method = 'xc becke88 perdew86\n vdw 3'):
	f = open(File_Name,'w')
	f.write("start %s\ncharge %d\ngeometry\n" % (Name,charge))
	for atom in Test_Molecule.Atom_List:
		f.write(' %s\t%f\t%f\t%f\n' % (atom.Element,atom.Position[0],atom.Position[1],atom.Position[2]))
	f.write("end\nbasis\n * library %s\nend\n" % (Basis))
	if Method == "rimp2":
		f.write("basis \"ri-mp2 basis\"\n * library %s\nend\n" % (Aux_Basis))
	if Method == "dft":
		f.write("dft\n %s\nend\n" % DFT_Method)
	if Implicit_Solvent_Dielectric != 0.0:
		f.write("cosmo\n dielec %f\nend\n" % (Implicit_Solvent_Dielectric))
	f.write("task %s energy" % Method)
	if ChelpG:
		f.write("\ntask esp")
	f.close()
	return

def Write_NWChem_Optimize_Geometry(File_Name,Test_Molecule,Name,Method = "dft",charge = 0,Basis = "def2-SVP",Aux_Basis = "cc-pVTZ",Implicit_Solvent_Dielectric = 0.0):
	return

def Write_LAMMPS_Minimize(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nread_data %s\n\ntimestep 1.0\nminimize 1.0e-8 1.0e-8 1000 100000\nwrite_data %s_Optimized_Geometry.data\nrun 0" % (Coul_Cutoff,Data_File,Name))
	f.close()

def Write_LAMMPS_Nonbonded_Intel(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\treal\natom_style\tfull\npair_style\tlj/charmm/coul/charmm 9.0 10.0 %.1f %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nnewton off\nread_data %s\n\ntimestep 2.0\ngroup frozen molecule == 1\ngroup solvent molecule != 1\nthermo 1000\nthermo_style custom step temp density press etotal pe ke\nneighbor 2.0 bin\nrun_style verlet\nfix 1 solvent npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0\ndump 1 all custom 1000 %s.lammpstrj id type mol x y z ix iy iz\nvariable a loop 100000\nlabel loop\nrun 20\nwrite_data Nonbonded_%s_*.data\nnext a\njump SELF loop" % (Coul_Cutoff,Coul_Cutoff + .2,Data_File,Name,Name))
	f.close()

def Write_LAMMPS_Nonbonded_Single_Point_Energy_Multi(File_Name,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nnewton off\n\n\ntimestep 2.0\ngroup frozen molecule == 1\ngroup solvent molecule != 1\n\nthermo 1000\ncompute nonbond frozen pe/atom pair\ncompute totalnonbond frozen reduce sum c_nonbond\nthermo_style custom step temp density press etotal pe ke c_nonbond c_totalnonbond\nneighbor 2.0 bin\nrun_style verlet\nfix 1 solvent npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0\nvariable a loop 1000000\nvariable b equal 20\nlabel loop\nread_data Nonbonded_%s_\"$b\".data\nrun 0\nvariable b equal $b + 20\nnext a\njump SELF loop" % (Coul_Cutoff,Name))
	f.close()	

def Write_LAMMPS_Averaged_Energy_No_Nonbonded(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nread_data %s.data\n\ntimestep 1.0\ncompute bonded all pe bond angle dihedral improper\nthermo 10\nthermo_style custom step temp density press etotal pe ke c_bonded\nneighbor 2.0 bin\nrun_style verlet\nfix 1 solvent npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0\ndump 1 all custom 1000 %s.lammpstrj id type mol x y z ix iy iz\nrun 2000000\nwrite_data Final_%s.data" % (Coul_Cutoff,Name,Name,Name))
	f.close()

def Write_LAMMPS_Single_Point_Energy(File_Name,Data_File,Name,Coul_Cutoff,dielectric=1.0,Plumed=""):
	f = open(File_Name,'w')
	f.write("units\t\treal\natom_style\tfull\npair_style\tlj/cut/coul/cut 10.0 %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\ndielectric %.1f\nread_data %s.data\n\ntimestep 1.0\n\nthermo 10\nthermo_style custom step temp density press evdwl ecoul ebond eangle edihed eimp etotal pe\nneighbor 2.0 bin\nrun_style verlet\n" % (Coul_Cutoff,dielectric,Name))
	if Plumed != "":
		f.write("fix MetaD all plumed plumedfile %s outfile %s.out\n" % (Plumed,Name))
	f.write("run 0")
	f.close()

def Write_LAMMPS_Test_Condensed():
	return

def Write_LAMMPS_Minimize(File_Name,Data_File,Name,Coul_Cutoff):
	f = open(File_Name,'w')
	f.write("units\t\treal\natom_style\tfull\npair_style\tlj/charmm/coul/charmm 9.0 10.0 %.1f %.1f\nbond_style\tharmonic\nangle_style\tharmonic\ndihedral_style\topls\nspecial_bonds\tlj/coul 0 0 0.5\nimproper_style\tcvff\npair_modify mix geometric\nnewton off\nread_data %s\n\ntimestep 2.0\ngroup frozen molecule == 1\ngroup solvent molecule != 1\nminimize 1e-4 1e-6 100 1000\nwrite_data %s_Minimized.data" % (Coul_Cutoff,Coul_Cutoff+.2,Data_File,Name))

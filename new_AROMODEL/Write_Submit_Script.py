#! usr/bin/python
import math

def Write_Qchem(f,nproc,In_File,Name):
	#f.write('export QCMACHINEFILE=`generate_pbs_nodefile`\ncat $QCMACHINEFILE|uniq >.tmp\nmv .tmp $QCMACHINEFILE\nmodule load qchem\nexport QCSCRATCH=/scratch/$USER/$SLURM_JOBID\nexport QCMACHINEFILE=`generate_pbs_nodefile`\nqchem -nt %d %s %s.out' % (nproc,In_File,Name))
	f.write('\nmodule load qchem\n\nqchem -nt %d %s %s.out' % (nproc,In_File,Name))

def Write_Qchem_Run_Only(f,nproc,In_File,Name):
	f.write('\nqchem -nt %d %s %s.out' % (nproc,In_File,Name))

def Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path):
	f.write('export PATH=%s:%s/bin:$PATH\nexport RSH_COMMAND="/usr/bin/ssh -x"\nexport LD_LIBRARY_PATH=%s/lib\n%s/orca %s >> %s.out\n' % (OMP_Path,Executable_Path,OMP_Path,Executable_Path,In_File,Name))

def Write_Orca_Run_Only(f,nproc,In_File,Name,Executable_Path,OMP_Path):
	f.write('orca %s >> %s.out\n' % (In_File,Name))

def Write_LAMMPS(f,nproc,In_File,Name):
	f.write('module load lammps/2018.12.12-knl\n\nsrun --cpu-bind=cores -n %d lmp_cori -pk intel 0 omp 4 -sf intel -in %s -log log.%s' % (nproc,In_File,Name))

def Write_LAMMPS_KNL(f,nproc,In_File,Name):
	f.write('export KMP_BLOCKTIME=0\n\nexport OMP_PROC_BIND=true\nexport OMP_PLACES=threads\nexport OMP_NUM_THREADS=4\nsource /opt/intel/parallel_studio_xe_2019.3.062/psxevars.sh\n\nmodule load lammps/2018.12.12-knl\n\nsrun --cpu-bind=cores -n %d lmp_cori -pk intel 0 omp 4 -sf intel -in %s -log log.%s\ns' % (nproc,In_File,Name))

def Write_LAMMPS_KNL_Run_Only(f,nproc,In_File,Name):
	f.write('\nsrun --cpu-bind=cores -n %d lmp_cori -pk intel 0 omp 4 -sf intel -in %s -log log.%s' % (nproc,In_File,Name))

def Write_LAMMPS_Run_Only(f,nproc,In_File,Name):
	f.write('srun --cpu-bind-cores -n %d lmp_cori -in %s -log log.%s' % (nproc,In_File,Name))

def Write_NWChem(f,nproc,In_File,constraint='cori'):
	if constraint != 'knl':
		f.write('module load nwchem\nsrun --cpu-bind=cores -np %d nwchem %s > %s' % (nproc,In_File,Out_File))
	else:
		f.write('module load nwchem\nsrun --cpu-bind=cores -S 4 -n %d nwchem %s > %s' % (nproc,In_File,Out_File))
	return

def Write_SLURM(File_Name,In_File,Name,nproc,Cluster_Location,Job_Type,queue="shared",proc_per_node=32,account = "m3047",walltime = 2,Executable_Path = "",OMP_Path = "",constraint='knl',nodes=1):
	f = open(File_Name,'w')
	"""nodes = math.floor(nproc/proc_per_node)
	if nodes == 0:
		nodes = 1
		queue = "shared"
		ppn = nproc
	else:
		queue = "compute"
		ppn = proc_per_node"""
	if constraint == "knl":
		spec_core = "#SBATCH --core-spec=4\n"
	else:
		spec_core = ""
	Out_File = In_File.split('.')[0] + ".out"
	if walltime == 0:
		f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t 0:30:00\n#SBATCH --constraint=%s\ncd %s\n' % (Name,Name,queue,nodes,proc_per_node,account,constraint,Cluster_Location))
	else:
		f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t %d:00:00\n#SBATCH --constraint=%s\ncd %s\n' % (Name,Name,queue,nodes,proc_per_node,account,walltime,constraint,Cluster_Location))
	if Job_Type == "QChem":	
		if Executable_Path == "":
			Write_Qchem(f,nproc,In_File,Name)

	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path)

	if Job_Type == "LAMMPS" and constraint == "knl":
		if Executable_Path == "":
			Write_LAMMPS_KNL(f,nproc,In_File,Name)

	elif Job_Type == "LAMMPS":
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File,Name)

	if Job_Type == "NWChem":
		if Executable_Path == "":
			Write_NWChem(f,nproc,In_File,constraint = constraint)
			
	f.close()

def Write_SLURM_Batch(File_Name,In_File_List,Name,nproc,Cluster_Location,Job_Type,queue="premium",proc_per_node=32,account = "m3047",walltime = 2,Executable_Path = "",OMP_Path = "",constraint = 'haswell',nodes = 1):
	f = open(File_Name,'w')
	"""nodes = math.floor(nproc/proc_per_node)
	if nodes == 0:
		nodes = 1
		queue = "shared"
		ppn = nproc
	else:
		queue = "compute"
		ppn = proc_per_node"""
	if constraint == "knl":
		spec_core = "#SBATCH --core-spec=4\n"
	else:
		spec_core = ""
	if walltime ==0:
		f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t 0:30:00\n#SBATCH --constraint=%s\n%scd %s\n' % (Name,Name,queue,nodes,proc_per_node,account,constraint,spec_core,Cluster_Location))
	else:
		f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t %d:00:00\n#SBATCH --constraint=%s\n%scd %s\n' % (Name,Name,queue,nodes,proc_per_node,account,walltime,constraint,spec_core,Cluster_Location))
	if Job_Type == "QChem":	
		if Executable_Path == "":
			Out_File = In_File_List[0].split('.')[0]
			Write_Qchem(f,nproc,In_File_List[0],Out_File)
			for In_File in In_File_List[1:]:
				Out_File = In_File.split('.')[0]
				Write_Qchem_Run_Only(f,nproc,In_File,Out_File)


	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		for In_File in In_File_List:
			Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path)

	if Job_Type == "LAMMPS" and constraint == "knl":
		if Executable_Path == "":
			Write_LAMMPS_KNL(f,nproc,In_File,Name)

	elif Job_Type == "LAMMPS":
		if Executable_Path == "":
			for In_File in In_File_List:
				Write_LAMMPS(f,nproc,In_File,Name)

	if Job_Type == "NWChem":
		if Executable_Path == "":
			for In_File in In_File_List:
				Write_NWChem(f,nproc,In_File,constraint = constraint)
				
	f.close()

def Write_TORQUE(File_Name,In_File,Name,nproc,Cluster_Location,Job_Type,proc_per_node=28,walltime = 2,Executable_Path = "",OMP_Path = "",queue = "condo"):
	f = open(File_Name,'w')
	nodes = math.floor(nproc/proc_per_node)
	if nodes == 0:
		nodes = 1
		ppn = nproc
	else:
		ppn = proc_per_node
	f.write('#!/bin/bash\n#PBS -N %s\n#PBS -o %s\n#PBS -q %s\n#PBS -l nodes=%d:ppn=%d\n#PBS -l walltime=%d:00:00\n\ncd %s\n\n' % (Name,Name,queue,nodes,ppn,walltime,Cluster_Location))
	if Job_Type == "QChem":	
		if Executable_Path == "":
			Write_Qchem(f,nproc,In_File,Name)


	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		Write_Orca(f,nproc,In_File,Name,Executable_Path,OMP_Path)

	elif Job_Type == "LAMMPS":
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File,Name)
			
	f.close()
import math
from string import Template

def Write_Qchem(f,In_File,Name):
	#f.write('export QCMACHINEFILE=`generate_pbs_nodefile`\ncat $QCMACHINEFILE|uniq >.tmp\nmv .tmp $QCMACHINEFILE\nmodule load qchem\nexport QCSCRATCH=/scratch/$USER/$SLURM_JOBID\nexport QCMACHINEFILE=`generate_pbs_nodefile`\nqchem -nt %d %s %s.out' % (nproc,In_File,Name))
	with open('./Templates/submit_qchem.template') as template_file:
		t = Template(template_file.read().replace('\r\n','\n'))
		f.write(t.substitute({'name':Name, 'in_file' : In_File}))

def Write_Qchem_Run_Only(f,In_File,Name):
	f.write('\nqchem %s >> %s.out' % (In_File,Name))

def Write_Orca(f,In_File,Name):
	with open('./Templates/submit_orca.template') as template_file:
		t = Template(template_file.read().replace('\r\n','\n'))
		f.write(t.substitute({'name':Name, 'in_file' : In_File}))

def Write_Orca_Run_Only(f,In_File,Name):
	f.write('\norca %s >> %s.out' % (In_File,Name))

def Write_LAMMPS(f,nproc,In_File,Name):
	f.write('module load lammps/2018.12.12-knl\n\nsrun --cpu-bind=cores -n %d lmp_cori -pk intel 0 omp 4 -sf intel -in %s -log log.%s' % (nproc,In_File,Name))

def Write_LAMMPS_KNL(f,nproc,In_File,Name):
	f.write('export KMP_BLOCKTIME=0\n\nexport OMP_PROC_BIND=true\nexport OMP_PLACES=threads\nexport OMP_NUM_THREADS=4\nsource /opt/intel/parallel_studio_xe_2019.3.062/psxevars.sh\n\nmodule load lammps/2018.12.12-knl\n\nsrun --cpu-bind=cores -n %d lmp_cori -pk intel 0 omp 4 -sf intel -in %s -log log.%s\ns' % (nproc,In_File,Name))

def Write_LAMMPS_KNL_Run_Only(f,nproc,In_File,Name):
	f.write('\nsrun --cpu-bind=cores -n %d lmp_cori -pk intel 0 omp 4 -sf intel -in %s -log log.%s' % (nproc,In_File,Name))

def Write_LAMMPS_Run_Only(f,nproc,In_File,Name):
	f.write('srun --cpu-bind-cores -n %d lmp_cori -in %s -log log.%s' % (nproc,In_File,Name))

# def Write_NWChem(f,nproc,In_File,constraint='cori'):
# 	if constraint != 'knl':
# 		f.write('module load nwchem\nsrun --cpu-bind=cores -np %d nwchem %s > %s' % (nproc,In_File,Out_File))
# 	else:
# 		f.write('module load nwchem\nsrun --cpu-bind=cores -S 4 -n %d nwchem %s > %s' % (nproc,In_File,Out_File))
# 	return
def Write_SLURM(File_Name,In_File,Job_Name,nproc,Cluster_Location,Job_Type,run_time = 60,tasks_per_node=1):	
	""" 
	A function for creating submit script to be used by system with SLURM scheduler

	File_Name: string, name of the submit script to be generated; 
	In_File: string, input file to be used by simulation TODO: is this actually true
	Job_Name: string, name of the job
	nproc: idfk TODO
	Cluster_Location: 
	Job_Type: string, 
	tasks_per_node: no clue TODO
	run_time: int, max time in min for the job
	
	"""
	hours = int(run_time) // 60
	minutes = int(run_time) % 60

	f = open(File_Name,'w')
	with open('./Templates/submit_slurm.template') as template_file:
		t = Template(template_file.read().replace('\r\n','\n')) #the replace part 'should' turn file from dos to unix
		f.write(t.substitute({'hour':hours,'min':minutes, 'dir' : Cluster_Location, 'tasks_per_node' : tasks_per_node, 'name':Job_Name}))
		f.write('\n')#just in case the template doesn't have this at the end
	#to be deleted
	# nodes = math.floor(nproc/tasks_per_node)
	# if nodes == 0:
	# 	nodes = 1
	# 	queue = "shared"
	# 	ppn = nproc
	# else:
	# 	queue = "compute"
	# 	ppn = tasks_per_node
	# if constraint == "knl":
	# 	spec_core = "#SBATCH --core-spec=4\n"
	# else:
	# 	spec_core = ""
	# Out_File = In_File.split('.')[0] + ".out"
	# if walltime == 0:
	# 	f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH \
	# 	--nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t 0:30:00\n#SBATCH --constraint=%s\ncd %s\n' % (Name,Name,queue,nodes,tasks_per_node,account,constraint,Cluster_Location))
	# else:
	# 	f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t %d:00:00\n#SBATCH --constraint=%s\ncd %s\n' % (Name,Name,queue,nodes,tasks_per_node,account,walltime,constraint,Cluster_Location))
	if Job_Type == "QChem":	
		Write_Qchem(f,In_File,Job_Name)

	if Job_Type == "Orca":
		Write_Orca(f,In_File,Job_Name)

	elif Job_Type == "LAMMPS":
		raise Exception("not implemented yet")
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File,Name)

	# if Job_Type == "NWChem":
	# 	if Executable_Path == "":
	# 		Write_NWChem(f,nproc,In_File,constraint = constraint)
			
	f.close()


def Write_SLURM_Batch(File_Name,In_File_List,Name,Cluster_Location,Job_Type,run_time = 180, tasks_per_node=64):

	print("Writing slurm batch: File_Name: %s,\n Name: %s,\n In_File_List: %s,\n Job_Type: %s"%(File_Name,Name,In_File_List,Job_Type))

	hours = int(run_time) // 60
	minutes = int(run_time) % 60

	f = open(File_Name,'w')
	with open('./Templates/submit_slurm.template') as template_file:
		t = Template(template_file.read().replace('\r\n','\n')) #the replace part 'should' turn file from dos to unix
		f.write(t.substitute({'hour':hours,'min':minutes, 'dir' : Cluster_Location, 'tasks_per_node' : tasks_per_node, 'name':Name}))

	# nodes = math.floor(nproc/tasks_per_node)
	# if nodes == 0:
	# 	nodes = 1
	# 	queue = "shared"
	# 	ppn = nproc
	# else:
	# 	queue = "compute"
	# 	ppn = tasks_per_node
	# if constraint == "knl":
	# 	spec_core = "#SBATCH --core-spec=4\n"
	# else:
	# 	spec_core = ""
	# if walltime ==0:
	# 	f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t 0:30:00\n#SBATCH --constraint=%s\n%scd %s\n' % (Name,Name,queue,nodes,tasks_per_node,account,constraint,spec_core,Cluster_Location))
	# else:
	# 	f.write('#!/bin/bash\n#SBATCH --job-name="%s"\n#SBATCH --output=%s\n#SBATCH --qos=%s\n#SBATCH --nodes=%d\n#SBATCH --ntasks-per-node=%d\n#SBATCH -A %s\n#SBATCH --export=ALL\n#SBATCH -t %d:00:00\n#SBATCH --constraint=%s\n%scd %s\n' % (Name,Name,queue,nodes,tasks_per_node,account,walltime,constraint,spec_core,Cluster_Location))
	if Job_Type == "QChem":	
		Out_File = In_File_List[0].split('.')[0]
		Write_Qchem(f,In_File_List[0],Out_File)
		for In_File in In_File_List[1:]:
			Out_File = In_File.split('.')[0]
			Write_Qchem_Run_Only(f,In_File,Out_File)


	if Job_Type == "Orca":
		Out_File = In_File_List[0].split('.')[0]
		Write_Orca(f,In_File_List[0],Out_File)
		for In_File in In_File_List[1:]:
			Out_File = In_File.split('.')[0]
			Write_Orca_Run_Only(f,In_File,Out_File)

	# if Job_Type == "LAMMPS" and constraint == "knl":
	# 	raise Exception("LAMMPS sucks")
	# 	if Executable_Path == "":
	# 		Write_LAMMPS_KNL(f,nproc,In_File,Name)

	elif Job_Type == "LAMMPS":
		raise Exception("LAMMPS sucks")
		if Executable_Path == "":
			for In_File in In_File_List:
				Write_LAMMPS(f,nproc,In_File,Name)

	f.close()

def Write_TORQUE(File_Name,In_File,Job_Name,nproc,Cluster_Location,Job_Type,tasks_per_node=28,walltime = 2,Executable_Path = "",OMP_Path = "",queue = "condo"):	
#def Write_TORQUE(config_dict, nproc, tasks_per_node=28, walltime=2, queue="condo"):

	# Job_Type = config_dict['Job_Type']
	# Job_Name = config_dict['Job_Name']
	# Executable_Path = config_dict['Executable_Path']
	# OMP_Path = config_dict['OMP_Path']
	# In_File = config_dict['In_File']
	# Cluster_Location = config_dict['Cluster_Location']

	f = open(File_Name,'w')
	print("I just opened %s in Write_Submit_script" %(File_Name))
	nodes = math.floor(nproc/tasks_per_node)
	if nodes == 0:
		nodes = 1
		ppn = nproc
	else:
		ppn = tasks_per_node
	
	f.write('#!/bin/bash\n#PBS -N %s\n#PBS -o %s\n#PBS -q %s\n#PBS -l nodes=%d:ppn=%d\n#PBS -l walltime=%d:00:00\n\ncd %s\n\n' % (Job_Name,Job_Name,queue,nodes,ppn,walltime,Cluster_Location))
	
	if Job_Type == "QChem":	
		raise Exception("not implemented yet")
		if Executable_Path == "":
			Write_Qchem(f,nproc,In_File,Job_Name)


	if Job_Type == "Orca":
		if Executable_Path == "":
			raise Exception("Executable_Path required for Orca")
		if OMP_Path == "":	
			raise Exception("OMP_Path required for Orca")
		Write_Orca(f,nproc,In_File,Job_Name)

	elif Job_Type == "LAMMPS":
		raise Exception("not implemented yet")
		if Executable_Path == "":
			Write_LAMMPS(f,nproc,In_File,Job_Name)
			
	f.close()
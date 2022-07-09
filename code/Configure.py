#! usr/bin/python

"""
    This file contains paths to different executables 
    and template files, so that this software can be 
    configured easily on different machines.
"""


# Local Paths
Aromodel_Path = "C:\\Users\\Leon\\Personal Document\\College class material\\2022 - 2.5 Summer\\aromodel\code"
Template_Path = Aromodel_Path + "/Templates/"


# Remote Paths
Comet_Login = "andrewk@comet.sdsc.edu"
NERSC_Login = "andrewk@edison.nersc.gov"
Comet_Path = "/oasis/scratch/comet/andrewk/temp_project/Aromodel/%s" # % Directory_Name
NERSC_Path = "/scratch2/scratchdirs/andrewk/"
Orca_Path = "/oasis/scratch/comet/cjpais/temp_project/programs/orca_3_0_3_linux_x86-64/orca"

local_dict = {
    "User_Path" : "C:\\Users\\Leon",
    "Aromodel_Path" : "C:\\Users\\Leon\\Personal Document\\College class material\\2022 - 2.5 Summer\\aromodel\code/",
    "Template_Path" : "C:\\Users\\Leon\\Personal Document\\College class material\\2022 - 2.5 Summer\\aromodel\code/Templates/", # Aromodel_Path + "/Templates/"
    "Lammps_Path" : "lmp_stable"
}

orca_dict = { # for aromodel_lib.py , ring.py,
    "Cluster_Login" : "theleonzhang@login.expanse.sdsc.edu",
    "Base_Cluster_Location" : '/expanse/lustre/scratch/theleonzhang/temp_project',
    "Scheduler_Type" : "SLURM",
    "End_Condition" : "Opt_Orca",
    "Shared_File_Location" : "../Shared_Files_Dihedral_Parameterization", #assumes the code is running from aromodel/code
    #TODO add these to Write_Orca(), and modify the template
    "module_name" : "orca/4.2.1", 
    "dependency_name" : "cpu/0.15.4 gcc/9.2.0 openmpi/3.1.6"
}

lammps_dict = { # for aromodel_lib.py, conjugated_polymer.py, ring.py
    "Cluster_Login" : "theleonzhang@login.expanse.sdsc.edu",
    "Base_Cluster_Location" : '/expanse/lustre/scratch/theleonzhang/temp_project',
    "Scheduler_Type" : "SLURM",
    "End_Condition" : "SPE_QChem",
    "Shared_File_Location" :  "../Shared_Files_Dihedral_Parameterization"
}

qchem_dict = { # for aromodel_lib.py, ring.py, was openmp_dict
    "Job_Type" : "Orca",
    "Folder_Name" : "Interring_Bonds",
    "Cluster_Login" :"theleonzhang@login.expanse.sdsc.edu",
    "Base_Cluster_Location" : '/expanse/lustre/scratch/theleonzhang/temp_project',
    "Cluster_Location": "/expanse/lustre/scratch/theleonzhang/temp_project/Interring_Bonds",
    "Scheduler_Type" : "SLURM",
    "End_Condition" : "SPE_QChem",
    "Shared_File_Location" : "../Shared_Files_Dihedral_Parameterization",
    #"Executable_Location" : "/home/andrewk/orca_4_2_0_linux_x86-64_openmpi314",
    #"OpenMP_Location" : "/home/andrewk/openmpi-3.1.4", #TODO: delete?
    "module_name" : "qchem/5.4 ",
    "dependency_name" : "cpu/0.15.4  gcc/10.2.0  mvapich2/2.3.6"
    
}

scratch_dir = "/scratch/andrewk" #Conjugated_polymer.py
molecule_orca = "/Users/andrewkleinschmidt/Library/Orca/orca"


c2c = "scp %s " + Comet_Login + ":" + Comet_Path # % (File, Directory_Name)
c2l = "scp " + Comet_Login + ":" + Comet_Path + "/%s ./" # % (Directory_Name, File)
SBATCH = "sbatch " + Comet_Path + "/%s" # %( Directory_Name, Submit_Script)


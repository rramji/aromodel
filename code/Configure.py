#! usr/bin/python

"""
    This file contains paths to different executables 
    and template files, so that this software can be 
    configured easily on different machines.
"""


# Local Paths
Aromodel_Path = "/Users/andrewkleinschmidt/AROMODEL/AROMODEL/"
Template_Path = Aromodel_Path + "/Templates/"


# Remote Paths
Comet_Login = "andrewk@comet.sdsc.edu"
NERSC_Login = "andrewk@edison.nersc.gov"
Comet_Path = "/oasis/scratch/comet/andrewk/temp_project/Aromodel/%s" # % Directory_Name
NERSC_Path = "/scratch2/scratchdirs/andrewk/"
Orca_Path = "/oasis/scratch/comet/cjpais/temp_project/programs/orca_3_0_3_linux_x86-64/orca"

local_dict = {
    "User_Path" : "/Users/andrewkleinschmidt",
    "Aromodel_Path" : "/Users/andrewkleinschmidt/AROMODEL/AROMODEL/",
    "Template_Path" : "/Users/andrewkleinschmidt/AROMODEL/AROMODEL/Templates/" # Aromodel_Path + "/Templates/"
}

cluster_dict = { # for aromodel_lib.py , ring.py
    "Cluster_Login" : "andrewk@cori.nersc.gov",
    "Base_Cluster_Location" : '/global/cscratch1/sd/andrewk',
    "Scheduler_Type" : "SLURM",
    "End_Condition" : "SPE_QChem",
    "Shared_File_Location" : "/Users/andrewkleinschmidt/Shared_Files_Dihedral_Parameterization"
}

lammps_dict = { # for aromodel_lib.py, conjugated_polymer.py, ring.py
    "Cluster_Login" : "andrewk@tscc-login.sdsc.edu",
    "Base_Cluster_Location" : '/oasis/tscc/scratch/andrewk',
    "Scheduler_Type" : "TORQUE",
    "End_Condition" : "SPE_QChem",
    "Shared_File_Location" : "/Users/andrewkleinschmidt/Shared_Files_Dihedral_Parameterization"
}

openmp_dict = { # for aromodel_lib.py, ring.py
    "Job_Type" : "Orca",
    "Folder_Name" : "Interring_Bonds",
    "Cluster_Login" :lammps_dict["Cluster_Login"],
    "Base_Cluster_Location" : lammps_dict["Base_Cluster_Location"],
    "Cluster_Location":lammps_dict["Base_Cluster_Location"]+"/Interring_Bonds",
    "Scheduler_Type" : "TORQUE",
    "End_Condition" : "Opt_Orca",
    "Executable_Location" : "/home/andrewk/orca_4_2_0_linux_x86-64_openmpi314",
    "OpenMP_Location" : "/home/andrewk/openmpi-3.1.4"
}

scratch_dir = "/scratch/andrewk" #Conjugated_polymer.py
molecule_orca = "/Users/andrewkleinschmidt/Library/Orca/orca"


c2c = "scp %s " + Comet_Login + ":" + Comet_Path # % (File, Directory_Name)
c2l = "scp " + Comet_Login + ":" + Comet_Path + "/%s ./" # % (Directory_Name, File)
SBATCH = "sbatch " + Comet_Path + "/%s" # %( Directory_Name, Submit_Script)


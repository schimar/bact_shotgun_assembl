#! /bin/bash 

# properties = {properties}

module load Anaconda3/2021.04/miniconda-base-2021.04

source $UIBK_CONDA_PROFILE
conda activate smk6 

{exec_job} 



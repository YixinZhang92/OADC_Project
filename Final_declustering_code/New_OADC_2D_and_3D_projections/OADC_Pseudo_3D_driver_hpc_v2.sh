#!/bin/sh
#BATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=02-00:00:00
#SBATCH --account=egdaublab
#SBATCH --job-name=simul1

# to run, use sbatch -n 1 -p computeq OADC_Pseudo_3D_driver_hpc.sh

module load matlab/R2018a

matlab -nodisplay < OADC_Pseudo_3D_driver_v2.m


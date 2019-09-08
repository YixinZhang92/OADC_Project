#!/bin/sh
#BATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01-00:00:00
#SBATCH --account=egdaublab
#SBATCH --job-name=testing_matlab_slurm

# to run, use sbatch -n 1 -p computeq test_matlab_slurm.sh

module load matlab/R2018a

matlab -nodisplay < test_matlab_for_slurm.m


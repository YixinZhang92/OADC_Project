#!/bin/sh
#SBATCH --partition computeq
#SBATCH --cpus-per-task 1
#SBATCH --array 0-9
#SBATCH --time=02-00:00:00
#SBATCH --account=egdaublab

# to run, use sbatch -n 1 -p computeq OADC_Pseudo_3D_driver_hpc.sh

module load matlab/R2018a

function echoMe {
    echo "matlab -nodisplay < OADC_Pseudo_3D_driver('Simul.real.err1_2.incr10.no$SLURM_ARRAY_TASK_ID')" 
    sleep 10
    exit 0
}

echoMe

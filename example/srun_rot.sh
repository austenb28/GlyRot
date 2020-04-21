#!/bin/bash
#SBATCH --job-name=GlyRot
#SBATCH --time=0-05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --partition=high
#SBATCH --exclusive

module load gromacs/5.1

GlyRot.py -gro BChEG_GMX.gro -top BChEG_GMX.top -gex gmx_mpi -cex srun
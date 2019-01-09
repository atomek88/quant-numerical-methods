#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
# should be ok for threads per task this instead of openmp
#SBATCH --exclusive

# set the job name to
#SBATCH --job-name=mpi_nbody

# send output to hello-hybrid.out
#SBATCH --output=mpi2.out

# load the openmpi default module

module load mvapich2

# set OMP_NUM_THREADS to the number of --cpus-per-task we asked for
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


mpirun ./mpi 1000
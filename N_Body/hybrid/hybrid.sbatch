#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
# should be ok for threads per task this instead of openmp
#SBATCH --exclusive

# set the job name to
#SBATCH --job-name=mpi_hybrid_nbody

# send output to hello-hybrid.out
#SBATCH --output=hybrid4.out

# load the openmpi default module

module load mvapich2

# set OMP_NUM_THREADS to the number of --cpus-per-task we asked for
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


mpirun ./hybrid 1000
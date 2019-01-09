#!/ bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

# set the job name to
#SBATCH --job-name=mpi

# send output to hello-hybrid.out
#SBATCH --output=julia


# load the openmpi default module

module load mvapich2

# Run the process with mpirun. Notice -n is not required. mpirun will
# automatically figure out how many processes to run from the slurm options
mpirun ./julia 12000 0
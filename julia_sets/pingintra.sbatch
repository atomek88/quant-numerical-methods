#!/ bin/bash
#SBATCH --time=00:05:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive

# set the job name to
#SBATCH --job-name=mpi_ping_intra


#SBATCH --output=ping_intra


# load the openmpi default module

module load mvapich2


# Run the process with mpirun. Notice -n is not required. mpirun will
# automatically figure out how many processes to run from the slurm options
mpirun ./ping
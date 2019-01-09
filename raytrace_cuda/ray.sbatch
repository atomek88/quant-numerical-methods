# This script will request one GPU device and 1 CPU core

#SBATCH --job-name=gpuSbatch
#SBATCH --output=gpurays.out
#SBATCH --error=gpuSbatch.err
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load cuda

./rays 256 100000000

BATCH --job-name=Job_2d_small
#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --partition=small
#SBATCH --output=slurm-%A_%a.out


module load devtoolset/gcc-11 intel/2022
srun --unbuffered ./cfdft-unopt.x

echo Job is donedone already

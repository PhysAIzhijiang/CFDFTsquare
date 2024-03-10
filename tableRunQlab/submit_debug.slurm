#!/bin/bash -l
#SBATCH --job-name=Job_2d
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=1
#SBATCH --time=10-00:00:00
#SBATCH --partition=debug
#SBATCH --output=slurm-%A_%a.out


module load devtoolset/gcc-11 intel/2022
./cfdft-unopt.x

echo Job is donedone already

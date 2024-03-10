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
./cfdft-unopt.x ../tri_lead231230/SeaFillingBGt0r00.csv -2.5 -1.0

echo Job is donedone already

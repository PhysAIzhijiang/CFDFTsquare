#!/bin/bash
echo copy all script files to current folder
cp ../* .
echo copy done, now you can leave
echo prepare fresh cfdft-unopt.x
module load devtoolset/gcc-11 intel/2022
make cfdft-unopt
echo excutable compipled
sbatch submit_mini.slurm
echo job submitted
squeue
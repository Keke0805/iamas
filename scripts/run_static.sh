#!/bin/bash
#SBATCH -p batch
#SBATCH -J static 
#SBATCH -o static.log
#SBATCH -e static.err
#SBATCH -N 1 -n 1
#SBATCH --ntasks-per-node=1
#####SBATCH --nodelist=node08
date
#####setenv LD_LIBRARY_PATH ${NETCDF}/lib:$LD_LIBRARY_PATH
mpirun -n 1 ./init_atmosphere
date
#if set config_native_gwd_static = true, we need very large memory, so we can only use half cpu per node instead of 18 per node

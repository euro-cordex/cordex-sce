#!/bin/bash

#SBATCH -A esp
#SBATCH -J CAPE_EUR_ERA
#SBATCH --mail-user=rhg11c@my.fsu.edu
#SBATCH --mail-type=ALL
## This is the maximum time allowed.
#SBATCH -t 24:00:00
## The partition name
#SBATCH -p esp1
## This means a total of processors
#SBATCH -N 5 --ntasks-per-node=20
#SBATCH -o /home/netapp-clima/users/rglazer/Euro-11/cape_logs/cape_mpi%j.out
#SBATCH -e /home/netapp-clima/users/rglazer/Euro-11/cape_logs/cape_mpi%j.err


year=$1
domain=$2 #CAM-22 AUS-22 SAM-22
model=$3 #ECMWF-ERAINT MPI-M-MPI-ESM-MR 
scenario=$4 #evaluation historical rcp26 rcp85



set -e
{
mpirun ./cape_mpi_paralell_yr.py $year $domain $model $scenario
}

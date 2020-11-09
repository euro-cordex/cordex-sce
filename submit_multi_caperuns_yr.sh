#!/bin/bash


year=1980 #this is the year it will start on
endyear=1980 #last year you want to do, inclusive
domain="EUR-11"
model="ECMWF-ERAINT" #ECMWF-ERAINT MPI-M-MPI-ESM-MR MOHC-HadGEM2-ES NOAA-GFDL-GFDL-ESM2M
scenario="evaluation" #evaluation historical rcp26 rcp85


jobId=$( sbatch submit_cape_runs_yr.sh ${year} ${domain} ${model} ${scenario} | cut -f 4 -d " ")


lastjobid=${jobId}
echo ${lastjobid}
let year++

#### have these jobs depend on the last job from the previous year ####
while [ $year -le $endyear ]; do


jobId=$( sbatch -d afterok:${lastjobid} submit_cape_runs_yr.sh ${year} ${domain} ${model} ${scenario}  | cut -f 4 -d " ")


lastjobid=${jobId}
echo ${lastjobid}
let year++
done


############################
SCRIPTS:
############################

LI_mpi_paralell.py - python script that will calculate LI given input files of: 

-surface pressure [ps] 
-temperature at 500mb [ta500] 
-temperature and specific humidity at (925, 850, 700) mb [ta925, ta850, ta700, hus925, hus850, hus700] -landmask [sftlf]

The LI is calculated by lifting parcels from the 925, 850 and 700 mb levels. The parcel temperature at 500 mb is then compared to the environmental temperature at 500 mb. The most unstable LI value out of the 3 lifting levels is then selected. If the surface is above one of the lifting levels then the LI calculation there won't take place. If no calculation is available (the surface is above 700mb) then the LI is set to missing. Calculations of the parcel profile are done using the MetPy module.

############

cape_mpi_paralell.py - cape script that will calculate CAPE given input files of:

-surface specific humidity [huss]
-surface air temperature [tas]
-surface pressure [ps]
-temperature and specific humidity at (925, 850, 700, 600, 500, 400, 300, 200) mb 
[ta925, ta850, ta700, ta600, ta500, ta400, ta300, ta200, 
hus925, hus850, hus700, hus600, hus500, hus400, hus300, hus200]
-landmask [sftlf]

The CAPE is calculated from a mixed parcel in the first 100mb of the profile that is given. The script also checks if the surface is above any of the lower levels and then removes any levels below the surface and placing the surface at the bottom of the revised profile. Calculations of the mixed parcel, parcel profile and CAPE are done using the MetPy module.

#############

submit_multi_runs.sh - along with "submit_runs" will submit MPI jobs for every year in a range of years given by the user (1 job per year). It is necessary to set the parameters at the beginning of the script according to the output from your simulations.

HOW TO RUN IT:
1) make sure you've given execution permissions to it
2) set "year" and "endyear" according to the years you want to run
3) set "domain", "model", and "scenario" according the simulations you are running and how your files
are named
4) make sure the name of the "submit_runs" scipt is correct 
5) execute ./

##############

submit_runs.sh - this is the job script from which "submit_multi_runs" will submit jobs based on the years you specifiy. This script you will probably need to modifiy according to your HPC setup.

##############

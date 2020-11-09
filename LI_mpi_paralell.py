#!/usr/bin/env python
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


import numpy as np
from netCDF4 import Dataset, num2date, date2num
import numpy.ma as ma
import subprocess
import metpy.calc as mpcalc
from metpy.plots import SkewT
from datetime import datetime,date,timedelta
from metpy.units import units
import matplotlib.pyplot as plt
import itertools
import time
import sys

#
# python script to parallelize the CAPE calculation
# Requires two external libraries in python, f90nml and mpi4py
#

try:
    from mpi4py import MPI
except:
    print("Please install mpi4py: https://pypi.python.org/pypi/mpi4py/2.0.0")
    sys.exit(1)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

############### Functions #################
#will enumerate a masked array skipping masked values
def maenumerate(marr):
    mask = ~marr.mask.ravel()
    for i, m in itertools.zip_longest(np.ndenumerate(marr), mask):
        if m: yield i

def pprint(args):
   if rank == 0:
      print(args)
###########################################
pprint('# of procs = '+str(size))

try:
    year = sys.argv[1]
    domain = sys.argv[2] #CAM-22 AUS-22 SAM-22
    model = sys.argv[3] #ECMWF-ERAINT MPI-M-MPI-ESM-MR MOHC-HadGEM2-ES NOAA-GFDL-GFDL-ESM2M
    scenario = sys.argv[4] #evaluation historical rcp26 rcp85
except:
    usage("Missing Argument")

pprint('getting basic definitions to all procs')
pprint('year is '+str(year))
pprint('domain is '+str(domain))
pprint('model is '+str(model))
pprint('The scenario is '+str(scenario))
# need to open one file to get some definitions for all procs

year2=int(year)+1
LI_lev = np.array([50000])
proftot = np.array([92500,85000,70000,60000,50000,40000,30000,20000])
proflevs = np.array([925,850,700])
z_levs=np.shape(proflevs)
totz_levs=np.shape(proftot)

tpar500=np.full((totz_levs[0],z_levs[0]),999)

dir_to_files = './input'
qas_path = dir_to_files+'/huss_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_3hr_'+str(year)+'01010300-'+str(year2)+'01010000.nc'

nc_id_qas=Dataset(qas_path, 'r')
lats=nc_id_qas.variables['lat'][:]
lons=nc_id_qas.variables['lon'][:]
dims = np.shape(lats)

## some of the data are 3hrly, get 18z (Europe) only from the files ##
t=nc_id_qas.variables['time'][:]
date_time=num2date(t,units='hours since 1949-12-01 00:00:00',calendar='gregorian')
zero_z_indexes3hr = []
n=0
for date_time in date_time: 
    if date_time.hour == 18:
       zero_z_indexes3hr.append(n)
    n=n+1
days=np.shape(zero_z_indexes3hr)

pprint('done getting dimensions and basic definitions')

## open most of the data needed to send to other procs on the RANK 0 proc ##
if rank == 0:
    ps_path = dir_to_files+'/ps_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_3hr_'+str(year)+'01010300-'+str(year2)+'01010000.nc'
    
    #get landmask file and data
    landmask_path = dir_to_files+'/sftlf_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_fx.nc'
    nc_id_landmask=Dataset(landmask_path, 'r')
    landmask_b=nc_id_landmask.variables['sftlf'][:]
    print('open and retrieve surface variable data')
    nc_id_ps=Dataset(ps_path, 'r')
    ps=nc_id_ps.variables['ps'][zero_z_indexes3hr,:,:]
    ps=np.asarray(ps)
    
    print('get all indexes of 6hrly 0z times')
    #get all the indexes where it is 00z from the 6hrly files
    ta_path = dir_to_files+'/ta500_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_6hr_'+str(year)+'01010600-'+str(year2)+'01010000.nc'
    nc_id_ta=Dataset(ta_path, 'r')
    t=nc_id_ta.variables['time'][:]
    date_time=num2date(t,units='hours since 1949-12-01 00:00:00',calendar='gregorian')
    zero_z_indexes6hr = []
    n=0
    for date_time in date_time: #if you take the date_time array not selecting anything then it will act like a numpy array and you cant get the hour just by doing '.hour' but if you select an individual index from date_time then you can extract the hour with '.hour'. This for loop is just iterating the date_time array itself - its just looping over itself    
       if date_time.hour == 18:
          zero_z_indexes6hr.append(n)
       n=n+1
       
       
    print('getting qa and ta profiles from input data')
    ta500=np.zeros((days[0],dims[0],dims[1]))
    ta=np.zeros((days[0],dims[0],dims[1],z_levs[0]))
    qa=np.zeros((days[0],dims[0],dims[1],z_levs[0]))
    dewtemp=np.zeros((days[0],dims[0],dims[1],z_levs[0]))
    
    temp_500=nc_id_ta.variables['ta500'][zero_z_indexes6hr,:,:]
    ta500=np.asarray(temp_500)
    
    # Get t and q profiles only at 00z times
    for i in range(0,(z_levs[0]-1)+1,1):
       ta_path = dir_to_files+'/ta'+str(proflevs[i])+'_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_6hr_'+str(year)+'01010600-'+str(year2)+'01010000.nc'
       qa_path = dir_to_files+'/hus'+str(proflevs[i])+'_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_6hr_'+str(year)+'01010600-'+str(year2)+'01010000.nc'

       nc_id_ta=Dataset(ta_path, 'r')
       temp_air=nc_id_ta.variables['ta'+str(proflevs[i])][zero_z_indexes6hr,:,:]
       temp_air=np.asarray(temp_air)

       nc_id_qa=Dataset(qa_path, 'r')
       q_air=nc_id_qa.variables['hus'+str(proflevs[i])][zero_z_indexes6hr,:,:]
       q_air=np.asarray(q_air)
       
       ta[:,:,:,i]=temp_air[:,:,:]
       qa[:,:,:,i]=q_air[:,:,:]
    
    p = np.broadcast_to(proflevs, (days[0],dims[0],dims[1],z_levs[0]))
    
    print('calculating dewpoint from q')
    #calculate dewpoint from the qa values
    dewtemp=mpcalc.dewpoint_from_specific_humidity(qa,ta*units.kelvin,(p*100)*units.Pa)


pprint('reached the mpi part')

comm.Barrier()



##########################


values = None
val1 = np.zeros((dims[0],dims[1]))


if rank == 0:
   values = np.zeros((days[0],dims[0],dims[1]))


ntasks=days[0] #each proc will take 1 day from the year
nloop = ntasks / size # you need to loop again if you haven't completed the year
remaining = ntasks % size # this is the modulus, for if the ntasks is not split evenly with the # of procs


#report the time before the loop and after
proc = subprocess.Popen(['date'], stdout=subprocess.PIPE, shell=True)
(out,err) = proc.communicate()
timedate=out.decode('utf-8')
pprint( 'Started calculating CAPE at '+ str(timedate) )

# buffers for sending data
buff = np.zeros([dims[0],dims[1],z_levs[0]])
buffp = np.zeros([dims[0],dims[1]])
bufft = np.zeros([dims[0],dims[1],z_levs[0]])
bufft500 = np.zeros([dims[0],dims[1]])
buffq = np.zeros([dims[0],dims[1],z_levs[0]])
buf2 = np.zeros([dims[0],dims[1]])
p_input = np.zeros([dims[0],dims[1]])
t_input = np.zeros([dims[0],dims[1],z_levs[0]])
t500_input = np.zeros([dims[0],dims[1]])
q_input = np.zeros([dims[0],dims[1],z_levs[0]])
landmask = np.zeros([dims[0],dims[1]])

ncount = 1
startrank=0

### BEGINNING OF MAIN LOOP ###

for timestep in range(rank,ntasks,size):

   #### To begin you need to send data set up on proc 0
   #### to the rest of the procs - each proc needs a 
   #### specific day of data. this requires using buffered arrays
   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
         if ip == 0:
             buffp[:,:] = ps[ip+startrank,:,:]
         else:
             buf2[:,:] = ps[ip+startrank,:,:]   
             comm.Send(buf2,dest=ip,tag=ip)
   else:
       comm.Recv(buffp,source=0,tag=rank)

   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
         if ip == 0:
             bufft500[:,:] = ta500[ip+startrank,:,:]
         else:
             buf2[:,:] = ta500[ip+startrank,:,:]
             comm.Send(buf2,dest=ip,tag=ip)
   else:
       comm.Recv(bufft500,source=0,tag=rank)

   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
         if ip == 0:
             bufft[:,:,:] = ta[ip+startrank,:,:,:]
         else:
             buff[:,:,:] = ta[ip+startrank,:,:,:]
             comm.Send(buff,dest=ip,tag=ip)
   else:
       comm.Recv(bufft,source=0,tag=rank)
#
   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
          if ip == 0:
            buffq[:,:,:] = dewtemp[ip+startrank,:,:,:]
          else:
             buff[:,:,:] = dewtemp[ip+startrank,:,:,:]
             comm.Send(buff,dest=ip,tag=ip)
   else:
       comm.Recv(buffq,source=0,tag=rank)
   
   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
         if ip == 0:
             landmask[:,:] = landmask_b[:,:]
         else:
             buf2[:,:] = landmask_b[:,:]
             comm.Send(buf2,dest=ip,tag=ip)            
   else:
       comm.Recv(landmask,source=0,tag=rank)
   landmask_m = ma.masked_equal(landmask,0)
   p_input = buffp #*units.Pa
   t_input = bufft*units.kelvin
   t500_input = bufft500 #*units.kelvin
   q_input = buffq*units.degC
   
   # COMPUTATION - 1 day per proc
   for i, val in maenumerate(landmask_m): #this will work, it skips any points 
                                              #that are masked and i is the index 
                                              #of non-masked points
      t0 = time.time()
      pprint( i )

   
      t1 = time.time()
      
      #calculate the parcel profile from 925, 850 and 700 if available   proftot*units.Pa
           
      # 925 hPa
      if p_input[i[0],i[1]] > proftot[0]:
         tpar500[:,0] = mpcalc.parcel_profile(proftot*units.Pa,t_input[i[0],i[1],0].to('degC'),q_input[i[0],i[1],0])
      else:
         continue
	 
      # 850 hPa
      if p_input[i[0],i[1]] > proftot[1]:
         tpar500[:,1] = mpcalc.parcel_profile(proftot*units.Pa,t_input[i[0],i[1],1].to('degC'),q_input[i[0],i[1],1])
      else:
         continue
	 
      # 700 hPa
      if p_input[i[0],i[1]] > proftot[2]:
         tpar500[:,2] = mpcalc.parcel_profile(proftot*units.Pa,t_input[i[0],i[1],2].to('degC'),q_input[i[0],i[1],2])
      else:
         continue

      t2 = time.time()
      pprint((t2-t1,'time doing parcel profiles'))
      
      li_vals = t500_input[i[0],i[1]] - tpar500[4,:]
      
      pprint(np.min(li_vals))
      val1[ i[0],i[1] ] = np.min(li_vals)
      
      t6 = time.time()
      pprint((t6-t0,'time in this loop iteration'))



     ## Gather LI values onto the rank 0 processor - send rank 0 the LI ##
     ##              values from all the other processors              ##
   if rank != 0:
      comm.Send(val1,dest=0,tag=rank)
   if rank == 0:
      values[startrank,:,:] = val1[:,:] #the LI on the 0 processor goes into the LI array in the correct spot
      nrecv = size
      if ncount > nloop:
         nrecv = remaining
      for ip in range(1,nrecv,1): #start at 1 because the 0 processor doesn't need to send to itself
         comm.Recv(values[ip+startrank,:,:],source=ip,tag=ip)

   
   startrank = startrank+size
   if ncount > nloop:
      pprint('# of days finished = '+str(int(startrank)-(int(size)-int(remaining))))
   else:
      pprint('# of days finished = '+str(startrank))

   ncount = ncount + 1
   ### END OF MAIN LOOP ###
   
comm.Barrier()
print('completed cape mpi and data sent')

### Set up output file ###

if rank == 0:
   lindex = values[:,:,:]
   print(np.shape(lindex),np.min(lindex),np.max(lindex),np.count_nonzero(lindex))

   
   dir_to_capeout='./output/'
   outfile='LI_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_'+str(year)+'0101-'+str(year2)+'0101.nc'
   out_nf = dir_to_capeout+outfile
   print('output file is '+str(outfile))
   out_id = Dataset(out_nf,'w',format='NETCDF4_CLASSIC')
   #create dims
   lat = out_id.createDimension('lat',dims[0])
   lon = out_id.createDimension('lon',dims[1])
   time = out_id.createDimension('time',None)
   
   #Create coordinate variables for 3-dimensions
   times = out_id.createVariable('time',np.float64, ('time',))
   latitudes = out_id.createVariable('latitude',np.float64, ('lat','lon'))
   longitudes = out_id.createVariable('longitude',np.float64,('lat','lon'))
   li_var = out_id.createVariable('li',np.float64,('time','lat','lon'))
   
   #write data to your variables
   latitudes[:,:] = lats
   longitudes[:,:] = lons
   li_var[:,:,:] = lindex
   
   #set global attributes
   import time
   out_id.description = 'lifted-index '
   out_id.history = 'created ' +time.ctime(time.time())
   out_id.source = 'netCDF4 python' 
   #set variable attributes
   latitudes.units = 'degree_north'
   latitudes.axis = 'X'
   longitudes.units = 'degree_east'
   longitudes.axis = 'Y'
   times.units = 'hours since 1949-12-01 00:00:00'
   times.calendar = 'gregorian'
   
   li_var.units = 'kelvin'
   li_var.long_name = 'most unstable lifted-index at 500hPa'
   
   dates=[]
   #### if we are dealing with monthly files - you need to
   #### specify the year AND month to get the right dates
   for n in range(cape_var.shape[0]):
      dates.append(datetime(int(year),1,1,0) + n * timedelta(hours=24))
   times[:] = date2num(dates,units=times.units, calendar=times.calendar)
   print('time values (in units %s): ' % times.units + '\n', times[:])

   dates = num2date(times[:], units=times.units, calendar=times.calendar)
   print('dates corresponding to time values:\n', dates) 
         
   for varname in out_id.variables.keys():
      vars = out_id.variables[varname]
      print(varname, vars.dtype, vars.dimensions, vars.shape, np.mean(vars))
   out_id.close() #and the file is written!


#report the time before the loop and after
proc = subprocess.Popen(['date'], stdout=subprocess.PIPE, shell=True)
(out,err) = proc.communicate()
timedate=out.decode('utf-8')
pprint( 'Done calculating CAPE at '+ str(timedate) )

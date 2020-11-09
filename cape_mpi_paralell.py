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

def pprint(arg):
   if rank == 0:
      print(arg)
###########################################

pprint('# of procs = '+str(size))

try:
    year = sys.argv[1]
    domain = sys.argv[2] #EUR-11
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

plevs = np.array([925,850,700,600,500,400,300,200])
z_levs=np.shape(plevs)
ex_hour=18 #this is the solar maximum time that you will calculate
           #instability at
dir_to_files = './input'
qas_path = dir_to_files+'/huss_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_3hr_'+str(year)+'01010300-'+str(year2)+'01010000.nc'

nc_id_qas=Dataset(qas_path, 'r')
lats=nc_id_qas.variables['lat'][:]
lons=nc_id_qas.variables['lon'][:]
dims = np.shape(lats)

t=nc_id_qas.variables['time'][:]
date_time=num2date(t,units='hours since 1949-12-01 00:00:00',calendar='gregorian')
zero_z_indexes3hr = []
n=0
   
for date_time in date_time: #if you take the date_time array not selecting anything then it will act like a numpy array and you cant get the hour just by doing '.hour' but if you select an individual index from date_time then you can extract the hour with '.hour'. This for loop is just iterating the date_time array itself - its just looping over itself
    if date_time.hour == ex_hour: ### 0z for CAM and SAM but 6z for AUS, 12z for WAS ####
       zero_z_indexes3hr.append(n)
    n=n+1
days=np.shape(zero_z_indexes3hr)

pprint('done getting dimensions and basic definitions')

if rank == 0:
    ps_path = dir_to_files+'/ps_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_3hr_'+str(year)+'01010300-'+str(year2)+'01010000.nc'
    tas_path = dir_to_files+'/tas_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_3hr_'+str(year)+'01010300-'+str(year2)+'01010000.nc'
    
    
    #get landmask file and data
    landmask_path = dir_to_files+'/sftlf_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_fx.nc'
    nc_id_landmask=Dataset(landmask_path, 'r')
    landmask_b=nc_id_landmask.variables['sftlf'][:]
    
    print('open and retrieve surface variable data')
    dewtemp_s=np.zeros((days[0],dims[0],dims[1]))

    # Get surface data only at 00z times
    nc_id_tas=Dataset(tas_path, 'r')
    tas=nc_id_tas.variables['tas'][zero_z_indexes3hr,:,:]
    tas=np.asarray(tas)
    qas=nc_id_qas.variables['huss'][zero_z_indexes3hr,:,:]
    qas=np.asarray(qas)
    nc_id_ps=Dataset(ps_path, 'r')
    ps=nc_id_ps.variables['ps'][zero_z_indexes3hr,:,:]
    ps=np.asarray(ps)

    
    #calculate sfc dewpoint from the surface t and q 
    dewtemp_s[:,:,:]=mpcalc.dewpoint_from_specific_humidity(qas[:,:,:],tas[:,:,:]*units.kelvin,ps[:,:,:]*units.Pa)
    print('get all indexes of 6hrly 0z times')
    
    #get all the indexes where it is 00z from the 6hrly files
    ta_path = dir_to_files+'/ta925_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_6hr_'+str(year)+'01010600-'+str(year2)+'01010000.nc'
    nc_id_ta=Dataset(ta_path, 'r')
    
    t=nc_id_ta.variables['time'][:]
    date_time=num2date(t,units='hours since 1949-12-01 00:00:00',calendar='gregorian')
    zero_z_indexes6hr = []
    n=0
    for date_time in date_time: #if you take the date_time array not selecting anything then it will act like a numpy array and you cant get the hour just by doing '.hour' but if you select an individual index from date_time then you can extract the hour with '.hour'. This for loop is just iterating the date_time array itself - its just looping over itself    
       if date_time.hour == ex_hour: #### 0z for CAM and SAM but 6z for AUS, 12z for WAS #####
          zero_z_indexes6hr.append(n)
       n=n+1
    print('getting qa and ta profiles from input data')
    ta=np.zeros((days[0],dims[0],dims[1],z_levs[0]))
    qa=np.zeros((days[0],dims[0],dims[1],z_levs[0]))
    dewtemp=np.zeros((days[0],dims[0],dims[1],z_levs[0]))
    
    # Get t and q profiles only at 00z times
    for i in range(0,(z_levs[0]-1)+1,1):
       ta_path = dir_to_files+'/ta'+str(plevs[i])+'_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_6hr_'+str(year)+'01010600-'+str(year2)+'01010000.nc'
       qa_path = dir_to_files+'/hus'+str(plevs[i])+'_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_r1i1p1_ICTP-RegCM4-6_v1_6hr_'+str(year)+'01010600-'+str(year2)+'01010000.nc'

       nc_id_ta=Dataset(ta_path, 'r')
       temp_air=nc_id_ta.variables['ta'+str(plevs[i])][zero_z_indexes6hr,:,:]
       temp_air=np.asarray(temp_air)

       nc_id_qa=Dataset(qa_path, 'r')
       q_air=nc_id_qa.variables['hus'+str(plevs[i])][zero_z_indexes6hr,:,:]
       q_air=np.asarray(q_air)
       
       ta[:,:,:,i]=temp_air[:,:,:]
       qa[:,:,:,i]=q_air[:,:,:]
    
    #make numpy array of plevs so you can use it in the dewpoint calc
    p=np.broadcast_to(plevs,(days[0],dims[0],dims[1],z_levs[0]))
    print('calculating dewpoint from q')
    #calculate dewpoint from the qa values
    dewtemp=mpcalc.dewpoint_from_specific_humidity(qa,ta*units.kelvin,(p*100)*units.Pa)
    
    #calculate mixed layer parcel and then calculate the parcel profile using that mixed layer parcel
    # first lets create the total profile by adding in the sfc arrays to the atm profile of t and q
    
    bottom_indexes=np.zeros((days[0],dims[0],dims[1]))
    
    tot_t=np.concatenate((np.expand_dims(tas,axis=3),ta),axis=3)
    tot_dew=np.concatenate((np.expand_dims(dewtemp_s,axis=3),dewtemp),axis=3)
    tot_p=np.concatenate((np.expand_dims(ps,axis=3),p*100),axis=3)
    
    #we need adjust the profile so that if the surface pressure is lower than 925 we remove those levels in between
    #whichever index is the surface pressue is given by bottom_indexes - initially it is all zeros assuming the surface
    #is the very first one at the bottom 
    for i in range(1,z_levs[0]+1,1):
       tot_p[:,:,:,i]=np.where(tot_p[:,:,:,i] > tot_p[:,:,:,0],tot_p[:,:,:,0],tot_p[:,:,:,i])
       x,y,z = np.where(tot_p[:,:,:,i] == tot_p[:,:,:,0])
       bottom_indexes[x,y,z]=int(i)
    

pprint('reached the mpi part')

comm.Barrier()


#####################################


values = None
values1 = None
val1 = np.zeros((dims[0],dims[1]))
val2 = np.zeros((dims[0],dims[1]))

if rank == 0:
   values = np.zeros((days[0],dims[0],dims[1]))
   values1 = np.zeros((days[0],dims[0],dims[1]))

ntasks=days[0] #each proc will take 1 day from the year
nloop = ntasks / size # you need to loop again if you haven't completed the year
remaining = ntasks % size # this is the modulus, for if the ntasks is not split evenly with the # of procs


#report the time before the loop and after
proc = subprocess.Popen(['date'], stdout=subprocess.PIPE, shell=True)
(out,err) = proc.communicate()
timedate=out.decode('utf-8')
pprint( 'Started calculating CAPE at '+ str(timedate) )

# buffers for sending data
buff = np.zeros([dims[0],dims[1],z_levs[0]+1])
buffp = np.zeros([dims[0],dims[1],z_levs[0]+1])
bufft = np.zeros([dims[0],dims[1],z_levs[0]+1])
buffq = np.zeros([dims[0],dims[1],z_levs[0]+1])
buf2 = np.zeros([dims[0],dims[1]])
p_input = np.zeros([dims[0],dims[1],z_levs[0]+1])
t_input = np.zeros([dims[0],dims[1],z_levs[0]+1])
q_input = np.zeros([dims[0],dims[1],z_levs[0]+1])
bot_index = np.zeros([dims[0],dims[1]])
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
             buffp[:,:,:] = tot_p[ip+startrank,:,:,:]
         else:
             buff[:,:,:] = tot_p[ip+startrank,:,:,:]
             comm.Send(buff,dest=ip,tag=ip)
   else:
       print(np.shape(buffp),type(buffp))
       comm.Recv(buffp,source=0,tag=rank)

   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
         if ip == 0:
             bufft[:,:,:] = tot_t[ip+startrank,:,:,:]
         else:
             buff[:,:,:] = tot_t[ip+startrank,:,:,:]
             comm.Send(buff,dest=ip,tag=ip)
   else:
       comm.Recv(bufft,source=0,tag=rank)

   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
          if ip == 0:
            buffq[:,:,:] = tot_dew[ip+startrank,:,:,:]
          else:
             buff[:,:,:] = tot_dew[ip+startrank,:,:,:]
             comm.Send(buff,dest=ip,tag=ip)
   else:
       comm.Recv(buffq,source=0,tag=rank)

   if rank == 0:
      nsend = size
      if ncount > nloop:
         nsend = remaining
      for ip in range(nsend):
         if ip == 0:
             bot_index[:,:] = bottom_indexes[ip+startrank,:,:]
         else:
             buf2 = bottom_indexes[ip+startrank,:,:]
             comm.Send(buf2,dest=ip,tag=ip)
   else:
       comm.Recv(bot_index,source=0,tag=rank)
   
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
   p_input = buffp*units.Pa
   t_input = bufft*units.kelvin
   q_input = buffq*units.degC
   # CAPE COMPUTATION - 1 day per proc
   for i, val in maenumerate(landmask_m): #this will work, it skips any points 
                                              #that are masked and i is the index 
                                              #of non-masked points
      ##t0 = time.time()
      pprint( i )

   
      ##t1 = time.time()
      
      # calculate the mixed parcel using the profiles corrected for lvls below the ground
      p,tmp,dew = mpcalc.mixed_parcel(p_input[i[0],i[1],int(bot_index[i[0],i[1]]):],t_input[i[0],i[1],
         int(bot_index[i[0],i[1]]):].to('degC'),
         q_input[i[0],i[1],int(bot_index[i[0],i[1]]):])
	 
      ##t2 = time.time()
      #pprint(t2-t1,'time calculating mixed parcel')
      
      # below is necessary for the parcel profile lcl function - an improved function for parcel profile that gives 
      # you back the lcl point so that a bug doesnt happen
      p_input[i[0],i[1],int(bot_index[i[0],i[1]])]=p # this will put the pressure value of the
                                                               # mixed parcel at the bottom (surface) of the profile 
      t_input[i[0],i[1],int(bot_index[i[0],i[1]])]=tmp.to('kelvin') # this will put the temp 
                                                                              #of the mixed parcel at the
                                                                              # bottom (surface) of the profile
      q_input[i[0],i[1],int(bot_index[i[0],i[1]])]=dew # this will put the dewpt of 
                                                                   #the mixed parcel at the bottom 
                                                                   #(surface) of the profile 

      #calculate parcel profile - lift mixed layer parcel to the lcl where it become saturated and then lift at
      #moist adiabatic lapse rate - what is returned from this function is the parcel temp profile in the atm
      ##t3 = time.time()
      p_pp,amb_tmp,amb_dpt,pp=mpcalc.parcel_profile_with_lcl(p_input[i[0],i[1],
         int(bot_index[i[0],i[1]]):],t_input[i[0],i[1],
         int(bot_index[i[0],i[1]]):],q_input[i[0],i[1],
         int(bot_index[i[0],i[1]]):])
      ##t4 = time.time()
      #pprint(t4-t3,'time of parcel profile calculation')
      try:
         cape,cin = mpcalc.cape_cin(p_pp,amb_tmp,amb_dpt.to('kelvin'),pp)
      # remember that if you are not on rank 0 but you encounter one of these errors, it will not print it out
      # because you have pprint here - the error might be happening on other processors but you wont know
      except IndexError:
         pprint('an index error has occurred, cape calc is bugged here') 
         continue
      except RuntimeError:
         pprint('The cape/cin calculation failed to converge - this point should be set to nan')
         continue
      except:
         pprint('something else has gone wrong')
         continue
      ##t5 = time.time()
      #pprint(t5-t4,'time of cape calculation')

      pprint(cape)
      
      val1[ i[0],i[1] ] = cape.magnitude
      val2[ i[0],i[1] ] = cin.magnitude

      ##t6 = time.time()
      #pprint(t6-t0,'time in this loop iteration')



     ## Gather cape values onto the rank 0 processor - send rank 0 the cape ##
     ##              values from all the other processors                   ##
   if rank != 0:
      comm.Send(val1,dest=0,tag=rank)
   if rank == 0:
      values[startrank,:,:] = val1[:,:] #the cape on the 0 processor goes into the cape array in the correct spot
      nrecv = size
      if ncount > nloop:
         nrecv = remaining
      for ip in range(1,nrecv,1): #start at 1 because the 0 processor doesn't need to send to itself
         comm.Recv(values[ip+startrank,:,:],source=ip,tag=ip)


   if rank != 0:
      comm.Send(val2,dest=0,tag=rank))
   if rank == 0:
      values1[startrank,:,:] = val2[:,:]
      nrecv = size
      if ncount > nloop:
         nrecv = remaining
      for ip in range(1,nrecv,1):
         comm.Recv(values1[ip+startrank,:,:],source=ip,tag=ip)
	 
	 
   startrank = startrank+size
   if ncount > nloop:
      pprint('# of days finished = '+str(int(startrank)-(int(size)-int(remaining))))
   else:
      pprint('# of days finished = '+str(startrank))

   ncount = ncount + 1
   ### END MAIN LOOP ###
comm.Barrier()
print('completed cape mpi and data sent')

### Set up output file ###
if rank == 0:
   cape = values[:,:,:]
   cin = values1[:,:,:]
   print(np.shape(cape),np.min(cape),np.max(cape),np.count_nonzero(cape))
   print(np.shape(cin),np.min(cin),np.max(cin),np.count_nonzero(cin))
   
   dir_to_capeout='./output/'
   outfile='cape_'+str(domain)+'_'+str(model)+'_'+str(scenario)+'_'+str(year)+'0101-'+str(year2)+'0101.nc'
   out_nf = dir_to_capeout+outfile
   print('output file is '+str(outfile))
   out_id = Dataset(out_nf,'w',format='NETCDF4_CLASSIC')
   #creat dims
   lat = out_id.createDimension('lat',dims[0])
   lon = out_id.createDimension('lon',dims[1])
   time = out_id.createDimension('time',None)
   
   #Create coordinate variables for 3-dimensions
   times = out_id.createVariable('time',np.float64, ('time',))
   latitudes = out_id.createVariable('latitude',np.float64, ('lat','lon'))
   longitudes = out_id.createVariable('longitude',np.float64,('lat','lon'))
   cape_var = out_id.createVariable('cape',np.float64,('time','lat','lon'))
   cin_var = out_id.createVariable('cin',np.float64,('time','lat','lon'))
   
   #write data to your variables
   latitudes[:,:] = lats
   longitudes[:,:] = lons
   cape_var[:,:,:] = cape
   cin_var[:,:,:] =  cin
   
   #set global attributes
   import time
   out_id.description = 'cape and cin '
   out_id.history = 'created ' +time.ctime(time.time())
   out_id.source = 'netCDF4 python' 
   #set variable attributes
   latitudes.units = 'degree_north'
   latitudes.axis = 'X'
   longitudes.units = 'degree_east'
   longitudes.axis = 'Y'
   times.units = 'hours since 1949-12-01 00:00:00'
   times.calendar = 'gregorian'
   
   cape_var.units = 'joules/kg'
   cape_var.long_name = 'Convective Available Potential Energy'
   cin_var.units = 'joules/kg'
   cin_var.long_name = 'Convective Inhibition'
   
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

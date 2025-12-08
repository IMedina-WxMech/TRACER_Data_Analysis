#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 10:32:17 2023

@author: isaac
"""

import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
from glob import glob
from siphon.catalog import TDSCatalog


### Paths ### 
Read_Path = '/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houdlfptM1.b1/*'
#Save_Path = '/Users/isaac/Desktop/Composite/'
Save_Path = '/Users/isaac/Desktop/DLFP/'
Prof_Top = 200
Hourly_Files = sorted(glob(Read_Path))

### See how many days we have ###
Dates = []
Data_Points =[]
Heights=[]
Points =0
for DL in range (len(Hourly_Files)):
    try:
        DF = Dataset(Hourly_Files[DL])
        print(f"This file was good: {Hourly_Files[DL]}")
    except OSError:
        print(f"This file was bad: {Hourly_Files[DL]}")
        continue
    Date  = str(Hourly_Files[DL][-19:-11])
    Day = str(Hourly_Files[DL][-13:-11])

    if DL == (len(Hourly_Files)-1):
        Dates.append(int(Date))

        Heights.append(Prof_Top)
        Points = Points + DF.dimensions['time'].size
        Data_Points.append(Points)
        Points = 0 
        print(' Last Day Read:',Date)
    else:
        if int(Hourly_Files[DL][-13:-11]) == int(Hourly_Files[DL+1][-13:-11]):
            Points = Points +  DF.dimensions['time'].size
            continue
        else:
            Dates.append(int(Date))
            Heights.append(Prof_Top)
            Points = Points +  DF.dimensions['time'].size
            Data_Points.append(Points)
            Points = 0 
            print('Day Read:',Date)
   
print('Done Reading in Days! There is',len(Dates),'days available.')
print('--')
print('Splitting Files')
### Create some Data frames ###
### Only making what I need: Range, time_offset, intensity, backscatter, and radial_velocity ###

for DAY in range (len(Dates)):
    counter = 0
    position =[]
    Range = np.full((Heights[DAY]),np.nan, np.double)
    Time = np.full((Data_Points[DAY]),np.nan, np.double)
    SNR = np.full((Data_Points[DAY],Heights[DAY]),np.nan, np.double)
    W = np.full((Data_Points[DAY],Heights[DAY]),np.nan, np.double)
    ABS = np.full((Data_Points[DAY],Heights[DAY]),np.nan, np.double)
    
    Data_Date = Dates[DAY]
### Break down the data and put it into the Dataframes   
    for DL in range (len(Hourly_Files)):
        ## Seperate by day, based on whats in the Dates list ##
        if int(Hourly_Files[DL][-19:-11]) == Data_Date:
           Temp_DF = Dataset(Hourly_Files[DL])
           Temp_Range = Temp_DF['range'][:Prof_Top].data
           Temp_TO = Temp_DF['time_offset'][:].data
           Temp_SNR = Temp_DF['intensity'][:,:Prof_Top].data
           Temp_W = Temp_DF['radial_velocity'][:,:Prof_Top].data
           Temp_ABS = Temp_DF['attenuated_backscatter'][:,:Prof_Top].data
           Time_points = len(Temp_TO)
           TLat = Temp_DF['lat'][0]
           TLon = Temp_DF ['lon'][0]
           
           
           if counter == 0:
               Range[:] = Temp_Range
               Time [0:Time_points] = Temp_TO
               SNR[0:Time_points,:] = Temp_SNR
               W[0:Time_points,:] = Temp_W
               ABS[0:Time_points,:]= Temp_ABS
           else: 
               Time [counter:counter+Time_points] = Temp_TO
               SNR[counter:counter+Time_points,:] = Temp_SNR
               W[counter:counter+Time_points,:] = Temp_W
               ABS[counter:counter+Time_points,:]= Temp_ABS
               
           counter = counter + Time_points
           
        else:
            counter =0   
    
    print("Saving File...")
### After theyre seperated, save the new file ###
    output_file = nc.Dataset(Save_Path+str(Dates[DAY])+'_ARM_DLFP.nc', 'w', clobber=True, format='NETCDF3_64BIT')
                
                
    # global attributes
    output_file.title = 'Combined files from the ARM DLFP'
    output_file.author = 'Isaac Medina'
    output_file.contact = 'isaac.j.medina@ou.edu'
    output_file.reference = 'Honestly if the data is saved weird Im sorry'
    output_file.site = 'ARM Peckham'
    output_file.basetime = str(Dates[DAY])
    output_file.platform = 'Lidar FP'
    #output_file.latitude = data['dlat'][int(site)]
    #output_file.longitude = data['dlon'][int(site)]


    # define dimensions       (name,value)
    output_file.createDimension('time', len(Time)) #time dimension, 10min
    output_file.createDimension('height', len(Range)) #height variable 
    output_file.createDimension('pos',2)


    #hourly output
    # create a variable file.createVariable(name, precision, dimensions) = values (usually some array)    
    output_file.createVariable('time_offset','f8',('time'))[:] = Time
    setattr(output_file.variables['time_offset'],'units','Hours')
    setattr(output_file.variables['time_offset'],'description','time since 0z')

    output_file.createVariable('height','f8',('height'))[:] = Range
    setattr(output_file.variables['height'],'units','m')
    setattr(output_file.variables['height'],'description','height in meters')

    output_file.createVariable('intensity','f8',('time','height'))[:] = SNR
    setattr(output_file.variables['intensity'],'units','unitless')
    setattr(output_file.variables['intensity'],'description','Intensity (signal to noise ratio + 1)')

    output_file.createVariable('radial_velocity','f8',('time','height'))[:] = W
    setattr(output_file.variables['radial_velocity'],'units','m/s')
    setattr(output_file.variables['radial_velocity'],'description','Radial velocity (Pretty much vertical velocity I think)')
    
    output_file.createVariable('backscatter','f8',('time','height'))[:] = ABS
    setattr(output_file.variables['backscatter'],'units','1/(m sr)')
    setattr(output_file.variables['backscatter'],'description','Attenuated backscatter')
    
    output_file.createVariable('pos','f8',('pos'))[:] = [TLat,TLon]
    setattr(output_file.variables['pos'],'units','lat /lon')
    setattr(output_file.variables['pos'],'description','Lattitude and Longitude')

    # close it up
    output_file.close()
    print("File written: "+str(Dates[DAY])+'_ARM_DLFP.nc')
    
    print(' --- ' ,str(Dates[DAY]),'COMPLETE ---')  
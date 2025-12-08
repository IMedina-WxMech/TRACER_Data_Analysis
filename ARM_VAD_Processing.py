#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 16:06:56 2022

@author: joshua.gebauer
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import glob
from datetime import datetime, timedelta
import datetime as dt
import Lidar_functions
from netCDF4 import Dataset
from siphon.catalog import TDSCatalog


#%%

files = sorted(glob.glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houdlppiM1.b1/*'))
data = Dataset(files[0])
print(data)
#%%

startdate = datetime(2022,6,1)
enddate = datetime(2022,8,1)

curr = startdate

while (curr <= enddate):
    
    print('Creating VAD file for ' + curr.strftime("%Y%m%d"))
    file_path = '/Users/isaac/Desktop/DLVAD/'

    uu = []
    vv = []
    ww = []
    speed = []
    wdir = []
    res = []
    cor = []
    u_err = []
    v_err = []
    w_err = []
    times = []
    hour = []
    temp_intensity =[]

   
    #files = files + sorted(glob.glob('/Users/joshua.gebauer/mnt/joshua.gebauer/ARM_data/lidar' + '/' + 'sgpdlppi*' + curr.strftime("%Y%m%d") + '*.cdf'))
   # files = (TDSCatalog("https://archive.arm.gov/orders/catalog/orders/medinai2/240625/catalog.html?ticket=ST-15000-iHRl-uVPjcOPqXcaEDVuxLK457wsso")).datasets
    files =sorted(glob.glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houdlppiM1.b1/*'))
    for i in (range(len(files))):
        if int(str(files[i])[-19:-11]) == int(str(curr)[0:4]+str(curr)[5:7]+str(curr)[8:10]):
            try:
                fid = Dataset(files[i])
            except OSError:
                print(f"This file was bad: {files[i]}")
                continue
            
            bt = fid.variables['base_time'][0]
            to = fid.variables['time_offset'][:].data
            for t in range (0,len(fid.variables['intensity'][:])):
                intense = fid.variables['intensity'][t,:].data
                temp_intensity.append(intense)
            rngx = fid.variables['range'][:].data/1000.
            azx = fid.variables['azimuth'][:].data
            elx = fid.variables['elevation'][:].data[0]
            vrx = fid.variables['radial_velocity'][:,:].data
            snrx = 10*np.log10(fid.variables['intensity'][:,:].data - 1)
            
            fid.close()
            
            secs = bt+np.nanmean(to)
            time = datetime.utcfromtimestamp(secs)
            
            foo = np.where((snrx < -20) | (snrx > -5))
            
            vrx[foo] = np.nan
            
            vad = Lidar_functions.ARM_VAD(vrx,rngx,elx,azx)
            
            uu.append(np.copy(vad.u[0]))
            vv.append(np.copy(vad.v[0]))
            ww.append(np.copy(vad.w[0]))
            speed.append(np.copy(vad.speed[0]))
            wdir.append(np.copy(vad.wdir[0]))
            res.append(np.copy(vad.residual[0]))
            cor.append(np.copy(vad.correlation[0]))
            u_err.append(np.copy(vad.du[0]))
            v_err.append(np.copy(vad.dv[0]))
            w_err.append(np.copy(vad.dw[0]))
            times.append(secs)
            
            hour.append(time.hour + time.minute/60. + time.second/3600.)

    nc_file = Dataset(file_path + 'CLAMPS_VAD_' + curr.strftime("%Y%m%d") + '.nc' ,'w',format='NETCDF4')
    nc_file.createDimension('time',None)
    nc_file.createDimension('height',len(vad.z))
        
    base_time = nc_file.createVariable('base_time','i')
    base_time[:] = times[0]
    base_time.string = 'Model start time'
    base_time.long_name = 'Base time in model time'
    base_time.units =  'seconds since start of model'
    base_time.ancillary_variables = 'time_offset'
    
    time_offset = nc_file.createVariable('time_offset','d','time')
    time_offset[:] = np.array(times) - times[0]
    time_offset.long_name = 'Time offset from base_time'
    time_offset.units = 'seconds' 
    time_offset.ancillary_variables = "base_time"
    
    height = nc_file.createVariable('height','f','height')
    height[:] = vad.z[:]
    height.long_name = 'Height above ground level'
    height.units = 'm'
    height.standard_name = 'height'
        
    elevation_angle = nc_file.createVariable('elevation_angle','f','time')
    elevation_angle[:] = vad.el
    elevation_angle.long_name = 'Beam elevation angle'
    elevation_angle.units = 'degree'
    
    u = nc_file.createVariable('u','f',('time','height'))
    u[:,:] = np.array(uu)
    u.long_name = 'Eastward component of wind vector'
    u.units = 'm/s'
    
    u_error = nc_file.createVariable('u_error','f',('time','height'))
    u_error[:,:] = np.array(u_err)
    u_error.long_name = 'Estimated error in eastward component of wind vector'
    u_error.units = 'm/s'
    
    v = nc_file.createVariable('v','f',('time','height'))
    v[:,:] = np.array(vv)
    v.long_name = 'Northward component of wind vector'
    v.units = 'm/s'
    
    v_error = nc_file.createVariable('v_error','f',('time','height'))
    v_error[:,:] = np.array(v_err)
    v_error.long_name = 'Estimated error in northward component of wind vector'
    v_error.units = 'm/s'
    
    w = nc_file.createVariable('w','f',('time','height'))
    w[:,:] = np.array(ww)
    w.long_name = 'Vertical component of wind vector'
    w.units = 'm/s'
    
    w_error = nc_file.createVariable('w_error','f',('time','height'))
    w_error[:,:] = np.array(w_err)
    w_error.long_name = 'Estimated error in vertical component of wind vector'
    w_error.units = 'm/s'
    
    wind_speed = nc_file.createVariable('wind_speed','f',('time','height'))
    wind_speed[:,:] = np.array(speed)
    wind_speed.long_name = 'Wind speed'
    wind_speed.units = 'm/s'
    
    wind_direction = nc_file.createVariable('wind_direction','f',('time','height'))
    wind_direction[:,:] = np.array(wdir)
    wind_direction.long_name = 'Wind direction'
    wind_direction.units = 'degree'
        
    residual = nc_file.createVariable('residual','f',('time','height'))
    residual[:,:] = np.array(res)
    residual.long_name = 'Fit residual'
    residual.units = 'm/s'
    
    correlation = nc_file.createVariable('correlation','f',('time','height'))
    correlation[:,:] = np.array(cor)
    correlation.long_name = 'Fit correlation coefficient'
    correlation.units = 'unitless'
        
    intensity = nc_file.createVariable('intensity','f',('time','height'))
    intensity[:,:] = np.array(temp_intensity)
    intensity.long_name = 'SNR Stuff'
    intensity.units = 'IDK'
    
    nc_file.history = 'created on ' + dt.datetime.utcnow().strftime('%Y/%m/%d %H:%M:%S UTC')
    
    nc_file.close()
        
    curr = curr + timedelta(days=1)
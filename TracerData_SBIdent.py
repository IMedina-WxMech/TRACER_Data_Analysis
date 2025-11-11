#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 10:57:53 2025

@author: isaac
"""

# =============================================================================
"""
Section I: Package Read in and Defined Functions
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['mathtext.fontset'] = 'stix'  
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Cambria Math'] + plt.rcParams['font.serif']
from netCDF4 import Dataset
from glob import glob
from scipy import signal
## Only used to get the date for graphing. This can also be done with hardcoding if this dependency is an issue
import datetime 



def butter_lowpass_filter(data, cutoff, sample_rate, order=5):
    """
    Applies a Butterworth low-pass filter to the input data.

    Parameters:
    - data: The input signal to be filtered.
    - cutoff: The cutoff frequency of the filter in Hz.
    - sample rate: The sampling rate of the signal in Hz.
    - order: The order of the Butterworth filter.

    Returns:
    - filtered_data: The filtered data.
    """
    nyquist = 0.5 * sample_rate
    normal_cutoff = cutoff / nyquist
    sos = signal.butter(order, normal_cutoff, btype='low', analog=False, output='sos')
    filtered_data = signal.sosfiltfilt(sos, data)
    return filtered_data


def windowed_average(data, time, window_size):
    """
    Perform a moving average
    
    Args:
      data: 1-d dataset
      time: 1-d tine array(for graphing purposes)
      window_size: size of averaging window (careful of time units)
      
    Returns:
    	result: averaged data
    """
    result = []
    time_array = []
    for i in range(0,int(len(data)/ window_size)+1):
        if i<int(len(data)/ window_size):
            window = data[int(i*window_size):int((i+1)*window_size)]
        if i == int(len(data)/ window_size):
            window = data[int(i*window_size):]
        average = sum(window) / window_size
        time_array.append(time[i*window_size]+time[int(window_size/2)])
        result.append(average)
    return result, time_array

    
def First_Derivative (U,Nx,dx):
    """
    Perform the first derivative of a given 1d-dataset using  2nd-order Central Differencing
    
    Args:
      U: Dataset
      Nx: size of dataset
      dx = spacing of dataset (here dt)
      
    Returns:
    	U_x: First x-derivative 
    """
    U_x = np.zeros((Nx)) 
    for x_int in range (1,Nx-1):
        U_x[x_int] =( U[x_int+1]-U[x_int-1]) / (2*dx)
    
    return (U_x)

#%%
# =============================================================================
"""
Section II: Data pathing/ read in, Pre-Processing, and Defined Constants 
All paths are defined with local 
All variables are named based on what physical variable they relate to such that:
Temperature: Temp
Water Vapor: WV
Wind Direction: WDir
Aerosol Data: Aero
"""
## Defined Constants, primarily for plotting

Plot_Show = True
Plot_Save = False 
## Wind data also has a height component
## Here we define what height we want (in m)
## Min height = 90m, Max height = 4325m
Wind_Height = 250

## There is a package that does this, but I will do it manually for now
st_h =  6 #sunrise hour
st_m = 20 #sunrise minute
sunrise_time = ((6+5)*3600)+(st_m*60) ## sunrise time in seconds

# If you do not want to loop through all data available you can hard code a singal file here
# This is days after 6-1 (so June 1st would be 0, July 1st would be 30)
Day_Index = 21

## Data pathing

Data_Path = '/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/'

## Mulitple paths for AOD data depending on which one you want to use

#AOD_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houaodbe5chM1.c1/*'))
AOD_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houaoscpcfM1.b1/*'))
#AOD_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houaoscpcuM1.b1/*'))

# Ultra high was planned to be used to identify what aersols were present to attempt to correlate aero data with 
# traditional energy generation. However, due to time constraints this was not feasible. These files are NOT on github
#AOD_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houaosuhsasM1.b1/*'))


MWR_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houmwrlosM1.b1/*'))
Lidar_Data_path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houdlprofwind4newsM1.c1/*'))

Num_files = len(Lidar_Data_path) ## Number of Files for looping over entire dataset
## This loop is commented out but can be uncommented and all following lines can be indented (command + ]) to process ALL files at once
#for Day_Index in range (0,Num_files):
    
## Data Read in 
Aero_Data = Dataset(AOD_Data_Path[Day_Index])
MWR_Data = Dataset(MWR_Data_Path[Day_Index])
Lidar_Data = Dataset(Lidar_Data_path[Day_Index])

st_index_MWR = np.where(np.abs(MWR_Data['time'][:]-sunrise_time) == np.min( np.abs(MWR_Data['time'][:]-sunrise_time)   ))[0][0]
st_index_Lidar = np.where(np.abs(Lidar_Data['time'][:]-sunrise_time) == np.min( np.abs(Lidar_Data['time'][:]-sunrise_time)   ))[0][0]

height_target_array = np.abs(Wind_Height- Lidar_Data['height'][:])
z_wind = np.where (height_target_array == np.min(height_target_array) )[0][0]

## Date gathered through datetime 
Date = (datetime.datetime.fromtimestamp(Aero_Data['time'][-1]+Aero_Data['base_time'][-1]).strftime('%Y-%m-%d %H:%M:%S'))[0:10]   

## Filtering read in data with lowpass filter function
Aero_Filt = butter_lowpass_filter((Aero_Data['concentration'][:]),1/700,7)
Temp_Filt = butter_lowpass_filter((MWR_Data['tkair'][:]),1/50,2)
WV_Filt = butter_lowpass_filter((MWR_Data['vap'][:]),1/50,5)
WDir_Filt = butter_lowpass_filter((Lidar_Data['wind_direction'][:,z_wind]),1/10,1)

#%%
# =============================================================================
"""
Section III: Data Processing

"""
## Finding gradients in the filtered data 
Grad_Aero = First_Derivative(Aero_Filt, len(Aero_Filt), Aero_Data['time'][1]/3600-Aero_Data['time'][0]/3600)
Grad_Temp = First_Derivative(Temp_Filt, len(Temp_Filt), MWR_Data['time'][1]/3600-MWR_Data['time'][0]/3600)
Grad_WV = First_Derivative(WV_Filt, len(WV_Filt), MWR_Data['time'][1]/3600-MWR_Data['time'][0]/3600)
Grad_WDir = First_Derivative(WDir_Filt, len(WDir_Filt), Lidar_Data['time'][1]-Lidar_Data['time'][0])

## Finding the max/min (depending on variable) of the gradient to find SB passage time  
Aero_Max_index = np.where(Grad_Aero==np.max(Grad_Aero))[0][0]
WV_Max_index = np.where(Grad_WV[st_index_MWR:]==np.max(Grad_WV[st_index_MWR:]))[0][0]+st_index_MWR #Only daytime 
Temp_Max_index = np.where(Grad_Temp==np.max(Grad_Temp))[0][0]
WDir_Max_index = np.where(Grad_WDir[st_index_Lidar:]==np.min(Grad_WDir[st_index_Lidar:]))[0][0]+st_index_Lidar

##Time array in dataset is in seconds, so 3600 = 1 hour
window = int(3600)
Average_Aero,Avg_Aero_T = windowed_average(Aero_Filt,Aero_Data['time'][:]/3600,window)


#%%
# =============================================================================
"""
Section IV: Graphing
Plot saving and show variables are in Section II

"""
## Wind Speed Figure

fig, ax = plt.subplots(figsize=(30,20))
rc('font',weight='normal',size=40) 
ax.plot(Lidar_Data['time'][:]/3600,Lidar_Data['wind_direction'][:,z_wind],lw=5,color='k',label='raw data')
ax.plot(Lidar_Data['time'][:]/3600,WDir_Filt,lw=7,color='g',label='Filtered data (Lowpass)')
ax.plot(Lidar_Data['time'][:1]/3600,(Grad_WDir[:1]),lw=5,color='g',ls='--',label=r'$\frac{d}{dt}$(Filtered data)')
plt.legend(loc='upper left')
ax.set_ylim(0,360)
ax.tick_params(labelsize=40)
ax.set_xlabel('Time - UTC [hour]',fontsize=40)
ax.set_ylabel('Wind Direction [degrees]',fontsize=40)
ax1 = ax.twinx() 
ax1.plot(Lidar_Data['time'][:]/3600,(Grad_WDir),lw=5,color='g',ls='--',label=r'$\frac{d}{dt}$ Filtered data')
ax1.plot(Lidar_Data['time'][WDir_Max_index]/3600,np.min(Grad_WDir),lw=5,marker='o',markersize=20,color='r')
ax1.set_ylim(-1,1)
ax1.set_ylabel(r'$d~/~dt$ (Wind Direction) [degrees/second]',fontsize=40)
Title = "Wind Direction and Direction Shift Rate: " + Date
ax.set_title(Title)
if Plot_Show == True:
    plt.show()
else:
    plt.close()

## Water Vapor Figure

fig, ax = plt.subplots(figsize=(30,20))
rc('font',weight='normal',size=40) 
ax.plot(MWR_Data['time'][:]/3600,MWR_Data['vap'][:],lw=3,color='k',label='raw data')
#ax.plot(MWR_Data['time'][:]/3600,WV_Filt,lw=5,color='b',label='Filtered data (Lowpass)')
ax.plot(MWR_Data['time'][:1]/3600,(Grad_WV[:1]+np.nanmean(MWR_Data['vap'][:])),lw=3,color='b',ls='--',label=r'$\frac{d}{dt}$ Filtered data')
plt.legend(loc='upper left')
ax.set_xlabel('Time - UTC [hour]',fontsize=40)
ax.set_ylabel('Water Vapor LoS Depth [cm]',fontsize=40)
ax1 = ax.twinx() 
ax1.plot(MWR_Data['time'][:]/3600,(Grad_WV),lw=3,color='b',ls='--',label=r'$\frac{d}{dt}$ Filtered data')
ax1.plot(MWR_Data['time'][WV_Max_index]/3600,np.max(Grad_WV),marker='o',markersize=20,color='r')
ax1.set_ylabel(r'$d~/~dt$ (Water Vapor LoS Depth) [cm/second]',fontsize=40)
Title = "Water Vapor LoS Depth and Gradients Rate: " + Date

ax.set_title(Title)
if Plot_Show == True:
    plt.show()
else:
    plt.close()


## Temperature Figure


fig, ax = plt.subplots(figsize=(30,20))
rc('font',weight='normal',size=40) 
ax.plot(MWR_Data['time'][:]/3600,MWR_Data['tkair'][:],lw=3,color='k',label='raw data')
ax.plot(MWR_Data['time'][:]/3600,Temp_Filt,lw=5,color='r',label='Filtered data (Lowpass)')
ax.plot(MWR_Data['time'][:1]/3600,(Grad_Temp[:1]+np.nanmean(MWR_Data['tkair'][:])),lw=3,color='r',ls='--',label=r'$\frac{d}{dt}$ Filtered data')## This is here for legend purposes
plt.legend(loc='upper left')
ax1 = ax.twinx() 
ax1.plot(MWR_Data['time'][:]/3600,(Grad_Temp),lw=3,color='r',ls='--',label=r'$\frac{d}{dt}$ Filtered data')
ax1.plot(MWR_Data['time'][Temp_Max_index]/3600,np.max(Grad_Temp),marker='o',markersize=20,color='r')

if Plot_Show == True:
    plt.show()
else:
    plt.close()

## Aerosol Figure

fig, ax = plt.subplots(figsize=(30,20))
rc('font',weight='normal',size=40) 
ax.plot(Aero_Data['time'][:]/3600,Aero_Data['concentration'][:],lw=3,color='k',label='raw data')
ax.plot(Aero_Data['time'][:]/3600,Aero_Filt,lw=5,color='orange',label='Filtered data (Lowpass)')
ax.plot(Avg_Aero_T,Average_Aero,lw=5,color='orange',ls=':',label='1-hour Moving Average')
ax.plot(Aero_Data['time'][:1]/3600,(Grad_Aero[:1]),lw=3,color='orange',ls='--',label=r'$\frac{d}{dt}$ (Filtered data)')## This is here for legend purposes
ax.axvline(MWR_Data['time'][WV_Max_index]/3600,lw=4,color='lightsteelblue',ls='--',zorder=0)
ax.axvline(Lidar_Data['time'][WDir_Max_index]/3600,lw=4,color='darkseagreen',ls='--',zorder=0)
ax.text(Lidar_Data['time'][WDir_Max_index]/3600,2,'Sb at %1.4fh'%(MWR_Data['time'][WV_Max_index]/3600))
ax.set_ylim(0,4*1e4)
#ax.set_ylim()
plt.legend(loc='upper left')
ax1 = ax.twinx() 
ax1.plot(Aero_Data['time'][:]/3600,(Grad_Aero),lw=3,color='orange',ls='--',label=r'$\frac{d Filtered data}{dt}$ ')
ax1.plot(Aero_Data['time'][Aero_Max_index]/3600,np.max(Grad_Aero),marker='o',markersize=20,color='r')
ax1.plot(Aero_Data['time'][Aero_Max_index]/3600,np.max(Grad_Aero),marker='o',markersize=20,color='r')
#plt.legend(loc='upper right')
ax.set_title(Date)
if Plot_Show == True:
    plt.show()
else:
    plt.show()

## Agregate Figure

#plt.figure(figsize=(30,20))
fig, ax4 = plt.subplots(figsize=(35,25))
rc('font',weight='normal',size=40)  

#ax4 = ax1.twinx() 
ax4.plot(Aero_Data['time'][:]/3600,Aero_Filt,color='orange',lw=2) #total_N_conc
#ax4.set_ylim(0.0,100.0)
ax4.set_ylabel('Aerosol Concentration',fontsize=40)  
ax4.yaxis.label.set_color('k')
ax4.axvline(MWR_Data['time'][WV_Max_index]/3600,lw=4,color='lightsteelblue',ls='--',zorder=0)
ax4.axvline(Lidar_Data['time'][WDir_Max_index]/3600,lw=4,color='darkseagreen',ls='--',zorder=0)

ax1 = ax4.twinx() 
ax1.plot(MWR_Data['time'][:]/3600,Temp_Filt,color='red',lw=2)
ax1.set_ylabel('Temperature')
ax1.yaxis.label.set_color('red')

ax1.set_ylim(290,315)

ax2 = ax4.twinx() 
ax2.plot(MWR_Data['time'][:]/3600,WV_Filt,lw=2,color='royalblue')
ax2.set_ylabel('Water Vapor')
ax2.yaxis.label.set_color('royalblue')
ax2.spines.right.set_position(("axes", 1.1))
ax2.set_ylim(2.5,5.5)

ax3 = ax4.twinx() 
ax3.plot(Lidar_Data['time'][:]/3600,WDir_Filt,lw=3,color='g')
ax3.set_ylim(0.0,360.0)

ax3.spines.right.set_position(("axes", 1.2))
ax3.set_ylabel('Wind Direction')
ax3.yaxis.label.set_color('g')



plt.title(Date)
if Plot_Show == True:
    plt.show()
else:
    plt.close()






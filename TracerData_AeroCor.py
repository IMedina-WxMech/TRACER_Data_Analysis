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

def SB_Passage_Time (WV,Wdir,WV_Time,WDir_Time):
    
    ddt_WV = First_Derivative(WV, len(WV), WV_Time[1]-WV_Time[0])
    ddt_WDir =  First_Derivative(Wdir, len(Wdir), WDir_Time[1]-WDir_Time[0])
    
    WV_index = np.where(ddt_WV[st_index_mwr:]==np.max(ddt_WV[st_index_mwr:]))[0][0]+st_index_mwr #Only daytime 
    WDir_index = np.where(ddt_WDir[st_index_lidar:]==np.min(ddt_WDir[st_index_lidar:]))[0][0]+st_index_lidar
    SB_time = np.nanmean(((WV_Time[WV_index]/3600)+(WDir_Time[WDir_index]/3600))/2)*3600
    #print(WV_Time[WV_index]/3600,WDir_Time[WDir_index]/3600)
    return SB_time

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

Plot_Show = False
Plot_Save = False 
## Wind data also has a height component
## Here we define what height we want (in m)
## Min height = 90m, Max height = 4325m
Wind_Height = 250

## There is a package that does this, but I will do it manually for now
st_h =  6 #sunrise hour
st_m = 25 #sunrise minute
sunrise_time = ((6+5)*3600)+(st_m*60) ## sunrise time in seconds


## Data pathing

Main_Path = '/Users/isaac/Desktop/Energy_Meteorology/Coding/SB_Data/*'
sub_folders = sorted( glob(Main_Path))

Aggregate_Aero_Data = np.zeros((3,3,2)) # number of classifications,time periods, (ave,st.error)

# Subfolders are our classification identifiers: Clear, SB, consecutive
for i in range (0,len(sub_folders)):
    #%%
    
    Data_Path = sub_folders[i] 
    print("Working on:" ,Data_Path)
    
    ## Mulitple paths for AOD data depending on which one you want to use
    #AOD_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houaodbe5chM1.c1/*'))
    AOD_Data_Path = sorted( glob(Data_Path+'/houaoscpcfM1.b1/*'))
    #AOD_Data_Path = sorted( glob('/Users/isaac/Desktop/Energy_Meteorology/Coding/Data/houaoscpcuM1.b1/*'))
    
    MWR_Data_Path = sorted( glob(Data_Path+'/houmwrlosM1.b1/*'))
    Lidar_Data_path = sorted( glob(Data_Path+'/houdlprofwind4newsM1.c1/*'))
    
    Num_files = len(Lidar_Data_path) ## Number of Files for looping over entire dataset
    Aero_averages = np.zeros((3,(Num_files) ))
    for Day_Index in range (0,Num_files):
        
        ## Data Read in 
        Aero_Data = Dataset(AOD_Data_Path[Day_Index])
        MWR_Data = Dataset(MWR_Data_Path[Day_Index])
        Lidar_Data = Dataset(Lidar_Data_path[Day_Index])
        
        st_index_mwr = np.where(np.abs(MWR_Data['time'][:]-sunrise_time) == np.min( np.abs(MWR_Data['time'][:]-sunrise_time)   ))[0][0]
        st_index_lidar = np.where(np.abs(Lidar_Data['time'][:]-sunrise_time) == np.min( np.abs(Lidar_Data['time'][:]-sunrise_time)   ))[0][0]
        
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
        Grad_Aero = First_Derivative(Aero_Filt, len(Aero_Filt), Aero_Data['time'][1]-Aero_Data['time'][0])
        ## Finding the max/min (depending on variable) of the gradient to find SB passage time  
        Aero_Max_index = np.where(Grad_Aero==np.max(Grad_Aero))[0][0]
        ##Time array in dataset is in seconds, so 3600 = 1 hour
        window = int(3600)
        Average_Aero,Avg_Aero_T = windowed_average(Aero_Filt,Aero_Data['time'][:]/3600,window)
        ##Time of SB passage (for SB passage days)
        if i ==0:
            SB_t = 17.5*3600 #Average SB time
        else:
            SB_t = SB_Passage_Time(WV_Filt,WDir_Filt,MWR_Data['time'],Lidar_Data['time'])
        
        sb_index = np.where(np.abs(Aero_Data['time'][:]-SB_t) == np.min( np.abs(Aero_Data['time'][:]-SB_t)   ))[0][0]
        st_aero_index = np.where(np.abs(Aero_Data['time'][:]-MWR_Data['time'][st_index_mwr]) == np.min( np.abs(Aero_Data['time'][:]-MWR_Data['time'][st_index_mwr])   ))[0][0]
        #print(SB_t/3600,Aero_Data['time'][sb_index]/3600)
        
        Ave_areo_Night = np.nanmean(Aero_Filt[0:st_aero_index])
        Ave_aero_Morning = np.nanmean(Aero_Filt[st_aero_index:sb_index])
        Ave_aero_postSB = np.nanmean(Aero_Filt[sb_index:])
        #print(Ave_areo_Night,Ave_aero_Morning,Ave_aero_postSB)
        
        Aero_averages[0,Day_Index] = Ave_areo_Night
        Aero_averages[1,Day_Index] = Ave_aero_Morning
        Aero_averages[2,Day_Index] = Ave_aero_postSB
    
    Aggregate_Aero_Data[i,0,0] = np.nanmean(Aero_averages[0,:])
    Aggregate_Aero_Data[i,0,1] = np.std(Aero_averages[0,:]) / np.sqrt(Num_files)
    
    Aggregate_Aero_Data[i,1,0] = np.nanmean(Aero_averages[1,:])
    Aggregate_Aero_Data[i,1,1] = np.std(Aero_averages[1,:]) / np.sqrt(Num_files)
        
    Aggregate_Aero_Data[i,2,0] = np.nanmean(Aero_averages[2,:])
    Aggregate_Aero_Data[i,2,1] = np.std(Aero_averages[2,:]) / np.sqrt(Num_files)
    
   
#%%

# =============================================================================
"""
Section IV: Data Graphing

"""
plt.figure(figsize=(15,10))
rc('font',weight='normal',size=20)
plt.grid(zorder=10)
plt.bar(0.0,Aggregate_Aero_Data[0,0,0],width=0.5,color='#BBBBBB',label='Clear')
plt.bar(0.5,Aggregate_Aero_Data[1,0,0],width=0.5,color='#66CCEE',label='Sea-Breeze')
plt.bar(1.0,Aggregate_Aero_Data[2,0,0],width=0.5,color='#4477AA',label='Consecutive')

plt.errorbar(0.0,Aggregate_Aero_Data[0,0,0],yerr=Aggregate_Aero_Data[0,0,1],ecolor='k',capsize=5)
plt.errorbar(0.5,Aggregate_Aero_Data[1,0,0],yerr=Aggregate_Aero_Data[1,0,1],ecolor='k',capsize=5)
plt.errorbar(1.0,Aggregate_Aero_Data[2,0,0],yerr=Aggregate_Aero_Data[2,0,1],ecolor='k',capsize=5)



plt.bar(2.0,Aggregate_Aero_Data[0,1,0],width=0.5,color='#BBBBBB')
plt.bar(2.5,Aggregate_Aero_Data[1,1,0],width=0.5,color='#66CCEE')
plt.bar(3.0,Aggregate_Aero_Data[2,1,0],width=0.5,color='#4477AA')

plt.errorbar(2.0,Aggregate_Aero_Data[0,1,0],yerr=Aggregate_Aero_Data[0,1,1],ecolor='k',capsize=5)
plt.errorbar(2.5,Aggregate_Aero_Data[1,1,0],yerr=Aggregate_Aero_Data[1,1,1],ecolor='k',capsize=5)
plt.errorbar(3.0,Aggregate_Aero_Data[2,1,0],yerr=Aggregate_Aero_Data[2,1,1],ecolor='k',capsize=5)




plt.bar(4.0,Aggregate_Aero_Data[0,2,0],width=0.5,color='#BBBBBB')
plt.bar(4.5,Aggregate_Aero_Data[1,2,0],width=0.5,color='#66CCEE')
plt.bar(5.0,Aggregate_Aero_Data[2,2,0],width=0.5,color='#4477AA')

plt.errorbar(4.0,Aggregate_Aero_Data[0,2,0],yerr=Aggregate_Aero_Data[0,2,1],ecolor='k',capsize=5)
plt.errorbar(4.5,Aggregate_Aero_Data[1,2,0],yerr=Aggregate_Aero_Data[1,2,1],ecolor='k',capsize=5)
plt.errorbar(5.0,Aggregate_Aero_Data[2,2,0],yerr=Aggregate_Aero_Data[2,2,1],ecolor='k',capsize=5)




labels= ['Pre-Dawn','Pre-SB','Post-SB']
tick = [0.5,2.5,4.5]
plt.xticks(tick,labels)
plt.ylabel(r'Aerosol Content [$1/cm^3$]')#,fontsize=40)
plt.legend(loc='upper left')
plt.xlim(-0.5,5.5)
plt.title("Average Aerosol Content per Time of Day and Case Type")
plt.show()





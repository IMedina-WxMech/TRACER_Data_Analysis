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
      time: 1-d time array(for graphing purposes)
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

def hour_average(data, time):
    """
    Perform a hourly average
    
    Args:
      data: 1-d dataset
      time: 1-d time array(for graphing purposes)
      
    Returns:
    	result: averaged data
    """
    result = []
    for h in range (0,24):
        
        hour =h
        hour_next = h+1
        temp = []
        for t_int in range (0,len(time)):
            t = time[t_int]
            if hour<t<hour_next:
                temp.append(data[t_int])
        
        result.append(np.average(temp))
    return result

    
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
    SB_time = ((np.nanmean(((WV_Time[WV_index]/3600)+(WDir_Time[WDir_index]/3600))/2)-.25)*3600)
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
Ag_Hour_Aero_Data  = np.zeros((3,24,2)) # number of classifications,time periods, (ave,st.error)
Ag_Hour_BL_Data  = np.zeros((3,24,2))
Ag_Hour_WW_Data  = np.zeros((3,24,2))
Sb_Time_ave = []

# Subfolders are our classification identifiers: Clear, SB, consecutive
for i in range (0,len(sub_folders)):

    
    Data_Path = sub_folders[i] 
    print("Working on:" ,Data_Path)
    
    ## Mulitple paths for AOD data depending on which one you want to use
    Dist_Data_Path = sorted( glob(Data_Path+'/houaosuhsasM1.b1/*'))
    #AOD_Data_Path = sorted( glob(Data_Path+'/houaoscpcfM1.b1/*'))
    AOD_Data_Path = sorted( glob(Data_Path+'/houaoscpcuM1.b1/*'))
    
    MWR_Data_Path = sorted( glob(Data_Path+'/houmwrlosM1.b1/*'))
    Lidar_Data_path = sorted( glob(Data_Path+'/houdlprofwind4newsM1.c1/*'))
    
    ABL_Data_Path = sorted( glob(Data_Path+'/BLData/*'))

    Num_files = len(Lidar_Data_path) ## Number of Files for looping over entire dataset
    Aero_averages = np.zeros((3,(Num_files) ))
    Aero_errors = np.zeros((3,(Num_files) ))
    Hour_ave = np.zeros((24,Num_files))
    BL_Hour_ave = np.zeros((24,Num_files))
    WW_Hour_ave = np.zeros((24,Num_files))
    SB_ave = []
    
    for Day_Index in range (0,Num_files):
        
        ## Data Read in 
        Aero_Data = Dataset(AOD_Data_Path[Day_Index])
        Dist_Data = Dataset(Dist_Data_Path[Day_Index])
        MWR_Data = Dataset(MWR_Data_Path[Day_Index])
        Lidar_Data = Dataset(Lidar_Data_path[Day_Index])
        BL_Data = Dataset(ABL_Data_Path[Day_Index])
        
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
        #Average_Aero = hour_average(Aero_Filt,Aero_Data['time'][:]/3600)
        Average_Aero,Avg_Aero_T = windowed_average(Aero_Filt,Aero_Data['time'][:]/3600,window)
        Hour_ave[:,Day_Index] = Average_Aero
        BL_Hour_ave[:,Day_Index] = BL_Data['pblh_1hr'][:]
        ##Time of SB passage (for SB passage days)
        if i ==0:
            SB_t = 17.5*3600 #Average SB time
        else:
            SB_t = SB_Passage_Time(WV_Filt,WDir_Filt,MWR_Data['time'],Lidar_Data['time'])
        SB_ave.append(SB_t/3600)
        sb_index = np.where(np.abs(Aero_Data['time'][:]-SB_t) == np.min( np.abs(Aero_Data['time'][:]-SB_t)   ))[0][0]
        st_aero_index = np.where(np.abs(Aero_Data['time'][:]-MWR_Data['time'][st_index_mwr]) == np.min( np.abs(Aero_Data['time'][:]-MWR_Data['time'][st_index_mwr])   ))[0][0]
        
        Ave_areo_Night = np.nanmean(Aero_Filt[0:st_aero_index])
        Ave_aero_Morning = np.nanmean(Aero_Filt[st_aero_index:sb_index])
        Ave_aero_postSB = np.nanmean(Aero_Filt[sb_index:])
       
        Aero_averages[0,Day_Index] = Ave_areo_Night
        Aero_averages[1,Day_Index] = Ave_aero_Morning
        Aero_averages[2,Day_Index] = Ave_aero_postSB
         
        Aero_errors[0,Day_Index] = np.std(Aero_Filt[0:st_aero_index])/np.sqrt(st_aero_index)
        Aero_errors[1,Day_Index] = np.std(Aero_Filt[st_aero_index:sb_index])/np.sqrt(sb_index-st_aero_index)
        Aero_errors[2,Day_Index] = np.std(Aero_Filt[sb_index:])/np.sqrt(len(Aero_Filt[:])-sb_index)
        
        ## ======================
        ## Here is where we find the veritcal integral of vertical velocity variance
        temp_WW_array = []
        for t_int in range (0,len(BL_Data['t_sig'][:])):
            time  = BL_Data['t_sig'][t_int]*60
            t_index = np.where(np.abs(BL_Data['t'][:]-time) == np.min( np.abs(BL_Data['t'][:]-time)   ))[0][0]
            PBLH = BL_Data['pblh'][t_index]
            try:
                bl_index = np.where(np.abs(BL_Data['z'][:]-PBLH) == np.min( np.abs(BL_Data['z'][:]-PBLH)   ))[0][0]
               # print(bl_index)
                ww_sum = np.nansum(BL_Data['sig_w'][t_int,:bl_index])
                temp_WW_array.append(ww_sum)
            except:
                temp_WW_array.append(0)     
        Average_WW = hour_average(temp_WW_array,BL_Data['t_sig'][:]/60)#,window)       
        WW_Hour_ave[:,Day_Index] = Average_WW
        ## ======================
    
    ## Now we save all that data to arrays outside the loop
    Aggregate_Aero_Data[i,0,0] = np.nanmean(Aero_averages[0,:])
    Aggregate_Aero_Data[i,0,1] = np.std(Aero_averages[0,:]) / np.sqrt(Num_files)
    
    Aggregate_Aero_Data[i,1,0] = np.nanmean(Aero_averages[1,:])
    Aggregate_Aero_Data[i,1,1] = np.std(Aero_averages[1,:]) / np.sqrt(Num_files)
        
    Aggregate_Aero_Data[i,2,0] = np.nanmean(Aero_averages[2,:])
    Aggregate_Aero_Data[i,2,1] = np.std(Aero_averages[2,:]) / np.sqrt(Num_files)
    Sb_Time_ave.append(np.nanmean(SB_ave))
    for h in range (0,24):
        Ag_Hour_Aero_Data[i,h,0] = np.nanmean(Hour_ave[h,:])
        Ag_Hour_Aero_Data[i,h,1] = np.std(Hour_ave[h,:]) / np.sqrt(Num_files)
        
        Ag_Hour_BL_Data[i,h,0] = np.nanmean(BL_Hour_ave[h,:])
        Ag_Hour_BL_Data[i,h,1] = np.std(BL_Hour_ave[h,:]) / np.sqrt(Num_files)
        
        Ag_Hour_WW_Data[i,h,0] = np.nanmean(WW_Hour_ave[h,:])
        Ag_Hour_WW_Data[i,h,1] = np.std(WW_Hour_ave[h,:]) / np.sqrt(Num_files)
        
    # =============================================================================
    """
    Section IV: Data Graphing
    
    """
    
    ## Plots for consecutive case days:
    if i ==2:
        
        ##Line plot
        cmap = plt.get_cmap('Oranges', 7)   
        fig, ax = plt.subplots(figsize=(30,20))
        rc('font',weight='normal',size=35)
        plt.grid(zorder=10)  
        plt.axvline(np.nanmean(SB_ave[:6]),color='grey',lw=3,ls='--',zorder=10)
        plt.text(np.nanmean(SB_ave[:])+0.25,1000,r"Average SB Passage",fontsize=35)
        plt.text(np.nanmean(SB_ave[:])+0.27,100,r"Time: %1.2f Hours"%(np.nanmean(SB_ave[:])),fontsize=35)
        plt.axvline(MWR_Data['time'][st_index_mwr]/3600,color='lightcoral',lw=3,ls='--',zorder=10)
        plt.text((MWR_Data['time'][st_index_mwr]/3600)+0.25,1000,r"Sunrise Time:",fontsize=35)
        plt.text((MWR_Data['time'][st_index_mwr]/3600)+0.27,100,r"%1.2f Hours"%(MWR_Data['time'][st_index_mwr]/3600),fontsize=35)
        
        
        for j in range (0,6):
            plt.plot(Hour_ave[:,j],label='SB Day %1.0f; Average:%1.3f'%(j+1,np.nanmean(Hour_ave[:,j])),color=cmap(j+1),lw=5)
      
        tick = np.arange(0,24,1)
        plt.xticks(tick,tick,fontsize=35)       
        plt.ylabel(r'Aerosol Content [$1/cm^3$]',fontsize=35)        
        plt.legend(loc='upper left')      
        plt.xlabel(r'Hour [UTC]',fontsize=35)      
        plt.xlim(-0.25,23.75)       
        plt.ylim(0,30000)
        plt.title("Aerosol Content Per Hour for July 3-8")
        plt.show()
       
        ##Bar plot
        fig, ax = plt.subplots(figsize=(30,20))
        rc('font',weight='normal',size=35)
        
        plt.bar(0.0 +(1.75*0),Aero_averages[0,0],width=0.25,color=cmap(1),label='SB Day 1')
        plt.bar(0.25+(1.75*0),Aero_averages[0,1],width=0.25,color=cmap(2),label='SB Day 2')
        plt.bar(0.5 +(1.75*0),Aero_averages[0,2],width=0.25,color=cmap(3),label='SB Day 3')
        plt.bar(0.75+(1.75*0),Aero_averages[0,3],width=0.25,color=cmap(4),label='SB Day 4')
        plt.bar(1.0 +(1.75*0),Aero_averages[0,4],width=0.25,color=cmap(5),label='SB Day 5')
        plt.bar(1.25+(1.75*0),Aero_averages[0,5],width=0.25,color=cmap(6),label='SB Day 6')
               
        for j in range (0,3):
            for index in range (0,6):
                plt.bar((0.25*index) +(1.75*j),Aero_averages[j,index],width=0.25,color=cmap(index+1))#,label='SB Day 1')
                #plt.errorbar((0.25*index) +(1.75*j),Aero_averages[j,index],yerr=np.std(Aero_Filt[j,:]) / np.sqrt(Num_files),ecolor='#BBBBBB',capsize=6,lw=3)
             
        plt.legend(loc='upper left',fontsize=35)
        labels= ['Pre-Dawn','Pre-SB','Post-SB']
        tick = [0.6,2.3,4.1]
        plt.xticks(tick,labels,fontsize=42)
        
        plt.ylabel(r'Aerosol Content [$1/cm^3$]',fontsize=35)
      
        #plt.xlim(-0.25,2.75)
        plt.ylim(0,20000)
      
        plt.title("Aerosol Content Per Time of Day for July 3-8",fontsize=45)
       
        plt.show()
        
        
        
        
#%%
# =============================================================================
"""
Section IV: Data Graphing (Cont.)

"""

## Primary Bar plot for aerosol concentrations
plt.figure(figsize=(15,10))
rc('font',weight='normal',size=20)
#plt.grid(zorder=10)
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
plt.xticks(tick,labels,fontsize=28)
plt.ylabel(r'Aerosol Content [$1/cm^3$]')#,fontsize=40)
plt.legend(loc='upper left')
plt.xlim(-0.5,5.5)
plt.ylim(0,20000)
plt.title("Average Aerosol Content per Time of Day and Case Type",fontsize=28)
plt.show()




#%%
## Primary line plot for aerosol

fig, ax = plt.subplots(figsize=(30,20))
rc('font',weight='normal',size=35)
plt.grid(zorder=10)


plt.axvline(Sb_Time_ave[1],color='grey',lw=3,ls='--',zorder=10)
plt.text(Sb_Time_ave[1]+0.25,2200,r"Average Singular")
plt.text(Sb_Time_ave[1]+0.25,1600,r"SB Time: %1.2f Hours"%(Sb_Time_ave[1]))
plt.axvline(Sb_Time_ave[2],color='grey',lw=3,ls='--',zorder=10)
plt.text(Sb_Time_ave[2]+0.25,800,r"Average Consecutive")
plt.text(Sb_Time_ave[2]+0.25,200,r"SB Time: %1.2f Hours"%(Sb_Time_ave[2]))


plt.plot(Ag_Hour_Aero_Data[0,:,0],'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
plt.plot(Ag_Hour_Aero_Data[1,:,0],'o',ls='-',color='#66CCEE',lw=6,label='Sea-Breeze',markersize=10)
plt.plot(Ag_Hour_Aero_Data[2,:,0],'o',ls='-',color='#4477AA',lw=6,label='Consecutive',markersize=10)
plt.axvline(MWR_Data['time'][st_index_mwr]/3600,color='lightcoral',lw=3,ls='--',zorder=10)

plt.text((MWR_Data['time'][st_index_mwr]/3600)+0.25,1000,r"Sunrise Time:")#,fontsize=35)
plt.text((MWR_Data['time'][st_index_mwr]/3600)+0.27,100,r"%1.2f Hours"%(MWR_Data['time'][st_index_mwr]/3600))#,fontsize=35)

for h in range (0,24):
    plt.errorbar(h,Ag_Hour_Aero_Data[1,h,0],yerr=Ag_Hour_Aero_Data[1,h,1],ecolor='#66CCEE',capsize=6,lw=3)
    plt.errorbar(h,Ag_Hour_Aero_Data[0,h,0],yerr=Ag_Hour_Aero_Data[0,h,1],ecolor='#BBBBBB',capsize=6,lw=3)
    plt.errorbar(h,Ag_Hour_Aero_Data[2,h,0],yerr=Ag_Hour_Aero_Data[2,h,1],ecolor='#4477AA',capsize=6,lw=3)


#labels= ['Pre-Dawn','Pre-SB','Post-SB']
tick = np.arange(0,24,1)
plt.xticks(tick,tick)
plt.ylabel(r'Aerosol Content [$1/cm^3$]',fontsize=40)
plt.legend(loc='upper left',fontsize=40)
plt.xlabel(r'Hour[UTC]')

plt.xlim(-0.25,23.75)
plt.ylim(0,20000)


#f.title("Average Aerosol Content per Hour of Day and Case Type")
#plt.tight_layout(pad=1.4,w_pad=0.5)
plt.title("Average Aerosol Content per Hour of Day and Case Type",fontsize=50)
plt.show()
#%%
## Seperated line plot for aerosol concentrations

#f, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(25,10),layout="constrained")
f, (ax1, ax2) = plt.subplots(2,1, sharex=True,figsize=(15,20),layout="constrained")
#plt.figure()
rc('font',weight='normal',size=24)
ax1.grid(zorder=10)
ax2.grid(zorder=10)


ax1.axvline(Sb_Time_ave[1],color='grey',lw=3,ls='--',zorder=10)
ax1.text(Sb_Time_ave[1]+0.25,800,r"Average SB Passage")
ax1.text(Sb_Time_ave[1]+0.25,200,r"Time: %1.2f Hours"%(Sb_Time_ave[1]))
ax2.axvline(Sb_Time_ave[2],color='grey',lw=3,ls='--',zorder=10)
ax2.text(Sb_Time_ave[2]+0.25,800,r"Average SB Passage")
ax2.text(Sb_Time_ave[2]+0.25,200,r"Time: %1.2f Hours"%(Sb_Time_ave[2]))


ax1.plot(Ag_Hour_Aero_Data[0,:,0],'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
ax2.plot(Ag_Hour_Aero_Data[0,:,0],'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
ax1.plot(Ag_Hour_Aero_Data[1,:,0],'o',ls='-',color='#66CCEE',lw=6,label='Sea-Breeze',markersize=10)
ax2.plot(Ag_Hour_Aero_Data[2,:,0],'o',ls='-',color='#4477AA',lw=6,label='Consecutive',markersize=10)
ax1.axvline(MWR_Data['time'][st_index_mwr]/3600,color='lightcoral',lw=3,ls='--',zorder=10)

ax1.text((MWR_Data['time'][st_index_mwr]/3600)+0.25,1000,r"Sunrise Time:")#,fontsize=35)
ax1.text((MWR_Data['time'][st_index_mwr]/3600)+0.27,100,r"%1.2f Hours"%(MWR_Data['time'][st_index_mwr]/3600))#,fontsize=35)
ax2.axvline(MWR_Data['time'][st_index_mwr]/3600,color='lightcoral',lw=3,ls='--',zorder=10)
ax2.text((MWR_Data['time'][st_index_mwr]/3600)+0.25,1000,r"Sunrise Time:")#,fontsize=35)
ax2.text((MWR_Data['time'][st_index_mwr]/3600)+0.27,100,r"%1.2f Hours"%(MWR_Data['time'][st_index_mwr]/3600))#,fontsize=35)

for h in range (0,24):
    ax1.errorbar(h,Ag_Hour_Aero_Data[0,h,0],yerr=Ag_Hour_Aero_Data[0,h,1],ecolor='#BBBBBB',capsize=6,lw=3)
    ax1.errorbar(h,Ag_Hour_Aero_Data[1,h,0],yerr=Ag_Hour_Aero_Data[1,h,1],ecolor='#66CCEE',capsize=6,lw=3)
    ax2.errorbar(h,Ag_Hour_Aero_Data[0,h,0],yerr=Ag_Hour_Aero_Data[0,h,1],ecolor='#BBBBBB',capsize=6,lw=3)
    ax2.errorbar(h,Ag_Hour_Aero_Data[2,h,0],yerr=Ag_Hour_Aero_Data[2,h,1],ecolor='#4477AA',capsize=6,lw=3)


#labels= ['Pre-Dawn','Pre-SB','Post-SB']
tick = np.arange(0,24,1)
ax1.set_xticks(tick,tick)
ax2.set_xticks(tick,tick)
ax1.set_ylabel(r'Aerosol Content [$1/cm^3$]',fontsize=30)
ax2.set_ylabel(r'Aerosol Content [$1/cm^3$]',fontsize=30)
#ax1.set_ylabel(r'Aerosol Content [$1/cm^3$]')#,fontsize=40)
ax1.legend(loc='upper left',fontsize=25)
ax2.legend(loc='upper left',fontsize=25)

#ax1.set_xlabel(r'Hour [UTC]')
ax2.set_xlabel(r'Hour[UTC]')

ax1.set_xlim(-0.25,23.75)
ax2.set_xlim(-0.25,23.75)
ax1.set_ylim(0,20000)
ax2.set_ylim(0,20000)


#f.title("Average Aerosol Content per Hour of Day and Case Type")
#plt.tight_layout(pad=1.4,w_pad=0.5)
f.suptitle("Average Aerosol Content per Hour of Day and Case Type",fontsize=30)
plt.show()



#%%
## Primary line plot for ABL structure

#f, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=(25,10),layout="constrained")
f, (ax1, ax2) = plt.subplots(2,1, sharex=True,figsize=(15,22),layout="constrained")
#plt.figure()
rc('font',weight='normal',size=22)
ax1.grid(zorder=10)
ax2.grid(zorder=10)


ax1.axvline(Sb_Time_ave[1],color='grey',lw=3,ls='--',zorder=10)
ax1.text(Sb_Time_ave[1]+0.25,360,r"Average Singular")
ax1.text(Sb_Time_ave[1]+0.25,280,r"SB Time: %1.2f Hours"%(Sb_Time_ave[1]))
ax1.axvline(Sb_Time_ave[2],color='grey',lw=3,ls='-.',zorder=10)
ax1.text(Sb_Time_ave[2]+0.25,180,r"Average Consecutive")
ax1.text(Sb_Time_ave[2]+0.25,110,r"SB Time: %1.2f Hours"%(Sb_Time_ave[2]))

ax2.axvline(Sb_Time_ave[1],color='grey',lw=3,ls='--',zorder=10)
ax2.text(Sb_Time_ave[1]+0.25,9,r"Average Singular")
ax2.text(Sb_Time_ave[1]+0.25,7,r"SB Time: %1.2f Hours"%(Sb_Time_ave[1]))
ax2.axvline(Sb_Time_ave[2],color='grey',lw=3,ls='-.',zorder=10)
ax2.text(Sb_Time_ave[2]+0.25,4,r"Average Consecutive")
ax2.text(Sb_Time_ave[2]+0.25,2,r"SB Time: %1.2f Hours"%(Sb_Time_ave[2]))





ax1.axvline(MWR_Data['time'][st_index_mwr]/3600,color='lightcoral',lw=3,ls='--',zorder=10)
ax1.text((MWR_Data['time'][st_index_mwr]/3600)+0.25,180,r"Sunrise Time:")#,fontsize=35)
ax1.text((MWR_Data['time'][st_index_mwr]/3600)+0.27,110,r"%1.2f Hours"%(MWR_Data['time'][st_index_mwr]/3600))#,fontsize=35)
ax2.axvline(MWR_Data['time'][st_index_mwr]/3600,color='lightcoral',lw=3,ls='--',zorder=10)
ax2.text((MWR_Data['time'][st_index_mwr]/3600)+0.25,4,r"Sunrise Time:")#,fontsize=35)
ax2.text((MWR_Data['time'][st_index_mwr]/3600)+0.27,2,r"%1.2f Hours"%(MWR_Data['time'][st_index_mwr]/3600))#,fontsize=35)



#ax1.plot(Ag_Hour_BL_Data[0,:,0]*1000,'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
ax1.plot(Ag_Hour_BL_Data[0,:,0]*1000,'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
ax1.plot(Ag_Hour_BL_Data[1,:,0]*1000,'o',ls='-',color='#66CCEE',lw=6,label='Sea-Breeze',markersize=10)
ax1.plot(Ag_Hour_BL_Data[2,:,0]*1000,'o',ls='-',color='#4477AA',lw=6,label='Consecutive',markersize=10)

for h in range (0,24):
    ax1.errorbar(h,Ag_Hour_BL_Data[0,h,0]*1000,yerr=Ag_Hour_BL_Data[0,h,1]*1000,ecolor='#BBBBBB',capsize=6,lw=3)
    ax1.errorbar(h,Ag_Hour_BL_Data[1,h,0]*1000,yerr=Ag_Hour_BL_Data[1,h,1]*1000,ecolor='#66CCEE',capsize=6,lw=3)
    ax1.errorbar(h,Ag_Hour_BL_Data[0,h,0]*1000,yerr=Ag_Hour_BL_Data[0,h,1]*1000,ecolor='#BBBBBB',capsize=6,lw=3)
    ax1.errorbar(h,Ag_Hour_BL_Data[2,h,0]*1000,yerr=Ag_Hour_BL_Data[2,h,1]*1000,ecolor='#4477AA',capsize=6,lw=3)

#ax1.plot(Ag_Hour_WW_Data[0,:,0],'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
ax2.plot(Ag_Hour_WW_Data[0,:,0],'o',ls='-',color='#BBBBBB',lw=6,label='Clear',markersize=10)
ax2.plot(Ag_Hour_WW_Data[1,:,0],'o',ls='-',color='#66CCEE',lw=6,label='Sea-Breeze',markersize=10)
ax2.plot(Ag_Hour_WW_Data[2,:,0],'o',ls='-',color='#4477AA',lw=6,label='Consecutive',markersize=10)


for h in range (0,24):
    ax2.errorbar(h,Ag_Hour_WW_Data[0,h,0],yerr=Ag_Hour_WW_Data[0,h,1],ecolor='#BBBBBB',capsize=6,lw=3)
    ax2.errorbar(h,Ag_Hour_WW_Data[1,h,0],yerr=Ag_Hour_WW_Data[1,h,1],ecolor='#66CCEE',capsize=6,lw=3)
    ax2.errorbar(h,Ag_Hour_WW_Data[0,h,0],yerr=Ag_Hour_WW_Data[0,h,1],ecolor='#BBBBBB',capsize=6,lw=3)
    ax2.errorbar(h,Ag_Hour_WW_Data[2,h,0],yerr=Ag_Hour_WW_Data[2,h,1],ecolor='#4477AA',capsize=6,lw=3)


#labels= ['Pre-Dawn','Pre-SB','Post-SB']
tick = np.arange(0,24,1)
ax1.set_xticks(tick,tick)
ax2.set_xticks(tick,tick)
ax1.set_ylabel(r'Boundary Layer Height Estimate [m]',fontsize=27)
ax2.set_ylabel(r'Integrated Vertical Velocity Variance [$m^2/s$]',fontsize=27)    
#ax1.set_ylabel(r'Aerosol Content [$1/cm^3$]')#,fontsize=40)
ax1.legend(loc='upper left',fontsize=25)
ax2.legend(loc='upper left',fontsize=25)

#ax1.set_xlabel(r'Hour [UTC]')
ax2.set_xlabel(r'Hour[UTC]')

ax1.set_xlim(-0.25,23.75)
ax2.set_xlim(-0.25,23.75)
ax1.set_ylim(0,2000)
ax2.set_ylim(0,70)
ax1.set_title(r"Average Boundary Layer Height per Hour of Day and Case Type",fontsize=27)

#f.title("Average Aerosol Content per Hour of Day and Case Type")
#plt.tight_layout(pad=1.4,w_pad=0.5)
ax2.set_title("Integrated Vertical Velocity Variance per Hour of Day",fontsize=30)
plt.show()











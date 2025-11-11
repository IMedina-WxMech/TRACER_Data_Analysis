# TRACER_Data_Analysis
This repository is for the class project for Energy Meteorology (EN.530.637) class at Johns Hopkins University. The scope of this code and project is to analyze data from the Tracing Aerosol Convection Interactions Experiment (TRACER) to correlate sea breeze passage with aresol content. This consists of two different code files: TracerData_SBIdent and TracerData_AeroCor. These are used for identifing sea-breezes in the data, and correlating seabreeze passage to aersol content respectively. 

## Data
The data files was accessed through the ARM data viewer: https://www.arm.gov/data 

All the data for this project was collected in 2022 from the ARM AMF1 site in La Porte, TX during the TRACER campaign. Due to time constraints and data availability, the data was only analyzed for June 1st to July 31st. Zip files of these data are uploaded to this repository. The following datastreams were accessed:

* AOSCPCF: _AOS Fine Condensation Particle Counter, 10 nm detection limit_
* AOSCPCU: _AOS Ultrafine Condensation Particle Counter, 2.5 nm detection limit_
* DLPROFWIND4NEWS: _Doppler lidar wind VAP_
* MWRLOS: _Microwave Water Radiometer (MWR) water liq. & vapor along line of sight (LOS) path_

## SB_Identification: Functions and Scope
This script is used to identify when a sea-breeze was present in the data using gradients in wind speeds and water vapor (further explained in the attatched brief). The script features three defined functions. These are used for designing and applying a lowpass filter, performing a windowed average, and performing a first derivate with a central differencing scheme. Each defined function gives a detailed description of inputs and outputs in the script itself. 

This script is broken down into 4 main sections:
* Defined functions and packages
* Data pathing, read in, and pre processing
* Data processing
* Data Graphing

## Aero_Correlation: Functions and Scope
Once cases are identified as either a baseline(no SB),SB, or consecutive SB case the correlation code can be run. This code reads in segregated data files (in individual folders based on cassification) and read in. For the SB and consecutive cases the SB passage time is calculated as a function seperately for data anlysis. Here it loops over classifications (1-3), then over case files. From here the mean aerosol content and standard error for each time period and classification are saved outside the loop and graphed. 

This script is broken down into 4 main sections:
* Defined functions and packages
* Data pathing, read in, and pre processing
* Data processing
* Data Graphing

## Dependencies 
The python codes require the following packages (all are generally available in standard python distributions):
* Numpy
* Matplotlib
* netCDF4
* glob
* datetime
* scipy

## Authors
Isaac Medina, Johns Hopkins University, imedina2@jhu.edu

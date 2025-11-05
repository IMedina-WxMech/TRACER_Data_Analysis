# TRACER_Data_Analysis
This repository is for the class project for Energy Meteorology (EN.530.637) class at Johns Hopkins University. The scope of this code and project is to analyze data from the Tracing Aerosol Convection Interactions Experiment (TRACER) to correlate sea breeze passage with aresol content. 

## Data
The data files was accessed through the ARM data viewer: https://www.arm.gov/data 

All the data for this project was collected in 2022 from the ARM AMF1 site in La Porte, TX during the TRACER campaign. Due to time constraints and data availability, the data was only analyzed for June 1st to July 31st. The following datastreams were accessed:

* AOSCPCF: _AOS Fine Condensation Particle Counter, 10 nm detection limit_
* AOSCPCU: _AOS Ultrafine Condensation Particle Counter, 2.5 nm detection limit_
* AOSUHSAS: _AOS Ultrahigh Sensitivity Aerosol Spectrometer_
* DLPROFWIND4NEWS: _Doppler lidar wind VAP_
* MWRLOS: _Microwave Water Radiometer (MWR) water liq. & vapor along line of sight (LOS) path_

## Functions and Scope
The script features three defined functions. These are used for designing and applying a lowpass filter, performing a windowed average, and performing a first derivate with a central differencing scheme. Futher description of the design of these functions and thier implimentation can be found in the attatched report. Each defined function gives a detailed description of inputs and outputs in the script itself. 

The script is broken down into 4 main sections:
* Defined functions and packages
* Data pathing, read in, and pre processing
* Data processing
* Data Graphing

## Dependencies 
The python code requires the following packages (all are generally available in standard python distributions):
* Numpy
* Matplotlib
* netCDF4
* glob
* datetime
* scipy

## Authors
Isaac Medina, Johns Hopkins University, imedina2@jhu.edu

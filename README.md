# TRACER_Data_Analysis
This repository is for the class project for Energy Meteorology (EN.530.637) class at Johns Hopkins University. The scope of this code and project is to analyze data from the Tracing Aerosol Convection Interactions Experiment (TRACER) to correlate sea breeze passage with aresol content. The primary consists of two different code files: TracerData_SBIdent and TracerData_AeroCor. These are used for identifing sea-breezes in the data, and correlating seabreeze passage to aersol content respectively. To extend our work to the structure of the atmospheric boundary layer (ABL) we also make use of Smith and Carlin's (2024) fuzzy logic algorithm for ABL height identification. A modified version of this code is here, with a link to the original and a citation below. 

## Data
The data files was accessed through the ARM data viewer: https://www.arm.gov/data 

All the data for this project was collected in 2022 from the ARM AMF1 site in La Porte, TX during the TRACER campaign. Due to time constraints and data availability, the data was only analyzed for June 1st to July 31st. Zip files of these data are uploaded to this repository. The following datastreams were accessed:

For Aerosol work:
* AOSCPCF: _AOS Fine Condensation Particle Counter, 10 nm detection limit_
* AOSCPCU: _AOS Ultrafine Condensation Particle Counter, 2.5 nm detection limit_
* DLPROFWIND4NEWS: _Doppler lidar wind VAP_
* MWRLOS: _Microwave Water Radiometer (MWR) water liq. & vapor along line of sight (LOS) path_

For fuzzy logic:
* TROPOE: _Tropospheric Optimal Estimation Retrieval_ (This data is too large to upload to github and must be accessed via the following link: https://adc.arm.gov/discovery/#/results/id::houtropoeM1.c1_waterVapor_sfcstate_tropoe_derivmod?dataLevel=c1&showDetails=true )
* DLFPT: _Doppler Lidar fixed pointing mode_ (This data is too large to upload to github and must be accessed via the following link: https://adc.arm.gov/discovery/#/results/id::houdlfptM1.b1_attenuated_backscatter_macro_dl_cloud?dataLevel=b1&showDetails=true )
* DLPPI: _Dopper Lidar Plan Position Indicator_ (This data is too large to upload to github and must be accessed via the following link: https://adc.arm.gov/discovery/#/results/id::houdlppiM1.b1_attenuated_backscatter_macro_dl_cloud?dataLevel=b1&showDetails=true )

Both DFLPT and DLPPI need to be processed before they can be used in the fuzzy logic code. Processing for the DFLP can be found in the DLFP_Combine file and processing for the PPI can be found in the ARM_VAD_Processing file. Both of these files read in the data from ARM and restructure it for easy use in the fuzzy logic code.

However, if you do not want to process these files, the outputted data from the fuzzy logic code with boundary layer heights and vertical velocity variances is available here in teh zipped BLData file.

## SB_Identification: Functions and Scope
This script is used to identify when a sea-breeze was present in the data using gradients in wind speeds and water vapor (further explained in the attatched report). The script features three defined functions. These are used for designing and applying a lowpass filter, performing a windowed average, and performing a first derivate with a central differencing scheme. Each defined function gives a detailed description of inputs and outputs in the script itself. 

This script is broken down into 4 main sections:
* Defined functions and packages
* Data pathing, read in, and pre processing
* Data processing
* Data Graphing

## Aero_Correlation: Functions and Scope
Once cases are identified as either a baseline(no SB),SB, or consecutive SB case the correlation code can be run. This code reads in segregated data files (in individual folders based on cassification) and read in. For the SB and consecutive cases the SB passage time is calculated as a function seperately for data anlysis. This process is the same as the previous python script, but more concise. The script loops over classifications (1-3), then over associated case files. From here the mean aerosol content and standard error for each time period and classification are saved outside the loop and graphed. We also read in ABL height estimations and vertical velocity variances from the fuzzy logic code file. Hourly average ABL height and hourly averaged vertically integrated velocity variances are saved outside the loop and graphed.

This script is broken down into 4 main sections:
* Defined functions and packages
* Data pathing, read in, and pre processing
* Data processing
* Data Graphing

##ARM_fuzzy_main
This code is a modified version of the fuzzy logic algorithm for finding atmospheric boundary layer (ABL) height developed by Elizabeth Smith and Jacob Carlin. The orignal fuzzy-logic code can be found on Smith's GitHub: https://github.com/eeeeelizzzzz/bliss-fl. The associated publication is Smith, E. N., and J. T. Carlin, 2024: A multi-instrument fuzzy logic boundary-layer-top detection algorithm. Atmospheric Measurement Techniques, 17, 4087–4107, https://doi.org/10.5194/amt-17-4087-2024.

##ARM_VAD_Processing
This file restrucutres the raw lidar PPI data into VADs for fuzzy logic

##DLFP_Combine
This file restructures the multiple daily ARM files into one file per day for fuzzy logic.


## Dependencies 
The python codes require the following packages (all are generally available in standard python distributions):
* Numpy
* Matplotlib
* netCDF4
* glob
* datetime
* scipy
_only for fuzzy logic:_
* pandas
* metpy
* suntime
* cmocean

## Authors
Isaac Medina, Johns Hopkins University, imedina2@jhu.edu

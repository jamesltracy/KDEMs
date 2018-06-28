## This code assembles three calibrated Training KDE models for each year and creates a minimum consensus Training Subset Ensemble TSE) model
## NOTE: The geoprocessing environment must have the output coordinates, processing extent, snap raster and cell size set to match the human
## population density raster "pop10kmn3"

# make arcmap python commands and spatial analyst (sa) operations available
import arcpy
from arcpy.sa import *

import os

# Overwrite pre-existing files
arcpy.env.overwriteOutput = True

################################################  
# set working directory
arcpy.env.workspace = "G:/KDEM-master/KDEMVignetteData"

# read arcmap datasets into memory
list = arcpy.ListDatasets("*")
list.sort()
print list

# Set input and output directories
InDirectPre = "G:/KDEM-master/KDEMVignetteData/MonRstKDEMCNAA/" 
OutDirect = "G:/KDEM-master/KDEMVignetteData/"

NObskfoldgrppL = ['A', 'B', 'C']
#YearL = range(2005,2017,1)
YearL = ['02_16']
#YearL = range(2011,2017,1)
LoopCountList = range(0,3)

for Year in YearL:
    #Year=2005
    KDEMRasterL = range(0,3) # Create list for KDE rasters 
    for LoopCount in LoopCountList:
        #LoopCount = 0
        InDirect = InDirectPre + "KDEMProj_" + str(Year) + "_TrainSet" + NObskfoldgrppL[LoopCount] + "_ThresholdCalib" + str(LoopCount + 1) + NObskfoldgrppL[LoopCount] + "/"
        # Read in calibrated Raster
        KDEMRasterL[LoopCount] = Raster(InDirect + "MonRstKDEM_" + str(LoopCount + 1) + "ThreshCal_TrainSet" + NObskfoldgrppL[LoopCount] + "Cal.tif")
    # Sum raster in list
    KDEMRasterSumTSE = sum(KDEMRasterL)
    # Find minimum consensus of raster sum
    KDEMRasterConTSE = Con(KDEMRasterSumTSE >= 1,1,0)
    # Save raster
    KDEMRasterConTSE.save(OutDirect + "KDEM" + str(Year) + "TSEMN")


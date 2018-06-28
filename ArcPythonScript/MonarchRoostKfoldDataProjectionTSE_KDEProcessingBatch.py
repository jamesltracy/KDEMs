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
arcpy.env.workspace = "G:/MonarchRoost"

# read arcmap datasets into memory
list = arcpy.ListDatasets("*")
list.sort()
print list

# Set input and output directories
InDirectPre = "G:/MonarchRoost/MonRstKDEMCNAA/" 
OutDirect = "G:/MonarchRoost/"

NObskfoldgrppL = ['A', 'B', 'C']
#YearL = range(2005,2017,1)
YearL = ['02_16']
#YearL = range(2011,2017,1)
LoopCountList = range(0,3)

for Year in YearL:
    #Year=2005
    KDERasterL = range(0,3) # Create list for KDE rasters 
    for LoopCount in LoopCountList:
        #LoopCount = 0
        InDirect = InDirectPre + "KDEMProj_" + str(Year) + "_TrainSet" + NObskfoldgrppL[LoopCount] + "_ThresholdCalib" + str(LoopCount + 1) + NObskfoldgrppL[LoopCount] + "/"
        # Read in calibrated Raster
        KDERasterL[LoopCount] = Raster(InDirect + "MonRstKDEM_" + str(LoopCount + 1) + "ThreshCal_TrainSet" + NObskfoldgrppL[LoopCount] + "Cal.tif")
    # Sum raster in list
    KDERasterSumTSE = sum(KDERasterL)
    # Find minimum consensus of raster sum
    KDERasterConTSE = Con(KDERasterSumTSE >= 1,1,0)
    # Save raster
    KDERasterConTSE.save(OutDirect + "KDE" + str(Year) + "TSEMN")


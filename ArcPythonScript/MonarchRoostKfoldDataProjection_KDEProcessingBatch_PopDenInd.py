## This creates three training KDEs for each year in North America Albers Equal Area Conic projection
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
InDirect = "G:/KDEM-master/KDEMVignetteData"
Species = "MonRst"
OutDirect = InDirect + "/" + Species + "Spat10kMask"

# Read in background raster and create equivalent raster with value "1"
bckgrndrst=Raster(Species + "bck")
bckgrndrst1 = Int(Divide(Plus(bckgrndrst, 1),Plus(bckgrndrst, 1))) 
bckgrndrst1.save(InDirect + "/"  + Species + "bck1")


NObskfoldgrppL = ['A', 'B', 'C']
#YearL = range(2005,2017,1)
YearL = ['02_16']
#YearL = range(2011,2017,1)
LoopCountList = range(0,3)

for Year in YearL:
    # Set shapefile prefix
    ShpName = "MonRst" + str(Year) + "_PresProjTrain"
    for LoopCount in LoopCountList:
        #LoopCount = 1

        # Calculate kernel density raster using shapefile with values 1 to 10 weighted for human population
        # using default bandwidth calculation
        arcpy.gp.KernelDensity_sa(ShpName + NObskfoldgrppL[LoopCount] + ".shp", "",
                                  InDirect + "/" + "rstkde" + str(Year) + NObskfoldgrppL[LoopCount] + "N",
                                  "9552.93165844373", "", "SQUARE_KILOMETERS", "EXPECTED_COUNTS", "GEODESIC")

        # Multiply kernel density raster by background raster
        outras = Int(Times(Raster(InDirect + "/" + "rstkde" + str(Year) + NObskfoldgrppL[LoopCount] + "N"), Times(bckgrndrst1,10000)))
        outras.save(OutDirect + "/"  + "rstkde" + str(Year) + NObskfoldgrppL[LoopCount])



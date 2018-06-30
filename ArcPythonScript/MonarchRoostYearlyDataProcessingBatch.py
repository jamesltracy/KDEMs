## This code takes the 2002-2016 monarch roost data in North America Albers Equal Area Conic, breaks into separate years for 2005-2016, and obtains
## the Human population density index value for each point from the "rstdenpop" shapefile
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
InDirect = "G:/MonarchRoost"
OutDirect = "G:/MonarchRoost"

# Set shapefile prefix
ShpName = "MonRst"

# Read in "env10k" raster and multiply by 10000
env10kx10k=Times(Raster("env10k"), 10000)

YearL = range(2005,2017,1)

for Year in YearL:
    #Year = 2005
    InPtShp = InDirect + "/MonarchRoosts_2002_2016EastAlbers.shp"

    # Make shapefile a feature so that it is selectable
    arcpy.MakeFeatureLayer_management(InPtShp, "MonRst.lyr")

    # Select year from monarch roost shapefile
    arcpy.SelectLayerByAttribute_management("MonRst.lyr",
                                           "NEW_SELECTION", "Year = %s" % Year)

    # Write the selected features of layer to a new shapefile
    arcpy.CopyFeatures_management("MonRst.lyr", OutDirect + "/MonRst" + str(Year) + ".shp")

    # Convert to 10k resolution raster matching background evaluation extent raster
    arcpy.PointToRaster_conversion(in_features="MonRst" + str(Year) + ".shp", value_field="FID",
                                   out_rasterdataset=OutDirect + "/monrst" + str(Year), cell_assignment="MOST_FREQUENT",
                                   priority_field="NONE", cellsize="9552.93165844373")

    # Convert raster back to point shapefile aligned with points for "rstdenpopall.shp"
    arcpy.RasterToPoint_conversion(in_raster="MonRst" + str(Year), out_point_features=OutDirect + "/monrst" + str(Year) + "bck.shp",
                                   raster_field="VALUE")

    # Intersect point shapefile for yearly roost to "rstdenpopall.shp" to obtain human population density index value for kernel density
    arcpy.Intersect_analysis(in_features=["monrst" + str(Year) + "bck.shp", "rstdenpopall.shp"], out_feature_class=OutDirect + "/MonRst" + str(Year) + "denpop.shp",
                         join_attributes="ALL", cluster_tolerance="-1 Unknown", output_type="INPUT")

    # Add X Y coordinates to the shapefile
    arcpy.AddXY_management(in_features="MonRst" + str(Year) + "denpop.shp")


######
## In R, partition yearly presence data into three groups and create training subset ensemble

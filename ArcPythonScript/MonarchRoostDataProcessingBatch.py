## This code takes raw WGS84 monarch roost point shapefile, projects to North America Albers Equal Area Conic, assigns a field of Value = 1,
## converts to a raster, divides by human population raster to output roost density by population index
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
OutDirect = "G:/KDEM-master/KDEMVignetteData"

# Project shapefile from WGS84 to North America Equal Area Albers
arcpy.Project_management(in_dataset=InDirect + "/MonarchRoosts_2002_2016East.shp", out_dataset=OutDirect + "/MonarchRoosts_2002_2016EastAlbers.shp",
                         out_coor_system="PROJCS['North_America_Albers_Equal_Area_Conic',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-96.0],PARAMETER['Standard_Parallel_1',20.0],PARAMETER['Standard_Parallel_2',60.0],PARAMETER['Latitude_Of_Origin',40.0],UNIT['Meter',1.0]]", transform_method="WGS_1984_(ITRF00)_To_NAD_1983", in_coor_system="GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],VERTCS['Unknown VCS',VDATUM['Unknown'],PARAMETER['Vertical_Shift',0.0],PARAMETER['Direction',1.0],UNIT['Meter',1.0]]",
                         preserve_shape="NO_PRESERVE_SHAPE", max_deviation="")

# Add field of "Value" to shapefile before converting to raster
arcpy.AddField_management(OutDirect + "/MonarchRoosts_2002_2016EastAlbers.shp",\
                          "Value", "SHORT", 3, "", "", "", "NULLABLE", "REQUIRED", "")

# Give "Id" field a value of 1 before converting to raster
cursor =  arcpy.UpdateCursor(OutDirect + "/MonarchRoosts_2002_2016EastAlbers.shp")  
      
for row in cursor:  
    row.setValue("Value", 1)
    cursor.updateRow(row)
    
del cursor, row   

# Convert shapefile to raster using MEAN option to only count 1 roost per cell to spatially thin to 10 km
arcpy.PointToRaster_conversion(OutDirect + "/MonarchRoosts_2002_2016EastAlbers.shp", "Value", 
                                 OutDirect + "/rstden10kall", 
                                 "MEAN", "", "")

# Divide number of roosts by human population density over a 3 x 3 cell radius at 10 km resolution
outras = (Raster(OutDirect + "/rstden10kall"))/(Raster((InDirect + "/pop10kmn3")) + 0.00001)

# Limit maximum value to 10 and minimum value to one
outras2 = Con(outras > 10, 10, outras)
outras3 = Con(outras2 < 1, 1, outras2)
outras3.save(OutDirect + "/rstdnpall")

# Convert above raster to point shapefile
arcpy.RasterToPoint_conversion(in_raster=Raster(OutDirect + "/rstdnpall"),
                               out_point_features=OutDirect + "/rstdenpopall.shp",
                               raster_field="Value")

# Convert original raster with values all one to point shapefile
arcpy.RasterToPoint_conversion(in_raster=Raster(OutDirect + "/rstden10kall"),
                               out_point_features=OutDirect + "/rst10kall.shp",
                               raster_field="Value")

# Calculate kernel density raster using shapefile with values 1 to 10 weighted for human population
#arcpy.gp.KernelDensity_sa("rstdenpopall", "GRID_CODE",
#                          OutDirect + "/rstdnpallkd",
#                          "9552.93165844373", "", "SQUARE_KILOMETERS", "EXPECTED_COUNTS", "GEODESIC")

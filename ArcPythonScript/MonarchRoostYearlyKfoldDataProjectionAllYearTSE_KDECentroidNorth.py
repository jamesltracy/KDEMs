## Obtains the centroid and average eastern and western longitudes for the Training Subset Ensemble monarch roost KDE model (KDEM) over all years
## for the 27 to 37N portion of the central flyway, calculating the average width of the KDE and the north-south and east-west distances
## of the training KDE centroids to the training subset ensemble centroid using data from all years. Data are collected in a csv file.
## NOTE: The geoprocessing environment must have the output coordinates, processing extent, snap raster and cell size set to match the human
## population density raster "pop10kmn3"

# make arcmap python commands and spatial analyst (sa) operations available
import arcpy
import numpy
from arcpy.sa import *
from dbfpy import dbf

import os

## Python function to convert dbf to csv from https://geonet.esri.com/thread/110894
from os import path as p  
  
def TableToCSV(fc,CSVFile):  
      
    fields = [f.name for f in arcpy.ListFields(fc) if f.type <> 'Geometry']  
    with open(CSVFile, 'w') as f:  
        f.write(','.join(fields)+'\n') #csv headers  
        with arcpy.da.SearchCursor(fc, fields) as cursor:  
            for row in cursor:  
                f.write(','.join([str(r) for r in row])+'\n')  
    print 'Created %s Successfully' %p.basename(CSVFile)  

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
InDirectPre = "G:/MonarchRoost/"
OutDirect = "G:/MonarchRoost/"

# Read in calibrated Raster
KDEMRaster = Raster(OutDirect + "kdem02_16tse")

# Set Null zero values for KDE model
KDEMRaster1 = SetNull(KDEMRaster<1,1)

# Convert KDE model to polygon
arcpy.RasterToPolygon_conversion(in_raster=KDEMRaster1, out_polygon_features=OutDirect + "KDEM1poly.shp",
                                 simplify="NO_SIMPLIFY", raster_field="VALUE")

# Intersect 37 N and above Central Flyway area and KDEM polygon
arcpy.Intersect_analysis(in_features=[OutDirect + "KDEM1poly.shp", OutDirect + "CentralFlywayRoughNorthAlb.shp"],
                         out_feature_class=OutDirect + "KDEM1EvalAreapoly.shp",
                         join_attributes="ALL", cluster_tolerance="-1 Unknown", output_type="INPUT")

### Create point shapefile of centroid for above intersect polygon and retrieve x y coordinates
# Convert above intersect polygon to raster
arcpy.PolygonToRaster_conversion(in_features=OutDirect + "KDEM1EvalAreapoly.shp", value_field="GRIDCODE",
                                 out_rasterdataset="KDEMRaster2",
                                 cell_assignment="CELL_CENTER", priority_field="NONE", cellsize="9552.93165844373")

# Find centroid cell of raster
arcpy.gp.ZonalGeometry_sa(OutDirect + "KDEMRaster2", "VALUE",
                          OutDirect + "KDEMRasterc", "CENTROID", "9552.93165844373")

# Convert centroid raster to point shapefile
arcpy.RasterToPoint_conversion(in_raster=OutDirect + "KDEMRasterc",
                               out_point_features=OutDirect + "KDEM02_16TSECentroidNorth.shp",
                               raster_field="VALUE")

# Find x y coordinates of centroid
# First add x y coordinates to centroid point shapefile
arcpy.AddXY_management(in_features=OutDirect + "KDEM02_16TSECentroidNorth.shp")



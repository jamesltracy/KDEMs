## This code generates 10 km resolution background and pseudoabsence rasters
## to use in R script for generating background and pseudoabsence points
## NOTE: The geoprocessing environment must have the output coordinates, processing extent, snap raster and cell size set to match the human
## population density raster "pop10kmn3"

# make arcmap python commands and spatial analyst (sa) operations available
import arcpy
from arcpy.sa import *

# Overwrite pre-existing files
arcpy.env.overwriteOutput = True

import os
###########################################################################################
# Set names of working directory, species shapefile, evaluation area shapefile, numbers of desired pseudoabasence
# and background points, and various buffer sizes
OutDirect = "G:/MonarchRoost/" # Name of working and output directory
SpecName= "MonRst"
SpecFileName = "rstdenpopall"  # Name of species point shapefile (not thinned)
PseudoabsenceBuffDist = "100" # In Kilometers. Pseudoabsence buffer distance from Presence points.
ExtentBuffer = "500 kilometers"

# Set work directory
arcpy.env.workspace = OutDirect

# Create a land mass raster with value 1 called "env10k"
outras = Int(Divide((Raster("pop10kmn3")+1),(Raster("pop10kmn3")+1)))
outras.save(OutDirect + "env10k")

# Create convex hull polygon encompassing unthinned presence points
arcpy.MinimumBoundingGeometry_management(in_features=SpecFileName,
                                         out_feature_class=OutDirect + SpecFileName + "ConvexHull.shp",
                                         geometry_type="CONVEX_HULL", group_option="ALL",
                                         group_field="", mbg_fields_option="NO_MBG_FIELDS")

# Buffer around above convex hull polygon by ExtentBuffer
arcpy.Buffer_analysis(in_features=OutDirect + SpecFileName + "ConvexHull.shp",
                      out_feature_class=OutDirect + SpecName + "Extent.shp",
                      buffer_distance_or_field=ExtentBuffer, line_side="FULL",
                      line_end_type="ROUND", dissolve_option="NONE", dissolve_field="", method="PLANAR")

# Convert extent polygon shapefile to raster for further processing in R where converted to points, randomly sampled to 20,000 and thinned to 10km
arcpy.PolygonToRaster_conversion(in_features=OutDirect + SpecName + "Extent.shp", value_field="FID",
                                 out_rasterdataset=OutDirect + SpecName + "bck1",
                                 cell_assignment="CELL_CENTER",
                                 priority_field="NONE")

# Multiply above raster by env10k to restrict to land mass
outras = Times(Raster(OutDirect + SpecName + "bck1"), Raster(OutDirect + "env10k"))
outras.save(OutDirect + SpecName + "bck")

# Create PseudoabsenceBuffDist km buffer shapefile around unthinned presence points (may take about 10 minutes)
arcpy.Buffer_analysis(in_features=OutDirect + SpecFileName + ".shp",
                      out_feature_class=OutDirect + SpecFileName + "_" + PseudoabsenceBuffDist + "kmbuf.shp",
                      buffer_distance_or_field=PseudoabsenceBuffDist + " Kilometers",
                      line_side="FULL", line_end_type="ROUND",
                      dissolve_option="ALL",
                      dissolve_field="",
                      method="GEODESIC")

# Use the buffer as a "cut out" to remove buffered areas from evaluation area shapefile
arcpy.Erase_analysis(in_features=OutDirect +  SpecName + "Extent.shp", erase_features=SpecFileName + "_" + PseudoabsenceBuffDist + "kmbuf",
                     out_feature_class=OutDirect + SpecName + "psaarea.shp",
                     cluster_tolerance="")

# Convert extent with cut out polygon shapefile to raster for further processing in R where converted to points, randomly sampled to 20,000 and thinned to 10km
arcpy.PolygonToRaster_conversion(in_features=OutDirect + SpecName + "psaarea.shp", value_field="FID",
                                 out_rasterdataset=OutDirect + SpecName + "psr1",
                                 cell_assignment="CELL_CENTER",
                                 priority_field="NONE")

# Multiply above raster by env10k to restrict to land mass
outras = Times(Raster(OutDirect + SpecName + "psr1"), Raster(OutDirect + "env10k"))
outras.save(OutDirect + SpecName + "psr")

## Use above created background and pseudoabsence rasters in R script for generating background and pseudoabsence points

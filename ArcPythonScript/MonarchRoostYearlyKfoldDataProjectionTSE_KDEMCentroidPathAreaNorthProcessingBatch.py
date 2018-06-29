## This code assembles three monarch roost Training KDEM models for each year and obtains the centroid and average eastern and western longitudes
## for the 27 to 37N portion of the central flyway, calculating the average width of the KDEM and the north-south and east-west distances
## of the training KDEM centroids to the training subset ensemble centroid using data from all years. Data are collected in a csv file.
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
arcpy.env.workspace = "G:/KDEM-master/KDEMVignetteData"

# read arcmap datasets into memory
list = arcpy.ListDatasets("*")
list.sort()
print list

# Set input and output directories
InDirectPre = "G:/KDEM-master/KDEMVignetteData/"
OutDirect = "G:/KDEM-master/KDEMVignetteData/"
#
NObskfoldgrppL = ['A', 'B', 'C']
YearL1 = range(2005,2017,1)
YearL = ['02_16'] + YearL1
#YearL = range(2014,2017,1)
#
LoopCountList = range(0,3)

#Retrieve x y coordinates of overall centroid computed from Training Subset Ensemble of 2002-2016 roost data

# First retrieve X coordinate
# Create SearchCursor
rows = arcpy.SearchCursor(OutDirect + "KDEM02_16TSECentroidNorth.dbf")
# Get value of YCent
field = arcpy.ListFields(OutDirect + "KDEM02_16TSECentroidNorth.dbf")[4]
name = field.name
for row in rows:
    XOCent = row.getValue(name)

# Second retrieve Y coordinate
# Create SearchCursor
rows = arcpy.SearchCursor(OutDirect + "KDEM02_16TSECentroidNorth.dbf")
# Get value of YCent
field = arcpy.ListFields(OutDirect + "KDEM02_16TSECentroidNorth.dbf")[5]
name = field.name
for row in rows:
    YOCent = row.getValue(name)

###

# Create DBF to hold output data of distances for annual training KDEM width and centroid north-south and east-west shifts
OutputTable = dbf.Dbf(OutDirect + "KDEMWidth_CentroidShift_AreaDistancesNorth.dbf", new=True)

# Create fields in table
OutputTable.addField(
    ("KDEMYear", "C", 15),
    ("TrainSet", "C", 15),
    ("AreaSqKm", "N", 15, 6),
    ("WidthDist", "N", 15, 6),
    ("NSCentDist", "N", 15, 6),
    ("EWCentDIst", "N", 15, 6),
)
OutputTable.close()
#print OutputTable

###############################################################################################################################
for Year in YearL:
    #Year=2008
    for LoopCount in LoopCountList:
        #LoopCount = 0
        InDirect = InDirectPre + "KDEMProj_" + str(Year) + "_TrainSet" + NObskfoldgrppL[LoopCount] + "_ThresholdCalib" + str(LoopCount + 1) + NObskfoldgrppL[LoopCount] + "/"
        # Read in calibrated Raster
        KDEMRaster = Raster(InDirect + "MonRstKDEM_" + str(LoopCount + 1) + "ThreshCal_TrainSet" + NObskfoldgrppL[LoopCount] + "Cal.tif")

        # Set Null zero values for KDEM model
        KDEMRaster1 = SetNull(KDEMRaster<1,1)

        # Convert KDEM model to polygon
        arcpy.RasterToPolygon_conversion(in_raster=KDEMRaster1, out_polygon_features=OutDirect + "KDEM1poly.shp",
                                         simplify="NO_SIMPLIFY", raster_field="VALUE")

        # Intersect 27 to 37 N Central Flyway area and KDEM polygon
        arcpy.Intersect_analysis(in_features=[OutDirect + "KDEM1poly.shp", OutDirect + "CentralFlywayRoughNorthAlb.shp"],
                                 out_feature_class=OutDirect + "KDEM1EvalAreapolyNorth.shp",
                                 join_attributes="ALL", cluster_tolerance="-1 Unknown", output_type="INPUT")


        # Convert above intersect polygon to raster
        arcpy.PolygonToRaster_conversion(in_features=OutDirect + "KDEM1EvalAreapolyNorth.shp", value_field="GRIDCODE",
                                         out_rasterdataset="KDEMRaster2",
                                         cell_assignment="CELL_CENTER", priority_field="NONE", cellsize="9552.93165844373")

        ############################################
        ## Calculate mean width of KDEM between 27N and 37N using sums of raster cells across rows

        # Convert Raster to numpy array
        arr = arcpy.RasterToNumPyArray(OutDirect + "KDEMRaster2",nodata_to_value=0)

        # Obtain sum of all rows of array from raster
        arrSum = arr.sum(1)

        # Keep only sum values not equal to zero
        arrSum2 = arrSum[arrSum != 0]

        # Find mean of all nonzero row values
        MeanarrSum2 = sum(arrSum2)/len(arrSum2)

        # Find width in cell by meters
        cellSize = Raster(OutDirect + "KDEMRaster2").meanCellWidth

        # Find mean width of rows in kilometers
        MeanWidthKm = (MeanarrSum2 * cellSize)/1000000

        ###################
        ## Calculate area using above array values

        # Find sum of all nonzero row values
        SumarrSum2 = sum(arrSum2)

        # Find area of raster
        AreaSqKm = (SumarrSum2 * cellSize)/1000

        ################################################
        ### Create point shapefile of centroid for above intersect raster and retrieve x y coordinates

        # Find centroid cell of raster
        arcpy.gp.ZonalGeometry_sa(OutDirect + "KDEMRaster2", "VALUE",
                                  OutDirect + "KDEMRasterc", "CENTROID", "9552.93165844373")

        # Convert centroid raster to point shapefile
        arcpy.RasterToPoint_conversion(in_raster=OutDirect + "KDEMRasterc",
                                       out_point_features=OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidNorth.shp",
                                       raster_field="VALUE")

        # Find x y coordinates of centroid
        # First add x y coordinates to centroid point shapefile
        arcpy.AddXY_management(in_features=OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidNorth.shp")

        ## Retrieve y coordinate from dbf file of shapefile

        # First retrieve X coordinate
        # Create SearchCursor
        rows = arcpy.SearchCursor(OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidNorth.dbf")
        # Get value of YCent
        field = arcpy.ListFields(OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidNorth.dbf")[4]
        name = field.name
        for row in rows:
            XCent = row.getValue(name)

        # Second retrieve Y coordinate
        # Create SearchCursor
        rows = arcpy.SearchCursor(OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidNorth.dbf")
        # Get value of YCent
        field = arcpy.ListFields(OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidNorth.dbf")[5]
        name = field.name
        for row in rows:
            YCent = row.getValue(name)

        ###############
        ## Calculate distance from current centroid to overall year combined centroid to find north-south and east-west shift

        # Align current centroid on same latitude (YOCent) as overall centroid to determine east-west shift

        # First create new dbf table
        table = dbf.Dbf(OutDirect + "KDEMCentroidLatitudeEWNorth.dbf", new=True)

        # Create fields in table
        table.addField(
            ("NAME", "C", 15),
            ("X", "N", 15, 6),
            ("Y", "N", 15, 6),
        )
        #print table

        ## fill DBF records with x of current KDEM centroid and y of overall KDEM centroid
        for name, x, y in (
            ("CentEW", XCent, YOCent),
        ):
            rec = table.newRecord()
            rec["NAME"] = name
            rec["X"] = x
            rec["Y"] = y
            rec.store()
        table.close()

        ## read DBF and print records
        ##table = dbf.Dbf(OutDirect + "KDEMCentroidLatitudeEW.dbf")
        ##for rec in table:
        ##    print rec
        ##print

        # Convert east/west point dbf file to x y feature
        arcpy.MakeXYEventLayer_management(table=OutDirect + "KDEMCentroidLatitudeEWNorth.dbf",
                                          in_x_field="X", in_y_field="Y", out_layer="KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidLatitudeEWPointNorth",
                                          spatial_reference="PROJCS['North_America_Albers_Equal_Area_Conic',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-96.0],PARAMETER['Standard_Parallel_1',20.0],PARAMETER['Standard_Parallel_2',60.0],PARAMETER['Latitude_Of_Origin',40.0],UNIT['Meter',1.0]];-16688100 -9068200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision", in_z_field="")

        # Save layer to shapefile
        arcpy.FeatureClassToShapefile_conversion("KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidLatitudeEWPointNorth", OutDirect)


        # Align current centroid on same longitude (XOCent) as overall centroid to determine north-south shift

        # First create DBF to hold point coordinates
        table = dbf.Dbf(OutDirect + "KDEMCentroidLongitudeNSNorth.dbf", new=True)

        # Create fields in table
        table.addField(
            ("NAME", "C", 15),
            ("X", "N", 15, 6),
            ("Y", "N", 15, 6),
        )
        #print table

        ## fill DBF records with x of overall KDEM centroid and y of current KDEM centroid
        for name, x, y in (
            ("CentNS", XOCent, YCent),
        ):
            rec = table.newRecord()
            rec["NAME"] = name
            rec["X"] = x
            rec["Y"] = y
            rec.store()
        table.close()

        ## read DBF and print records
        ##table = dbf.Dbf(OutDirect + "KDEMCentroidLongitudeNS.dbf")
        ##for rec in table:
        ##    print rec
        ##print

        # Convert north/south point dbf file to x y feature
        arcpy.MakeXYEventLayer_management(table=OutDirect + "KDEMCentroidLongitudeNSNorth.dbf",
                                          in_x_field="X", in_y_field="Y", out_layer="KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidLongitudeNSPointNorth",
                                          spatial_reference="PROJCS['North_America_Albers_Equal_Area_Conic',GEOGCS['GCS_North_American_1983',DATUM['D_North_American_1983',SPHEROID['GRS_1980',6378137.0,298.257222101]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Albers'],PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],PARAMETER['Central_Meridian',-96.0],PARAMETER['Standard_Parallel_1',20.0],PARAMETER['Standard_Parallel_2',60.0],PARAMETER['Latitude_Of_Origin',40.0],UNIT['Meter',1.0]];-16688100 -9068200 10000;-100000 10000;-100000 10000;0.001;0.001;0.001;IsHighPrecision", in_z_field="")

        # Save layer to shapefile
        arcpy.FeatureClassToShapefile_conversion("KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidLongitudeNSPointNorth", OutDirect)

        #######

        ## Calculate distance between "CentNS" and overall centroid and "CentEW" and overall centroid

        # First calculate east to west shift
        arcpy.PointDistance_analysis(in_features=OutDirect + "KDEM02_16TSECentroidNorth.shp",
                                     near_features=OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidLatitudeEWPointNorth.shp",
                                     out_table=OutDirect + "KDEMCentroidEasttoWestDistanceNorth.dbf", search_radius="")

        # Extract centroid east to west shift distance from dbf file
        # Create SearchCursor
        rows = arcpy.SearchCursor(OutDirect + "KDEMCentroidEasttoWestDistanceNorth.dbf")
        # Get value of distance
        field = arcpy.ListFields(OutDirect + "KDEMCentroidEasttoWestDistanceNorth.dbf")[3]
        name = field.name
        for row in rows:
            CentroidEWShift = row.getValue(name)

        ##

        # Then calculate north to south shift
        arcpy.PointDistance_analysis(in_features=OutDirect + "KDEM02_16TSECentroidNorth.shp",
                                     near_features=OutDirect + "KDEM" + str(Year) + NObskfoldgrppL[LoopCount] + "CentroidLongitudeNSPointNorth.shp",
                                     out_table=OutDirect + "KDEMCentroidNorthtoSouthDistanceNorth.dbf", search_radius="")

        # Extract centroid east to west shift distance from dbf file
        # Create SearchCursor
        rows = arcpy.SearchCursor(OutDirect + "KDEMCentroidNorthtoSouthDistanceNorth.dbf")
        # Get value of distance
        field = arcpy.ListFields(OutDirect + "KDEMCentroidNorthtoSouthDistanceNorth.dbf")[3]
        name = field.name
        for row in rows:
            CentroidNSShift = row.getValue(name)

        ## Convert distances from meters to kilometers
        CentroidNSShiftkm = CentroidNSShift/1000
        CentroidEWShiftkm = CentroidEWShift/1000

        # Adjust CentroidNSShiftkm to negative if YCent < YOCent
        if YCent < YOCent:
            CentroidNSShiftkm = CentroidNSShiftkm * -1

        # Adjust CentroidEWShiftkm to negative if XCent < XOCent
        if XCent < XOCent:
            CentroidEWShiftkm = CentroidEWShiftkm * -1

        #########################################################################
        ## Save KDEM width and centroid displacement output by appending to DBF file

        # Open existing DBF
        table = dbf.Dbf(OutDirect + "KDEMWidth_CentroidShift_AreaDistancesNorth.dbf")
        # Add new records
        for year, trainset, areasqkm, meanwidthkm, centroidnsshift, centroidewshift in (
            (Year, NObskfoldgrppL[LoopCount], AreaSqKm, MeanWidthKm, CentroidNSShiftkm, CentroidEWShiftkm),
        ):
            rec = table.newRecord()
            rec["KDEMYear"] = year
            rec["TrainSet"] = trainset
            rec["AreaSqKm"] = areasqkm
            rec["WidthDist"] = meanwidthkm
            rec["NSCentDist"] = centroidnsshift
            rec["EWCentDist"] = centroidewshift    
            rec.store()
        table.close()

############################################################################################

# Convert above OutPut table to CSV file
TableToCSV(OutDirect + "KDEMWidth_CentroidShift_AreaDistancesNorth.dbf", OutDirect + "KDEMWidth_CentroidShift_AreaDistancesNorth.csv")

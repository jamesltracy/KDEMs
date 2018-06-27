##############################################################################
# This R program generates Kernel Density Estimate Models (KDEMs)
# in conjunction with KDE surface rasters produced using Arc Python
# script with ArcGIS (referenced in comments) for migration modeling
# described in Tracy et al. (in prep.)
#
# Tracy JL, Kantola T, Baum KA, Coulson RN in prep.) Modeling fall migration
#   pathways and spatially identifying potential migratory hazards for the
#   eastern monarch butterfly
#
# The program specifically processes species occurrence data (such as for
# migratory movement) to generate, calibrate and evaluate spatially interpolative KDEMs
##############################################################################
#
##################################################################################
# NOTE: Some sections are commented since pre-processed presence, pseudoabsence,
# and background data is provided in the KDEMVignetteData folder
# In addition, several k-fold partition schemes are commented since they are also provided
# in KDEMVignetteData folder
# Modify the three input directories to match where you put input files on your system
# 1) For R Functions: C:/Users/James/Documents/R/win-library/
# 2) For csv files with occurrence data and k-fold partition schemes
#    G:/MonarchRoost
#################################################################################
#
##############################################################################
# This section loads libraries, sets working directory, and establishes source
# file locations and output file naming conventions
##############################################################################
# Specify  directory for saved functions
FunctDirect <- "C:/Users/James/Documents/R/win-library/"
setwd(FunctDirect)
# Read in User Defined Functions
source("SpatialFilter_Function.R")
source("Round2_Function.R")
source("KDEMSubset_GridTrainTestEvalAIC_Calib_Function.R")
source("EvalStatVars_Summary_Function.R")
source("GamesHowell.Test.Padj_Function.R")
#startSocketServer(port=8889)
#
# Load needed packages of raster, rgdal, dismo, rjava, and maptools (printouts not shown)
library(dismo)
library(maptools)
library(rgdal)
library(sp)
library(PresenceAbsence)
library(raster)
library(car)
library(caret)
library(ENMeval)
library(devtools)
#install_github('johnbaums/rmaxent')  # for installing rmaxent
#library(rmaxent)
#library(rasterVis)
#library(viridis)
# Packages used in functions:
library(plyr)
# Create definition for a geographical projection
#crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
#crs.geo <- CRS("+init=epsg:26912")  # geographical, datum NAD83 12N

# Specify  directory for general data and shapefiles
InDirect <- "G:/MonarchRoost"
setwd(InDirect)
## Specify extent and various descriptors for maxent run
Extent <- "CNA"
PresLimit <- 10000 # Limited number of presence points through random sample after spatial thinning
SpatFiltBuff <- 10 # units in km
SpecFileName <- "MonRst"
Species <- paste0("MonRst")
SpecShort <- paste0("MonR")
SpeciesGen <- "Monarch"
Species2 <- "MnRst"
MaxBackgrndPseudoabsLimit <- 10000 # specify limit for background and pseudoabsence data
# Specify resolution of rasters in square km
ResolutionEnv <- 100.0  ## in units of square kilometers
Climate <- "Current"
Tag <- "A"
LongSpecies <- "MonarchRoosts"
ExtentRaster <- "monrstbck" # name of background evaluation extent raster
# Specify units of buffers
Units <- "km"
# Specify buffers for spatial thinning and pseudoabsenct points from presence points
SpatFiltBuffs <- "10km"
# Specify if spatial filtering
SpatFilt <- "Filtered" # there is spatial filtering
PsAbsBuff <- 100  # in km
PsAbsBuffs <- "100km"
Conditions <- paste0(PsAbsBuff, "mPsAbsBuff_", SpatFiltBuff, "mSpatFiltBuff")
Conditions2 <- paste("Only Native Data for Model Input")
# Set number of kfold partitions
NObsJ <- 3
# Set projection
CRS.WGS84 <- CRS("+init=epsg:4326")
#
ModelName <- paste("KDEM")
###############
## Create directory for storing KDE rasters
GridDirect <- paste0(InDirect, "/", Species, "Spat10kMask")
dir.create(GridDirect)
################
#### IMPORTANT- Specify output names for files and directories
# Identifying part of output rasters
## Specify run id consisting of date
Run <- 1
Set <- Sys.Date()
output1 <- paste0(Species, ModelName, Extent, Tag)
# Identifying part of output directory
# Create directory for storing envelope score related grids
dir.create(paste0(InDirect, "/", output1))
OutDirect <- paste0(InDirect, "/", output1)
#
dir.create(paste0(OutDirect))
OutDirectpsa <- paste0(OutDirect)
#
#########################################################################################
### In ArcGIS, use Python script MonRoostDataProcessingBatch.py
### to generate human population density weighted roost points at
### 10 km resolution
########################################################################################

##############################################################################
## Read in 10 km gridded presence points for eastern population
setwd(InDirect)
PresThinShp <- readOGR(".", paste0("rstdenpopall"))  # NOTE: readShapePoly was giving error, so used
summary(PresThinShp)
#Set projection
CRS.NAAlbersEqualAreaConic <- CRS(proj4string(PresThinShp))
head(PresThinShp)
#str(PresThinShpa)
AdminShp <- readOGR(".", paste0("AdminGlobalAlbers"))
#
dev.new()
plot(AdminShp)
plot(PresThinShp, add=TRUE, col="blue")
#
##############################################################################
###################################################################
### In ArcGIS, use Python script MonarchRoostBackgroundExtentRaster_PseudoabsenceRasterGeneration.py
### to generate a background raster and pseudoabsence raster for
### deriving Background and Pseudoabsence points
########################################################################################
### Generate and spatially filter pseudoabsence points
###
## Read in above created raster in ArcGIS for pseudoabsence points
PSARaster <- raster(paste0(SpecFileName, "psr"))
#plot(PSARaster)
## Convert  raster to spatial point data frames, restrict to 20,000 max
system.time(psaraw.spdf <- rasterToPoints(PSARaster, spatial=TRUE))
nrow(psaraw.spdf)
head(psaraw.spdf)
tail(psaraw.spdf)
########
## Randomly thin points to 10,000
psarawPresLimit.spdf <- psaraw.spdf[sample(1:length(psaraw.spdf), size=10000), ]
head(psarawPresLimit.spdf)
nrow(psarawPresLimit.spdf)
### Spatially thin presence points to minimum distance of SpatFiltBuff (generally 10 km) (takes 9 minutes with 10,000 points)
#system.time(psathin.spdf <- SpatialFilter(psarawPresLimit.spdf , dist=SpatFiltBuff, mapUnits=F))
#nrow(psathin.spdf)
#head(psathin.spdf)
psathin.spdf <- psarawPresLimit.spdf
#plot(psathin.spdf, add=TRUE, col='red')
# Convert to data frame
psathin.df <- data.frame(psathin.spdf)
head(psathin.df)
# Reformat
psathin.df <- psathin.df[,c(2,3,1)]
nrow(psathin.df)
write.table(psathin.df, file=paste0(SpecFileName, "psa.csv"))
## Convert to shapefile for display in ArcGIS
coordinates(psathin.df)=~x+y
proj4string(psathin.df)=CRS.NAAlbersEqualAreaConic
# Convert first to spatial points data frame
psathin.spdf <- SpatialPointsDataFrame(psathin.df, data.frame(id=1:length(psathin.df)))
head(psathin.spdf)
#plot(psathin.spdf, add=TRUE, col='green')
# Write to shapefile
writeOGR(psathin.spdf, InDirect, paste0(SpecFileName, "psa"), driver = "ESRI Shapefile", overwrite=TRUE)
##############################
## Read in evaluation area raster created in ArcGIS for generating background points
bkgrndRaster <- raster(paste0(SpecFileName, "bck"))
#plot(bkgrndRaster, add=TRUE)
## Convert background raster to points
system.time(bkgrndraw.spdf <- rasterToPoints(bkgrndRaster, spatial=TRUE))
nrow(bkgrndraw.spdf)
head(bkgrndraw.spdf)
tail(bkgrndraw.spdf)
########
## Randomly thin points to 10,000
bkgrndrawPresLimit.spdf <- bkgrndraw.spdf[sample(1:length(bkgrndraw.spdf), size=10000), ]
head(bkgrndrawPresLimit.spdf)
nrow(bkgrndrawPresLimit.spdf)
### Spatially thin presence points to minimum distance of SpatFiltBuff (generally 10 km) (takes 9 minutes with 10,000 points)
#system.time(bkgrndthin.spdf <- SpatialFilter(bkgrndrawPresLimit.spdf , dist=SpatFiltBuff, mapUnits=F))
#nrow(bkgrndthin.spdf)
#head(bkgrndthin.spdf)
bkgrndthin.spdf <- bkgrndrawPresLimit.spdf
#plot(bkgrndthin.spdf, add=TRUE, col='blue')
# Convert to data frame
bkgrndthin.df <- data.frame(bkgrndthin.spdf)
head(bkgrndthin.df)
# Reformat
bkgrndthin.df <- bkgrndthin.df[,c(2,3,1)]
nrow(bkgrndthin.df)
write.table(bkgrndthin.df, file=paste0(SpecFileName, "bkgrnd.csv"))
## Convert to shapefile for display in ArcGIS
coordinates(bkgrndthin.df)=~x+y
proj4string(bkgrndthin.df)=CRS.NAAlbersEqualAreaConic
# Convert first to spatial points data frame
bkgrndthin.spdf <- SpatialPointsDataFrame(bkgrndthin.df, data.frame(id=1:length(bkgrndthin.df)))
head(bkgrndthin.spdf)
# Write to shapefile
writeOGR(bkgrndthin.spdf, InDirect, paste0(SpecFileName, "bkgrnd"), driver = "ESRI Shapefile")
#


################################################################################
### Create Kernel Density Estimate surfaces for 2/3 training data for projections and calibration
########################################################################
### Create and save NObsJ usually three kfold partition scheme for multiple rounds of test with presence
### pseudoabsence and background points using all data
### NOTE: Only create and save partition scheme once in order to reuse
###### Read in presence data for kfold assignment
setwd(InDirect)
PresThinShp <- readOGR(".", paste0("rstdenpopall"))  # NOTE: readShapePoly was giving error, so used
summary(PresThinShp)
#Set projection
CRS.NAAlbersEqualAreaConic <- CRS(proj4string(PresThinShp))
head(PresThinShp)
nrow(PresThinShp)
#plot(PresThinShp, add=TRUE, col='blue')
#
CRS.In <- CRS.NAAlbersEqualAreaConic
#
# Read PresThinShp to spatial points data frame
PresThin.spdf <- readShapePoints("rstdenpopall.shp", proj4string=CRS.NAAlbersEqualAreaConic)
head(PresThin.spdf)
# Convert to data frame
PresThin.df <- data.frame(PresThin.spdf)
head(PresThin.df)
# Reformat
PresThin.df <- PresThin.df[,c(3,4,2,1)]
nrow(PresThin.df)
head(PresThin.df)
##
## Read in pseudoabsence points
library(dismo)
#
################################################################################
### Create NObs kfold group for presence points and use to define training KDEs
### NOTE: Only create NObs kfold groups once
NObskfoldgrpp <- kfold(PresThin.df, NObsJ)
#########
## For the three sets of training NObskfoldgrpp sets above, convert kfold presence training data of “rstdenpopall.shp” to
## shapefiles below and process with python "MonarchRoostKFoldDataProcessingBatch.py" to create
## kernel density estimate surfaces of chosen bandwidth for the training data
##
NObskfoldgrppL <- c("A", "B", "C")
GroupFirst <- c(1,2,3)
GroupSecond <- c(2,3,1)
count = 0
for(Group in NObskfoldgrppL) {
  count = count + 1
  PresThinProjTrain.df1 <- PresThin.df[NObskfoldgrpp==GroupFirst[count],]
  PresThinProjTrain.df2 <- PresThin.df[NObskfoldgrpp==GroupSecond[count],]
  PresThinProjTrain.df <- rbind(PresThinProjTrain.df1, PresThinProjTrain.df2)
  nrow(PresThinProjTrain.df)
  head(PresThinProjTrain.df)
  xy <- PresThinProjTrain.df[,1:2]
  PresThinProjTrain.spdf <- SpatialPointsDataFrame(coords = xy, data = PresThinProjTrain.df,
                                 proj4string = CRS.NAAlbersEqualAreaConic)
  head(PresThinProjTrain.spdf)
  nrow(PresThinProjTrain.spdf)
  #
  # Write to shapefile
  writeOGR(PresThinProjTrain.spdf, InDirect, paste0(SpecFileName, "PresProjTrain", Group), driver = "ESRI Shapefile", overwrite=TRUE)
  # Save initial kfold partition for presence data
  write.table(data.frame(NObskfoldgrpp), file=paste0(Species, "ProjPresenceDat", NObsJ, "Kfoldi", NObskfoldgrppL[count], ".csv"), sep=",")
}


#########################################################################################

#####################################################################################################
### Loop through individual or combined years for roost data and create training subset ensembles of
### roost kde models
#####################################################################################################
#
# Specify final model variable
FinalModelVariables <- 1
#FinalModelVariables <- 3
NObsJ <- 3
NumberModSets <- 1
NumProjections <- 3  # Specify number of models to project
#NumberModSets <- 30
#NumProjections <- 15
################################################################################
### Create Kernel Density Estimate surfaces for 2/3 training data for projections and calibration
########################################################################
### Create and save NObsJ usually three kfold partition scheme for multiple rounds of test with presence
### pseudoabsence and background points using all data
### NOTE: Only create and save partition scheme once in order to reuse
###### Read in presence data for kfold assignment
#YearL <- seq(2005,2016,1)
YearL <- c("02_16")
for(Year in YearL) {
  # Year=2005
  setwd(InDirect)
  PresThinShp <- readOGR(".", paste0("MonRst", Year, "denpop"))  # NOTE: readShapePoly was giving error, so used
  summary(PresThinShp)
  #Set projection
  CRS.NAAlbersEqualAreaConic <- CRS(proj4string(PresThinShp))
  head(PresThinShp)
  nrow(PresThinShp)
  #plot(PresThinShp, col='blue')
  #
  CRS.In <- CRS.NAAlbersEqualAreaConic
  #
  # Read PresThinShp to spatial points data frame
  PresThin.spdf <- readShapePoints(paste0("MonRst", Year, "denpop.shp"), proj4string=CRS.NAAlbersEqualAreaConic)
  head(PresThin.spdf)
  # Convert to data frame
  PresThin.df <- data.frame(PresThin.spdf)
  head(PresThin.df)
  ## Thin out and reorder columns
  if(Year=="02_16") {
    PresThin.df <- PresThin.df[c(3,4,1)]
    colnames(PresThin.df) <- c("x", "y", "denpop")
  } else {
    PresThin.df <- PresThin.df[c(9,10,6)]
    colnames(PresThin.df) <- c("x", "y", "denpop")
  }
  nrow(PresThin.df)
  head(PresThin.df)
  ##
  library(dismo)
  #
  ################################################################################
  ### Create NObs kfold group for presence points and use to define training KDEs
  ### NOTE: Only create NObs kfold groups once
  NObskfoldgrpp <- kfold(PresThin.df, NObsJ)
  # Save initial kfold partition for presence data
  write.table(data.frame(NObskfoldgrpp), file=paste0(Species, Year, "_ProjPresenceDat", NObsJ, "Kfold", ".csv"), sep=",")
  #########
  ## For the three sets of training NObskfoldgrpp sets above, convert kfold presence training data to
  ## shapefiles below and process with python "MonarchRoostKFoldDataProcessingBatch.py" to create
  ## kernel density estimate surfaces of chosen bandwidth for the training data
  ##
  NObskfoldgrppL <- c("A", "B", "C")
  GroupFirst <- c(1,2,3)
  GroupSecond <- c(2,3,1)
  count = 0
  for(Group in NObskfoldgrppL) {
    count = count + 1
    PresThinProjTrain.df1 <- PresThin.df[NObskfoldgrpp==GroupFirst[count],]
    PresThinProjTrain.df2 <- PresThin.df[NObskfoldgrpp==GroupSecond[count],]
    PresThinProjTrain.df <- rbind(PresThinProjTrain.df1, PresThinProjTrain.df2)
    nrow(PresThinProjTrain.df)
    head(PresThinProjTrain.df)
    xy <- PresThinProjTrain.df[,1:2]
    PresThinProjTrain.spdf <- SpatialPointsDataFrame(coords = xy, data = PresThinProjTrain.df,
                                   proj4string = CRS.NAAlbersEqualAreaConic)
    head(PresThinProjTrain.spdf)
    nrow(PresThinProjTrain.spdf)
    #
    # Write to shapefile
    writeOGR(PresThinProjTrain.spdf, InDirect, paste0(SpecFileName, Year, "_PresProjTrain", Group), driver = "ESRI Shapefile", overwrite=TRUE)
  }
}
################################################################################


############################################################################
### Process above kfold training data with "MonarchRoostKfoldDataProjection_KDEProcessingBatchPopDenInd.py"
### to create three training KDEs
#### NOTE  Default bandwidth was 234,698 m for "PROJRSTKDE1BW" (8/2/17)
###################################################################################


###################################################################################
### Generate KDE models (KDEMs) and obtain evaluation statistics from each training
### KDE
#################################################################################
# First read in presence and background data points
setwd(InDirect)
PresThinPntShp <- readOGR(".", paste0("rstdenpopall"))
nrow(PresThinPntShp)
#Set projection
CRS.NAAlbersEqualAreaConic <- CRS(proj4string(PresThinPntShp))
TotPres <- nrow(PresThinPntShp)
#
PseudoabsencePntShp <- readOGR(".", paste0(SpecFileName, "psa"))  # NOTE: readShapePoly was giving error, so used readOGR instead
#str(PseudoabsencePntShp)
#
BackgroundPntShp <- readOGR(".", paste0(SpecFileName, "bkgrnd"))
################################################################################
### Create NObs kfold group for pseudoabsence points for testing KDEMs
### NOTE: Only create NObs kfold groups once
NObskfoldgrpa <- kfold(PseudoabsencePntShp, NObsJ)
head(PseudoabsencePntShp)
if(nrow(PseudoabsencePntShp)>MaxBackgrndPseudoabsLimit) {
  NObskfoldgrpa <- NObskfoldgrpa[1:MaxBackgrndPseudoabsLimit ]
}
length(NObskfoldgrpa)
# Save initial kfold partition for presence data
write.table(data.frame(NObskfoldgrpa), file=paste0(Species, "ProjPseudoabsenceDat", NObsJ, "Kfold.csv"), sep=",")
##########################################################################################
### Loop through three training KDE models for each year and extract presence,
### background and pseudoabsence data for model calibration and calibrate models
### while calculating statistics of AUC, AICc etc.
###
setwd(InDirect)
#YearL <- seq(2005,2016,1)
YearL <- c("02_16")
#YearL <- c(2005)
NumberModSets <- 1
NObskfoldgrppL <- c("A", "B", "C")
TestObskfoldgrp <- c(3,1,2)
BigLoopCount <- seq(1,3,1)
ModelType <- "KDE"
BandwidthType <- "" #1 means 110% of default bandwidth
CRS.In <- CRS.NAAlbersEqualAreaConic
Count <- 0
KDEMKeepEvalStatsL <- list()
##################################################################################
t1 <- Sys.time()
for(Year in YearL) {
# Year="02_16"
  for(BigLoop in BigLoopCount) {
    #BigLoop=1
    #Count=1
    ## Read in KDE Training Model
    setwd(GridDirect)
    SetNameF <- paste0("TrainSet", NObskfoldgrppL[BigLoop])
    KDETrainModel <- raster(paste0("rst", BandwidthType, "kde", Year, NObskfoldgrppL[BigLoop]))
    #plot(KDETrainModel)
    #
    #########################################################
    ### For three different NObskfoldgrpp groups use the different
    ### group training KDEs to generate different models of a
    ### Training Subset Ensemble
    #########################################################
    ## Include training KDE raster with original environmental rasters
    #############################################################
    ##########################
    ## Loop through sets 1 through 10
    #ModelTypeNames <- c(paste0("Random", FullNObs))
    Count <- Count + 1
    #
    ######################################################################################
    ## Run Threshold Calibration Projection for Top Selected Variable Model for Given Variable Set (j)
    ######################################################################################
    ##
    #################
    ## NOTE: DO NOT RESORT BY THIS AUCfinaltest- THIS WILL INFLATE AUC- AUCfiltertest already used to sort
    #################
    VariableSubset <- data.frame(paste0("rst", BandwidthType, "kde", Year, NObskfoldgrppL[BigLoop]), stringsAsFactors=FALSE)
    colnames(VariableSubset) <- "VarNames"
    #str(VariableSubset)
    VarNames <- unlist(strsplit(as.character(VariableSubset[,1]),"-"))
    SubsetVariableNumber <- length(VarNames)
    Rep <- BigLoop
    PredictorsProj <- stack(KDETrainModel)
    # Add KDE to VariableNames
    VariableNames <- VariableSubset
    ######## Read back in kfold partition scheme for year ##################
    setwd(InDirect)
    NObskfoldgrpp <- as.vector(unlist(read.csv(paste0(Species, Year, "_ProjPresenceDat", NObsJ, "Kfold", ".csv"))))  # Read in kfolgrpp with NAs omitted
    length(NObskfoldgrpp)
    ####
    kfoldgrppin <- NObskfoldgrpp
    kfoldgrpain <- NObskfoldgrpa
    #
    ########################################################
    Run=TestObskfoldgrp[BigLoop]
    ##
    VariableNamesIn <- VariableNames
    #
    ######################
    #### If producing KDE models from single KDE variable:
    setwd(OutDirectpsa)
    output3 <- paste0("BW", BandwidthType, "_KDEMProj_", Year, "_Subsetof", NumberModSets, "_", SetNameF, "_", SubsetVariableNumber, "Vars_ThresholdCalib", BigLoop, NObskfoldgrppL[BigLoop])
    #output3 <- paste0("MaxentProj_", ModelType, "_Subsetof", NumberModSets, "_", SetNameF, SubsetVariableNumber, "_ThresholdCalib", Loop, "Minus", 14)
    dir.create(paste0(OutDirectpsa, "/", output3))
    OutDirectSub <- paste0(OutDirectpsa, "/", output3)
    OutDirectIn <- OutDirectSub
    ##
    OutGridID <- paste0(ModelType, "_", BigLoop, "ThreshCal_", SetNameF)
    # Takes 19 minutes for 25 variable grid
    #
    SetNameIn <- Year
    PredictorsIn <- PredictorsProj
    #
    Time <- system.time(KDEMKeepEvalStats.df <- KDEMSubset.GridTrainTestEvalAIC_Calib(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetNameIn, VariableSubset, PresThinPntShp, PseudoabsencePntShp, BackgroundPntShp, kfoldgrppin, kfoldgrpain, CRS.In, OutDirectIn, Output=TRUE))
    #
    head(KDEMKeepEvalStats.df)
    # Save results
    setwd(OutDirectIn)
    write.table(KDEMKeepEvalStats.df , file=paste0(Species, "BW", BandwidthType, "_KDEM", "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SubsetVariableNumber, "_", BigLoop, NObskfoldgrppL[BigLoop], ".csv"), sep=",", col.names=NA)
    #
    #KDEMod <- readRDS("KDEModel.rds") # Reads in stored KDE model from above function
    ## Add ModelType to output
    KDEMKeepEvalStats.df$ModelType <- ModelType
    ncol(KDEMKeepEvalStats.df)
    # Re-order columns
    KDEMKeepEvalStats.df <- KDEMKeepEvalStats.df[,c(1:2, 16, 3:15)]
    # Replace Run with Rep
    KDEMKeepEvalStats.df$Rep <- Rep
    ncol(KDEMKeepEvalStats.df)
    KDEMKeepEvalStats.df <- KDEMKeepEvalStats.df[,c(1:6, 17, 8:16)]
    KDEMKeepEvalStatsL[[Count]] <-  KDEMKeepEvalStats.df
  }
}
###
t2 <- Sys.time()
#############################################################################
difftime(t2,t1, units = "mins")
###
KDEMKeepEvalStatsAll.df <- do.call(rbind, KDEMKeepEvalStatsL)
## Sort by DataType
KDEMKeepEvalStatsAll.dfs  <- arrange(KDEMKeepEvalStatsAll.df, DataType)
## Save data
setwd(OutDirectpsa)
write.table(KDEMKeepEvalStatsAll.dfs, file=paste0(Species, "BW", BandwidthType, "_KDEMbyYears", "_TrainVsTestStatsThreshCal_of", NumberModSets, "_", SubsetVariableNumber, "Vars_", NumProjections, "Sets.csv"), sep=",", col.names=NA)






############################################################################################
### NOTE: Code below is for processing and analyzing annual KDEMs
#############################################################################################

################################################################################
### Use python scripts listed below to get data on displacement of annual centroids
### of KDEs between 27N and 37N in comparison to Training Subset Ensemble centroid of all years
### Also obtain average width in km of this area of KDEs
###  MonarchRoostYearlyKfoldDataProjectionAllYearTSE_KDECentroid.py
###  MonarchRoostYearlyKfoldDataProjectionTSE_KDECentroidPathProcessingBatch.py
################################################################################
##
## Read in output from above python code
KDEDirect <- "G:/MonarchRoost/MonRstYearlyKDE"
setwd(KDEDirect)
#
KDEDimensions.df <- data.frame(read.csv("KDEWidth_CentroidShift_Distances.csv"), stringsAsFactors=FALSE)

## Read back in KDE stats data
setwd(OutDirectpsa)
#
###Calculate Mean and Standard Deviation Values
EvaluationStats <- KDEDimensions.df
OutName <- paste0(Species, "_AnnualKDEDimensionsSummary.csv")
SortGroups <- c("KDEYEAR")
## All statistics, and only statistics, should be at and after the column specified below
StatVarFirstColumn <- 3
OutDirectIn <- OutDirectpsa
EvaluationStatsSummary.df <- EvalStatVars.Summary(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn)
head(EvaluationStatsSummary.df)
##
EvalStatsOut <- EvaluationStatsSummary.df
#
write.table(EvalStatsOut, file=paste0(Species, "_AnnualKDEDimensionsSummaryTable.csv"), sep=",", col.names=NA)
##















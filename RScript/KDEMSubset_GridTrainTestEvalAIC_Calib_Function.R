KDEMSubset.GridTrainTestEvalAIC_Calib <-
function(Species, VariableNamesIn, PredictorsIn, SubsetVariableNumber, TotPres, Run, SetName, VariableSubset, PresThinPntShp, PseudoabsencePntShp, BackgroundPntShp, kfoldgrppin, kfoldgrpain, CRS.In, OutDirectIn..., CatVarsPrefix, Output) {
  #SetName=SetNameIn
  #NOTE: Function parameters after the "...", such as TestDataType, have to be set with an equal sign in the function call, such as ScoreType="mean"
  setwd(OutDirectIn)
  library(dismo)
  library(raster)
  library(ENMeval)
  #
  # Make sure VariableNames is a data frame
  VariableNamesIn <- data.frame(VariableNamesIn, stringsAsFactors=FALSE)
  #
  RowNames.df <- data.frame(rownames(VariableSubset), stringsAsFactors=FALSE)  # Save row names to use in output
  SubsetSize.mat <- matrix(c("Singlets", "Doublets", "Triplets", "Quartets", "Quintets", "Sextets", "Septets", "Octets", "Nonets",
  "Dectets", "Undectets", "Duodectets","Tredectets", "Quattuordectets", "Quindectets", "Sexdectets", "Septendectets", "Octodectets", "Novemdectets",
  "Vigetets", "Unvigetets", "Duovigetets", "Trevigetets", "Quattuorvigetets", "Quinvigetets", "Sexvigetets", "Septenvigetets", "Octovigetet",
  "Novemvigetets", "Trigetets", "Untrigetets", "Duotrigetets", "Tretrigetets", "Quottuortrigetets", "Quintrigetets",
  "Sextrigetets", "Septentrigetets", "Octotrigetets", "Novemtrigetets", "Quadragetets", "Unquadragetets", "Duoquadragetets", "Trequadragetets",
  "Quattuorquadragetets", "Quinquadragetets", "Sexquadragetets", "Octoquadragetets", "Octoquadragetets", "Novemquadragetets", "Quinquagetets",
  "Unquinquagetets", "Duoquinquagetets", "Trequinguagetets", "Quattuorquinquagetets", "Quinquinquagetets",
  "Sexquinquagetets", "Septenquinquagetets", "Octoquinquagetets", "Novemquinquagetets", "Sexagetets"), ncol=1, nrow=60, byrow=TRUE, dimnames=list(c
   (seq(1:60)), c("Subset")))
  SubsetSize.df <- as.data.frame(SubsetSize.mat, stringsAsFactors=FALSE)
  Subset <- SubsetSize.df[SubsetVariableNumber,]
  if(is.na(Subset)) {
    Subset <- ""
  }
  TotVars <- nrow(VariableNamesIn)
  #
  if(missing(Output)) { Output=TRUE }
  #
  ModelType <- paste("KDEM", Subset, sep="")
  Model <- paste("KDEM for Training Subset")
  tail(VariableSubset)
  #VariableSubsets[2998,]
  ################
  #i=1
  ## Use function to round up from .5 from http://stackoverflow.com/questions/12688717/round-up-from-5-in-r
  Round2 <- function(x, n) {
    posneg = sign(x)
    z = abs(x)*10^n
    z = z + 0.5
    z = trunc(z)
    z = z/10^n
    z*posneg
  }
  #
  # Split VarNames of VariableSubsets into separate variables
  #str(VariableSubsets)
  VariableNamesSel <- c(unlist(VariableSubset[1,1]))
  # Check if dash used to separate variables
  if(grepl("-", VariableNamesSel)==TRUE) {
    VarNames <- unlist(strsplit(VariableNamesSel, "-"))
  } else {
    VarNames <- unlist(VariableNamesSel)
  }
  SubsetVarNum <- length(VarNames)
  ##
  ############ Convert KDE Environmental Layer to Niche Model
  ############ by normalizing from 0 to 1 with floating value
  ## Identify KDE raster
  KDE1 <- subset(PredictorsIn, VarNames)
  #names(PredictorsIn)
  ## Find minimum value for KDE1 raster
  MinVal <- minValue(KDE1)
  ## Subtract minimum value from all raster values
  KDE2 <- KDE1 - MinVal
  ## Find maximum value for KDE2 raster
  MaxVal <- maxValue(KDE2)
  KDEM.score <- KDE2/MaxVal
  #plot(KDE)
  # Multiply continuous model grid by 1000 and convert to integer
  KDEM.score1 <- calc(KDEM.score, function(x) as.integer(x * 1000) )
  # Save grid
  setwd(OutDirectIn)
  writeRaster(KDEM.score1, paste0(Species, "_", OutGridID, SubsetVariableNumber), format = "GTiff", overwrite=TRUE)
  #
  #####################################
  ## Calculate AICc_bg with point values
  # Designate number of parameters in KDE model
  nparams <- length(VarNames)
  ## From ENMeval Package documentation: AICc is the Akaike Information Criterion corrected for small
  ## sample sizes calculated as: (2 * K - 2 * logLikelihood) + (2 * K) * (K + 1)=(n - K - 1)
  ## where K is the number of parameters in the model (i.e., number of non-zero parameters in Maxent
  ## lambda file) and n is the number of occurrence localities.
  ## The logLikelihood is sum(log(vals/total))
  ## vals is vector of Maxent raw values at occurence localities
  ## total is the sum of Maxent raw values across the entire study area
  ##
  ### Create a version of KDE raster where all values sum to one like in MaxEnt raw raster
  ##
  ## Sum all values in the normalized raster
  SumVal <- cellStats(KDEM.score, stat='sum', na.rm=TRUE)
  ## Divide all values in raster by SumVal
  SumValDivFun <- function(x) {x/SumVal}
  KDEM.raw <- calc(KDEM.score, SumValDivFun)
  # Replace values of 0 with 0.000001
  KDEM.raw[KDEM.raw == 0] <- 0.000001
  # Check
  #SumCheck <- cellStats(KDEM.raw, stat='sum', na.rm=TRUE)
  ### Extract KDEprob values for presence training data
  ## Find Training Presence Values from original KDE raster from PresenceDat.df
  # Specify model training data
  PresThinPntShpK <- PresThinPntShp[kfoldgrppin == Run, ]
  head(PresThinPntShpK)
  ncol(PresThinPntShpK)
  nrow(PresThinPntShpK)
  # Convert point data.frame to SpatialPointsDataFrame
  ## Extract values of presence test points from KDE model grid
  system.time(prestestraw.df <- data.frame(extract(KDEM.raw, PresThinPntShpK)))
  colnames(prestestraw.df) <- "KDEMRaw"
  head(prestestraw.df)
  # Omit any rows with NA values
  prestestraw.df <- na.omit(prestestraw.df)
  # Keep presence values as vector
  vals <- prestestraw.df[,1]
  head(vals)
  n <- length(vals) # number of occurrence localities
  # total is sum of values across entire study area, includes background and occurrence localities
  ## Obtain background values for KDEM.raw
  ## Extract values of presence test points from KDE model grid
  system.time(backgroundraw.df <- data.frame(extract(KDEM.raw, BackgroundPntShp)))
  colnames(backgroundraw.df) <- "KDEMRaw"
  head(backgroundraw.df)
  # Omit any rows with NA values
  backgroundraw.df <- na.omit(backgroundraw.df)
  # Keep values as vector
  backgroundvals <-  backgroundraw.df[,1]
  #
  # Calculate sum of all values
  totalocc <- sum(vals) # sum from occurrence localities
  totalbg <- sum(backgroundvals)  # sum from background localities
  total <- totalocc + totalbg  # grand total sum
  #
  logLikelihood <- sum(log(vals/total))
  K <- nparams
  #K=5
  AICc_bg <- (2*K - 2*logLikelihood) + (2*K)*(K+1)/(n-K-1)
  NumDVars <- K
  ###
  ####################################################
  #################################################
  ### Calculate AICc with ENMeval package
  # Obtain number of parameters in Maxent model
  nparam <- K
  # Derive occurence data
  ## Extract coordinates for later AICc calculation
  proj4string(PresThinPntShpK)=CRS.In
  #
  KDEMPresTrainData.df <- coordinates(PresThinPntShpK)
  head(KDEMPresTrainData.df)
  KDEMPresTrainData.csv <- write.csv(KDEMPresTrainData.df, file="KDEMPresTrainData.csv", row.names=FALSE)
  occ <- read.table("KDEMPresTrainData.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
  #
  KDEMAICcResults.df <- calc.aicc(nparam, occ, KDEM.raw)
  AICc <- KDEMAICcResults.df$AICc
  ####################################################
  ########### Evaluate model using training and test data
  KDEMSubsetEvalL <- list()
  TestDataTypes <- c("Train", "Test")
  for(k in 1:2) {
    #k=2
    TestDataType <- TestDataTypes[k]
    if(k==1) {
    # Specify model testing data
      KDEMPresTestDataShp <- PresThinPntShp[kfoldgrppin != Run, ]
      head(KDEMPresTestDataShp)
      nrow(KDEMPresTestDataShp)
      KDEMAbsTestDataShp <- PseudoabsencePntShp[kfoldgrpain != Run, ]
      head(KDEMAbsTestDataShp)
      nrow(KDEMAbsTestDataShp)
      tail(KDEMAbsTestDataShp)
    } else {
      KDEMPresTestDataShp <- PresThinPntShp[kfoldgrppin == Run, ]
      head(KDEMPresTestDataShp)
      nrow(KDEMPresTestDataShp)
      tail(KDEMPresTestDataShp)
      KDEMAbsTestDataShp <- PseudoabsencePntShp[kfoldgrpain == Run, ]
      head(KDEMAbsTestDataShp)
      tail(KDEMAbsTestDataShp)
      nrow(KDEMAbsTestDataShp)
    }
    ########################
    ### Query projection grid to get values for test presence and absence points
    ## First presence points
    setwd(OutDirectIn)
    head(KDEMPresTestDataShp)
    # First, define xy coordinates
    xy <- coordinates(KDEMPresTestDataShp)
    ## Extract values of presence test points from KDE model grid
    system.time(prespred.df <- data.frame(extract(KDEM.score, KDEMPresTestDataShp)))
    colnames(prespred.df) <- "KDEMScore"
    # Rejoin coordinates to prespred.df and save as shapefile for checking
    prespred.spdf <- SpatialPointsDataFrame(coords=xy, data=prespred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(prespred.spdf, ".", paste0(Species, "Presence", TestDataType, "KDEMVals"), driver="ESRI Shapefile", overwrite=TRUE)
    # Omit any rows with NA values
    prespred.df <- na.omit(prespred.df)
    ## Then absence points
    head(KDEMAbsTestDataShp)
    # First, define xy coordinates
    xy <- coordinates(KDEMAbsTestDataShp)
    ## Extract values of presence test points from KDE model grid
    system.time(abspred.df <- data.frame(extract(KDEM.score, KDEMAbsTestDataShp)))
    colnames(abspred.df) <- "KDEMScore"
    # Rejoin coordinates to abspred.df and save as shapefile for checking
    abspred.spdf <- SpatialPointsDataFrame(coords=xy, data=abspred.df, proj4string=CRS.In)
    # Write shapefile including the correct projection
    writeOGR(abspred.spdf, ".", paste0(Species, "Pseudoabsence", TestDataType, "KDEMVals"), driver="ESRI Shapefile", overwrite=TRUE)
    # Omit any rows with NA values
    abspred.df <- na.omit(abspred.df)
    #############################################################################################
    # This section evaluates the KDE model using the PresenceAbsence package
    #############################################################################################
    #### Create a dataset with model predictions for presence and absence points
    # Use extracted prediction values for presence and absence points for each of three
    # models previously calculated in loop
    #
    library(gtools)
    ## Create directory of output for class pair run
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    # For presence data, assign a column of "1" under the name "OBSERVED" to indicate presence data
    # and assign the model results a name "KDEMScoreN" where N is the name of the rep
    # Also assign a column "id" for the row numbers to use in merging the data frames later
    presa.df <- data.frame(c(rep(1, nrow(prespred.df))))
    names(prespred.df) <- c("KDEM")
    names(presa.df) <- c("OBSERVED")
    pres.df <- data.frame(cbind(id=1:nrow(presa.df), presa.df, prespred.df))
    nrow(pres.df)
    # Repeat above process with absence data, but assign "OBSERVED" a value of 0
    absa.df <- data.frame(c(rep(0, nrow(abspred.df))))
    names(abspred.df) <- c("KDEM")
    names(absa.df) <- c("OBSERVED")
    abs.df <- data.frame(cbind(id=1:nrow(absa.df), absa.df, abspred.df))
    # For each model output, merge presence and absence data using "id' column as guide when all=TRUE
    # NOTE: PresenceAbsence package cannot handle several models at one time if the sample sizes differ
    # so have to analyze each model output separately
    presabspred <- rbind(pres.df, abs.df)
    tail(presabspred)
    head(presabspred)
    # Drop the id column used in merging for each dataset
    presabspred$id <- NULL
    # Make a column of data with the species name with same number of rows as data from each model
    SPECIES <- data.frame(c(rep(Species, nrow(presabspred))))
    names(SPECIES) <- c("SPECIES")
    # Make final dataset SPDATA by putting together SPECIES with extracted environmental data.
    SPDATA <- data.frame(SPECIES, presabspred)
    head(SPDATA)
    #SPDATA[100:160,]
    ################################################################################
    ### Run this block of code to evaluate model results with PresenceAbsence package
    ################################################################################
    library(PresenceAbsence)
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    #starttime <- Sys.time()
    #### FOR OLD WORLD DATA EVALUATION STATISTICS
    ### Define variables for later use.
    accurun <- list()
    accusum <- matrix(data=NA, ncol=9, nrow=1, byrow=TRUE, dimnames=list(NULL, c("MaxTSS", "Specificity_maxTSS", "Sensitivity_maxTSS", "AUC", "MaxKappa", "ThresholdMaxTSS", "AICc_bg", "NumDVars", "AICc")))
    species <- as.character(unique(SPDATA$SPECIES))
    model.names <- as.character(names(SPDATA)[-c(1, 2)])
    N.models <- ncol(SPDATA) - 2
    N.sp <- length(species)
    N.obs <- length(SPDATA$SPECIES[SPDATA$SPECIES == species[1]])
    Obs.prev <- table(SPDATA$SPECIES, SPDATA$OBSERVED)[, 2]/N.obs
    Obs.prev <- Round2(Obs.prev, 2)
    ### Mainly just run this code
    graphics.off()
    sp <- 1
    # Read in dataset for loop
    DATA <- SPDATA[SPDATA$SPECIES == species[sp], ]
    head(DATA)
    #
    # To assess accuracy per threshold, use limited threshold available for
    #  model based upon number of environmental layers in model
    # ("NumGrids")
    #NumGrids <- max(40, SubsetVarNum)
    #PossThresholds <- seq(1/NumGrids,1,length=NumGrids)
    PossThresholds <- 100
    #accu <- data.frame(presence.absence.accuracy(SPDATA, which.model = 1, threshold = PossThresholds, st.dev=FALSE))
    # accu <- presence.absence.accuracy(DATA, which.model = 1, threshold = c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.975, 0.98, 0.99, 0.999999))
    #accu <- presence.absence.accuracy(DATA, which.model = 1, threshold = 100, st.dev=FALSE)
    # print(paste("Species:", species[sp], "Model:", model.names))  not used
    accu <- data.frame(presence.absence.accuracy(DATA, which.model = 1, threshold = 100, st.dev=FALSE))
    # print(paste("Species:", species[sp], "Model:", model.names))  not used
    head(accu)
    maxSSS <-  data.frame(accu$sensitivity + accu$specificity)
    names(maxSSS) <- c("maxSSS")
    head(maxSSS)
    TSS <-  data.frame(accu$sensitivity + accu$specificity - 1)
    names(TSS) <- c("TSS")
    accurun <- data.frame(accu, maxSSS, TSS)
    head(accurun)
    accurun$Conditions <- paste("KDEM", Subset, sep="")
    maxKappa <- max(accurun$Kappa)
    maxTSS <- max(accurun$TSS)
    AUC <- max(accurun$AUC)
    # Find and average thresholds at TSS = maxTSS. In the case of tied optimal
    # thresholds, we select the mean threshold producing maximum TSS following
    # (Freeman and Moisen 2008). But, in the case of discrete thresholds as found
    # in envelope models, if the mean optimal threshold does not represent an
    # actual discrete threshold, we select the nearest discrete threshold to the
    # mean among 3 or more thresholds, or the smaller of two adjacent discrete thresholds.
    ThresholdsMaxTSS <- accurun$threshold[which(accurun$TSS == maxTSS)]
    ThresholdMaxTSS <- mean(ThresholdsMaxTSS)
    #ThresholdMaxTSSM <- mean(ThresholdsMaxTSS)
    ## Following commented code for envelope score
    #if (length(ThresholdsMaxTSS) < 3) {
    #ThresholdMaxTSS <- min(ThresholdsMaxTSS)
    #} else {  ThresholdMaxTSS <- PossThresholds[which(abs(PossThresholds - ThresholdMaxTSSM)== min(abs(PossThresholds - ThresholdMaxTSSM)))]
    #}
    # Calculate specificity and sensitivity at maxTSS
    Specificity_maxTSS <- accurun$specificity[which(accurun$TSS == maxTSS)]
    Specificity_maxTSS <- mean(Specificity_maxTSS)
    Sensitivity_maxTSS <- accurun$sensitivity[which(accurun$TSS == maxTSS)]
    Sensitivity_maxTSS <- mean(Sensitivity_maxTSS)
    #CheckTSS <- Specificity_maxTSS + Sensitivity_maxTSS - 1 # should equal maxTSS
    accusum[1,1] <- max(maxTSS, 0.0001)
    accusum[1,2] <- max(Specificity_maxTSS, 0.0001)
    accusum[1,3] <- max(Sensitivity_maxTSS, 0.0001)
    accusum[1,4] <- max(AUC, 0.0001)
    accusum[1,5] <- max(maxKappa, 0.0001)
    accusum[1,6] <- max(mean(ThresholdMaxTSS), 0.0001)
    accusum[1,7] <- AICc_bg
    accusum[1,8] <- NumDVars
    accusum[1,9] <- AICc
    # Save Threshold value and multiply by 1000 for grid calibration
    ThresholdK <-  (max(ThresholdMaxTSS, 0.0001))*1000
    ###############################
    #endtime <- Sys.time()
    #durtime <- endtime - starttime
    # Save evaluation statistics to .csv file
    #setwd(paste("C:/Users/JLTracy/Documents/R/win-library/3.0/10minClimIntEvalSTB1000", "/", output2, sep=""))
    accusum.df <- data.frame(accusum)
    accusum.df$Model <- ModelType
    #
    # Save variable matrix coordinates and Old World evaluation statistics to a vector
    KDEMEvalStats.df <- data.frame(t(as.matrix(c(VariableNamesSel, SubsetVariableNumber, TestDataType, TotPres, SetName, Run, accusum.df$MaxTSS, accusum.df$MaxKappa, accusum.df$AUC, accusum.df$Specificity_maxTSS, accusum.df$Sensitivity_maxTSS, accusum.df$ThresholdMaxTSS, accusum.df$AICc_bg, accusum.df$NumDVars, accusum.df$AICc))))
    #
    setwd(OutDirectIn)
    # Save threshold verus evaluation statitistic data for run
    write.table(SPDATA, file=paste(Species, TestDataType, "SPDATA_Run", Run, ".csv", sep=""), sep=",")
    write.table(accurun, file=paste(Species, TestDataType, "StatsVsThresholds_Run", Run, ".csv", sep=""), sep=",", col.names=NA)
    #print(paste("Species:", species[sp], "Model:", model.names))
    #print(accurun)
    bmp(paste(Species, TestDataType,"TSSvsThreshold_Run", Run, ".bmp", sep=""))
    #dev.new()
    plot(accurun$threshold, accurun$TSS, type="l")
    dev.off()
    # Plot a ROC plot
    # save following plot of grid as .bmp file
    bmp(paste(Species, TestDataType, "KDEMROC_Run", Run, ".bmp", sep=""))
    #dev.new()
    auc.roc.plot(DATA, color = TRUE, legend.cex = 1.2, main = "")
    dev.off()
    # Calculate optimal thresholds by various methods and save output to text file
    #outthresholds <- capture.output(optimal.thresholds(DATA, opt.methods = 1:12, threshold=10001))
    #out1 <- paste(outthresholds)
    #cat(out1, file=paste(Species, "KDEResults", TestDataType, "Thresholds.txt", sep=""), sep="\n", append=TRUE)
    # Can calculate the confusion matrix for a given threshold (not necessary)
    confmatrix <- cmx(DATA, threshold = ThresholdMaxTSS, na.rm = FALSE)
    # Save confusion matrix at threshold of maxTSS
    write.table(confmatrix, file=paste(Species, "KDEMResults", TestDataType, "ConfusionMatrixatMaxTSSThreshold_Run", Run, ".csv", sep=""), sep=",", col.names=c("Obs Pres","Obs Abs"), row.names=c("Pred Pres", "Pred Abs"))
    #Save evaluation stastics to file
    write.table(accusum.df, file=paste(Species, "KDEMResults", TestDataType, "Stats_Run", Run, ".csv", sep=""), sep=",", col.names=NA)
    graphics.off()
    # Save output evaluation statistics and summary statistics to text file
    # Output codes of variables used
    outA <- paste("\n\nSpecies:", Species, "\nVariable Subset:", VariableNamesSel, "\nModel:", Model, "\nScoring Alogorithm: KDE")
    cat(outA, file=paste(Species, "KDEMResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #NotUsedVars <- capture.output(print(DeleteBands))
    #outB <- paste("\nVariables Not Used: Code-", NotUsedVars, "; Variable-", GridNamesDrop)
    #cat(outB, file=paste(Species, InfClass, "KDEResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #UsedVars <- capture.output(print(SRandNumbers))
    #outB <- paste("\n  Code Nos. for Variables Used:", UsedVars)
    #cat(outB, file=paste("KDEMResults", output2, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #out1 <- paste("\nVariable Names")
    #cat(out1, file=paste(Species, "KDEMResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    #EnvVariables <- as.matrix(GridNamesKeepSel)
    #out2 <- capture.output(print(EnvVariables))
    #cat(out2, file=paste(Species, "KDEMResults", GridNamesKeepSelT, "Summary.txt", sep=""), sep="\n", append=TRUE)
    out3 <- paste("\n  Evaluation Statistics per Run")
    cat(out3, file=paste(Species, "KDEMResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    out4 <- capture.output(print(accusum.df))
    cat(out4, file=paste(Species, "KDEMResults", TestDataType, "Summary_Run", Run, ".txt", sep=""), sep="\n", append=TRUE)
    #
    #############
    KDEMSubsetEvalStats.df1 <- KDEMEvalStats.df
    head(KDEMSubsetEvalStats.df1)
    tail(KDEMSubsetEvalStats.df1)
    nrow(KDEMSubsetEvalStats.df1)
    #str(ESWrapperOutput.df)
    colnames(KDEMSubsetEvalStats.df1) <- c("VarNames", "SubsetVariableNumber", "DataType", "TotalPresPnts", "SetName", "Run", "TSS", "Kappa", "AUC", "Spec", "Sens", "ThreshMxTSS", "AICc_bg", "NumDVars", "AICc")
    rownames(KDEMSubsetEvalStats.df1) <- c(seq(1:nrow(KDEMEvalStats.df)))
    KDEMSubsetEvalStats.df1 <- as.data.frame(KDEMSubsetEvalStats.df1, stingsAsFactors=FALSE)
    head(KDEMSubsetEvalStats.df1)
    ncol(KDEMSubsetEvalStats.df1)
    # Save output matrix
    # Convert sixth through 15th columns from character to numeric
    KDEMSubsetEvalStats.df1[,6:15] <- sapply(KDEMSubsetEvalStats.df1[,6:15], function(x) as.numeric(as.character(x)))
    KDEMSubsetEvalL[[k]] <- KDEMSubsetEvalStats.df1
    #Calculate difference between training and test statistics (overfitting)
    if(k==2) {
      KDEMSubsetEvalL[[3]] <- KDEMSubsetEvalL[[2]]
      KDEMSubsetEvalL[[3]][,7:15] <- KDEMSubsetEvalL[[1]][,7:15] - KDEMSubsetEvalL[[2]][,7:15]
      KDEMSubsetEvalL[[3]][,3] <- "Diff"
    }
  }
  #####################################################################################
  KDEMSubsetEvalStats.df <- do.call(rbind, KDEMSubsetEvalL)
  #
  ################################################################################
  ### Calibrate model using using above calculated ThresholdK
  ### at maximum TSS from above evaluation 
  KDEM.scorecal <- calc(KDEM.score1, function(x) ifelse(x < ThresholdK, 0, 1) )
  #plot(KDE.scorecal)
  # SubsetVariableNumber <- 19
  #OutGridID <- "Full"
  writeRaster(KDEM.scorecal, paste0(Species, OutGridID, SubsetVariableNumber, "Cal"), format = "GTiff", overwrite=TRUE)
  ##############
  #
  if(Output==TRUE) {
    setwd(OutDirectIn)
    Sets <- nrow(KDEMSubsetEvalStats.df)
    write.table(KDEMSubsetEvalStats.df, file=paste(Species, "KDEM_TrainTest_", TotVars, "_TotVars_", SubsetVariableNumber, "Var_", Sets, Subset, "_", SetName, Run, ".csv", sep=""), sep=",", col.names=NA)
  }
  #
  gc()
  return(KDEMSubsetEvalStats.df)
}

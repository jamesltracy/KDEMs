EvalStatVars.Summary <-
function(EvaluationStats, SortGroups, StatVarFirstColumn, OutName, OutDirectIn..., OrigFormat) {
library(plyr)
library(dplyr)
  if(missing(OrigFormat)) { OrigFormat=FALSE }
  #  
  StatVars <- colnames(EvaluationStats)[StatVarFirstColumn:ncol(EvaluationStats)]
  rownames(EvaluationStats) <- seq(1,nrow(EvaluationStats),1)
  head(EvaluationStats)
  Freq <- count_(EvaluationStats, vars=SortGroups)
  if(is.factor(EvaluationStats$Run)) {
    EvaluationStats$Run2 <- as.character(EvaluationStats$Run)
    EvaluationStats$Run2 <- as.numeric(with(EvaluationStats, substring(EvaluationStats$Run2, nchar(EvaluationStats$Run2))))
  } else {
    EvaluationStats$Run2 <- EvaluationStats$Run
  }
  if("Run2" %in% colnames(EvaluationStats)) {
    MaxRun <- max(EvaluationStats$Run2)
  } else {
    MaxRun <- 0
  }
  SortGroupsLength <- length(SortGroups)
  # Set default value for plot options if left out of function call
  #OrigFormat <- TRUE
  # Calcluate statistics on mean and standard deviation 
  EvalStatVarSummL <- list()
  for(i in 1:2) {
    #i=1
    if(i==1) {
      EvalStatVarSumm <- ddply(EvaluationStats, SortGroups, numcolwise(mean))
      Stat <- EvalStatVarSumm[,1:2]
      colnames(Stat) <- c("Statistic", "N")
      Stat$Statistic <- "Mean"
      Stat$N <- Freq$n
    } else {
      EvalStatVarSumm <- ddply(EvaluationStats, SortGroups, numcolwise(sd))
      Stat <- EvalStatVarSumm[,1:2]
      colnames(Stat) <- c("Statistic", "N")
      Stat$Statistic <- "SD"
      Stat$N <- Freq$n
      rownames(Stat) <- rownames(EvalStatVarSumm)
    }
    FirstColNames <- data.frame(c(colnames(EvalStatVarSumm)[1:(StatVarFirstColumn-1)], "Run2"), stringsAsFactors=FALSE)
    # Identify column names not in sort groups or statistics that should be dropped
    DropColNames <- FirstColNames[which(!FirstColNames[,1] %in% SortGroups),]
    EvalStatVarSumm <- EvalStatVarSumm[,!(names(EvalStatVarSumm) %in% DropColNames)]
    if(nrow(EvalStatVarSumm)==1) {
      rownames(EvalStatVarSumm) <- i
    } else {
      rownames(EvalStatVarSumm) <- seq(1, (nrow(EvalStatVarSumm)),1)
    }
    EvalStatVarSumm.df1 <- cbind(EvalStatVarSumm[1:(length(SortGroups))], Stat, EvalStatVarSumm[(length(SortGroups)+1):ncol(EvalStatVarSumm)])
    head(EvalStatVarSumm.df1)
    EvalStatVarSummL[[i]] <- EvalStatVarSumm.df1
  }
  #
  EvalStatVarSumm.df <- do.call(rbind, EvalStatVarSummL)
  rownames(EvalStatVarSumm.df)
  # Write output
  setwd(OutDirectIn)
  write.table(EvalStatVarSumm.df, file=OutName, sep=",", col.names=NA)
  #
  return(EvalStatVarSumm.df)
}
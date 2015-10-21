rm(list=ls())

##the name of this script
script.name=parent.frame(2)$ofile
## the location (absolute path) to this script
script.location=normalizePath(dirname(parent.frame(2)$ofile))

library(MASS)

AFNI_R_DIR=Sys.getenv("AFNI_R_DIR", unset=NA)

## use the functions for loading and saving briks from the AFNI
## distribution as they can cleanly handle floats/doubles
if ( ! is.na(AFNI_R_DIR) ) {
  source(file.path(AFNI_R_DIR, "AFNIio.R"))
} else {
  stop("Couldn't find AFNI_R_DIR in environment. This points to the location from which to load functions for reading and writing AFNI BRIKS. Stopping!")
}

os=system('uname -s', intern=T)
if (os == "Darwin") {
  Group.results.dir=file.path("/Volumes/opt/mriAnalyses/TMARC/Group.results.MNI.13Apr2012/")
} else {
  Group.results.dir=file.path("/data/TMARC/Group.results.MNI.13Apr2012/")
}

ncpus=1
## flag to indincate whether a progress bar should be printed or one
## line per slice if no progress bar is to be used
useProgressBar=TRUE

## should the AFNI read/write routines produce verobse output
verbose=FALSE

## this indicates whether bootstrapping should be performed
doBootstrapping=FALSE

bootstrappingBrikLabelSuffxes=c("bias", "t.value", "ciLower", "ciUpper")
####################################################################################################
## function definitions

## this function takes the row and column names from an rlm
## coefficient matrix and combines these (along row) to get a list of
## labels for later use as labels to the 3drefit command so the briks
## in the bucket are correctly labeled. Spaces are replaced with . and
## ( or ) are deleted.
makeBrikLabels <- function (inRlmCoef, inBoot=FALSE) {
  rns=rownames(inRlmCoef)
  cns=colnames(inRlmCoef)
  
  labels=c()
  
  for (i in 1:length(rns) ) {
    for (j in 1:length(cns) ) {
      name=gsub("..", ".", gsub(" ", ".", paste(gsub("[(]|[)]", "", rns[i]), cns[j], sep=".")), fixed=TRUE)
      labels=c(labels, name)
    }
  }
  if (inBoot) {
    for (i in 1:length(rns) ) {
      for (j in 1:length(bootstrappingBrikLabelSuffxes)) {
        name=gsub("..", ".", gsub(" ", ".", paste(gsub("[(]|[)]", "", rns[i]), bootstrappingBrikLabelSuffxes[j], sep=".")), fixed=TRUE)
        labels=c(labels, name)
      }
    }
  }

  outList=list("numberOfLabels"=length(labels),
    "bootstrapLabelsStartAt"=ifelse(inBoot, (length(rns)*length(cns)) + 1, NA),
    "labels"=labels)

  return (outList)
}

## make a list of the AFNI brik indices that correspond to the t-stat
## from the matrix of coefficients produced by rlm. This take into
## account that AFN inumbers from 0 not 1
makeAfniTtestBrikIds <- function(inRlmCoef, inBoot=FALSE) {
  nrows=length(rownames(inRlmCoef))
  ncols=length(colnames(inRlmCoef))
  
  indices=seq(3, nrows*ncols, by=3)

  if (inBoot) {
    indices=c(indices, seq(indices[length(indices)]+2, (nrows*ncols) + (nrows*length(bootstrappingBrikLabelSuffxes)), by=4))
  }
  
  return (indices)
}

## given a list of AFNI brik indices and a degrees of freedom (inDf)
## make an appropriate list of statpar arguments to add to a 3drefit
## command line
makeAfniStatparArguments <- function(inDf, inAfniFtestBrikIds) {
  return( paste(sapply(inAfniFtestBrikIds, function(x) { sprintf("-substatpar %d fitt %d", x[1], inDf) }), collapse=" ") )
}

cleanRegressionVariable <- function(inName) {
  ## Replace 2 or more consequtive . with one . and remove any
  ## trailing . from the name using the gsub function
  return (gsub("\\.{2,}", ".", gsub("\\.$", "", inName)))
}

## this is a very trimmed down version of the runRegression function
## below. It is only for use with bootstrapping
bootRegression <- function(inData, inIndices, inModelFormula, inMaxIt=50, inNumberOfStatsBriks) {
  outStats <- vector(mode="numeric", length=inNumberOfStatsBriks)
  if ( ! inherits(myrlm <- try(rlm(inModelFormula, data=inData[inIndices, ], maxit=inMaxIt), silent=FALSE),
                  "try-error") ) {
    cat ("length(coefficients(myrlm)) is ", coefficients(myrlm), "\n")
    return(coefficients(myrlm))
  } else {
    cat("Got an exception\n")
    return(outStats)
  }
}

runRegression <- function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
  outStats <- vector(mode="numeric", length=inNumberOfStatsBriks)
  if ( ! all(inData == 0 ) ) {

    ## if inData is all zero then we are in a portion of the masked
    ## out data and should therefore not perform the rlm on it just
    ## return a vector of zeros, if it's not then we get in this
    ## branch and can try the rlm
    
    inModel$mri<-inData
    
    myrlm <- rlm(inModelFormula, data = inModel, maxit=inMaxIt)
    
    outStats[1] = mean(inModel$mri)
    ## > coef(summary(model))
    ## Value Std. Error    t value
    ## (Intercept) 89.09177404  8.0221830 11.1056770
    ## Age          0.04087267  0.1146234  0.3565821
    ## Educ         0.71157766  0.5370092  1.3250753
    ## as.vector(t(coef(summary(model))))
    ## [1] 89.09177404  8.02218304 11.10567705  0.04087267  0.11462344  0.35658210  0.71157766  0.53700922  1.32507531
    ##outStats[2:length(outStats)]=as.vector(t(coef(summary(myrlm))))

    if (inBoot) {
      if (is.na(inBootstrapStatsStartAt)) {
        stop("***ERROR in runRegression: inBootstrapStatsStartAt has not been set. It is currently NA. Cannot continue. Stopping\n")
      }
      ## cat ("\nRLM Indices: ", 2:(inBootstrapStatsStartAt-1), "\n")
      ## cat ("coef(summary(myrlm)): ", as.vector(t(coef(summary(myrlm)))), "\n")
      ## cat ("Length of stats vector is", length(as.vector(t(coef(summary(myrlm))))), "\n")
      ## cat ("inBootstrapStatsStartAt is ", inBootstrapStatsStartAt, "\n")
      ## outStats[2:(inBootstrapStatsStartAt-1)]=as.vector(t(coef(summary(myrlm))))
      ## cat ("inside inBoot outStats is: ", outStats, "\n")
      ## cat ("number of stats briks should be: ", length(2:(inBootstrapStatsStartAt-1)), "\n")
      
      bootStats=boot(inModel, bootRegression, R=inR, inModelFormula=inModelFormula, inMaxIt=inMaxIt, inNumberOfStatsBriks=length(2:(inBootstrapStatsStartAt-1)))
      if (is.vector(bootStats)) {
        bootBias=apply(bootStats$t, 2, mean) - bootStats$t0
        ## these are vectors with as many columns as there are terms in
        ## the regression model. Dont forget to include the intercept
        ## when you're ttrying to mentalize this
        bootSdCoeff=apply(bootStats$t, 2, sd)
        bootTCoeff=bootStats$t0 / bootSdCoeff
        bootCi=matrix(0, nrow=length(bootTCoeff), ncol=2)
        outputBootstrapStats=c()
        for (termIndex in seq(1, length(bootTCoeff))) {
          ## the normal element of the CI value contains 3 elements: 1) the CI
          ## level (in this case 0.95), 2) the lower bound on the CI, 3) the
          ## upper bound on the CI
          ci=boot.ci(bootStats, conf = c(0.95), type = c("norm"), index = termIndex)$normal[2:3]
          bootCi[termIndex, ]=ci
          outputBootstrapStats=c(outputBootstrapStats, bootBias[termIndex], bootTCoeff[termIndex], bootCi[termIndex, ])
          ##cat ("outputBootstrapStats now: ", outputBootstrapStats, "\n")
        }
        outStats[inBootstrapStatsStartAt:inNumberOfStatsBriks]=outputBootstrapStats
      }
    } ## end of if (inBoot) {
    else {
      outStats[2:length(outStats)]=as.vector(t(coef(summary(myrlm))))
    }
  } ## end of if ( ! all(inData == 0 ) ) {
  ##cat ("outStats is: ", outStats, "\n")  
  return(outStats)
}

## enable debugging of runRlm
##debug(runRegression)
##trace("runRegression", quote(if(! all(inData == 0 ) ) browser()))
####################################################################################################

## this file stores the demographic/neuropsych/neuromed information
## for each subject. These two are later merged (the mgd variable)
## so that only some of the the info from demogrpahics (only from
## those subjects included in the subjectOrder file, which should
## match exactly the order of the inputBrik) is added to the data
## from subjectOrder
demographicsFilename=file.path(Group.results.dir, "overview2.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

regressionVariables=list(
  list(variable="BDI",                                             name="BDI (at time of scan)"),
  list(variable="BDI",                                             name="BDI (at time of scan)")  
  )

## extract the list of variable names from the regressionVariables list rvs=regression variables
rvs=unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))]
## select only the columns we want to perform regressions on from the hnrcData data frame
m=match(unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))], colnames(hnrcData))
## remove NAs caused by the BDI and POMS not coming from the hnrcData frame
m=m[!is.na(m)]
## m now contains only the column numbers of the neuropsych variables we'll regress the %cs against

group="mddOnly"

contrasts=c("emotiveVsNeutral", "happyVsSad", "happyVsFearful", "happyVsNeutral", "happyRemVsHappyNotrem",
  "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")
reml="REML"

for ( regressionVariableCount in 1:length(regressionVariables ) ) {

  ## this is the version of the regressionVariable id that will be
  ## used in the model matrix and its accompaning formula and the
  ## output brik filename
  rvName=cleanRegressionVariable( regressionVariables[[regressionVariableCount]]$variable)
  
  for (contrast in contrasts) { 
    cat("####################################################################################################\n")
    cat(sprintf("*** Running RLM for the %s contrast of the %s group\an", contrast, group))
    
    ## setup all the filenames
    inputBucketFilename=paste("pine.bucket", groups, contrast, reml, "masked+tlrc.HEAD", sep=".")
    outputBucketPrefix= paste("pine", groups, contrast, reml, rvName, "rlm.bucket", sep=".")

    ## this file stores the order of the subjects in each of the following BRIK files
    subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, contrast, reml, "csv", sep="."))
    
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)
    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

    inputBrikFilename=file.path(Group.results.dir, inputBucketFilename)
    
    cat("*** Reading", inputBrikFilename, "\n")
    inputBrik=read.AFNI(inputBrikFilename, verb=verbose)
    
    ## dim is the number of parameters from the rlm that will be stored as
    ## subbriks and exported back for use in AFNI
    dimX=inputBrik$dim[1]
    dimY=inputBrik$dim[2]
    dimZ=inputBrik$dim[3]
    
    ## this stored the total number of subbriks for all subjects for all types (i.e. wins, losses, and ties)
    numberOfBriks=inputBrik$dim[4]
    
    ## now merge the subject order and the demographics file so we have
    ## the table completed for each subbrik in in the bucket file
    
    ## the demographics file will contain only one row for each subject
    ## but the subjectOrder file contains one row for each subject for
    ## each period type and hence will be much longer, this allows us to
    ## merge the two tables so we can complete the table for use in
    ## creation of the SPM. mgd stands for merged
    
    mgd=cbind(subjectOrder,
      demographics[match(subjectOrder$subject, demographics$ID), m]
      )
    
    mrData = inputBrik$brk
    dim(mrData) = c(dimX, dimY, dimZ, numberOfBriks)
    
    ## set up the data frame with the dependant variables we will want
    ## to include in the regression formula
    
    model = data.frame(
      as.vector(mgd[, regressionVariables[[regressionVariableCount]]$variable])
      )
    colnames(model) = c(rvName)
    
    ## set up the formula for the model that we want to regress
    modelFormula  = as.formula(paste("mri ~", rvName, sep=" "))
    
    ## if we got this far we should be able to run a single voxel to work
    ## out the DOF for the various stats, generate labels and indices for
    ## use with the unlisted rlm in the runRlm function
    i = 21
    j = 21
    k = 32
    model$mri = mrData[i, j, k, ]
    if( inherits(tempRlm <- try(rlm(modelFormula, data = model),
                                silent=FALSE),
                 "try-error") ) {
      tempRlmCoef=0
      tempRlmDf=0
      stop("Got an exception trying to setup the tempRlmCoef and tempRlmDf variables. Cannot continue beyond this point. Stopping.")
      traceback()
    } else {
      s=summary(tempRlm)
      cat("*** rlm summary is\n")
      print(s)
      tempRlmCoef = coef(s)
      tempRlmDf=s$df
    }

    ## These are the labels that will be attached to the subbriks in the
    ## output stats bucket file. The labels will be dictated by the
    ## model formula and importantly the order in which the coefficients
    ## from the RLM coefficient matrix are concatenated
    ## the makeBrikLabels
    ## function does not include the mean of the fMRI at the start of
    ## the list of labels so include it here with the c()
    d=makeBrikLabels(tempRlmCoef, inBoot=doBootstrapping)
    outputStatsBrikLabels=c("Mean", d[["labels"]])
    
    ## we add 1 to both numberOfStatsBriks and numberOfStatsBriks
    ## becasue makeBrikLabels does not take into account the
    ## factthat we will add an additional subbrik (the mean) to
    ## the output stats. Hence the output of makeBrikLabels is
    ## always 1 too small
    
    ## the number of stats subbriks to write out. This is dictated by the
    ## model Formula, changes to it likely imply changes to this number
    ##numberOfStatsBriks=length(outputStatsBrikLabels)
    numberOfStatsBriks=d[["numberOfLabels"]] + 1

    bootstrapStatsStartAt=0
    if (doBootstrapping) {
      bootstrapStatsStartAt=d[["bootstrapLabelsStartAt"]] + 1
    }
    
    ##stop("Stopping")
    
    maxIter=100
    resamples=25
    Stats = array(0, c(dimX, dimY, dimZ, numberOfStatsBriks))
    cat(paste("Starting at", date(), "\n"))
    startTime=proc.time()
    if (useProgressBar) {
      pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
    }
    if (ncpus > 1 ) {
      ## multiple cpus
      library(snow)
      
      cluster = makeCluster(ncpus, type = "SOCK")
      clusterEvalQ(cluster, library(MASS));

      ##function (inData, inNumberOfStatsBriks, inModel, inModelFormula, inMaxIt=50, inBoot=FALSE, inR=25, inBootstrapStatsStartAt=NA) {
      for ( kk in 1:dimZ) {
        if (useProgressBar) {
          setTxtProgressBar(pb, kk)
        } else {
          cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
        }
        Stats[ , , kk, ] = aperm(parApply(cluster, mrData[ , , kk, ],  c(1, 2), runRegression,
               inNumberOfStatsBriks=numberOfStatsBriks, inModel=model, inModelFormula=modelFormula, inMaxIt=maxIter, inBoot=doBootstrapping, inR=resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
      }
      stopCluster(cluster)
    } else {
      ##for ( kk in 1:dimZ ) {              
      for ( kk in 32:32 ) {    
        ## single cpu
        if (useProgressBar) {
          setTxtProgressBar(pb, kk)
        } else {
          cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
        }
        Stats[ , , kk, ] = aperm(apply(mrData[ , , kk, ],  c(1, 2), runRegression,
               inNumberOfStatsBriks=numberOfStatsBriks, inModel=model, inModelFormula=modelFormula, inMaxIt=maxIter, inBoot=doBootstrapping, inR=resamples, inBootstrapStatsStartAt=bootstrapStatsStartAt), c(2, 3, 1))
      }
    }

    if (useProgressBar) {
      close(pb)
    }
    
    cat(paste("Ended at", date(), "\n"))
    cat("Time consumed\n")
    print(proc.time() - startTime)
    
    ##stop()
    
    rlmOutBrikFilename=paste(outputBucketPrefix, ".", format(Sys.time(), "%Y%m%d-%H%M%Z"), view.AFNI.name(inputBrikFilename), sep="")
    rlmOutBrikFqfn=file.path(Group.results.dir, rlmOutBrikFilename)
    hostname=system('hostname', intern=T)
    user=Sys.getenv("USER")
    cat("*** Writing bucket file ", rlmOutBrikFqfn, "\n")
    write.AFNI(rlmOutBrikFqfn,
               Stats, verb=verbose,
               ##label = baselineBrik$head$DATASET_NAME,
               label=outputStatsBrikLabels,
               note = paste(paste("[", user, "@", hostname, ": ",  date(), "]", sep=""), file.path(script.location, script.name)),
               origin = inputBrik$origin,
               delta = inputBrik$delta,
               orient= inputBrik$orient)

    statpar = "3drefit"
    
    if ( is.matrix(tempRlmCoef) ) {
      cat ("Making statpar arguments\n")
      
      afniTtestBrikIds=makeAfniTtestBrikIds(tempRlmCoef)
      statparArguments=makeAfniStatparArguments(tempRlmDf[2], afniTtestBrikIds)
      statpar = paste(statpar, statparArguments)    
    }
    statpar = paste("(cd",  Group.results.dir, ";", statpar, " -view tlrc -space MNI -addFDR -newid ", rlmOutBrikFilename, ")")
    cat(statpar, "\n")
    system(statpar)
  } ## end of  for (contrast in contrasts)
} ## end of for (term in regressionterms) {



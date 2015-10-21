rm(list=ls())
gc()
##the name of this script
script.name=parent.frame(2)$ofile
## the location (absolute path) to this script
script.location=normalizePath(dirname(parent.frame(2)$ofile))

library(nlme)
##library(parallel)
##library(contrast)

verbose=F

ncpus=1
##detach("package:parallel", unload=TRUE)

if ( Sys.info()["sysname"] == "Darwin" ) {
  Admin.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin"
  Config.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/config"
  Ppi.seeds.data.dir=file.path(Config.data.dir, "ppi_seeds")
  
  Group.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data"
  Group.results.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results/ppi"  
  ncpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
  cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else if ( Sys.info()["sysname"] == "Linux" ) {
  Sys.setenv(AFNI_R_DIR="/abin")
  stop("Cannot setup paths for this system")
  ## Group.data.dir="/mnt/nfs/yangdata/restingstate/data/Group.data"
  ## Group.results.dir="/mnt/nfs/yangdata/restingstate/data/Group.results"  
  ## ncpus=as.integer(system("cat /proc/cpuinfo |grep -c processor", intern=T))
  ## cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else {
    stop("Sorry can't setup directories on this architecture\n")
}

AFNI_R_DIR=Sys.getenv("AFNI_R_DIR", unset=NA)
## use the functions for loading and saving briks from the AFNI
## distribution as they can cleanly handle floats/doubles
if ( ! is.na(AFNI_R_DIR) ) {
  source(file.path(AFNI_R_DIR, "AFNIio.R"))
} else {
  Sys.setenv(PATH=system("bash -login -c 'echo $PATH'", intern=TRUE))
  source(file.path(Sys.getenv("MRI_SOFTWARE_ROOT"), "afni", "AFNIio.R"))
}

if ( ! file.exists(Group.data.dir) ) {
    stop(paste("Cannot continue:", Group.data.dir, "does not exist. Stopping"))
}

if ( ! file.exists(Group.results.dir) ) {
    if ( ! dir.create(Group.results.dir)) {
        stop(paste("Cannot continue:", Group.results.dir, "does not exist and could not create it. Stopping"))
    } else {
        cat(paste(Group.results.dir, "created."))
    }
}

## load the MASC scoring code
## source("scoreMasc.r")

## flag to indincate whether a progress bar should be printed or one line per slice if no progress bar is to be used
useProgressBar=TRUE

## this function takes the row and column names from an anova and
## combines them (along row) to get a list of labels for later use as
## labels to the 3drefit command so the briks in the bucket are
## correctly labeled
makeBrikLabels <- function (inAnova) {
  rns=rownames(inAnova)
  cns=colnames(inAnova)
  
  labels=c()
  
  for (i in 1:length(rns) ) {
    for (j in 1:length(cns) ) {
      name=paste(gsub("[(]|[)]", "", rns[i]), cns[j], sep=".")
      labels=c(labels, name)
    }
  }
  return (labels)
}

makeContrastBrikLabels <- function(inContrastsList) {
  return (unlist(lapply(inContrastsList, function(x) { paste(x$name, c("contrast", "t-value"), sep=".") } )))
}

## This function generates lists if indices that can be applied to an
## unlisted anova to extract the elements of said list in the correct
## order for writing to an AFNI BRIK. It works for odd and even
## numbers of rows.

makeAnovaIndices <- function (inAnova, inFirstTwoOnly=FALSE) {
  nrows=length(rownames(inAnova))
  ncols=length(colnames(inAnova))
  indices=c()
  for ( r in 1:nrows ) {
    colCount=0
    for ( cl in 1:ncols ) {
      index = (r + (cl)*nrows)-nrows
      if (inFirstTwoOnly) {
        if (colCount < 2) {
          indices=c(indices, index)
        }
      }
      else {
        indices=c(indices, index)
      }
      colCount=colCount+1
    }
  }
  return (indices)
}

makeAfniFtestBrikIds <- function(inAnova) {
  nrows=length(rownames(inAnova))
  ncols=length(colnames(inAnova))
  
  return(seq(2, nrows*ncols, by=4))
}

## this function generates a numeric list of BRIK IDs (to be used when
## setting statpars in the call to 3drefit) from an anova and the list
## of preplanned contrasts
makeAfniTtestBrikIds <- function(inAnova, inContrastsList)  {
  numberOfContrasts=length(inContrastsList)
  if ( numberOfContrasts > 0 ) {
    numberOfAnovaBriks=length(colnames(inAnova)) * length(rownames(inAnova))
    return(seq(numberOfAnovaBriks + 1, numberOfAnovaBriks+(numberOfContrasts*2), by=2))
  } else {
    return( c() )
  }
}

## makes a appropriate statpar command line argument for use with
## 3drefit. Requires an anova object and a list of BRIK IDs
## corresponding to the list of f-stat briks. The t-test brik IDs
## (from the contrasts) and the corresponding degrees of freedom are
## optional and default to NULL indicating that there are no contrats
## BRIKs in the bucket file saved by this script
makeAfniStatparArguments <- function(inAnova, inAfniFtestBrikIds, inAfniTtestBrikIds=NULL, inContrastsDf=NULL) {
  anovaIndices=makeAnovaIndices(inAnova, inFirstTwoOnly=TRUE)
  temp = unlist(inAnova)
  indices=cbind(inAfniFtestBrikIds, matrix(anovaIndices, ncol=2, byrow=TRUE))
  anovaStatpar=paste(apply(indices, 1, function(x) { sprintf("-substatpar %d fift %d %d", x[1], temp[x[2]], temp[x[3]]) }), collapse=" ")
  contrastsStatpar=""
  if (! is.null(inAfniTtestBrikIds) ) {
    indices=cbind(inAfniTtestBrikIds, inContrastsDf)
    contrastsStatpar=paste(apply(indices, 1, function(x) { sprintf ("-statpar %d fitt %d", x[1], x[2])} ), collapse=" ")
  }
  return(paste(anovaStatpar, contrastsStatpar))
}

####################################################################################################
## runLme function definition

runLme <- function (inData, inZ, inNumberOfOutputBriks, inContrasts, inAnovaIndices, inModel, inModelFormula, inRandomFormula) {
  ## cat("There will be " , inNumberOfOutputBriks, " stats briks\n")
  outStats <- vector(mode="numeric", length=inNumberOfOutputBriks)
  ##cat (inData, "\n")
  if ( ! all(inData == 0 ) ) {
    
    ## if inData is all zero then we are in a portion of the masked
    ## out data and should therefore not perform the lme on it just
    ## return a vector of zeros, if it's not then we get in this
    ## branch and can try the lme
    
    inModel$fmri<-inData
    
    ## If lme throws an exception it returns an object, assigned to
    ## mylme, which inherits from the class try-error. If mylme
    ## inherits from this class an error was generated by a particular
    ## voxel. It will default to having its corresponding voxels in
    ## the outStats array set to 0.
    
    if( inherits(
                 mylme <- try(lme(fixed=inModelFormula, random=inRandomFormula, data = inModel),
                              silent=FALSE),
                 "try-error") ) {
      temp <- 0
      cat (paste("Error on slice", inZ, "\n"))
    } else {
      temp <- as.vector(unlist(anova(mylme, type="marginal")))
    }
    
    if(length(temp) > 1) {
      numberOfContrasts=length(inContrasts)
      outStats[1:(inNumberOfOutputBriks-(numberOfContrasts*2))] = temp[c(inAnovaIndices)]
      
      if (numberOfContrasts > 0 ) {
        contrastsStartAt=(inNumberOfOutputBriks-(numberOfContrasts*2))+1
        for (i in seq(1, numberOfContrasts, by=1)) {
          con1=contrast(mylme, a =inContrasts[[i]]$a, b =inContrasts[[i]]$b, type="average")
          if (length(con1) > 1) {
            outStats[contrastsStartAt:contrastsStartAt+1] <- c(con1$Contrast, con1$testStat)
          }
          contrastsStartAt=contrastsStartAt+2
        }
      }
    }
  } ## else {
  ##   cat("Got all zeros :-(\n")
  ## }
  
  return(outStats)
}
## enable debugging of runLme
##debug(runLme)
##trace("runLme", quote(if(! all(inData == 0 ) ) browser()))
####################################################################################################

prePlannedContrasts=list()

## the formulae used for the fixed effects (modelFormula) and random effects (randomFormula)
modelFormula  = as.formula("fmri ~ group")

randomFormula = as.formula("random = ~ 1 | subject")

contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")

##contrasts=c("fearfulVsHappyExtractedRightSgAcc")

##contrasts=c("fearfulVsHappy", "fearfulVsNeutral") ##, "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")
##contrasts=c("happyVsNeutral", "happyVsSad", "neutralVsSad")

## a.priori.rois=c("HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala", "HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala", "anteriorInsula.3mm", "posteriorInsula.3mm", "sgacc.left.3mm", "sgacc.right.3mm")
## a.priori.rois=c("HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala", "HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala")


#a.priori.rois=c("apriori.rois.noVmpfc.3mm", "apriori.rois.withVmpfc.3mm")
a.priori.rois=c("apriori.rois.withVmpfc.3mm")


##a.priori.rois=c("sgacc.left.3mm", "sgacc.right.3mm")

## makeContrasts <- function ( inContrastVector ) {
##     suffixes=c()
##     for ( contrast in inContrastVector ) {
##         seedListFile=file.path(Ppi.seeds.data.dir, paste(contrast, "seeds.txt", sep="."))
##         numberOfSeeds=strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2]
##         stimuli=unlist(strsplit(contrast, "Vs", fixed=FALSE))
##         for ( ii in seq(1, numberOfSeeds) ) {
##             for ( stimulus in stimuli) {
##                 suffixes=c(suffixes, sprintf("roi%d.seed.%s.ROIx%s", ii, contrast, stimulus))
##             }
##         }
##     }

##     return(sufficies)
## }

makeContrastsFromBetweenGroupResults <- function ( inContrastVector ) {
    suffixes=c()
    stimulus="ROIinteraction"

    for ( contrast in inContrastVector ) {
        seedListFile=file.path(Ppi.seeds.data.dir, paste(contrast, "seeds.txt", sep="."))
        numberOfSeeds=strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2]
        for ( ii in seq(1, numberOfSeeds) ) {
            suffixes=c(suffixes, sprintf("roi%d.seed.%s.%s", ii, contrast, stimulus))
        }
    }

    return(suffixes)
}

makeContrastsFromBetweenGroupResultsWithAPrioriMasks <- function ( inContrastVector, inMasks ) {
    suffixes=c()
    stimulus="ROIinteraction"

    for (contrast in inContrastVector) {
        for (mask in inMasks) {
            seedListFile=file.path(Ppi.seeds.data.dir, paste(contrast, mask, "seeds.txt", sep="."))
            print(seedListFile)
            numberOfSeeds=strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2]
            for ( ii in seq(1, numberOfSeeds) ) {
                suffixes = c(suffixes, sprintf("roi%d.seed.%s%s.%s", ii, contrast, mask, stimulus))
            }
        }
    }
            
    return(suffixes)
}

makeContrastsFromBetweenAprioriRois <- function ( inContrastVector, inRois ) {

    contrastByRois=expand.grid(inContrastVector, inRois)
    suffixes = apply(contrastByRois, 1,
        function(x) {
            return(sprintf("%s.seed.%s.ROIinteraction", x[2], x[1]))
        }
        )
    
    return(suffixes)
}

##suffixes=makeContrastsFromBetweenAprioriRois(contrasts)

suffixes=makeContrastsFromBetweenGroupResultsWithAPrioriMasks(contrasts, a.priori.rois)

##suffixes=makeContrastsFromBetweenGroupResults(contrasts)

##stop("Check contrast suffixes\n")

numberOfContrasts=length(suffixes)
contrastCount=1

reml="REML"
groups="mddAndCtrl"

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(Admin.data.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

for (contrast in suffixes) {

  cat("####################################################################################################\n")
  cat(sprintf("*** Running LME for the %s contrast (%02d of %02d)\n", contrast, contrastCount, numberOfContrasts))

  ## setup all the filenames
  inputBucketFilename=paste("pine.bucket", groups, contrast, reml, "z-score", "masked+tlrc.HEAD", sep=".")
  outputBucketPrefix= paste("pine.lme.bucket", groups, contrast,  reml, sep=".")

  clusterLogFilename=paste("pine.cluster", groups, contrast, "log", sep=".")
  
  ## this file stores the order of the subjects in each of the following BRIK files
  subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, contrast, reml, "z-score", "csv", sep="."))

  inputBrikFilename=file.path(Group.data.dir, inputBucketFilename)

  cat("*** Reading", subjectOrderFilename, "\n")
  subjectOrder=read.csv(subjectOrderFilename, header=T)
  subjectOrder$subject=gsub("_A", "", subjectOrder$subject, fixed=TRUE)

  ## yay for subjects with 2 IDs :-(
  subjectOrder$subject[subjectOrder$subject=="300"]="169/300"
  subjectOrder$subject=as.factor(subjectOrder$subject)
  
  cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

  cat("*** Reading",  inputBrikFilename, "\n")
  inputBrik=read.AFNI(inputBrikFilename, verb=verbose)
  
  ## dim is the number of parameters from the lme that will be stored as
  ## subbriks and exported back for use in AFNI
  dimX=inputBrik$dim[1]
  dimY=inputBrik$dim[2]
  dimZ=inputBrik$dim[3]
  
  ## this stored the total number of subbriks for all subjects for all types
  numberOfBriks=inputBrik$dim[4]
  
  ## now merge the subject order and the demographics file so we have
  ## the table completed for each subbrik in in the bucket file
  
  ## the demographics file will contain only one row for each subject
  ## but the subjectOrder file contains one row for each subject for
  ## each period/block type in the experiment and hence will be much
  ## longer, this allows us to merge the two tables so we can complete
  ## the table for use in creation of the SPM. mgd stands for merged
  
  ##mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "MASC.total", "MASC.tscore")])
  mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI")])  
  ## ensure that subject is s factor
  mgd$subject=as.factor(mgd$subject)
  ## mgd$MASC.tscore=as.numeric(mgd$MASC.tscore)
  
  ## this complicated looking regexp stuff cleans up the years in DOB
  ## and MRI with 4 digits to be just 2 digits
  mgd$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", mgd$DOB)
  mgd$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", mgd$MRI)
  mgd$DOB=as.Date(mgd$DOB, "%m/%d/%y")
  mgd$MRI=as.Date(mgd$MRI, "%m/%d/%y")
  age.in.days=difftime(mgd$MRI, mgd$DOB, units="days")
  age.in.days=as.numeric(age.in.days)
  age.in.weeks=difftime(mgd$MRI, mgd$DOB, units="weeks")
  age.in.weeks=as.numeric(age.in.weeks)
  mgd$age.in.years=(age.in.weeks)/52

  ## print(mgd$MASC.total)
  
  ## mgd.dim=dim(mgd)
  ## for (r in seq(1, mgd.dim[1]) ) {
  ##     cat("##################################################\n")
  ##     subjectNumber=mgd[r, "subject"]
  ##     gender=mgd[r, "Gender"]
  ##     age=round(mgd[r, "age.in.years"], 0)
  ##     old.masc.tscore=mgd[r, "MASC.tscore"]
      
  ##     cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f,\n", r, subjectNumber, gender, age, mgd[r, "MASC.total"], old.masc.tscore))
      
  ##     new.masc.tscore=scoreMasc(gender, age, mgd[r, "MASC.total"])
  ##     if (is.na(new.masc.tscore) ) {
  ##         warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, mgd[r, "MASC.total"]))
  ##     }
      
  ##     mgd[r, "MASC.tscore"]=new.masc.tscore
      
  ##     cat (sprintf("Old MASC tscore=%0.0f, new MASC tscore=%0.0f\n", old.masc.tscore, new.masc.tscore))
  ## }
  
  mrData = inputBrik$brk
  dim(mrData) = c(dimX, dimY, dimZ, numberOfBriks)
  
  ## now multiply by the mask to ensure that only those voxels in the
  ## mask will be included in the lme computation, all others are set to
  ## 0. The check for whether a voxel is in the mask of not is to be
  ## found in runLme at if ( ! all(inData == 0 ) )
  ##mrData = array(apply(mrData, 4, function (x) x * maskBrik$brk[,,,1]), dim=c(dimX, dimY, dimZ, numberOfBriks))

  model = data.frame(
      "subject" = as.vector(mgd$subject),
      "group"   = as.vector(mgd$Grp)
    )

  contrs=list()
  
  nContrasts=length(contrs)
  
  ## if we got this far we should be able to run a single voxel to work
  ## out the DOF for the various stats, generate labels and indices for
  ## use with the unlisted anova in the runLme function
  i = 21
  j = 21
  k = 32
  model$fmri = mrData[i, j, k, ]
  if( inherits(tempLme <- try(lme(fixed=modelFormula, random=randomFormula, data = model),
                              silent=FALSE),
               "try-error") ) {
    tempAnova <- 0
    stop("Got an exception trying to setup the tempLme variable. Cannot continue beyond this point. Stopping.")
  } else {
    tempAnova <- anova(tempLme, type="marginal")
    cat("The temporary ANOVA is\n")
    print(tempAnova)
  }
  contr.df <- vector(mode="numeric", length=nContrasts)
  if (nContrasts > 0 ) {
    for (i in seq(1, nContrasts, by=1)) {
      con1=try(contrast(tempLme, a =contrs[[i]]$a, b =contrs[[i]]$b, type="average"))
      if (length(con1) > 1) {
        contr.df[i] = con1$df
      }
    }
  }
  
  
  ## These are the labels that will be attached to the subbriks in the
  ## output stats bucket file. The labels will be dictated by the model
  ## formula and importantly the order in which the coefficients form
  ## the ANOVA matrix are concatenated
  baseStatsBrikLabels=makeBrikLabels(anova(tempLme))
  contratStatsBrikLabels=makeContrastBrikLabels(contrs)
  outputStatsBrikLabels=c(baseStatsBrikLabels, contratStatsBrikLabels)
  
  ## total number of briks in the output file, including everything
  ## from the anovas and contrasts
  numberOfOutputBriks=length(outputStatsBrikLabels)
  anovaIndices=makeAnovaIndices(anova(tempLme))
  

  ##stop("Stopping. Everything ok so far?")

  
  Stats = array(0, c(dimX, dimY, dimZ, numberOfOutputBriks))
  cat(paste("Starting at", date(), "\n"))
  startTime=proc.time()
  if (ncpus > 1 ) {
    ## multiple cpus
    library(snow)
    
    ## outfile="" should result in slave output being printed on the
    ## controlling tty device, i.e. on the master
    cat("Making a cluster with" , ncpus, "cpus\n")
    cluster = makeCluster(ncpus, type = "SOCK", outfile=file.path(Group.results.dir, clusterLogFilename))
    clusterEvalQ(cluster, library(nlme));
    ##clusterEvalQ(cluster, library(contrast))
    ## may need to add contrasts here later
    if (useProgressBar) {
      pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
    }
    for ( kk in 1:dimZ ) {
      if (useProgressBar) {
        setTxtProgressBar(pb, kk)
      } else {
        cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
      }
      a=parApply(cluster, mrData[ , , kk, ],  c(1, 2), runLme, inZ = kk, inNumberOfOutputBriks=numberOfOutputBriks, inContrasts=contrs,
        inAnovaIndices=anovaIndices, inModel=model, inModelFormula=modelFormula, inRandomFormula=randomFormula)
      ##cat(dim(a), "\n")
      Stats[ , , kk, ] = aperm(a, c(2, 3, 1))
    }
    stopCluster(cluster)
  } else {
    if (useProgressBar) {
      pb <- txtProgressBar(min = 0, max = dimZ, style = 3)
    }
    for ( kk in 1:dimZ ) {    
      if (useProgressBar) {
        setTxtProgressBar(pb, kk)
      } else {
        cat(paste("Processing Z slice", kk, "started at" , date(), "\n"))
      }
      a=apply(mrData[ , , kk, ],  c(1, 2), runLme, inZ = kk, inNumberOfOutputBriks=numberOfOutputBriks, inContrasts=contrs,
        inAnovaIndices=anovaIndices, inModel=model, inModelFormula=modelFormula, inRandomFormula=randomFormula)
      ##cat(dim(a), "\n")
      Stats[ , , kk, ] = aperm(a, c(2, 3, 1))
    }
  } ## single cpu
  
  if (useProgressBar) {
    close(pb)
  }
  
  cat(paste("Ended at", date(), "\n"))
  cat("Time consumed\n")
  print(proc.time() - startTime)
  
  ##stop()
  
  lmeOutBrikFilename=paste(outputBucketPrefix, ".", format(Sys.time(), "%Y%m%d-%H%M%Z"), view.AFNI.name(inputBrikFilename), sep="")
  lmeOutBrikFqfn=file.path(Group.results.dir, lmeOutBrikFilename)
  hostname=system('hostname', intern=T)
  user=Sys.getenv("USER")
  cat("*** Writing bucket file ", lmeOutBrikFqfn, "\n")
  write.AFNI(lmeOutBrikFqfn,
             Stats, verb=verbose,
             ##label = baselineBrik$head$DATASET_NAME,
             label=outputStatsBrikLabels,
             note = paste(paste("[", user, "@", hostname, ": ",  date(), "]", sep=""), file.path(script.location, script.name)),
             origin = inputBrik$origin,
             delta = inputBrik$delta,
             orient= inputBrik$orient)
  
  
  statpar = "3drefit"
  if ( is.list(tempAnova) ) {
    cat ("Making statpar arguments\n")
    
    afniFtestBrikIds=makeAfniFtestBrikIds(tempAnova)
    afniTtestBrikIds=makeAfniTtestBrikIds(tempAnova, contrs)
    statparArguments=makeAfniStatparArguments(tempAnova, afniFtestBrikIds, afniTtestBrikIds, contr.df)
    statpar = paste(statpar, statparArguments)    
  }
  statpar = paste("(cd",  Group.results.dir, ";", statpar, " -view tlrc -space MNI -addFDR -newid ", lmeOutBrikFilename, ")")
  cat(statpar, "\n")
  system(statpar)
  contrastCount=contrastCount + 1
} ## end of  for (contrast in contrasts)

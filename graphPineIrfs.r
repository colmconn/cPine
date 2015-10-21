rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
## for the turnpoints function that IDs where turnpoint son a curve
library(pastecs)
## the following two packages are for parallelizing plyr
library(doMC)
library(parallel)
library(nlme)

####################################################################################################
### Function Definitions 
####################################################################################################

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

substituteShortLabels <- function(inLevel) {
  returnSubstitutedLabel = gsub("[0-9]+ ", "", gsub("Inf", "Inferior", gsub("Sup", "Superior", gsub("Gy", "Gyrus",
      gsub("^R", "Right", gsub("^L", "Left", inLevel, fixed=FALSE), fixed=FALSE), fixed=TRUE), fixed=TRUE), fixed=TRUE))
  
  return (returnSubstitutedLabel)
}

stderror <- function(x) sd(x)/sqrt(length(x))

## http://wiki.stdout.org/rcookbook/Graphs/Plotting%20means%20and%20error%20bars%20(ggplot2)/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    ## New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    ## This is does the summary; it's not easy to understand...
      datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c(N    = length2(xx[,col], na.rm=na.rm),
                     mean = mean (xx[,col], na.rm=na.rm),
                     median = median (xx[,col], na.rm=na.rm),
                     mad = mad (xx[,col], na.rm=na.rm),                       
                     sd   = sd   (xx[,col], na.rm=na.rm),
                     min  = min  (xx[,col], na.rm=na.rm),
                     max  = max  (xx[,col], na.rm=na.rm),                       
                     nacount  = sum  (is.na((xx[,col])))
                     )
                 },
                 measurevar,
                 na.rm
                 )

    ## Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  ## Calculate standard error of the mean

    ## Confidence interval multiplier for standard error
    ## Calculate t-statistic for confidence interval: 
    ## e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

stack <- function(){ 
  it <- list() 
  res <- list( 
      push=function(x){ 
          it[[length(it)+1]] <<- x 
      }, 
      pop=function(){ 
          val <- it[[length(it)]] 
          it <<- it[-length(it)] 
          return(val) 
      }, 
      value=function(){ 
          return(it) 
      } 
      ) 
  class(res) <- "stack" 
  res 
} 

print.stack <- function(x,...) { 
    print(x$value()) 
} 

push <- function(stack,obj) { 
    stack$push(obj) 
} 

pop <- function(stack) { 
    stack$pop() 
}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

  Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
    symbols   =  c("***", "**", "*", ".", " "))
  f=format(Signif)

  ## only return the first one as we're only interested in marking significant group effects
  return(f[which.pValues])
}


readSubjectLists <- function( inSubjects ) {
    subjects=do.call(rbind, lapply(subjectListFilenames,
        function(X) {
            data.frame(read.table(X, header=F, sep=""))
        }
        ))

    return (as.vector(subjects[, 1]))
}

buildRoiStatsFileNames <- function (inSubjects, inStimuli, inContrast ) {

    subjectsByStimuli = expand.grid(inSubjects, inStimuli)

    filenames=apply(subjectsByStimuli, 1,
        function(x) {
            sprintf("%s/%s/functional/roistats.%s.acquisition.%s.stimulus.%s.contrast.iresp.txt", data.root, x[1], x[1], x[2], inContrast)
        }
        )
    return(filenames)
}

readRoiStatsFiles <- function(inRoiStatsFiles, inStimuli) {

    roistats=do.call(rbind, lapply(inRoiStatsFiles,
        function(X) {
            data.frame(read.table(X, header=T, sep=""))
        }
        ))
    
    roistats$stimulus=as.factor(rep(capwords(inStimuli), each=dim(roistats)[1]/length(inStimuli)))
    ##roistats=roistats[, -1]
    
    return(roistats)
}


gammaVar <- function (x, inDelay, inK, inRise, inDecay) {
    retVal=rep(0.0,length(x))
    ## which indices of x have positive differences with inDelay, no
    ## negative logs
    idx=which (! x-inDelay < 0 )
    if ( inDecay >= 0 ) {
        retVal[idx]=inK * exp( log(x[idx]-inDelay) * inRise ) * exp(-(x[idx]-inDelay)/inDecay)
    }
    return(retVal)
}

gammaVarSSE <- function (par, x, y) {
    delay=par[1]
    k=par[2]
    rise=par[3]
    decay=par[4]

    yhat=gammaVar(x, delay, k, rise, decay)

    sse=sum((y - yhat)^2)

    return(sse)
}

optimizeModelParameters <- function ( inDf ) {
    if ( ! all(c("time", "value") %in% colnames(inDf) ) ) {
        stop ("*** optimizeModelParameters: The column names of inDf do not contain (time, value)\n")
    }
    
    model=try(
        optim(c(1, 1, 2, 1.5), gammaVarSSE,
              method="L-BFGS-B",
              lower=c(-Inf, -Inf, -Inf, 0),
              upper=c(30, 30, 30, 30),
              hessian=TRUE, x=inDf$time, y=inDf$value)
        ## model=powell(c(1, 1, 2, 1.5), gammaVarSSE,
        ##     x=inDf$time, y=inDf$hrf)
        ##        )
        )
    return(model)
}

## inListElementName is the name of a element from inListOfHrfModels
## to which the gamma variate model should be fit. inListElementName
## is if the form Group.stimulus.cluster
## 
fitModelsToNewData <- function ( inListElementName, inListOfHrfModels, inTimes, inGroupVars ) {

    numberOfNewTimepoints=length(inTimes)

    groupAndStimulusAndCluster=unlist(strsplit(inListElementName, ".", fixed=TRUE))
    group=groupAndStimulusAndCluster[1]
    stimulus=groupAndStimulusAndCluster[2]
    cluster=groupAndStimulusAndCluster[3]    
    
    model=inListOfHrfModels[[inListElementName]]

    if ( ! inherits (model, "try-error") ) {
        fit=gammaVar(inTimes, model$par[1], model$par[2], model$par[3], model$par[4])
    } else {
        warning(sprintf("Got a try-error for %s %s %s. Setting fitted model to 0 for each time point\n", group, stimulus, cluster))
        fit=rep(0, numberOfNewTimepoints)
    }
    
    ## construct the data frame to return
    predicted.df=data.frame(
        inTimes,
        rep(group, numberOfNewTimepoints),
        rep(stimulus, numberOfNewTimepoints),
        rep(cluster, numberOfNewTimepoints),
        fit
        )
    colnames(predicted.df)=c("time", inGroupVars, "value")
    return(predicted.df)
}

## find the most extreme values of the HRF. the hard work (actually
## finding the extrema) is done by the turnpoints function
findExtrema <- function ( inDf ) {
    retVal = c(time=NA, value=NA)
    tp=turnpoints(inDf$value)
    if ( length(tp$tppos) > 1 ) {
        print(tp$tppos)
        stop("Can't handle more than one turnpoint in the model HFR")
    } else {
        ## turnpoints will not throw errors but will return NAs chich
        ## can cause havock when using NA as a subscript into data
        ## frame
        if (is.na(tp$tppos)) {
            retVal = c(time=NA, value=NA)
        } else {
            retVal = c(time=inDf[tp$tppos, "time"], value=inDf[tp$tppos, "value"])
        }
    }
    return(retVal)
}
   

fitGammaVariateModels <- function ( inRoistats, inGroupVars ) {

    ##complexHrfModels=dlply( inRoistats, inGroupVars, optimizeModelParameters, .parallel=TRUE)
    complexHrfModels=dlply( inRoistats, inGroupVars, optimizeModelParameters)    
    
    newTimes=seq(0, max(inRoistats$time), by=0.1)
    fittedModels=do.call(rbind, lapply(names(complexHrfModels),
        fitModelsToNewData,
        ## extra arguments to fitModelsToNewData
        inListOfHrfModels=complexHrfModels, inTimes=newTimes, inGroupVars=inGroupVars))
    
    ## peaksDf=ddply(fittedModels, inGroupVars, findExtrema, .parallel=TRUE )
    peaksDf=ddply(fittedModels, inGroupVars, findExtrema )    

    return(list(complexHrfModels=complexHrfModels, fittedModels=fittedModels, peaks=peaksDf))
}

## Code snippet from 
## http://stackoverflow.com/questions/3472980/ggplot-how-to-change-facet-labels
## facet_labeller <- function(var, value){
##     value <- as.character(value)
##     if (var=="cluster") {
##         roiId=unlist(strsplit("Mean_4", "_", fixed=TRUE))[2]
##         value[value=="Female"] <- "Woman"
##         value[value=="Male"]   <- "Man"
##     }
##     return(value)
## }


createRoiStatsDataFrame <- function (inContrast) {

    roistatsFilenames=buildRoiStatsFileNames(subjectList, stimuli[[inContrast]], inContrast)
    roistatsFileExists=file.exists(roistatsFilenames)
    numberOfNonexistantRoistatsFiles=sum(!roistatsFileExists)
    
    if (numberOfNonexistantRoistatsFiles > 0 ) {
        cat(sprintf("*** The following %d files do not exist and will be removed from the ROI stats file list\n", numberOfNonexistantRoistatsFiles))
        print(roistatsFilenames[! roistatsFileExists])
    }
    roistatsFilenames=roistatsFilenames[roistatsFileExists]
    
    if ( length(roistatsFilenames) == 0 ) {
        warning("*** After removing nonexistant ROI stats files there are none left to process\n")
        retVal=NULL
    } else {
        roistats=readRoiStatsFiles(roistatsFilenames, stimuli[[inContrast]])
        ## change the file name to be the subject name
        roistats$File=as.vector(roistats$File)
        roistats$File=sub("(([0-9]*)_[ABCD])(.*)", "\\2", roistats$File)
        ## rename the File column to subject
        colnames(roistats)=c("subject", colnames(roistats)[-1])
        ## now fix up the mess with 169/300 having two subject IDs
        roistats[roistats$subject=="300", "subject"] = "169/300"
        ## convert back to a factor
        roistats$subject=as.factor(roistats$subject)
        
        roistats$Group=demographics[match(roistats$subject, demographics$ID), c("Grp")]
        roistats$Group=drop.levels(roistats$Group)
        ## rename the Sub.brick column to time
        roistats=rename(roistats, replace=c("Sub.brick"="time"))
        ## multiply by 2 to convert TRs stored in the time column to seconds
        roistats$time=roistats$time*2

        ## the EPI timeseries were grand mean scaled to 10,000 before
        ## deconvolution so the beta weights estimates for the IRF
        ## will be 100 times to big, so compensate for this by
        ## divinding by 100
        roistats[, grep ("Mean_", colnames(roistats))] = roistats[, grep ("Mean_", colnames(roistats))] / 100

        ## ### the following 3 lines of code delete all the Mean_ columns that are not Mean_4
        ## a=colnames(roistats)[grep("Mean", colnames(roistats))]
        ## drops=a[a!="Mean_4"]
        ## roistats=roistats[, !colnames(roistats) %in% drops]
        
        clusterCount=length(grep("Mean", colnames(roistats)))
        
        clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, inContrast, fLabel))
        cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
        ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
        clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))

        attr(roistats, "cluster.names") <- paste(seq(1, clusterCount), clusterWhereAmI)

        melted.roistats=melt(roistats,
            id.vars=c("subject", "Group", "time", "stimulus"),
            measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
            variable_name="cluster")
        ## melted.roistats=melt(roistats,
        ##     id.vars=c("subject", "Group", "time", "stimulus"),
        ##     measure.vars="Mean_4",
        ##     variable_name="cluster")
        
        melted.roistats$cluster=factor(melted.roistats$cluster,
            levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
            labels=paste(seq(1, clusterCount), clusterWhereAmI))
        ## melted.roistats$cluster=factor(melted.roistats$cluster,
        ##     levels=c("Mean_4"),
        ##     labels=paste("4", clusterWhereAmI))
        
        retVal=list(roistats=roistats, melted.roistats=melted.roistats)
    }
    return(retVal)
}

graphIndividualSubjectIrfs <- function (inContrast, inMeltedRoistats, inCreateIrfGraphs=FALSE ) {

    imageDirectory=file.path(Group.results.dir, inContrast, "subjects")
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }

    groupVars=c("subject", "stimulus", "cluster")
    complexAndFittedModels=fitGammaVariateModels(inMeltedRoistats, inGroupVars=groupVars)        

    ## a list containins one entry for each unique combination of Group,
    ## stimulus, and cluster. The names of teh list elements correspond to
    ## the combinaiton of these as well. Each entry consists of the output
    ## of the optim call to fit the model hemodynamic function.
    complexHrfModels=complexAndFittedModels[["complexHrfModels"]]

    ## a data frams much likle roistats in that it contains columns for
    ## Group, stimulus, and cluster. It alsso contains time and value
    ## columns which correspond to the a time (in seconds) and the value
    ## of the fitted function at that time
    fittedModels=complexAndFittedModels[["fittedModels"]]
    ## the peak points of the fitted models. This is a data frame with
    ## Group, stimulus, clluster columns in addition to ipeak, x, and y
    ## cordinated of the peaks of the HRF for that combination of Group,
    ## stimulus, and cluster
    modelPeaks=complexAndFittedModels[["peaks"]]

    ## extract all the model parameters (delay, k, rise, decay) from the compexHrfModels list
    modelAttributes=attr(complexHrfModels, "split_labels")
    ## modelParams=lapply(names(complexHrfModels),
    ##     function(x) {
    ##         complexHrfModels[[x]]$par
    ##     } )
    modelParams=lapply(names(complexHrfModels),
        function(x) {
            if ( ! inherits ( complexHrfModels[[x]], "try-error") ) { 
                complexHrfModels[[x]]$par
            } else {
                ## there were 4 parameters estimated by the call to
                ## optim, that's wehere the 4 below comes from
                rep(NA, 4)
            }
        } )

    fittedMaxima=ddply(fittedModels, groupVars, .fun = function(x) { mx=max(x[, 'value'], na.rm = TRUE); names(mx)=c("max"); mx } )
    fittedMinima=ddply(fittedModels, groupVars, .fun = function(x) { mn=min(x[, 'value'], na.rm = TRUE); names(mn)=c("min"); mn } )

    parameters.roistats=cbind(modelAttributes,
        matrix(unlist(modelParams), ncol=4, byrow=TRUE),
        fittedMinima$min,
        fittedMaxima$max,
        modelPeaks$time,
        modelPeaks$value
        )
    colnames(parameters.roistats)=c(colnames(modelAttributes),  "delay", "k", "rise", "decay", "minAmplitude", "maxAmplitude", "timeOfExtrema", "extrema")
    parametersCsvFile=file.path(imageDirectory, "modelParameters.csv")
    cat (sprintf("*** Writing model parameters to %s\n", parametersCsvFile))
    write.csv(parameters.roistats, parametersCsvFile, quote=F, row.names=FALSE)


    if (inCreateIrfGraphs) {
        x.axis="time"
        y.axis="value"
        stimulus="stimulus"
        
        for ( sj in levels(inMeltedRoistats$subject ) ) {
            for ( level in levels(inMeltedRoistats$cluster) ) {

                df = subset (inMeltedRoistats, subject %in% sj & cluster %in% level)
                mp = subset (modelPeaks,       subject %in% sj & cluster %in% level)
                fm = subset (fittedModels,     subject %in% sj & cluster %in% level)                

                ## cat(sprintf("Graphing sj = %s, level = %s\n", sj, level))
                
                graph=ggplot(df, aes_string(x=x.axis, y=y.axis, color=stimulus, group=stimulus)) +
                    geom_line(position=pd, linetype=2) +
                        ##geom_line(data=fm) +
                        ##geom_point(data=mp, position=pd, size=3, shape=21, fill="black") + # 21 is filled circle                                                    
                        scale_color_brewer(name="Group:", palette="Set1") +
                            geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
                                scale_x_continuous(
                                    breaks=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2),
                                    labels=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2)) +
                                        ggtitle(paste(sj, ": ", substituteShortLabels(level), sep="")) +
                                            ylab("Impulse Response Function") + xlab("Time (sec)") +
                                                my_theme #+ theme(legend.position="bottom")

                ## the gsub on the sj variable is to handle the
                ## 169/300 subject: the embedded / results in a file
                ## creation error as it is interpreted as a directory
                ## seperator
                imageFilename=file.path(imageDirectory, sprintf("%s.%s.irf.%s.%s.%s.pdf", gsub("/", "-", sj), gsub(" +", ".", level),  task, inContrast, fLabel))
                cat(paste("*** Creating individual image named", imageFilename, "\n"))
                ggsave(imageFilename, graph)
            } ## end of for ( level in levels(roistats.summary$cluster) ) {
        } ## end of for ( sj in levels(inMeltedRoistats$subject ) ) {
    } ## end of if (inCreateIrfGraphs) {

    return(parameters.roistats)
    
} ## end of graphIndividualSubjectIrfs
        


graphGroupAverageIrfs <- function (inContrast, inMeltedRoistats, inCreateSingleIrfPdfFile=FALSE, inCreateIndividualIrfImageFiles=FALSE, inIndividualTitles=TRUE) {

    roistats.summary=summarySE(inMeltedRoistats, measurevar="value", groupvars=c("Group", "time", "cluster", "stimulus"))
    
    imageDirectory=file.path(Group.results.dir, inContrast)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }
    
    ## now collapse across subject and fit the Gamma variate model
    ## to to the average HRF for each group
    ## groupAverageRoistats=ddply(roistats, .(time, stimulus, Group),
    ##     .fun=colwise(
    ##         .fun=function (xx) {
    ##             ## the EPI timeseries were grand mean scaled to 10,000 before
    ##             ## deconvolution so the beta weights estimates for the IRF
    ##             ## will be 100 times to big, so compensate for this by
    ##             ## dividing by 100
    ##             c(mean=mean(xx) / 100)
    ##         },
    ##         ## which columns to apply the function to are listed below
    ##         c(paste("Mean_", seq(1, clusterCount), sep=""))
    ##         )
    ##     )

    ## ## now melt it 
    ## melted.groupAverageRoistats = melt(groupAverageRoistats,
    ##     id.vars=c("Group", "time", "stimulus"),
    ##     ##measure.vars="Mean_4",
    ##     measure.vars=c(paste("Mean_", seq(1, clusterCount), sep="")),
    ##     variable_name="cluster")
    ## ## replace the Mean_? column names with the corresponding
    ## ## brain region labels
    ## melted.groupAverageRoistats$cluster=factor(melted.groupAverageRoistats$cluster,
    ##     levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
    ##     labels=paste(seq(1, clusterCount), clusterWhereAmI))

    groupVars=c("Group", "stimulus", "cluster")

    complexAndFittedModels=fitGammaVariateModels(roistats.summary, inGroupVars=groupVars)        

    ## a list containins one entry for each unique combination of Group,
    ## stimulus, and cluster. The names of teh list elements correspond to
    ## the combinaiton of these as well. Each entry consists of the output
    ## of the optim call to fit the model hemodynamic function.
    complexHrfModels=complexAndFittedModels[["complexHrfModels"]]

    ## a data frams much likle roistats in that it contains columns for
    ## Group, stimulus, and cluster. It alsso contains time and value
    ## columns which correspond to the a time (in seconds) and the value
    ## of the fitted function at that time
    fittedModels=complexAndFittedModels[["fittedModels"]]
    ## the peak points of the fitted models. This is a data frame with
    ## Group, stimulus, clluster columns in addition to ipeak, x, and y
    ## cordinated of the peaks of the HRF for that combination of Group,
    ## stimulus, and cluster
    modelPeaks=complexAndFittedModels[["peaks"]]

    ## extract all the model parameters (delay, k, rise, decay) from the compexHrfModels list
    modelAttributes=attr(complexHrfModels, "split_labels")
    ## modelParams=lapply(names(complexHrfModels),
    ##     function(x) {
    ##         complexHrfModels[[x]]$par
    ##     } )
    modelParams=lapply(names(complexHrfModels),
        function(x) {
            if ( ! inherits ( complexHrfModels[[x]], "try-error") ) { 
                complexHrfModels[[x]]$par
            } else {
                ## there were 4 parameters estimated by the call to
                ## optim, that's wehere the 4 below comes from
                rep(NA, 4)
            }
        } )

    fittedMaxima=ddply(fittedModels, groupVars, .fun = function(x) { mx=max(x[, 'value'], na.rm = TRUE); names(mx)=c("max"); mx } )
    fittedMinima=ddply(fittedModels, groupVars, .fun = function(x) { mn=min(x[, 'value'], na.rm = TRUE); names(mn)=c("min"); mn } )

    parameters.roistats=cbind(modelAttributes,
        matrix(unlist(modelParams), ncol=4, byrow=TRUE),
        fittedMinima$min,
        fittedMaxima$max,
        modelPeaks$time,
        modelPeaks$value
        )
    colnames(parameters.roistats)=c(colnames(modelAttributes),  "delay", "k", "rise", "decay", "minAmplitude", "maxAmplitude", "timeOfExtrema", "extrema")
    parametersCsvFile=file.path(imageDirectory, "modelParameters.csv")
    cat (sprintf("*** Writing model parameters to %s\n", parametersCsvFile))
    write.csv(parameters.roistats, parametersCsvFile, quote=F, row.names=FALSE)
    print(parameters.roistats)

    x.axis="time"
    y.axis="value"
    group="Group"

    if (inCreateSingleIrfPdfFile) {
        singlePdfFileName=file.path(imageDirectory, sprintf("allInOne.irf.%s.%s.%s.pdf", task, inContrast, fLabel))
        
        ## Error bars represent standard error of the mean
        graph=ggplot(roistats.summary, aes_string(x=x.axis, y=y.axis, color=group, group=group)) +
            geom_errorbar(aes(ymin=value-se, ymax=value+se), width=0.5, position=pd) +
                geom_line(position=pd, linetype=2) +
                    ##geom_line(data=fittedModels) +
                    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
                        facet_grid(cluster ~ stimulus, scales="free_y") +
                            scale_color_brewer(name="Group:", palette="Set1") +
                                ## geom_point(data=modelPeaks, position=pd, size=3, shape=21, fill="black") + # 21 is filled circle                        
                                scale_x_continuous(
                                    breaks=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2),
                                    labels=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2)) +
                                        ggtitle(graphTitles[[inContrast]]) + ylab("Impulse Response Function") + xlab("Time (sec)") +
                                            my_theme + theme(legend.position="bottom")
        
        cat(paste("*** Creating single PDF file named", singlePdfFileName, "\n"))
        pdf(singlePdfFileName, paper="letter")
        print(graph)
        dev.off()
    } ## end of if (inCreateSingleIrfPdfFile) {

    if (inCreateIndividualIrfImageFiles) {
        
        for ( level in levels(roistats.summary$cluster) ) {
            
            graph=ggplot(roistats.summary[roistats.summary$cluster %in% level, ], aes_string(x=x.axis, y=y.axis, color=group, group=group)) +
                geom_errorbar(aes(ymin=value-se, ymax=value+se), position=pd) +
                    geom_line(position=pd, linetype=2) +
                        ## geom_line(data=fittedModels[fittedModels$cluster %in% level, ]) +                            
                        ## geom_point(data=modelPeaks[modelPeaks$cluster %in% level, ], position=pd, size=3, shape=21, fill="black") + # 21 is filled circle                                                    
                        facet_wrap(~stimulus, scales="free_y") +
                            scale_color_brewer(name="Group:", palette="Set1") +
                                geom_point(position=pd, size=2, shape=21, fill="white") + # 21 is filled circle
                                    xlab("Time (sec)") +
                                        ## scale_x_continuous(
                                        ##     breaks=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2),
                                        ##     labels=c("0", "", "4", "", "8", "", "12", "")) +
                                        scale_x_continuous(
                                            breaks=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2),
                                            labels=seq(min(inMeltedRoistats$time), max(inMeltedRoistats$time), by=2)) +
                                                my_theme ##+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
            
            if (inIndividualTitles) {
                graph = graph + ggtitle(substituteShortLabels(level)) + 
                    ylab("Impulse Response Function")
            } else {
                graph = graph + ylab("Impulse Response\nFunction")
            }
            
            ##+ theme(legend.position="bottom")
            
            imageFilename=file.path(imageDirectory, sprintf("%s.irf.%s.%s.%s.pdf", gsub(" +", ".", level),  task, inContrast, fLabel))
            cat(paste("*** Creating individual image named", imageFilename, "\n"))
            ggsave(imageFilename, graph, width=5.5, height=3, units="in")
        } ## end of for ( level in levels(roistats.summary$cluster) ) {
    } ## end of if (inCreateIndividualIrfImageFiles) {

    return(parameters.roistats)
    
} ## end of graphGroupAverageIrfs

makeTableString <- function(inGroup, inMean, inSe, inMin, inMax, inNaCount, inMissingData=TRUE) {
  ##  st=paste(round(inMean, 1), " / ", round(inSe, 1),    

  st=paste(round(inMean, 1), " ± ", round(inSe, 1),
    " (", round(inMin, 1), "-", round(inMax, 1), ")", ifelse(inMissingData & inNaCount > 0, paste(" [", inNaCount, "]", sep=""), ""), sep="")
  return(st)
}

runModelStatistics <- function (inIndividualModels, inCreateBoxplots=FALSE, inAlpha=0.05) {

    corrected.p.value=inAlpha / ( length(levels(inIndividualModels$Group)) * length(levels(inIndividualModels$cluster)) )
    cat (sprintf("*** Bonferroni correct p value is %0.4f\n", round(corrected.p.value, 4)))

    parameters=c("delay", "k", "rise", "decay", "timeOfExtrema", "extrema")
    mystack <- stack()
    header="Parameter,Cluster,Stimulus,Control,MDD,DF,Stat.,pValue,Signif."

    for ( param in parameters ) {
        individualModels.summarySE=summarySE(individualModels, measurevar=param, groupvars=c("stimulus", "cluster", "Group"), na.rm=TRUE)
        for ( cl in levels(inIndividualModels$cluster) ) {
            for ( stim in levels(inIndividualModels$stimulus )) {

                ctrl.string=""
                mdd.string=""
                test=""
                sm.df=subset(individualModels.summarySE, stimulus==stim & cluster==cl)
                ctrl.string=makeTableString(sm.df[2, 1], inMean=sm.df[2, "median"],  sm.df[2, "mad"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
                mdd.string =makeTableString(sm.df[1, 1], inMean=sm.df[1, "median"],  sm.df[1, "mad"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
                
                cat ("####################################################################################################\n")
                cat (sprintf("Wilcox test for cluster: %s stimulus: %s parameter: %s\n", cl, stim, param ))
                
                test <- wilcox.test(inIndividualModels[inIndividualModels$Group=="NCL" & inIndividualModels$cluster==cl, param],
                                    inIndividualModels[inIndividualModels$Group=="MDD" & inIndividualModels$cluster==cl, param])
                
                if( inherits(test <- try(wilcox.test(inIndividualModels[inIndividualModels$Group=="NCL" & inIndividualModels$cluster==cl, param],
                                                     inIndividualModels[inIndividualModels$Group=="MDD" & inIndividualModels$cluster==cl, param]),
                                         silent=FALSE),
                             "try-error") ) {
                    test <- 0
                }
                
                var.statistic=""
                var.pvalue=""
                var.parameter=""
                var.significance=""
                
                if (is.list(test)) {
                    var.statistic=round(test$statistic, 2)
                    var.parameter=ifelse(is.null(test$parameter), "NA", round(test$parameter, 2))
                    var.pvalue=round(test$p.value, 2)
                    var.significance=make.significance.indications(test$p.value)
                } ## end of if (is.list(test)) {
                print (test)
                print(sm.df)
                st=paste(param, cl, stim, ctrl.string, mdd.string,
                    var.parameter, var.statistic, var.pvalue, var.significance, sep=",")
                push(mystack, st)
            } ## end of for ( stim inlevels(melted.individual.models$stimulus )) {
        } ## end of for ( cl in levels(melted.individual.models$cluster) ) {
    } ## end of for ( param in parameters ) {
    
    cat("################################################################################\n");
    cat("Summary statistics table\n")  
    l=mystack$value()
    cat(header, "\n")
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
} ## end of runModelStatistics <- function (inNndividualModels) {


graphParameterBoxplots <- function (inIndividualModels, inContrast) {

    imageDirectory=file.path(Group.results.dir, inContrast)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }

    parameters=c("delay", "k", "rise", "decay", "timeOfExtrema")
    melted.individual.models=melt(inIndividualModels,
        id.vars=c("subject", "Group", "cluster", "stimulus"),
        measure.vars=parameters,
        variable_name="parameter")

    for (cl in levels(melted.individual.models$cluster) ) {

        df.cl=subset ( melted.individual.models, cluster %in% cl )
        
        graph=ggplot(df.cl, aes(x=parameter, y=value)) +
            geom_boxplot() +
                scale_y_log10() +
                    facet_grid ( stimulus ~ Group ) +
                        ggtitle (substituteShortLabels(cl)) +
                            xlab("Parameter") +
                                ylab("Value") +
                                    my_theme + theme(legend.position="bottom")
        quartz()
        print(graph)
        stop("Hit the stop in graphParameterBoxplots")
        ## imageFilename=file.path(imageDirectory, sprintf("boxplot.irf.%s.%s.%s.pdf", gsub(" +", ".", level),  task, inContrast, fLabel))
        ## cat(paste("*** Creating individual boxplot image named", imageFilename, "\n"))
        ## ggsave(imageFilename, graph)
    } ## end of for (cl in levels(melted.individual.models$cluster) ) {
} ## end of graphParameterBoxplots <- function (inIndividualModels) {


computeTemporalDifferences <- function(inRoistats, inLag=1) {

    retDf=ddply (inRoistats, .(subject, stimulus, Group),
        .fun=colwise(
            .fun=function (xx) {
                diff(xx, inLag)
            },
            c(paste("Mean_", seq_along(grep("Mean_", colnames(inRoistats))), sep=""))
            )
        )
    attr(retDf, "cluster.names") = attr(inRoistats, "cluster.names")

    return(retDf)
}


####################################################################################################
### Main code 
####################################################################################################

groups="mddAndCtrl"
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
  "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
  "MI RL", "MI AP", "MI IS")

## the vector of stimuli that were presented to the subject
## stimuli=c("happy", "fearful", "sad", "neutral", "happyRem", "happyNotrem", "fearfulRem", "fearfulNotrem", "sadRem", "sadNotrem", "neutralRem", "neutralNotrem")
stimuli=list(
    "fearfulVsHappy"            = c("fearful",    "happy"),
    "fearfulVsNeutral"          = c("fearful",    "Neutral"),
    "fearfulVsSad"              = c("fearful",    "sad"),
    "happyVsNeutral"            = c("happy",      "neutral"),
    "happyVsSad"                = c("happy",      "sad"),
    "neutralVsSad"              = c("neutral",    "sad"),
    "allEmotiveVsNeutral"       = c("fearful",    "happy", "sad", "neutral"),
    "happyRemVsHappyNotrem"     = c("happyRem",   "happyNotrem"),
    "fearfulRemVsFearfulNotrem" = c("fearfulRem", "fearfulNotrem"),
    "neutralRemVsNeutralNotrem" = c("neutralRem", "neutralNotrem"),
    "sadRemVsSadNotrem"         = c("sadRem",     "sadNotRem"),
    "allRemVsAllNotrem"         = c("happyRem",   "happyNotrem", "fearfulRem", "fearfulNotrem", "sadRem", "sadNotrem", "neutralRem", "neutralNotrem")
    )

## contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral", "happyRemVsHappyNotrem",
##     "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")

contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")

##contrasts=c("allEmotiveVsNeutral")

graphTitles=list(
    "fearfulVsHappy" = "Fearful vs. Happy",
    "fearfulVsNeutral" = "Fearful vs. Neutral",
    "fearfulVsSad" = "Fearful vs. Sad",
    "happyVsNeutral" = "Happy vs. Neutral",
    "happyVsSad" = "Happy vs. Sad",
    "neutralVsSad" = "Neutral vs. Sad",
    "allEmotiveVsNeutral" = "( Fearful, Happy, & Sad ) vs. Neutral",
    "happyRemVsHappyNotrem" = "Happy Rememebred vs. Happy Not Remembered",
    "fearfulRemVsFearfulNotrem" = "Fearful Remembered vs. Fearful Not Remembered",
    "neutralRemVsNeutralNotrem" = "Neutral Remembered vs. Neutral Not Remembered",
    "sadRemVsSadNotrem" = "Sad Remembered vs. Sad Not Remembered",
    "allRemVsAllNotrem" = "All Remembered vs. All Not Remembered")

root="/Volumes/PROMISEPEGASUS/yangdata/cPine"
data.root=file.path(root, "data")
admin.dir=file.path(data.root, "admin/")
config.dir=file.path(data.root, "config/")    
Group.data.dir=file.path(data.root, "Group.data/")
Group.results.dir=file.path(data.root, "Group.results/")

subjectListFilenames=c(file.path(config.dir, "control.subjectList.txt"), file.path(config.dir, "mdd.subjectList.txt"))

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

subjectList=readSubjectLists(subjectListFilenames)

my.base.size=14   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        legend.position="bottom",        
        ##panel.grid.major = element_blank(),
        ##panel.grid.minor = element_blank(),
        ## axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))

pd <- position_dodge(.2) # move them .1 to the left and right

## find out how manu cores are available
nCores=parallel::detectCores()
## tell doMC to use all the cores detected on the machine
doMC::registerDoMC(nCores)
timePointSelectors=list(
    "all" = c(0, 2, 4, 6, 8, 10, 12, 14),
    "start" = c(0, 2, 4, 6),
    "end" = c(8, 10, 12, 14)
    )
random.formula = as.formula("random = ~ 1 | subject")

my.stack = stack()
timeCourseTxtFilename=file.path(Group.results.dir, "timeCourseAnalysis.txt")


cat(paste("Check", timeCourseTxtFilename, "for time course ANOVAs\n"))

sink(timeCourseTxtFilename)
push(my.stack, "contrast;cluster;selector;term;Statistic;p-value;significance")
for (my.contrast in contrasts ) {
    
    ##my.contrast="fearfulVsSad"
    ##push(my.stack, my.contrast)
    cat("####################################################################################################\n")
    cat(sprintf("*** Graphing ROIs for the %s contrast\n", toupper(my.contrast)))
    
    statsList=createRoiStatsDataFrame(my.contrast)
    roistats=statsList[["roistats"]]
    melted.roistats=statsList[["melted.roistats"]]

    if ( 1 == 1 ) {
        cat(" *** Fitting and graphing models of the group average HRF\n")
        groupAverageModels=graphGroupAverageIrfs(my.contrast, melted.roistats, TRUE, TRUE, FALSE)

        ##cat(" *** Fitting and graphing models of the individual HRFs\n")
        ##individualModels=graphIndividualSubjectIrfs(my.contrast, melted.roistats, TRUE)
    }
    
    ## the model parameters data frame comping out of
    ## graphIndividualSubjectIrfs will not have a group column so add on to
    ## it

    ## individualModels$Group=demographics[match(individualModels$subject, demographics$ID), c("Grp")]
    ## drop the minAmplitude and maxAmplitude columns
    ## individualModels=individualModels[, ! colnames(individualModels) %in% c("minAmplitude", "maxAmplitude", "extrema")]
    ## runModelStatistics(individualModels)

    ## roistats.td=computeTemporalDifferences(roistats)

    ## add back in a temporal dimension, there are 7 rows becasue in the
    ## subtraction accomplished by the call to computeTemporalDifferences
    ## we go from 8 down to 7
    ## roistats.td$time=rep(0:6, times=length(roistats.td$stimulus)/7)

    for ( ii in seq_along(grep("Mean_", colnames(roistats))) ) {
        column.name=c(paste("Mean_", ii, sep=""))
        cluster.names = attr(roistats, "cluster.names")
        cluster.count = length(grep("Mean_", colnames(roistats)))
        
        for (selectorName in names(timePointSelectors) ) {
            cat("####################################################################################################\n")
            cat(sprintf("Running LME  (%02d of %02d) for %s = %s for the %s selector\n", ii, cluster.count, column.name, cluster.names[ii], selectorName ))
            
            selector = timePointSelectors[[selectorName]]
            
            ## only keep the foirst 4 TRs
            roistats.short=roistats[roistats$time %in% selector, ]
            
            model.formula  = as.formula(sprintf("%s ~ Group * stimulus * time", column.name))
            model.lme = lme(fixed=model.formula, random=random.formula, data = roistats.short)
            model.anova = anova(model.lme, type="marginal")

            ## the [-1] excludes the "(intercept)" term from the model
            ##significantRows=0
            for ( rname in rownames(model.anova)[-1] ) {
                ##if ( model.anova[rname, "p-value"] < 0.05 ) {
                st = sprintf("%s;%s;%s;%s;F(%0.2f, %0.2f) = %0.2f;%0.2f;%s", my.contrast, cluster.names[ii], selectorName, rname,
                    model.anova[rname, "numDF"], model.anova[rname, "denDF"], model.anova[rname, "F-value"],  model.anova[rname, "p-value"],
                    make.significance.indications( model.anova[rname, "p-value"] ) )
                push(my.stack, st)
                ##significantRows=significantRows + 1
                ##}
            } ## for ( rname in rownames(model.anova) ) {

            ## if ( significantRows == 0 ) {
            ##     push(my.stack, "No significant terms")
            ## }
            ## print(summary(model.lme))
            print(model.anova)
        } ## end of for (selectorName in names(timePointSelectors) ) {
    } ## end of for ( ii in seq_along(grep("Mean_", colnames(roistats.short))) ) {
} ## end of for (my.contrast in contrasts ) {
sink()

timeCourseCsvFilename=file.path(Group.results.dir, "timeCourseAnalysis.csv")
cat(paste("Check", timeCourseCsvFilename, "for the CSV version of the time course ANOVAs\n"))
sink(timeCourseCsvFilename)
l=my.stack$value()
for (i in 1:length(l)) {
    cat (l[[i]], "\n")
}
sink()

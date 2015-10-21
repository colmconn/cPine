rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
##source('pcor.R')

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

print.stack <- function(x,...){ 
    print(x$value()) 
} 

push <- function(stack,obj){ 
    stack$push(obj) 
} 

pop <- function(stack){ 
    stack$pop() 
}

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

  Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
    symbols   =  c("***", "**", "*", ".", " "))
  f=format(Signif)

  ## only return the first one as we're only interested in marking significant group effects
  return(f[which.pValues])
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
                           c( N    = length2(xx[,col], na.rm=na.rm),
                              mean = mean   (xx[,col], na.rm=na.rm),
                              sd   = sd     (xx[,col], na.rm=na.rm)
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

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats, inRoistats.averageFvalue=NULL, inCom=TRUE) {
  hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
  ##print(hemisphere)
  if ( inCom ) {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
  } else {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
  }

  ##cat("Locations: Volume and coordinates\n")
  ##print(locations)
  agg=aggregate(inRoistats[, grep("Mean", colnames(inRoistats))], list(inRoistats$Group), mean)
  ##cat("agg: mean for each group in each ROI\n")  
  ##print(agg)
  mns=round(t(agg[, -1]), 2)
  ##cat("mns: transposed mean for each group in each ROI\n")    
  colnames(mns)=levels(agg[,1])
  ##print(mns)

  if (! is.null(inRoistats.averageFvalue) ) {
      pubTable=cbind(locations, round(t(roistats.averageFvalue), 2), mns)
  } else {
      pubTable=cbind(locations, mns)
  }
  ##pubTable=locations
  ##cat("pubTable\n")
  ##print(pubTable)
  if (! is.null(inRoistats.averageFvalue) ) {
      colnames(pubTable)=
          c("Structure", "Hemisphere", "Volume", "CM RL", "CM AP", "CM IS", "Average F value", colnames(mns))
  } else {
      colnames(pubTable)=c("Structure", colnames(pubTable)[-1])
  }
  rownames(pubTable)=NULL
  ##print(pubTable)
  return(pubTable)
}

loadPctChangeRoiStats <- function (inContrast) {
    ## stores the names of the percentage change file names to be read    
    fnames=c()
    ## stores the names of the stimuli that went into computing the
    ## contrast. this is used to add a column to the roistats
    ## data.frame that identifies each row as comming from a
    ## particular %cs score for the corresponding stimulus
    stims=c()
    for ( stimulus in stimuli[[inContrast]] )  {
        roistats.filename=
            file.path(Group.results.dir, sprintf("roiStats.pctChange.%s.stimulus.fwhm%0.1f.%s.%s.%s.%s.contrast.%s.txt", stimulus, usedFwhm, task, groups, contrast, fLabel, groups))
        
        if ( file.exists (roistats.filename) ) {
            fnames=c(fnames, roistats.filename)
            stims=c(stims, stimulus)
        } else {
            cat(sprintf("*** No such file : %s\n", roistats.filename))
        }
    }
    roistats=do.call(rbind, lapply(fnames,
        function(X) {
            data.frame(read.table(X, header=T, sep=""))
        }
        ))
    
    roistats$stimulus=as.factor(rep(capwords(stims), each=dim(roistats)[1]/length(stims)))
    roistats=roistats[, -1]
    
    return(roistats)
}

makeGraphTitle <- function ( inContrast, inSeedCount, inSeed ) {
    
    return( 
        sprintf("%s: %d. %s Seed", 
                paste (capwords(unlist(strsplit(contrast, "Vs", fixed=TRUE))), collapse=" Vs. "),
                inSeedCount,
                inSeed)
        )
}


########################################################################################################################################################################################################
### Main code below
########################################################################################################################################################################################################


os=system('uname -s', intern=T)
if (os == "Darwin") {
    admin.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin/")
    Config.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/config"
    Ppi.seeds.data.dir=file.path(Config.data.dir, "ppi_seeds")
 
    Group.data.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data/")
    Group.results.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results")
    Group.results.ppi.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results/ppi")    
    Subject.data=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/")
} else {
    Group.data.dir=file.path("/mnt/nfs/yangdata/restingstate/data/Group.data/")
    Group.results.ppi.dir=file.path("/mnt/nfs/yangdata/restingstate/data/Group.results/")
    Subject.data=file.path("/mnt/nfs/yangdata/restingstate/data/")
}

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

##contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral", "happyRemVsHappyNotrem",
##    "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")

contrasts=c("fearfulVsHappyExtractedRightSgAcc", "fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")

##contrasts=c("fearfulVsHappyExtractedRightSgAcc")
    
##a.priori.rois=c("HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala", "HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala", "anteriorInsula.3mm", "posteriorInsula.3mm", "sgacc.left.3mm", "sgacc.right.3mm")
##a.priori.rois=c("HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala", "HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala", "sgacc.left.3mm", "sgacc.right.3mm")

#a.priori.rois=c("sgacc.left.3mm", "sgacc.right.3mm")

##a.priori.rois=c("bilateralAmygdalaeAndSgaccRois")


contrastsToSuffices=list()
for ( contrast in contrasts ) {
    seedListFile=file.path(Ppi.seeds.data.dir, paste(contrast, "seeds.txt", sep="."))
    numberOfSeeds=strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2]

    contrastsToSuffices[[contrast]] = c(
                           sapply(seq(1, numberOfSeeds),
                                  function(ii) {
                                      sprintf("roi%d.seed.%s.%s", ii, contrast, "ROIinteraction")
                                  })
                           ,
                           if (exists("a.priori.rois") ) {
                               apply(expand.grid(contrast, a.priori.rois), 1,
                                     function(x) {
                                         return(sprintf("%s.seed.%s.ROIinteraction", x[2], x[1]))
                                     })
                           }
                           )
} ## end of for ( contrast in contrasts ) {

graphTitles=list(
    "fearfulVsHappyExtractedRightSgAcc" = "Fearful vs. Happy (sgACC)",
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

####################################################################################################
## These variables control what barcharts are created

createIndividualBargraphImageFiles=TRUE
createSingleBarchartPdfFile=TRUE

####################################################################################################

groups="mddAndCtrl"
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

publicationTableFilename=file.path(Group.results.ppi.dir, "publicationTable.csv")
if (file.exists(publicationTableFilename)) {
  file.remove(publicationTableFilename)
}
cat("Writing publication table to", publicationTableFilename, "\n")
outlier.matrix.row.index=1

my.base.size=12   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        #panel.grid.major = element_line(color="white"),
        ##panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))


## ajp_theme=
##     theme(base_size =  my.base.size,
##           ##legend.position="none",
##           panel.grid.major = element_blank(),
##           panel.grid.minor = element_blank(),
##           axis.title.x=element_blank(),
##           axis.title.x = element_text(size=my.base.size, vjust=0),
##           axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
##           plot.title=element_text(size=my.base.size*1.2, vjust=1))



contrastCount=1
for ( contrast in contrasts ) {

    cat("####################################################################################################\n")
    cat(sprintf("*** Graphing ROIs for the %s contrast\n", contrast))

    seed.clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, contrast, fLabel))
    cat("*** Reading seed cluster locations from", seed.clusterLocationsFilename, "\n")
    ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
    seed.clusterWhereAmI=gsub(" +", " ", scan(file=seed.clusterLocationsFilename, what='character', sep=','))

    seedCount=1
    for ( seed in contrastsToSuffices[[contrast]] ) {
        
        cat(sprintf("*** Graphing seed: %s\n", seed))
        
        roistats.filename=file.path(Group.results.ppi.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.txt", usedFwhm, task, groups, seed, fLabel))
        roistats.averageFvalue.filename=file.path(Group.results.ppi.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.averageFvalue.txt", usedFwhm, task, groups, seed, fLabel))        
        
        if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the seed in each ROI,
            ## you do not what to graph this
            roistats=read.table(roistats.filename, header=T, sep="")
            roistats.averageFvalue=read.table(roistats.averageFvalue.filename, header=T, sep="")
            
            ## dump the first column as it's only the file name
            roistats=roistats[, -1]
            roistats.averageFvalue=roistats.averageFvalue[, -1]
            clusterCount=length(grep("Mean", colnames(roistats)))
            if (clusterCount > 0 ) {
                cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                clustersFilename=sprintf("clust.fwhm%0.1f.%s.%s.%s.%s.txt", usedFwhm, task, groups, seed, fLabel)
                cat("*** Reading", file.path(Group.results.ppi.dir, clustersFilename), "\n")
                clusters=read.table(file.path(Group.results.ppi.dir, clustersFilename))
                colnames(clusters) = clust.header
                ## this command contains the locations, as text, of the clusters and is the output of a perl script
                
                clusterLocationsFilename=file.path(Group.results.ppi.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, seed, fLabel))
                cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
                ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
                clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))
                
                
                ## this file stores the order of the subjects in each of the following BRIK files
                subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, seed, "REML.z-score", "csv", sep="."))
                cat("*** Reading", subjectOrderFilename, "\n")
                subjectOrder=read.csv(subjectOrderFilename, header=T)
                subjectOrder$subject=gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE)
                subjectOrder[subjectOrder$subject=="300", "subject"] = "169/300"
                subjectOrder$subject=as.factor(subjectOrder$subject)
                cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))
                
                roistats$subject=as.factor(as.vector(subjectOrder$subject))
                roistats$Group=demographics[match(roistats$subject, demographics$ID), c("Grp")]
                roistats$Group=drop.levels(roistats$Group)


                ## drop this column as we don't need it
                roistats=roistats[, -1]
                roistats.averageFvalue=roistats.averageFvalue[, -1]                
                ## cat ("####################################################################################################\n")
                ## cat ("### roistats START\n")
                ## print(roistats)
                ## cat ("### roistats END\n")
                ## cat ("### roistats.averageFvalue START\n")
                ## print(roistats.averageFvalue)
                ## cat ("### roistats.averageFvalue END\n")
                ## cat ("####################################################################################################\n")
                if (seed == "allRemVsAllNotrem") {
                    ## roistats$memory=ifelse(grepl("Notrem", roistats$stimulus), "Forgotten", "Remembered")
                    ## melted.roistats=melt(roistats, id.vars=c("subject", "Group", "memory"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                } else if (seed == "allEmotiveVsNeutral") {
                    ## roistats$emotion=ifelse(grepl("Neutral", roistats$stimulus), "Neutral", "Emotive")              
                    ## melted.roistats=melt(roistats, id.vars=c("subject", "Group", "emotion"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                } else {
                    melted.roistats=melt(roistats, id.vars=c("subject", "Group"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                }
                
                melted.roistats$cluster=factor(melted.roistats$cluster,
                    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                    labels=paste(seq(1, clusterCount), clusterWhereAmI))
                ##labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))
                
                ##stop("Check the data before graphics\n")
                
                xlabel="Group"
                ylabel="Functional Connectivity (Z-score)"
                graphTitle=makeGraphTitle ( contrast, seedCount, seed.clusterWhereAmI[seedCount] )

### Now make the publication table

                cat("*** Attempting to make publication table\n")
                header=graphTitle
                cat(header)
                cat("\n")
                cat(header, file=publicationTableFilename, sep="\n", append=TRUE)
                publicationTable=makePublicationTable(clusterWhereAmI, clusters, roistats, roistats.averageFvalue, inCom=TRUE)
                print(publicationTable)

                ############# stop("Check the publication table before graphics\n")
                
                if ( contrastCount == 1 && seedCount == 1 ) {
                    write.table(publicationTable, file=publicationTableFilename, quote=F, col.names=TRUE, row.names=FALSE, sep=",", append=TRUE)
                } else {
                    write.table(publicationTable, file=publicationTableFilename, quote=F, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
                }
                cat("\n", file=publicationTableFilename, append=TRUE)

### ###################################################################################################################################################
                if (createSingleBarchartPdfFile) {
                    imageDirectory=file.path(Group.results.ppi.dir, contrast)
                    if ( ! file.exists(imageDirectory) ) {
                        dir.create(imageDirectory)
                    }
                    singlePdfFileName=file.path(imageDirectory, sprintf("allInOne.fwhm%0.1f.%s.%s.%s.pdf", usedFwhm, task, seed, fLabel))
                    cat ("Using all in one pdf file", singlePdfFileName, "\n")
                    
                    pdf(singlePdfFileName, paper="letter")
                    
                    if (seed == "allRemVsAllNotrem") {
                        roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster", "memory"))
                        x.axis="memory"
                        y.axis="value"
                        group="Group"
                    } else if (seed == "allEmotiveVsNeutral") {
                        roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster", "emotion"))                  
                        x.axis="emotion"
                        y.axis="value"
                        group="Group"
                    } else {
                        roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster"))                  
                        x.axis="Group"
                        y.axis="value"
                        group="Group"
                    }
                    
                    ## Error bars represent standard error of the mean
                    graph=ggplot(roistats.summary, aes_string(x=x.axis, y=y.axis, fill=group)) +
                        geom_bar(stat="identity", position=position_dodge()) +
                            geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
                                facet_wrap(~cluster, scales="fixed") +
                                    scale_fill_brewer(name="Group:", palette="Set1") +
                                        labs(y=ylabel) +
                                            ggtitle(graphTitle) +
                                                my_theme + theme(legend.position="bottom")
                    
                    print(graph)
                    
                    dev.off()
                } ## end of if (createSingleBarchartPdfFile)
                
### ###################################################################################################################################################
                
                
                if ( createIndividualBargraphImageFiles ) {
                    imageDirectory=file.path(Group.results.ppi.dir, contrast)
                    if ( ! file.exists(imageDirectory) ) {
                        dir.create(imageDirectory)
                    }

                    followup.stack <- stack()
                    push(followup.stack, "ROI,Var,stat,pvalue,signif,pairwise differences") 

                    for ( level in levels(roistats.summary$cluster) ) {
                        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%0.1f.%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, seed, fLabel))
                        cat(paste("*** Creating", imageFilename, "\n"))
                        
                        if (seed == "allRemVsAllNotrem") {
                            roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster", "memory"))
                            x.axis="memory"
                            y.axis="value"
                            group="Group"
                        } else if (seed == "allEmotiveVsNeutral") {
                            roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster", "emotion"))                  
                            x.axis="emotion"
                            y.axis="value"
                            group="Group"
                        } else {
                            roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster"))                  
                            x.axis="Group"
                            y.axis="value"
                            group="Group"
                        }
                        
                        ## Error bars represent standard error of the mean
                        graph=ggplot(roistats.summary[roistats.summary$cluster %in% level, ], aes_string(x=x.axis, y=y.axis, fill=group)) +
                            geom_bar(stat="identity", position=position_dodge()) +
                                geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
                                    labs(title = substituteShortLabels(level), y=ylabel) +
                                        scale_fill_brewer(name="Group:", palette="Set1") +
                                            my_theme + theme(legend.position="none")
                        
                        ggsave(imageFilename, graph, width=3, height=3, units="in")

####################################################################################################
### now do the follow tests
####################################################################################################
                        ## if ( ! contrast %in% c("allRemVsAllNotrem", "allEmotiveVsNeutral") ) { 
                        ##     cat (sprintf("ANOVA for the %s cluster\n", level))
                        ##     small.df=subset(melted.roistats, subset=cluster==level)
                        ##     my.aov = aov(value ~ Group * stimulus, data=small.df)
                        ##     print(summary(my.aov))
                        ##     print(model.tables(my.aov, effect="means"))
                        ##     for ( gg in levels(small.df$Group) ) {
                        ##         if (length(levels(small.df$stimulus)) > 2 ) {
                        ##             stop ("Got more than 2 levels for small.df$stimulus:\n")
                        ##         }
                                
                        ##         stimulus.level1=levels(small.df$stimulus)[1]
                        ##         stimulus.level2=levels(small.df$stimulus)[2]
                                
                        ##         cat(sprintf("*** t test for %s (stimulus.level1) vs %s  (stimulus.level2) in %s group\n", stimulus.level1, stimulus.level2, gg))
                        ##         smaller.df=subset(small.df, subset=Group==gg)
                        ##         my.ttest=t.test(smaller.df[smaller.df$stimulus==stimulus.level1, "value"], smaller.df[smaller.df$stimulus==stimulus.level2, "value"])
                        ##         print(my.ttest)
                        ##         csvLine=sprintf("%s,%s,t(%0.2f)=%0.2f,%0.3f,%s,%s",
                        ##             level,
                        ##             gg,
                        ##             my.ttest$parameter,
                        ##             my.ttest$statistic,
                        ##             my.ttest$p.value,
                        ##             make.significance.indications(my.ttest$p.value),
                        ##             ifelse(my.ttest$estimate[1] < my.ttest$estimate[2], sprintf("%s < %s", stimulus.level1, stimulus.level2),  sprintf("%s > %s", stimulus.level1, stimulus.level2)))
                        ##         push(followup.stack, csvLine) 
                                
                        ##     }
                        ##     for ( gg in levels(small.df$stimulus) ) {
                        ##         cat(sprintf("*** t test for MDD vs NCL in %s stimulus\n", gg))
                        ##         smaller.df=subset(small.df, subset=stimulus==gg)
                        ##         my.ttest=t.test(smaller.df[smaller.df$Group=="MDD", "value"], smaller.df[smaller.df$Group=="NCL", "value"])
                        ##         print(my.ttest)
                        ##         csvLine=sprintf("%s,%s,t(%0.2f)=%0.2f,%0.3f,%s,%s",
                        ##             level,
                        ##             gg,
                        ##             my.ttest$parameter,
                        ##             my.ttest$statistic,
                        ##             my.ttest$p.value,
                        ##             make.significance.indications(my.ttest$p.value),
                        ##             ifelse(my.ttest$estimate[1] < my.ttest$estimate[2], "MDD < NCL", "MDD > NCL")
                        ##             )
                        ##         push(followup.stack, csvLine)
                        ##     }
                        ## } ## end of if ( ! contrast %in% c("allRemVsAllNotrem", "allEmotiveVsNeutral") ) {   
                        
                    } ## end of for ( level in levels(roistats.summary$cluster) )

                    l=followup.stack$value()
                    for (i in 1:length(l)) {
                        cat (l[[i]], "\n")
                    }
                    cat("\n")
                    
                } ## end of if ( createIndividualBargraphImageFiles ) 
                
            } ## end of if (clusterCount > 0 )
        } else {
            cat("### ", roistats.filename, "does not exist. Skipping\n")
        } ## end of if(file.exists(roistats.filename))
        seedCount=seedCount+1
    } ##     for ( seed in contrastsToSuffices[[contrast]]
    contrastCount=contrastCount + 1
} ## end of for (contrast in contrasts)

graphics.off()

## now remove all the NAs, if there are any
## outlier.matrix=outlier.matrix[rowSums(is.na(outlier.matrix))!=3, ]
## outlier.df=as.data.frame(outlier.matrix)
## colnames(outlier.df)=c("contrast", "cluster", "subject")

## outlier.df.count=
##   ddply(outlier.df, .variables=c("contrast", "subject"),
##   .fun=
##         function(xx, col) {
##           N=nrow(xx)
##         }
##         )
## ## Rename the "V1" column    
## outlier.df.count <- rename(outlier.df.count, c("V1"="count"))
## old.colnames=colnames(outlier.df.count)
## outlier.df.count=cbind(outlier.df.count,
##   mgd[match(outlier.df.count$subject, mgd$subject), c("Group")]
##   )
## colnames(outlier.df.count) = c(old.colnames, "Group")

## outlier.df.count.by.group=ddply(outlier.df.count, .variables=c("subject", "Group"),
##   .fun=
##   function(xx, col) {
##     N=nrow(xx)
##   }
##   )
## ## Rename the "V1" column    
## outlier.df.count.by.group <- rename(outlier.df.count.by.group, c("V1"="count"))

## ## now try to count the number of volumes censored for each subject in the table
## countCensoredVolumes <- function (inSubjectNumber) {

##   tbl=scan(file.path(Subject.data, paste(inSubjectNumber, "_A", sep=""), "functional", paste("censor_", inSubjectNumber, "_A_combined_2.1D", sep="")), what=integer())
##   censoredCount=length(tbl)-sum(tbl)

##   return(censoredCount)
## }

## ##apply(outlier.df.count.by.group$subject, 
## outlier.df.count.by.group$motionAndOutlierCensorCount=unlist(lapply(outlier.df.count.by.group$subject, countCensoredVolumes))

## #now print it out
## print(with(outlier.df.count.by.group, outlier.df.count.by.group[order(count, Group, decreasing = TRUE), ]))
## print(with(outlier.df.count.by.group, outlier.df.count.by.group[order(motionAndOutlierCensorCount, count, decreasing = TRUE), ]))

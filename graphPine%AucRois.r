rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
##source('pcor.R')

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

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats, inCom=TRUE) {
  hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
  ##print(hemisphere)
  if ( inCom ) {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
  } else {
      locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
  }
      
  ##print(locations)
  ##agg=aggregate(inRoistats[, grep("Mean", colnames(inRoistats))], list(inRoistats$Group), mean)
  ##print(agg)
  ##mns=round(t(agg[, -1]), 2)
  ##print(mns)
  ##colnames(mns)=levels(agg[,1])
  ##pubTable=cbind(locations, mns)
  pubTable=locations  
  ##print(pubTable)
  colnames(pubTable)=c("Structure", colnames(pubTable)[-1])
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
##    "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")

## contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "happyRemVsHappyNotrem",
##    "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem")

contrasts=c("allEmotions")

## graphTitles=list(
##     "fearfulVsHappy" = "Fearful vs. Happy",
##     "fearfulVsNeutral" = "Fearful vs. Neutral",
##     "fearfulVsSad" = "Fearful vs. Sad",
##     "happyVsNeutral" = "Happy vs. Neutral",
##     "happyVsSad" = "Happy vs. Sad",
##     "neutralVsSad" = "Neutral vs. Sad",
##     "allEmotiveVsNeutral" = "( Fearful, Happy, & Sad ) vs. Neutral",
##     "happyRemVsHappyNotrem" = "Happy Rememebred vs. Happy Not Remembered",
##     "fearfulRemVsFearfulNotrem" = "Fearful Remembered vs. Fearful Not Remembered",
##     "neutralRemVsNeutralNotrem" = "Neutral Remembered vs. Neutral Not Remembered",
##     "sadRemVsSadNotrem" = "Sad Remembered vs. Sad Not Remembered",
##     "allRemVsAllNotrem" = "All Remembered vs. All Not Remembered")

graphTitles=list(
    "task.F-value" = "Task Effect",
    "group.F-value" = "Group Effect"
    )

os=system('uname -s', intern=T)
if (os == "Darwin") {
    admin.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin/")
    Group.data.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data/")
    Group.results.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results/")
    Subject.data=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/")
} else {
    Group.data.dir=file.path("/mnt/nfs/yangdata/restingstate/data/Group.data/")
    Group.results.dir=file.path("/mnt/nfs/yangdata/restingstate/data/Group.results/")
    Subject.data=file.path("/mnt/nfs/yangdata/restingstate/data/")
}

####################################################################################################
## These variables control what barcharts are created

createIndividualBargraphImageFiles=TRUE
createSingleBarchartPdfFile=TRUE

####################################################################################################

groups="mddAndCtrl"
fLabels=c("task.F-value", "group.F-value")
usedFwhm=4.2
task="pine"

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

publicationTableFilename=file.path(Group.results.dir, "publicationTable.csv")
if (file.exists(publicationTableFilename)) {
  file.remove(publicationTableFilename)
}

cat("Writing publication table to", publicationTableFilename, "\n")
contrastCount=1
outlier.matrix.row.index=1

## assume the worst case scenario that every subject for every ROI is
## an outlier and alloctae the matrix accordingly.  The 60 is hard
## coded as ideally it should be length(mgd$subject) but that data
## frame is not initialized until much later in the code
outlier.matrix=matrix(NA, nrow=length(contrasts) * 60, ncol=3)

my.base.size=12   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        ##panel.grid.major = element_blank(),
        ##panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))


for ( fLabel in fLabels) {
for (contrast in contrasts ) {

  cat("####################################################################################################\n")
  cat(sprintf("*** Graphing ROIs for the %s contrast\n", toupper(contrast)))

  roistats.filename=file.path(Group.results.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.%%auc.txt", usedFwhm, task, groups, contrast, fLabel))
  if(file.exists(roistats.filename)) {
      
      ## roistats contains the average from the contrast in each ROI,
      ## you do not what to graph this
      roistats=read.table(roistats.filename, header=T, sep="")
      ## dump the first column as it's only the file name
      roistats=roistats[, -1]
      clusterCount=length(grep("Mean", colnames(roistats)))
      if (clusterCount > 0 ) {
          cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
          
          roistats=read.table(roistats.filename, header=TRUE)
          ## print (pct.roistats)
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
          
          clustersFilename=sprintf("clust.fwhm%0.1f.%s.%s.%s.%s.%%auc.txt", usedFwhm, task, groups, contrast, fLabel)
          cat("*** Reading", file.path(Group.results.dir, clustersFilename), "\n")
          clusters=read.table(file.path(Group.results.dir, clustersFilename))
          colnames(clusters) = clust.header
          ## this command contains the locations, as text, of the clusters and is the output of a perl script

          clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.%%auc.csv", usedFwhm, task, groups, contrast, fLabel))
          cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
          ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
          clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))
          
          ## this file stores the order of the subjects in each of the following BRIK files
          subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, contrast, "%auc.csv", sep="."))
          cat("*** Reading", subjectOrderFilename, "\n")
          subjectOrder=read.csv(subjectOrderFilename, header=T)
          subjectOrder$subject=gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE)
          subjectOrder[subjectOrder$subject=="300", "subject"] = "169/300"
          subjectOrder$subject=as.factor(subjectOrder$subject)
          cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))
          
          roistats$subject=as.factor(as.vector(subjectOrder$subject))
          roistats$Group=demographics[match(roistats$subject, demographics$ID), c("Grp")]
          roistats$Group=drop.levels(roistats$Group)
          roistats$stimulus=rep(c("Happy", "Fearful", "Neutral", "Sad"), each=dim(roistats)[1]/4)

### Now make the publication table
          
          cat("*** Attempting to make publication table\n")
          header=sprintf("%s", paste (capwords(unlist(strsplit(contrast, "Vs", fixed=TRUE))), collapse=" Vs. "))
          cat(header)
          cat("\n")
          cat(header, file=publicationTableFilename, sep="\n", append=TRUE)

          publicationTable=makePublicationTable(clusterWhereAmI, clusters, roistats, inCom=TRUE)
          print(publicationTable)
          if ( contrastCount == 1 ) {
              write.table(publicationTable, file=publicationTableFilename, quote=F, col.names=TRUE, row.names=FALSE, sep=",", append=TRUE)
          } else {
              write.table(publicationTable, file=publicationTableFilename, quote=F, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
          }
          cat("\n", file=publicationTableFilename, append=TRUE)          

          ## drop this column as we don't need it
          roistats=roistats[, -1]
          ## cat ("####################################################################################################\n")
          ## cat ("### pct.roistats START\n")
          ## print(pct.roistats)
          ## cat ("### pct.roistats END\n")
          ## cat ("####################################################################################################\n")
          melted.roistats=melt(roistats, id.vars=c("subject", "Group", "stimulus"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
          
          melted.roistats$cluster=factor(melted.roistats$cluster,
              levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
              labels=paste(seq(1, clusterCount), clusterWhereAmI))
          ##labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))
          
          ## stop("Check the data before graphics\n")

          xlabel="Group"
          ylabel="Mean Percentage Change"
          ##graphTitle=sprintf("Controls vs MDD\n%s", toupper(contrast))

## ###################################################################################################################################################
          if (createSingleBarchartPdfFile) {
              singlePdfFileName=file.path(Group.results.dir, sprintf("allInOne.fwhm%0.1f.%s.%s.%s.pdf", usedFwhm, task, contrast, fLabel))
              cat ("Using all in one pdf file", singlePdfFileName, "\n")
              
              pdf(singlePdfFileName, paper="letter")

              if (fLabel == "task.F-value") {
                  roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("stimulus", "cluster"))
                  x.axis="stimulus"
                  y.axis="value"
                  scale.name="Stimulus:"
              } else if (fLabel == "group.F-value") {
                  roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster"))                  
                  x.axis="Group"
                  y.axis="value"
                  scale.name="Group:"
              }
              
              ## Error bars represent standard error of the mean
              graph=ggplot(roistats.summary, aes_string(x=x.axis, y=y.axis, fill=x.axis)) +
                  geom_bar(stat="identity", position=position_dodge()) +
                      geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
                          facet_wrap(~cluster, scales="free_y") +
                              scale_fill_brewer(name=scale.name, palette="Set1") +
                                  labs(title = graphTitles[[fLabel]], y=ylabel) +
                                      my_theme + theme(legend.position="bottom")
              
              print(graph)
              dev.off()
          } ## end of if (createSingleBarchartPdfFile)
          
          
## ###################################################################################################################################################

      
      if ( createIndividualBargraphImageFiles ) {
          imageDirectory=file.path(Group.results.dir, contrast)
          if ( ! file.exists(imageDirectory) ) {
              dir.create(imageDirectory)
          }
          
          for ( level in levels(roistats.summary$cluster) ) {
              imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%0.1f.%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, contrast, fLabel))
              cat(paste("*** Creating", imageFilename, "\n"))

              if (fLabel == "task.F-value") {
                  cat (sprintf("*** Pairwise t tests for the %s cluster\n", level))
                  print(pairwise.t.test(melted.roistats[melted.roistats$cluster==level, "value"],
                                        melted.roistats[melted.roistats$cluster==level, "stimulus"]))
                  
                  roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("stimulus", "cluster"))
                  x.axis="stimulus"
                  y.axis="value"
                  scale.name="Stimulus:"
               } else if (fLabel == "group.F-value") {
                  roistats.summary=summarySE(melted.roistats, measurevar="value", groupvars=c("Group", "cluster"))                  
                  x.axis="Group"
                  y.axis="value"
                  scale.name="Group:"
              }

              ## Error bars represent standard error of the mean
              graph=ggplot(roistats.summary[roistats.summary$cluster %in% level, ], aes_string(x=x.axis, y=y.axis, fill=x.axis)) +
                  geom_bar(stat="identity", position=position_dodge()) +
                      geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
                          labs(title = substituteShortLabels(level), y=ylabel) +
                              scale_fill_brewer(name=scale.name, palette="Set1") +
                                  my_theme + theme(legend.position="bottom")
              
              ggsave(imageFilename, graph)
              
          } ## end of for ( level in levels(roistats.summary$cluster) )
      } ## end of if ( createIndividualBargraphImageFiles ) 
          
      } ## end of if (clusterCount > 0 )
  } else {
      cat("### ", roistats.filename, "does not exist. Skipping\n")
  } ## end of if(file.exists(roistats.filename))
  contrastCount=contrastCount+1
} ## end of for (contrast in contrasts)
}
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

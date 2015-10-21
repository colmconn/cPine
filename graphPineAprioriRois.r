rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
library(grid)
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

substituteShortLabels <- function(inLevel, inGyrusLabelOnNewLine=FALSE) {
    returnSubstitutedLabel = gsub("[0-9]+ ", "", gsub("Inf", "Inferior", gsub("Sup", "Superior", 
        gsub("^R", "Right", gsub("^L", "Left", inLevel, fixed=FALSE), fixed=FALSE),  fixed=TRUE), fixed=TRUE))

    if (inGyrusLabelOnNewLine) {
        returnSubstitutedLabel=gsub(" Gy", "\nGyrus", returnSubstitutedLabel, fixed=TRUE)        
    } else {
        returnSubstitutedLabel=gsub("Gy", "Gyrus", returnSubstitutedLabel, fixed=TRUE)
    }
    
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

## the following very useful functions are from
## http://stackoverflow.com/questions/13297155/add-floating-axis-labels-in-facet-wrap-plot
facetAdjust <- function(x, pos = c("up", "down"))
{
  pos <- match.arg(pos)
  p <- ggplot_build(x)
  gtable <- ggplot_gtable(p); dev.off()
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  n <- space - panels
  if(panels != space){
    idx <- (space - ncol - n + 1):(space - ncol)
    gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  class(gtable) <- c("facetAdjust", "gtable", "ggplot"); gtable
}

print.facetAdjust <- function(x, newpage = is.null(vp), vp = NULL) {
  if(newpage)
    grid.newpage()
  if(is.null(vp)){
    grid.draw(x)
  } else {
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(x)
    upViewport()
  }
  invisible(x)
}

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats=NULL, inMeanFValues=NULL, inCom=TRUE) {
    hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
    ## cat("Hemisphere\n")
    ## print(hemisphere)
    if ( inCom ) {
        locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
    } else {
        locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
    }
    ## cat("Locations\n")
    ## print(locations)
    ##agg=aggregate(inRoistats[, grep("Mean", colnames(inRoistats))], list(inRoistats$Group), mean)
    ##print(agg)
    ##mns=round(t(agg[, -1]), 2)
    ##print(mns)
    ##colnames(mns)=levels(agg[,1])
    ##pubTable=cbind(locations, mns)
    pubTable=locations
    ## cat("pubTable\n")
    ## print(pubTable)
    ## print("pubTable colnames\n")
    ## print(colnames(pubTable))

    if (! is.null(inMeanFValues)) {
        pubTable = cbind(pubTable, t(inMeanFValues[, grep("Mean", colnames(inMeanFValues))]))
        ## cat("pubTable with mean f value\n")
        ## print(pubTable)
    }

    rownames(pubTable)=NULL
    ## print(pubTable)
    ## print(colnames(pubTable))
    ## print(colnames(pubTable)[-1])

    ## print(dim(pubTable))
    ## print(colnames(pubTable)[-c(1, dim(pubTable)[2])])

    if ( ! is.null(inMeanFValues) ) {
        ## print(dim(pubTable))
        ## print(length(c("Structure", colnames(pubTable)[-1], "Mean F-value")))


        ## print(c("Structure", colnames(pubTable)[-c(1, dim(pubTable)[2])], "Mean F-value"))

        colnames(pubTable)=c("Structure", colnames(pubTable)[-c(1, dim(pubTable)[2])], "Mean F-value")

    } else {
        colnames(pubTable)=c("Structure", colnames(pubTable)[-1])
    }
        
    ##cat("Final Table\n")
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
            file.path(Group.results.dir, sprintf("roiStats.pctChange.%s.stimulus.fwhm%0.1f.%s.%s.%s.%s.%s.contrast.%s.txt", stimulus, usedFwhm, task, groups, contrast, fLabel, mask, groups))
        
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

contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral", "happyRemVsHappyNotrem",
    "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")

## contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "happyRemVsHappyNotrem",
##    "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem")

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

os=system('uname -s', intern=T)
if (os == "Darwin") {
    admin.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin/")
    rois.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/rois/")
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
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))


my.base.size=14
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        legend.position="bottom",        
        ##panel.grid.major = element_blank(),
        ##panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))

##masks=c("HarvardOxford-sub-maxprob-thr25-3mm-bilateral-amygdala",  "sgacc.bilateral.3mm+tlrc.HEAD")

##masks=c("bilateralAmygdalaeAndSgaccRois")
##masks=c("apriori.rois.noVmpfc.3mm",  "apriori.rois.withVmpfc.3mm")
masks=c("apriori.rois.withVmpfc.3mm")
for ( mask in masks ) {

    publicationTableFilename=file.path(Group.results.dir, sprintf("publicationTable.%s.csv", mask))
    if (file.exists(publicationTableFilename)) {
        file.remove(publicationTableFilename)
    }
    cat("*** Writing publication table to", publicationTableFilename, "\n")
    contrastCount=1
    outlier.matrix.row.index=1

    ## assume the worst case scenario that every subject for every ROI is
    ## an outlier and alloctae the matrix accordingly.  The 60 is hard
    ## coded as ideally it should be length(mgd$subject) but that data
    ## frame is not initialized until much later in the code
    outlier.matrix=matrix(NA, nrow=length(contrasts) * 60, ncol=3)

    results.file=file.path(Group.results.dir, sprintf("%s.RoiGraphingAndFollowupResults.txt", mask))
    cat(sprintf("*** Check %s for errors and results of follow up tests\n", results.file))
    sink(results.file, append=FALSE)
    
    followup.stack <- stack()
    push(followup.stack, "Contrast,ROI,Var,stat,pvalue,signif,pairwise differences") 
    
    for ( contrast in contrasts ) {

        cat("####################################################################################################\n")
        cat(sprintf("*** Graphing ROIs for the %s contrast\n", toupper(contrast)))


        roistats.filename=file.path(Group.results.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.%s.txt", usedFwhm, task, groups, contrast, fLabel, mask))
        
        if(file.exists(roistats.filename)) {
            
            ## roistats contains teh avergae from the contrast in each ROI,
            ## you do not what to graph this
            roistats=read.table(roistats.filename, header=T, sep="")
            ## dump the first column as it's only the file name
            roistats=roistats[, -1]
            clusterCount=length(grep("Mean", colnames(roistats)))

            ## now load the roistats file with mean F value
            roistats.mean.f.value.filename=file.path(Group.results.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.%s.meanFValue.txt", usedFwhm, task, groups, contrast, fLabel, mask))
            roistats.mean.f.value=read.table(roistats.mean.f.value.filename, header=T, sep="")
            roistats.mean.f.value=roistats.mean.f.value[, -1]
            
            if (clusterCount > 0 ) {
                cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
                
                pct.roistats=loadPctChangeRoiStats(contrast)
                ## print (pct.roistats)
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order

                clustersFilename=sprintf("clust.fwhm%0.1f.%s.%s.%s.%s.%s.txt", usedFwhm, task, groups, contrast, fLabel, mask)
                cat("*** Reading", file.path(Group.results.dir, clustersFilename), "\n")
                clusters=read.table(file.path(Group.results.dir, clustersFilename))
                colnames(clusters) = clust.header
                ## this command contains the locations, as text, of the clusters and is the output of a perl script

                clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.%s.csv", usedFwhm, task, groups, contrast, fLabel, mask))        
                cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
                ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
                clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))

                ## Now make sure it says Amygdala not parahippocampas gyrus
                clusterWhereAmI=gsub("Parahippocampal Gy", "Amygdala", clusterWhereAmI, fixed=TRUE)

                ## Now make sure it says sgACC not anterior cingulate gyrus
                clusterWhereAmI=gsub("Anterior Cingulate", "sgACC", clusterWhereAmI, fixed=TRUE)

                ## this file stores the order of the subjects in each of the following BRIK files
                subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, contrast, "REML.csv", sep="."))
                cat("*** Reading", subjectOrderFilename, "\n")
                subjectOrder=read.csv(subjectOrderFilename, header=T)
                subjectOrder$subject=gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE)
                subjectOrder[subjectOrder$subject=="300", "subject"] = "169/300"
                subjectOrder$subject=as.factor(subjectOrder$subject)
                cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))
                
                pct.roistats$subject=as.factor(rep(as.vector(subjectOrder$subject), times=length(levels(pct.roistats$stimulus))))
                pct.roistats$Group=demographics[match(pct.roistats$subject, demographics$ID), c("Grp")]
                pct.roistats$Group=drop.levels(pct.roistats$Group)

                roistats$subject=as.factor(as.vector(subjectOrder$subject))
                roistats$Group=demographics[match(roistats$subject, demographics$ID), c("Grp")]
                roistats$Group=drop.levels(roistats$Group)

### Now make the publication table
                
                cat("*** Attempting to make publication table\n")
                header=sprintf("%s", paste (capwords(unlist(strsplit(contrast, "Vs", fixed=TRUE))), collapse=" Vs. "))
                cat(header)
                cat("\n")
                cat(header, file=publicationTableFilename, sep="\n", append=TRUE)

                publicationTable=makePublicationTable(clusterWhereAmI, clusters, roistats, roistats.mean.f.value,inCom=TRUE)
                print(publicationTable)
                if ( contrastCount == 1 ) {
                    write.table(publicationTable, file=publicationTableFilename, quote=F, col.names=TRUE, row.names=FALSE, sep=",", append=TRUE)
                } else {
                    write.table(publicationTable, file=publicationTableFilename, quote=F, col.names=FALSE, row.names=FALSE, sep=",", append=TRUE)
                }
                cat("\n", file=publicationTableFilename, append=TRUE)          
                
                ## drop this column as we don't need it
                pct.roistats=pct.roistats[, -1]
                ## cat ("####################################################################################################\n")
                ## cat ("### pct.roistats START\n")
                ## print(pct.roistats)
                ## cat ("### pct.roistats END\n")
                ## cat ("####################################################################################################\n")
                if (contrast == "allRemVsAllNotrem") {
                    pct.roistats$memory=ifelse(grepl("Notrem", pct.roistats$stimulus), "Forgotten", "Remembered")
                    melted.pct.roistats=melt(pct.roistats, id.vars=c("subject", "Group", "memory"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                } else if (contrast == "allEmotiveVsNeutral") {
                    pct.roistats$emotion=ifelse(grepl("Neutral", pct.roistats$stimulus), "Neutral", "Emotive")              
                    melted.pct.roistats=melt(pct.roistats, id.vars=c("subject", "Group", "emotion"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                } else {
                    melted.pct.roistats=melt(pct.roistats, id.vars=c("subject", "Group", "stimulus"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                }

                melted.pct.roistats$cluster=factor(melted.pct.roistats$cluster,
                    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                    labels=paste(seq(1, clusterCount), clusterWhereAmI))
                ##labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))

                ## stop("Check the data before graphics\n")

                xlabel="Group"
                ylabel="Mean Percentage\nSignal Change"
                ##ylabel="Mean % Change"                
                ##graphTitle=sprintf("Controls vs MDD\n%s", toupper(contrast))

                ## End of the book keeping code
                ## ###################################################################################################################################################
                if (createSingleBarchartPdfFile) {
                    singlePdfFileName=file.path(Group.results.dir, sprintf("allInOne.fwhm%0.1f.%s.%s.%s.%s.pdf", usedFwhm, task, contrast, fLabel, mask))
                    cat ("Using all in one pdf file", singlePdfFileName, "\n")
                    
                    ##pdf(singlePdfFileName, paper="letter")

                    if (contrast == "allRemVsAllNotrem") {
                        roistats.summary=summarySE(melted.pct.roistats, measurevar="value", groupvars=c("Group", "cluster", "memory"))
                        x.axis="memory"
                        y.axis="value"
                        group="Group"
                    } else if (contrast == "allEmotiveVsNeutral") {
                        roistats.summary=summarySE(melted.pct.roistats, measurevar="value", groupvars=c("Group", "cluster", "emotion"))                  
                        x.axis="emotion"
                        y.axis="value"
                        group="Group"
                    } else {
                        roistats.summary=summarySE(melted.pct.roistats, measurevar="value", groupvars=c("Group", "cluster", "stimulus"))                  
                        x.axis="stimulus"
                        y.axis="value"
                        group="Group"
                    }
                    
                    ## ## Error bars represent standard error of the mean
                    ## graph=ggplot(roistats.summary, aes_string(x=x.axis, y=y.axis, fill=group)) +
                    ##     geom_bar(stat="identity", position=position_dodge()) +
                    ##         geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
                    ##             facet_wrap(~cluster, scales="free_y") +
                    ##                 scale_fill_brewer(palette="Set1") +
                    ##                     labs(title = graphTitles[[contrast]], y=ylabel) +
                    ##                         my_theme + theme(legend.position="bottom")

                    ## START of code make scatter plots with error
                    ## bars and lines connecting the means of each of
                    ## the two stimulus types
                    my_dodge=position_dodge(.2)
                    melted.pct.roistat.means.df=ddply(melted.pct.roistats, .(Group, stimulus, cluster), summarise, mean=mean(value))
                    mn="mean"
                    graph=ggplot(melted.pct.roistats, aes_string(x=x.axis, y=y.axis, color=group, shape=group)) +
                        geom_jitter(alpha=0.4) +
                            geom_errorbar(data=roistats.summary, aes(ymin=value-se, ymax=value+se),  position=my_dodge) +
                                geom_line (aes_string(x=x.axis, y=mn, color=group, group=group), melted.pct.roistat.means.df, position=my_dodge) +    
                                    geom_point(aes_string(x=x.axis, y=mn, shape=group), melted.pct.roistat.means.df, position=my_dodge) +
                                        facet_wrap(~cluster, scales="free_y") +
                                            scale_color_brewer(palette="Set1") +
                                                labs(title = graphTitles[[contrast]], y=ylabel) +
                                                    my_theme + theme(legend.position="bottom")
                    graph=facetAdjust(graph)
                    ## END of code make scatter plots with error
                    ## bars and lines connecting the means of each of
                    ## the two stimulus types
                    
                    ##print(graph)
                    ggsave(singlePdfFileName, graph)
                    
                    dev.off()
                } ## end of if (createSingleBarchartPdfFile)
                
                
                ## ###################################################################################################################################################

                
                if ( createIndividualBargraphImageFiles ) {
                    imageDirectory=file.path(Group.results.dir, contrast)
                    if ( ! file.exists(imageDirectory) ) {
                        dir.create(imageDirectory)
                    }

                    
                    for ( level in levels(roistats.summary$cluster) ) {
                        imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%0.1f.%s.%s.%s.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, contrast, fLabel, mask))
                        cat(paste("*** Creating", imageFilename, "\n"))

                        if (contrast == "allRemVsAllNotrem") {
                            roistats.summary=summarySE(melted.pct.roistats, measurevar="value", groupvars=c("Group", "cluster", "memory"))
                            x.axis="memory"
                            y.axis="value"
                            group="Group"
                        } else if (contrast == "allEmotiveVsNeutral") {
                            roistats.summary=summarySE(melted.pct.roistats, measurevar="value", groupvars=c("Group", "cluster", "emotion"))                  
                            x.axis="emotion"
                            y.axis="value"
                            group="Group"
                        } else {
                            roistats.summary=summarySE(melted.pct.roistats, measurevar="value", groupvars=c("Group", "cluster", "stimulus"))                  
                            x.axis="stimulus"
                            y.axis="value"
                            group="Group"
                        }

                        ## Error bars represent standard error of the mean
                        graph=ggplot(roistats.summary[roistats.summary$cluster %in% level, ], aes_string(x=x.axis, y=y.axis, fill=group)) +
                            geom_bar(stat="identity", position=position_dodge()) +
                                geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
                                    ##labs(title = substituteShortLabels(level, inGyrusLabelOnNewLine=TRUE), y=ylabel) +
                                    labs(y=ylabel) +                                        
                                        scale_fill_brewer(name=paste(group, ":", sep=""), palette="Set1") +
                                            my_theme + theme(legend.position="none")
                        
                        ggsave(imageFilename, graph, width=3, height=3, units="in")


####################################################################################################
### now do the follow tests
####################################################################################################
                        if ( ! contrast %in% c("allRemVsAllNotrem", "allEmotiveVsNeutral") ) { 
                            cat (sprintf("ANOVA for the %s cluster\n", level))
                            small.df=subset(melted.pct.roistats, subset=cluster==level)
                            my.aov = aov(value ~ Group * stimulus, data=small.df)
                            print(summary(my.aov))
                            print(model.tables(my.aov, effect="means"))
                            for ( gg in levels(small.df$Group) ) {
                                if (length(levels(small.df$stimulus)) > 2 ) {
                                    stop ("Got more than 2 levels for small.df$stimulus:\n")
                                }

                                stimulus.level1=levels(small.df$stimulus)[1]
                                stimulus.level2=levels(small.df$stimulus)[2]
                                
                                cat(sprintf("*** t test for %s (stimulus.level1) vs %s  (stimulus.level2) in %s group\n", stimulus.level1, stimulus.level2, gg))
                                smaller.df=subset(small.df, subset=Group==gg)
                                my.ttest=t.test(smaller.df[smaller.df$stimulus==stimulus.level1, "value"], smaller.df[smaller.df$stimulus==stimulus.level2, "value"])
                                print(my.ttest)
                                csvLine=sprintf("%s,%s,%s,t(%0.2f)=%0.2f,%0.3f,%s,%s",
                                    contrast,
                                    level,
                                    gg,
                                    my.ttest$parameter,
                                    my.ttest$statistic,
                                    my.ttest$p.value,
                                    make.significance.indications(my.ttest$p.value),
                                    ifelse(my.ttest$estimate[1] < my.ttest$estimate[2], sprintf("%s < %s", stimulus.level1, stimulus.level2),  sprintf("%s > %s", stimulus.level1, stimulus.level2)))
                                push(followup.stack, csvLine) 
                                
                            }
                            for ( gg in levels(small.df$stimulus) ) {
                                cat(sprintf("*** t test for MDD vs NCL in %s stimulus\n", gg))
                                smaller.df=subset(small.df, subset=stimulus==gg)
                                my.ttest=t.test(smaller.df[smaller.df$Group=="MDD", "value"], smaller.df[smaller.df$Group=="NCL", "value"])
                                print(my.ttest)
                                csvLine=sprintf("%s,%s,%s,t(%0.2f)=%0.2f,%0.3f,%s,%s",
                                    contrast,
                                    level,
                                    gg,
                                    my.ttest$parameter,
                                    my.ttest$statistic,
                                    my.ttest$p.value,
                                    make.significance.indications(my.ttest$p.value),
                                    ifelse(my.ttest$estimate[1] < my.ttest$estimate[2], "MDD < NCL", "MDD > NCL")
                                    )
                                push(followup.stack, csvLine)
                            }
                        }   
                    } ## end of for ( level in levels(roistats.summary$cluster) )
                } ## end of if ( createIndividualBargraphImageFiles ) 
                
            } ## end of if (clusterCount > 0 )
        } else {
            cat("### ", roistats.filename, "does not exist. Skipping\n")
        } ## end of if(file.exists(roistats.filename))
        contrastCount=contrastCount+1
    } ## end of for (contrast in contrasts)
    l=followup.stack$value()
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
    cat("\n")
    
    sink()
} ## end of for (mask in masks) {
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
